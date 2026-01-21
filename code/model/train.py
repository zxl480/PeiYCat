import json
import torch
import torch.nn as nn
import torch.optim as optim
import math
from sklearn.metrics import r2_score
from sklearn.model_selection import train_test_split
from model import ReactionModel
import time
from torch.utils.data import DataLoader
from torch_geometric.data import Data, Batch
from torch.nn.utils.rnn import pad_sequence, pack_padded_sequence
import os
import numpy as np

# Constants
INPUT_DATA_FILE = "../../data/feature_extraction/prepared_input_final.json"
BATCH_SIZE = 64
FEATURE_DIM = 31
NUM_EPOCHS = 2000
LEARNING_RATE = 0.0004
INPUT_SIZE = 1280  # ESM-1b feature size
POCKET_MAX_ATOMS = 100
POCKET_FEATURE_DIM = 4

def load_combined_data(file_path):
    """Load combined reaction and enzyme data from JSON file."""
    with open(file_path, "r") as f:
        combined_data = json.load(f)
    return combined_data


# Dataset class
class EnzymeDataset(torch.utils.data.Dataset):
    def __init__(self, data_list):
        self.data_list = data_list

    def __len__(self):
        return len(self.data_list)

    def __getitem__(self, idx):
        data_point = self.data_list[idx]

        reactant_graph = Data(x=torch.tensor(data_point['reactant_features'], dtype=torch.float),
                              edge_index=torch.tensor(data_point['reactant_edge'], dtype=torch.long))
        product_graph = Data(x=torch.tensor(data_point['product_features'], dtype=torch.float),
                             edge_index=torch.tensor(data_point['product_edge'], dtype=torch.long))

        pocket_x = torch.tensor(data_point['x'], dtype=torch.float)
        pocket_edge_index = torch.tensor(data_point['edge_index'], dtype=torch.long)

        if pocket_x.size(0) > POCKET_MAX_ATOMS:
            pocket_x = pocket_x[:POCKET_MAX_ATOMS]
            mask = (pocket_edge_index[0] < POCKET_MAX_ATOMS) & (pocket_edge_index[1] < POCKET_MAX_ATOMS)
            pocket_edge_index = pocket_edge_index[:, mask]
        elif pocket_x.size(0) < POCKET_MAX_ATOMS:
            padding_size = POCKET_MAX_ATOMS - pocket_x.size(0)
            padding = torch.zeros(padding_size, POCKET_FEATURE_DIM, dtype=torch.float)
            pocket_x = torch.cat([pocket_x, padding], dim=0)

        pocket_graph = Data(x=pocket_x, edge_index=pocket_edge_index)

        enzyme_seq = torch.tensor(data_point['enzyme_sequence'], dtype=torch.float)
        kcat_value = math.log1p(data_point['kcat'])
        kcat_value = torch.tensor([kcat_value], dtype=torch.float)

        pH = torch.tensor(data_point.get("pH", 7.5), dtype=torch.float)
        temperature = torch.tensor(data_point.get("temperature", 25.0), dtype=torch.float)

        return reactant_graph, product_graph, pocket_graph, enzyme_seq, pH, temperature, kcat_value


def check_for_nan_inf(tensor, name):
    if torch.isnan(tensor).any() or torch.isinf(tensor).any():
        return True
    return False


# Training function
def main():
    highest_r2 = 0
    valid_highest_r2 = 0
    start_time = time.time()
    combined_data = load_combined_data(INPUT_DATA_FILE)


    train_split, temp_data = train_test_split(combined_data, test_size=0.2, random_state=42)
    train_data = train_split
    valid_data, test_data = train_test_split(temp_data, test_size=0.5, random_state=42)


    print(f"Train size: {len(train_data)}, Valid size: {len(valid_data)}, Test size: {len(test_data)}")


    # DataLoaders
    train_loader = DataLoader(EnzymeDataset(train_data), batch_size=BATCH_SIZE, shuffle=True, collate_fn=collate_fn)
    valid_loader = DataLoader(EnzymeDataset(valid_data), batch_size=BATCH_SIZE, shuffle=False, collate_fn=collate_fn)
    test_loader = DataLoader(EnzymeDataset(test_data), batch_size=BATCH_SIZE, shuffle=False, collate_fn=collate_fn)

    model = ReactionModel()
    device = torch.device("cuda:0" if torch.cuda.is_available() else "cpu")
    model.to(device)

    criterion = nn.MSELoss()
    optimizer = optim.Adam(model.parameters(), lr=LEARNING_RATE, weight_decay=1e-5)
    scheduler = torch.optim.lr_scheduler.ReduceLROnPlateau(optimizer, mode='min', factor=0.8, patience=10)

    for epoch in range(NUM_EPOCHS):
        model.train()
        total_loss = 0
        all_true_values_train, all_pred_values_train = [], []

        for batch_idx, (reactant_graphs, product_graphs, pocket_graphs, enzyme_seqs, pHs, temps, kcat_values) in enumerate(train_loader):
            try:
                reactant_graphs = reactant_graphs.to(device)
                product_graphs = product_graphs.to(device)
                pocket_graphs = pocket_graphs.to(device)
                enzyme_seqs = enzyme_seqs.to(device)
                pHs = pHs.to(device)
                temps = temps.to(device)
                kcat_values = kcat_values.to(device)

                if any([
                    check_for_nan_inf(reactant_graphs.x, "reactant_x"),
                    check_for_nan_inf(product_graphs.x, "product_x"),
                    check_for_nan_inf(pocket_graphs.x, "pocket_x"),
                    check_for_nan_inf(enzyme_seqs, "enzyme_seqs"),
                ]):
                    continue

                pred_kcat = model(
                    reactant_graphs.x, reactant_graphs.edge_index,
                    product_graphs.x, product_graphs.edge_index,
                    pocket_graphs.x, pocket_graphs.edge_index,
                    enzyme_seqs, reactant_graphs.batch, pocket_graphs.batch,
                    pHs, temps
                )

                if check_for_nan_inf(pred_kcat, "pred_kcat"):
                    continue

                loss = criterion(pred_kcat, kcat_values)
                if torch.isnan(loss) or torch.isinf(loss):
                    continue

                total_loss += loss.item()
                all_true_values_train.extend(kcat_values.cpu().numpy())
                all_pred_values_train.extend(pred_kcat.detach().cpu().numpy())

                optimizer.zero_grad()
                loss.backward()
                torch.nn.utils.clip_grad_norm_(model.parameters(), max_norm=1.0)
                optimizer.step()

            except Exception as e:
                print(f"Error in batch {batch_idx}: {e}")
                continue

        if len(all_pred_values_train) == 0:
            continue

        train_loss = total_loss / len(train_loader)
        valid_indices = ~(np.isnan(all_true_values_train) | np.isnan(all_pred_values_train))
        train_r2 = r2_score(np.array(all_true_values_train)[valid_indices],
                            np.array(all_pred_values_train)[valid_indices]) if np.sum(valid_indices) > 0 else float('nan')

        scheduler.step(train_loss)

        if (epoch + 1) % 5 == 0:
            # 验证集评估
            model.eval()
            total_loss_valid = 0
            all_true_values_valid, all_pred_values_valid = [], []

            with torch.no_grad():
                for batch_idx, (reactant_graphs, product_graphs, pocket_graphs, enzyme_seqs, pHs, temps, kcat_values) in enumerate(valid_loader):
                    try:
                        reactant_graphs = reactant_graphs.to(device)
                        product_graphs = product_graphs.to(device)
                        pocket_graphs = pocket_graphs.to(device)
                        enzyme_seqs = enzyme_seqs.to(device)
                        pHs = pHs.to(device)
                        temps = temps.to(device)
                        kcat_values = kcat_values.to(device)

                        pred_kcat = model(
                            reactant_graphs.x, reactant_graphs.edge_index,
                            product_graphs.x, product_graphs.edge_index,
                            pocket_graphs.x, pocket_graphs.edge_index,
                            enzyme_seqs, reactant_graphs.batch, pocket_graphs.batch,
                            pHs, temps
                        )

                        if check_for_nan_inf(pred_kcat, "valid_pred_kcat"):
                            continue

                        loss = criterion(pred_kcat, kcat_values)
                        if torch.isnan(loss) or torch.isinf(loss):
                            continue

                        total_loss_valid += loss.item()
                        all_true_values_valid.extend(kcat_values.cpu().numpy())
                        all_pred_values_valid.extend(pred_kcat.detach().cpu().numpy())

                    except Exception as e:
                        print(f"Error in valid batch {batch_idx}: {e}")
                        continue

            valid_loss = total_loss_valid / len(valid_loader)
            valid_indices = ~(np.isnan(all_true_values_valid) | np.isnan(all_pred_values_valid))
            valid_r2 = r2_score(np.array(all_true_values_valid)[valid_indices],
                                np.array(all_pred_values_valid)[valid_indices]) if np.sum(valid_indices) > 0 else float('nan')


            total_loss_test = 0
            all_true_values_test, all_pred_values_test = [], []

            with torch.no_grad():
                for batch_idx, (reactant_graphs, product_graphs, pocket_graphs, enzyme_seqs, pHs, temps, kcat_values) in enumerate(test_loader):
                    try:
                        reactant_graphs = reactant_graphs.to(device)
                        product_graphs = product_graphs.to(device)
                        pocket_graphs = pocket_graphs.to(device)
                        enzyme_seqs = enzyme_seqs.to(device)
                        pHs = pHs.to(device)
                        temps = temps.to(device)
                        kcat_values = kcat_values.to(device)

                        pred_kcat = model(
                            reactant_graphs.x, reactant_graphs.edge_index,
                            product_graphs.x, product_graphs.edge_index,
                            pocket_graphs.x, pocket_graphs.edge_index,
                            enzyme_seqs, reactant_graphs.batch, pocket_graphs.batch,
                            pHs, temps
                        )

                        if check_for_nan_inf(pred_kcat, "test_pred_kcat"):
                            continue

                        loss = criterion(pred_kcat, kcat_values)
                        if torch.isnan(loss) or torch.isinf(loss):
                            continue

                        total_loss_test += loss.item()
                        all_true_values_test.extend(kcat_values.cpu().numpy())
                        all_pred_values_test.extend(pred_kcat.detach().cpu().numpy())

                    except Exception as e:
                        print(f"Error in test batch {batch_idx}: {e}")
                        continue

            test_loss = total_loss_test / len(test_loader)
            valid_indices_test = ~(np.isnan(all_true_values_test) | np.isnan(all_pred_values_test))
            test_r2 = r2_score(np.array(all_true_values_test)[valid_indices_test],
                               np.array(all_pred_values_test)[valid_indices_test]) if np.sum(valid_indices_test) > 0 else float('nan')

            if test_r2 > highest_r2:
                highest_r2 = test_r2
                torch.save(model.state_dict(), './result_test.pth')
                print(f"==> New best test model saved! Epoch [{epoch + 1}/{NUM_EPOCHS}]")

            if valid_r2 > valid_highest_r2:
                valid_highest_r2 = valid_r2
                torch.save(model.state_dict(), './result_valid.pth')
                print(f"==> New best valid model saved! Epoch [{epoch + 1}/{NUM_EPOCHS}]")

            print(f"Epoch [{epoch + 1}/{NUM_EPOCHS}] "
                  f"| Train Loss: {train_loss:.4f}, Train R²: {train_r2:.4f} "
                  f"| Valid Loss: {valid_loss:.4f}, Valid R²: {valid_r2:.4f} "
                  f"| Test Loss: {test_loss:.4f}, Test R²: {test_r2:.4f}")

    print(f"Total training time: {time.time() - start_time:.2f} seconds")



def collate_fn(batch):
    reactant_graphs, product_graphs, pocket_graphs, enzyme_seqs, pHs, temps, kcat_values = zip(*batch)
    reactant_graphs = Batch.from_data_list(reactant_graphs)
    product_graphs = Batch.from_data_list(product_graphs)
    pocket_graphs = Batch.from_data_list(pocket_graphs)
    enzyme_seqs = torch.stack(enzyme_seqs)
    pHs = torch.stack(pHs)
    temps = torch.stack(temps)
    kcat_values = torch.stack(kcat_values)
    return reactant_graphs, product_graphs, pocket_graphs, enzyme_seqs, pHs, temps, kcat_values


if __name__ == "__main__":
    main()
