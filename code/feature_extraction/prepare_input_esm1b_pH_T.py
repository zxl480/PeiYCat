import json
import torch
from utility import one_hot_encode_sequence
from esm import pretrained

# Constants
REACTION_DATA_FILE = "../../data/feature_extraction/combined_reactions_amber.json"
ENZYME_DATA_FILE = "../../data/preprocess/final_data.json"
OUTPUT_DATA_FILE = "../../data/feature_extraction/prepared_input.json"


def load_reaction_data(file_path):
    """Load reaction data from JSON and extract necessary features."""
    with open(file_path, "r") as f:
        reactions = json.load(f)

    data_list = []
    for reaction in reactions:
        if "reactant_node_features" not in reaction or not isinstance(reaction["reactant_node_features"], list):
            print(f"Error: 'reactant_node_features' is missing or not a list in reaction {reaction}")
            continue
        try:
            sub_x = torch.tensor(reaction["reactant_node_features"], dtype=torch.float)
        except Exception as e:
            print(f"Error converting reactant_node_features to tensor: {e}")
            continue

        if "reactant_edge_index" in reaction:
            try:
                sub_edge = torch.tensor(reaction["reactant_edge_index"], dtype=torch.long)
            except Exception as e:
                print(f"Error converting reactant_edge_index to tensor: {e}")
                continue
        else:
            print("Error: 'reactant_edge_index' is missing in reaction.")
            continue

        try:
            prod_x = torch.tensor(reaction["product_node_features"], dtype=torch.float)
        except Exception as e:
            print(f"Error converting product_node_features to tensor: {e}")
            continue

        try:
            prod_edge = torch.tensor(reaction["product_edge_index"], dtype=torch.long)
        except Exception as e:
            print(f"Error converting product_edge_index to tensor: {e}")
            continue

        data_list.append((sub_x, sub_edge, prod_x, prod_edge, reaction["ReactionID"]))
    return data_list


# Load ESM-1b model
device = torch.device("cuda:0" if torch.cuda.is_available() else "cpu")
esm_model, alphabet = pretrained.esm1b_t33_650M_UR50S()
esm_batch_converter = alphabet.get_batch_converter()
esm_model = esm_model.to(device)
esm_model.eval()


def extract_esm1b_feature(sequence):
    """Extract ESM-1b feature representation of the enzyme sequence using GPU."""
    data = [("enzyme", sequence)]
    batch_labels, batch_strs, batch_tokens = esm_batch_converter(data)
    batch_tokens = batch_tokens.to(device)
    with torch.no_grad():
        results = esm_model(batch_tokens, repr_layers=[33], return_contacts=False)
    return results["representations"][33].mean(dim=1).squeeze()


def load_enzyme_data(file_path):
    """Load enzyme sequence, kcat, temperature, and pH data from the enzyme file."""
    with open(file_path, "r", encoding='utf-8-sig') as f:
        enzymes = json.load(f)

    enzyme_data = {}
    total_sequences = len(enzymes)
    processed_sequences = 0

    for reaction_id, enzyme in enzymes.items():
        if "mutant" in enzyme:
            enzyme_seq_str = enzyme["mutant"]
        elif "chain" in enzyme:
            enzyme_seq_str = enzyme["chain"]
        else:
            print(f"Error: Neither 'mutant' nor 'chain' found for reaction {reaction_id}")
            continue

        enzyme_seq = extract_esm1b_feature(enzyme_seq_str)
        if enzyme_seq is None:
            print("Warning: Failed to extract ESM-1b features!")
            continue

        try:
            kcat_value = float(enzyme["kcat"])
        except ValueError:
            print(f"Error: Unable to convert kcat value for reaction {reaction_id}")
            continue

        # Handle temperature and pH with defaults
        temperature_raw = enzyme.get("temperature", "RT")
        pH_raw = enzyme.get("pH", "NT")
        try:
            temperature_value = 25.0 if temperature_raw == "RT" else float(temperature_raw)
        except:
            print(f"Warning: Invalid temperature '{temperature_raw}' in reaction {reaction_id}, using 25.0")
            temperature_value = 25.0
        try:
            pH_value = 7.5 if pH_raw == "NT" else float(pH_raw)
        except:
            print(f"Warning: Invalid pH '{pH_raw}' in reaction {reaction_id}, using 7.5")
            pH_value = 7.5

        enzyme_data[reaction_id] = {
            "esm_feature": enzyme_seq,
            "kcat": kcat_value,
            "temperature": temperature_value,
            "pH": pH_value
        }

        processed_sequences += 1
        if processed_sequences % 10 == 0:
            print(f"Processed {processed_sequences}/{total_sequences} sequences")

    print(f"Finished processing all {processed_sequences} sequences.")
    return enzyme_data


def combine_data(reaction_data, enzyme_data):
    """Combine reaction data and enzyme data based on reaction_ID."""
    combined_data = []

    for sub_x, sub_edge, prod_x, prod_edge, reaction_id in reaction_data:
        if reaction_id in enzyme_data:
            enzyme_info = enzyme_data[reaction_id]
            combined_data.append({
                'reaction_id': reaction_id,
                'reactant_features': sub_x.tolist(),
                'reactant_edge': sub_edge.tolist(),
                'product_features': prod_x.tolist(),
                'product_edge': prod_edge.tolist(),
                'enzyme_sequence': enzyme_info["esm_feature"].tolist(),
                'kcat': enzyme_info["kcat"],
                'temperature': enzyme_info["temperature"],
                'pH': enzyme_info["pH"]
            })
        else:
            print(f"Warning: No enzyme data found for reaction {reaction_id}")

    return combined_data


def save_combined_data(combined_data):
    """Save the combined data to a JSON file."""
    with open(OUTPUT_DATA_FILE, "w") as f:
        json.dump(combined_data, f, indent=4)
    print(f"Combined data saved to {OUTPUT_DATA_FILE}")


def main():
    # Load reaction data and enzyme data
    reaction_data = load_reaction_data(REACTION_DATA_FILE)
    enzyme_data = load_enzyme_data(ENZYME_DATA_FILE)

    # Combine the reaction data with enzyme data
    combined_data = combine_data(reaction_data, enzyme_data)

    # Save the combined data to JSON file
    save_combined_data(combined_data)


if __name__ == "__main__":
    main()
