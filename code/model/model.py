# === model.py ===
import torch
import torch.nn as nn
import torch.nn.functional as F
from torch_geometric.nn import MessagePassing
from torch_geometric.utils import add_self_loops

class MPNNLayer(MessagePassing):
    def __init__(self, in_channels, out_channels):
        super(MPNNLayer, self).__init__(aggr='mean')
        self.lin = nn.Linear(in_channels, out_channels)

    def forward(self, x, edge_index):
        edge_index = add_self_loops(edge_index, num_nodes=x.size(0))[0]
        return self.propagate(edge_index, x=self.lin(x))

    def message(self, x_j):
        return x_j

# Multi-layer MPNN
class MPNN(nn.Module):
    def __init__(self, num_layers=5, in_dim=30, hidden_dim=30):
        super(MPNN, self).__init__()
        self.layers = nn.ModuleList([MPNNLayer(in_dim, hidden_dim) for _ in range(num_layers)])

    def forward(self, x, edge_index):
        for layer in self.layers:
            x = F.relu(layer(x, edge_index))
        return x

class ReactionModel(nn.Module):
    def __init__(self, num_layers=5, feature_dim=31, pocket_feature_dim=4, dropout_rate=0.2,
                 num_amino_acids=20, num_atom_types=4):
        super(ReactionModel, self).__init__()
        self.mpnn_sub_prod = MPNN(num_layers, feature_dim, feature_dim)

        # Pocket特征处理层
        # 类别特征使用embedding
        self.amino_acid_embedding = nn.Embedding(num_amino_acids, 8)  # 氨基酸类型 -> 8维
        self.atom_type_embedding = nn.Embedding(num_atom_types, 4)    # 原子类型 -> 4维

        # 连续特征使用Linear + LayerNorm
        self.distance_transform = nn.Sequential(
            nn.Linear(1, 4), nn.LayerNorm(4), nn.ReLU()
        )
        self.charge_transform = nn.Sequential(
            nn.Linear(1, 4), nn.LayerNorm(4), nn.ReLU()
        )

        # 计算pocket特征的总维度: 8 + 4 + 4 + 4 = 20
        pocket_embedded_dim = 20
        self.mpnn_pocket = MPNN(num_layers, pocket_embedded_dim, pocket_embedded_dim)

        # 温度和pH的嵌入层
        self.pH_layer = nn.Sequential(
            nn.Linear(1, 8), nn.LayerNorm(8), nn.ReLU())
        self.temp_layer = nn.Sequential(
            nn.Linear(1, 8), nn.LayerNorm(8), nn.ReLU())

        # 更新FNN的输入维度：加上温度和pH的嵌入维度（8+8）
        # 1280 (enzyme) + 1240 (sub) + 1240 (diff) + 100*20 (pocket) + 8 + 8 = 5776
        self.fnn = nn.Sequential(
            nn.Linear(5776, 10240), nn.ReLU(), nn.Dropout(dropout_rate),
            nn.Linear(10240, 5120), nn.ReLU(), nn.Dropout(dropout_rate),
            nn.Linear(5120, 2048), nn.ReLU(), nn.Dropout(dropout_rate),
            nn.Linear(2048, 1024), nn.ReLU(), nn.Dropout(dropout_rate),
            nn.Linear(1024, 128), nn.ReLU(), nn.Dropout(dropout_rate),
            nn.Linear(128, 1)
        )
        # self.fnn = nn.Sequential(
        #     nn.Linear(5776, 5120), nn.ReLU(), nn.Dropout(dropout_rate),
        #     nn.Linear(5120, 2048), nn.ReLU(), nn.Dropout(dropout_rate),
        #     nn.Linear(2048, 1024), nn.ReLU(), nn.Dropout(dropout_rate),
        #     nn.Linear(1024, 128), nn.ReLU(), nn.Dropout(dropout_rate),
        #     nn.Linear(128, 1)
        # )

    def forward(self, sub_x, sub_edge, prod_x, prod_edge, pocket_x, pocket_edge,
                enzyme_seq, sub_batch, pocket_batch, pH, temperature):
        # Apply MPNN to substrates and products
        sub_feat = self.mpnn_sub_prod(sub_x, sub_edge)
        prod_feat = self.mpnn_sub_prod(prod_x, prod_edge)

        # 处理pocket特征
        # pocket_x的形状: (total_pocket_atoms, 4)
        # 特征顺序: [氨基酸类型, 原子类型, 距离, 电荷]
        amino_acid_ids = pocket_x[:, 0].long()
        atom_type_ids = pocket_x[:, 1].long()
        distances = pocket_x[:, 2:3]  # 保持2D形状
        charges = pocket_x[:, 3:4]

        amino_acid_emb = self.amino_acid_embedding(amino_acid_ids)
        atom_type_emb = self.atom_type_embedding(atom_type_ids)
        distance_emb = self.distance_transform(distances)
        charge_emb = self.charge_transform(charges)

        pocket_embedded = torch.cat([amino_acid_emb, atom_type_emb, distance_emb, charge_emb], dim=1)
        pocket_feat = self.mpnn_pocket(pocket_embedded, pocket_edge)

        # Compute per-atom feature difference
        feature_diff = torch.norm(sub_feat - prod_feat, dim=1)
        feature_diff_real = sub_feat - prod_feat

        batch_size = sub_batch.max().item() + 1
        top_features_list, sub_features_list = [], []

        for i in range(batch_size):
            node_mask = (sub_batch == i)
            sub_feat_graph = sub_feat[node_mask]
            feature_diff_graph = feature_diff[node_mask]
            feature_diff_real_graph = feature_diff_real[node_mask]
            num_atoms = feature_diff_graph.size(0)
            num_top_atoms = min(40, num_atoms)
            top_atoms = torch.topk(feature_diff_graph, num_top_atoms).indices
            top_features = feature_diff_real_graph[top_atoms]
            sub_feat_top = sub_feat_graph[top_atoms]

            # Padding if the graph has fewer than 40 atoms
            if num_atoms < 40:
                pad = lambda x: torch.cat([x, torch.zeros(40 - num_atoms, x.size(1), device=x.device)], dim=0)
                top_features = pad(top_features)
                sub_feat_top = pad(sub_feat_top)

            top_features_list.append(top_features)
            sub_features_list.append(sub_feat_top)

        top_features_tensor = torch.stack(top_features_list, dim=0)
        sub_features_tensor = torch.stack(sub_features_list, dim=0)
        pocket_feature_flattened = pocket_feat.view(batch_size, -1)

        top_feature_flattened = top_features_tensor.view(batch_size, -1)
        sub_feature_flattened = sub_features_tensor.view(batch_size, -1)

        # 融合 pH 和 temperature 的嵌入
        pH_emb = self.pH_layer(pH.unsqueeze(1))
        temp_emb = self.temp_layer(temperature.unsqueeze(1))

        # 拼接所有特征：enzyme + substrate + difference + pocket + pH + temperature
        x = torch.cat([enzyme_seq, sub_feature_flattened, top_feature_flattened,
                       pocket_feature_flattened, pH_emb, temp_emb], dim=1)
        return self.fnn(x)
