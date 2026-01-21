import os
import json
import numpy as np
from rdkit import Chem
import deepchem as dc


MOL2_FOLDER = "../../data/feature_extraction/amber_result"
SMI_PATH = "../../data/feature_extraction/all_smile_molecules.smi"


def load_smiles_mol2_mapping(smi_path):
    mapping = {}
    with open(smi_path, "r") as f:
        for idx, line in enumerate(f, start=1):
            if line.strip():
                smiles = line.strip().split('\t')[0]
                mol2_path = os.path.join(MOL2_FOLDER, f"{idx}_bcc_noH.mol2")
                mapping[smiles] = mol2_path
    return mapping

def extract_charges_from_mol2(mol2_path):
    charges = []
    with open(mol2_path, "r") as f:
        lines = f.readlines()
    atom_section = False
    for line in lines:
        if line.startswith("@<TRIPOS>ATOM"):
            atom_section = True
            continue
        if line.startswith("@<TRIPOS>"):
            if atom_section:

                break
        if atom_section:
            tokens = line.strip().split()
            if len(tokens) >= 9:
                try:
                    charge = float(tokens[-1])
                    charges.append(charge)
                except ValueError:
                    # 防止解析失败
                    charges.append(0.0)
    return charges


def featurize_molecule(smiles, mol2_path=None):
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        raise ValueError(f"Invalid SMILES: {smiles}")


    predefined = {
        "C": [1.0] + [0.0]*29,
        "O": [0.0, 0.0, 1.0] + [0.0]*27,
        "N": [0.0, 1.0, 0.0] + [0.0]*27,
    }
    if smiles in predefined:
        node_feature = predefined[smiles] + [0.0]
        if len(node_feature) != 31:
            print(f"Warning: Predefined SMILES {smiles} has feature dimension {len(node_feature)}")
        return {
            "node_features": [node_feature],
            "edge_index": None
        }



    charged_smiles_map = {
        "[NH4+]": +1,
        "[SH-]": -1,
        "[S-2]": -2,
        "[F-]": -1,
        "[Hg+2]": +2
    }

    neutral_smiles = ["[S]", "Br", "[Fe]", "Cl", "S", "[Hg]"]

    if smiles in charged_smiles_map or smiles in neutral_smiles:
        charge = charged_smiles_map.get(smiles, 0.0)
        zero_feat = [0.0] * 30 + [charge]
        return {
            "node_features": [zero_feat],
            "edge_index": None
        }

    featurizer = dc.feat.MolGraphConvFeaturizer()
    features = featurizer.featurize([mol])[0]

    node_feats = features.node_features
    node_feats_list = []
    for nf in node_feats:
        node_feats_list.append(nf.tolist() if isinstance(nf, np.ndarray) else list(nf))


    if mol2_path and os.path.exists(mol2_path):
        charges = extract_charges_from_mol2(mol2_path)
        if len(charges) != len(node_feats_list):
            raise ValueError(f"Mismatch between mol2 atoms ({len(charges)}) and node features ({len(node_feats_list)}) for {smiles}")
    else:
        charges = [0.0] * len(node_feats_list)


    new_node_feats = []
    for i, nf in enumerate(node_feats_list):
        full_feat = nf + [charges[i]]
        if len(full_feat) != 31:
            print(f"Feature length error for atom {i} in {smiles}: expected 31, got {len(full_feat)}. Original: {nf}, charge: {charges[i]}")
        new_node_feats.append(full_feat)

    edge_index = features.edge_index
    edge_index_list = (
        edge_index.tolist() if isinstance(edge_index, np.ndarray)
        else list(edge_index) if isinstance(edge_index, (list, tuple))
        else None
    )

    return {
        "node_features": new_node_feats,
        "edge_index": edge_index_list
    }



def convert_ndarray_to_list(obj):
    if isinstance(obj, np.ndarray):
        return obj.tolist()
    elif isinstance(obj, dict):
        return {k: convert_ndarray_to_list(v) for k, v in obj.items()}
    elif isinstance(obj, list):
        return [convert_ndarray_to_list(i) for i in obj]
    else:
        return obj


smiles_to_mol2 = load_smiles_mol2_mapping(SMI_PATH)


def process_json_data(json_data):
    results = {}

    for reaction_id, value in json_data.items():
        substrates = value.get("reactants", [])
        products = value.get("products", [])
        ec_number = reaction_id

        if ec_number not in results:
            results[ec_number] = {
                "ReactionID": reaction_id,
                "Reactants": [],
                "Products": []
            }

        for reactant in substrates:
            smiles = reactant["smiles"]
            mol2_path = smiles_to_mol2.get(smiles, None)
            try:
                features = featurize_molecule(smiles, mol2_path)
                results[ec_number]["Reactants"].append({
                    "SMILES": smiles,
                    "Heavy Atoms": reactant["heavy_atoms"],
                    "Mapped Indices": reactant["mapped_indices"],
                    "node_features": features["node_features"],
                    "edge_index": features["edge_index"]
                })
                if len(features["node_features"]) != len(reactant["mapped_indices"]):
                    print(reaction_id, "Reactant nonmatch")
            except Exception as e:
                print(f"Reactant error ({reaction_id}, {smiles}): {e}")

        for product in products:
            smiles = product["smiles"]
            mol2_path = smiles_to_mol2.get(smiles, None)
            try:
                features = featurize_molecule(smiles, mol2_path)
                results[ec_number]["Products"].append({
                    "SMILES": smiles,
                    "Heavy Atoms": product["heavy_atoms"],
                    "Mapped Indices": product["mapped_indices"],
                    "node_features": features["node_features"],
                    "edge_index": features["edge_index"]
                })
                if len(features["node_features"]) != len(product["mapped_indices"]):
                    print(reaction_id, "Product nonmatch")
            except Exception as e:
                print(f"Product error ({reaction_id}, {smiles}): {e}")

    return list(results.values())

# 主流程
if __name__ == "__main__":
    with open("../../data/feature_extraction/reaction_indice.json", "r") as json_file:
        processed_reactions = json.load(json_file)

    processed_data = process_json_data(processed_reactions)

    processed_data = convert_ndarray_to_list(processed_data)

    with open("../../data/feature_extraction/reaction_feature_amber.json", "w") as output_file:
        json.dump(processed_data, output_file, indent=4)


