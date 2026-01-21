import json
import re
from localmapper import localmapper

# Load the JSON data
input_file = "../../data/preprocess/final_data.json"  # Replace with your JSON file path
output_file_json = "../../data/feature_extraction/reaction_mapper.json"

with open(input_file, 'r', encoding='utf-8-sig') as f:
    data = json.load(f)

# Initialize the atom mapper
mapper = localmapper()

# Prepare the output data
output_data = {}

# Total reactions to process
total_reactions = len(data)
processed_reactions = 0

# Regular expression pattern for removing specific ionic SMILES: Cl-, Fe3+, Mg2+, C-#N, Fe+2, Cu+2, Na+

#ion_pattern = r"Cl|\(Cl\)|\[Cl-\]|\[Fe\+3\]|\[Mg\+2\]|\[C-\]#N|\[Fe\+2\]|\[Cu\+2\]|\[Na\+\]|\[Br-\]|\(Br\)|\[I-\]|\(I\)|\[Mn\+2\]|\[Mn\+3\]|\[Na\+]|\[K\+]|\[Ca\+2\]|\[Zn\+2\]|\[Cu\+]|\[Al\+3\]|\[Li\+]|\[Ni\+2\]|\[Ni\]|\(Br\)|\[SH-\]|\[Cr\+3\]|\[Cr\+6\]|\[NH4\+]|\[Fe\]|\[SeH2\]|\[Se-2\]|\[Se\]|\[F-\]|\[Cl\]|\[Co\]|\[Co\+3\]|\[Hg\+2\]|\[Hg\]|\[Fe\+5\]|\[Fe\+4\]"

ion_pattern = r"Cl|\(Cl\)|\[Cl-\]|\[Fe\+3\]|\[Mg\+2\]|\[C-\]#N|\[Fe\+2\]|\[Cu\+2\]|\[Na\+\]|\[Br-\]|\(Br\)|\[I-\]|\(I\)|\[Mn\+2\]|\[Mn\+3\]|\[Na\+]|\[K\+]|\[Ca\+2\]|\[Zn\+2\]|\[Cu\+]|\[Al\+3\]|\[Li\+]|\[Ni\+2\]|\[Ni\]|\(Br\)|\[SH-\]|\[Cr\+3\]|\[Cr\+6\]|\[NH4\+]|\[Fe\]|\[SeH2\]|\[Se-2\]|\[F-\]|\[Cl\]|\[Co\]|\[Co\+3\]|\[Hg\+2\]|\[Hg\]|\[Fe\+5\]|\[Fe\+4\]"


def clean_smiles(smiles_list, pattern):
    """Clean a list of SMILES by removing specified patterns and collapsing extra dots."""
    # Remove specified ions and trim whitespace
    cleaned_list = [re.sub(pattern, "", smiles).strip() for smiles in smiles_list]
    # Filter out empty strings
    filtered_list = [s for s in cleaned_list if s]
    # Join the SMILES with dots and collapse multiple dots
    return re.sub(r"\.+", ".", ".".join(filtered_list)).strip(".")

# Process each reaction in the JSON
for key, details in data.items():
    processed_reactions += 1
    substrate_smiles = details.get("SubstrateSmiles", "").split(";")
    product_smiles = details.get("ProductSmiles", "").split(";")
    # Clean substrates and products
    substrates = clean_smiles(substrate_smiles, ion_pattern)
    products = clean_smiles(product_smiles, ion_pattern)
    reaction_smiles = f"{substrates}>>{products}"

    try:
        # Get atom map
        atom_map_result = mapper.get_atom_map(reaction_smiles)
        # Add to output data
        output_data[key] = {
            "Reaction": reaction_smiles,
            "Mapped Reaction": atom_map_result
        }
        print(f"[{processed_reactions}/{total_reactions}] Successfully processed reaction: {key}")
    except Exception as e:
        # Add error message for failed mapping
        output_data[key] = {
            "Reaction": reaction_smiles,
            "Mapped Reaction": f"Error: {str(e)}"
        }
        print(f"[{processed_reactions}/{total_reactions}] Failed to map reaction: {key} | Error: {str(e)}")

# Save the output data to a JSON file
with open(output_file_json, 'w') as jsonfile:
    json.dump(output_data, jsonfile, indent=4)


