import json
import re
from rdkit import Chem

# Function to extract heavy atoms (atoms other than Hydrogen) from a SMILES string
def extract_heavy_atoms(smiles):
    mol = Chem.MolFromSmiles(smiles)
    if mol:
        return [atom.GetSymbol() for atom in mol.GetAtoms() if atom.GetSymbol() != 'H']
    return []

# Function to extract mapped indices from mapped SMILES strings
def extract_mapped_indices(mapped_smiles):
    return [int(match[1]) - 1 for match in re.findall(r'\[([^:]+):(\d+)\]', mapped_smiles)]

# Function to process the reaction SMILES into substrates and products
def process_reaction(reaction_smiles, mapped_reaction):
    reactants_smiles, products_smiles = reaction_smiles.split(">>")
    reactant_smiles_list = reactants_smiles.split(".")
    product_smiles_list = products_smiles.split(".")
    reactant_heavy_atoms = [extract_heavy_atoms(reactant) for reactant in reactant_smiles_list]
    product_heavy_atoms = [extract_heavy_atoms(product) for product in product_smiles_list]
    #print(reactant_heavy_atoms)
    #print(reactant_smiles_list,mapped_reaction)
    mapped_reactant_atoms = [extract_mapped_indices(mapped_reaction.split(">>")[0].split(".")[i]) for i in range(len(reactant_smiles_list))]
    mapped_product_atoms = [extract_mapped_indices(mapped_reaction.split(">>")[1].split(".")[i]) for i in range(len(product_smiles_list))]

    # Ensure that mapped_reactant_atoms and mapped_product_atoms have the same total number of atoms and are not empty
    total_reactant_atoms = sum(len(reactant) for reactant in mapped_reactant_atoms)
    total_product_atoms = sum(len(product) for product in mapped_product_atoms)

    # Check if the lengths of the total atoms and heavy atoms match
    total_reactant_heavy_atoms = sum(len(heavy_atoms) for heavy_atoms in reactant_heavy_atoms)
    total_product_heavy_atoms = sum(len(heavy_atoms) for heavy_atoms in product_heavy_atoms)

    if total_reactant_atoms != total_product_atoms :
        print("nonmatch1")
        return None  # Skip this reaction if the total number of atoms or heavy atoms doesn't match
    elif  total_reactant_heavy_atoms != total_reactant_atoms :
        print("nonmatch2")
        print(reaction_id)
        print(total_reactant_heavy_atoms,total_reactant_atoms)
        print(mapped_reactant_atoms)
        return None  # Skip this reaction if the total number of atoms or heavy atoms doesn't match
    elif  total_product_heavy_atoms != total_product_atoms:
        print("nonmatch3")
        print(total_product_heavy_atoms,total_product_atoms)
        print(mapped_product_atoms)
        return None  # Skip this reaction if the total number of atoms or heavy atoms doesn't match

    elif not mapped_reactant_atoms or not mapped_product_atoms:
        print("empty")
        return None  # Skip if either mapped_reactant_atoms or mapped_product_atoms are empty
    else:
      return {
        "reactants": [
            {
                "smiles": reactant_smiles_list[i],
                "heavy_atoms": reactant_heavy_atoms[i],
                "mapped_indices": mapped_reactant_atoms[i],
            }
            for i in range(len(reactant_smiles_list))
        ],
        "products": [
            {
                "smiles": product_smiles_list[i],
                "heavy_atoms": product_heavy_atoms[i],
                "mapped_indices": mapped_product_atoms[i],
            }
            for i in range(len(product_smiles_list))
        ],
      }

# Load the JSON file with reactions
with open("../../data/feature_extraction/reaction_mapper.json", "r") as json_file:
    reactions_data = json.load(json_file)

# Process each reaction in the JSON file
processed_reactions = {}
for reaction_id, reaction_info in reactions_data.items():
    reaction_smiles = reaction_info["Reaction"]
    mapped_reaction = reaction_info.get("Mapped Reaction", "")  # Ensure it's a string
    
#    print(reaction_id)
    
#    if reaction_id == "67ba5b50fdc5cfa6ad7e5d8967f83539":
 #       print(reaction_smiles, mapped_reaction)

    # Skip if mapped_reaction contains "Error"
    if "Error" in mapped_reaction:
        print(f"Skipping reaction {reaction_id} due to error in Mapped Reaction.")
        continue

    processed_reaction = process_reaction(reaction_smiles, mapped_reaction)
    if processed_reaction:  # Only add to the result if the reaction is valid
        processed_reactions[reaction_id] = processed_reaction

# Save the processed reactions to a new JSON file
with open("../../data/feature_extraction/reaction_indice.json", "w") as output_file:
    json.dump(processed_reactions, output_file, indent=4)


