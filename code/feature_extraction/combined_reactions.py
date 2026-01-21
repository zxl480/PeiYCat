import json

# Function to read and load the processed data from the JSON file
def read_processed_data(file_path):
    with open(file_path, "r") as json_file:
        data = json.load(json_file)
    return data

# Function to convert edge indices to real atom indices based on mapped indices
def convert_edge_index_to_real_index(edge_index, mapped_indices):
    index_mapping = {i: mapped_indices[i] for i in range(len(mapped_indices))}
    real_edge_index = [[], []]
    for atom1, atom2 in zip(edge_index[0], edge_index[1]):
        try:
            real_edge_index[0].append(index_mapping[atom1])
            real_edge_index[1].append(index_mapping[atom2])
        except KeyError:
            print(f"Warning: Atom indices {atom1} or {atom2} not found in mapped indices!")
    return real_edge_index

# Function to reorder combined node features and mapped indices
def reorder_combined_features_and_indices(combined_features, combined_indices):
    # Ensure they are paired by their indices and then sorted
    paired_data = sorted(zip(combined_indices, combined_features), key=lambda x: x[0])
    # Unzip sorted pairs into separate lists
    sorted_indices, sorted_features = zip(*paired_data)
    return list(sorted_features), list(sorted_indices)

# Modify the combine_reaction_data_custom function to add the reorder step for combined node features and mapped indices
def combine_reaction_data_custom(reactants, products):
    def combine_smiles(smiles_list):
        return ".".join(smiles_list)

    def combine_mapped_incidence(mapped_incidence_list):
        return [item for sublist in mapped_incidence_list for item in sublist]

    def combine_node_features(features_list):
        return [feature for sublist in features_list for feature in sublist]

    def combine_edge_index(edge_index_list):
        combined_edge_index = [[], []]
        for edge_index in edge_index_list:
            if edge_index is not None:  # Skip None values
                combined_edge_index[0].extend(edge_index[0])
                combined_edge_index[1].extend(edge_index[1])
        return combined_edge_index
   
    # Combine Reactants
    reactant_smiles = combine_smiles([r['SMILES'] for r in reactants])
    reactant_mapped_incidence = combine_mapped_incidence([r['Mapped Indices'] for r in reactants])
    reactant_node_features = combine_node_features([r['Node Features'] for r in reactants])
    reactant_edge_index = combine_edge_index([r['Real Edge Index'] for r in reactants if r['Real Edge Index'] is not None])

    # Combine Products
    product_smiles = combine_smiles([p['SMILES'] for p in products])
    product_mapped_incidence = combine_mapped_incidence([p['Mapped Indices'] for p in products])
    product_node_features = combine_node_features([p['Node Features'] for p in products])
    product_edge_index = combine_edge_index([p['Real Edge Index'] for p in products if p['Real Edge Index'] is not None])

    # Reorder the combined node features and mapped indices for reactants and products
    reordered_reactant_node_features, reordered_reactant_mapped_indices = reorder_combined_features_and_indices(
        reactant_node_features, reactant_mapped_incidence
    )
    reordered_product_node_features, reordered_product_mapped_indices = reorder_combined_features_and_indices(
        product_node_features, product_mapped_incidence
    )
    return {
        'reactant_SMILES': reactant_smiles,
        'reactant_mapped_incidence': reordered_reactant_mapped_indices,
        'reactant_node_features': reordered_reactant_node_features,
        'reactant_edge_index': reactant_edge_index,
        'product_SMILES': product_smiles,
        'product_mapped_incidence': reordered_product_mapped_indices,
        'product_node_features': reordered_product_node_features,
        'product_edge_index': product_edge_index,
    }

# Main script to process data and combine reactions (no change here)
processed_data = read_processed_data("../../data/feature_extraction/reaction_feature_amber.json")
combined_reactions = []

for i, entry in enumerate(processed_data):
    reaction_data = {
        'ReactionID': entry['ReactionID'],
    }
    reactants = []
    for reactant in entry['Reactants']:
        # Check if edge_index is not None before converting it
        real_edge_index = None
        if reactant['edge_index'] is not None:
            real_edge_index = convert_edge_index_to_real_index(
                reactant['edge_index'], reactant['Mapped Indices']
            )

        reactants.append({
            'SMILES': reactant['SMILES'],
            'Mapped Indices': reactant['Mapped Indices'],
            'Node Features': reactant['node_features'],
            'Real Edge Index': real_edge_index,
        })

    products = []
    for product in entry['Products']:
        # Check if edge_index is not None before converting it
        real_edge_index = None
        if product['edge_index'] is not None:
            real_edge_index = convert_edge_index_to_real_index(
                product['edge_index'], product['Mapped Indices']
            )

        products.append({
            'SMILES': product['SMILES'],
            'Mapped Indices': product['Mapped Indices'],
            'Node Features':  product['node_features'],
            'Real Edge Index': real_edge_index,
        })

    combined_data = combine_reaction_data_custom(reactants, products)
    reaction_data.update(combined_data)
    combined_reactions.append(reaction_data)

# Save or print the results
with open("../../data/feature_extraction/combined_reactions_amberr.json", "w") as output_file:
    json.dump(combined_reactions, output_file, indent=4)



