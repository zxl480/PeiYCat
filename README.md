# PeiYCat

## Introduction
The PeiYcat is a Python package for predicting enzyme kcat values using a graph-based deep learning model. The repository provides a complete pipeline to encode protein sequences with ESM-1b embeddings, represent reactions using atom-mapped substrates and products, integrate pocket features, and train the ReactionModel for accurate kcat prediction.


## Required software
- Python v3.9
- pytorch-lightning v2.5.0.post0
- scikit-learn v1.6.0
- pandas v1.5.3
- SciPy v1.10.1
- NumPy v1.26.3
- Biopython v1.85
- RDKit 2024.9.4
- seaborn v0.13.2
- Matplotlib v3.9.2

## Data collection and preprocessing
### Data collection and preprocessing from the SABIO-RK database
- run sabio_data_to_json.py to extract reactions with complete stoichiometry and valid kinetic parameters and convert the data into JSON format
- run sabio_filter_Kcat.py to filter kcat values and unify all units to s⁻¹
- run sabio_data_break_uniprotId.py to resolve enzyme complexes and assign a single representative UniProt identifier
- run sabio_get_sub_pro_clean_ion.py to separate substrates and products, clean stoichiometric coefficients, and remove ubiquitous ions

### Data collection and preprocessing from the BRENDA database
- run brenda_ec_data_filter.py to parse and filter downloaded EC-specific CSV files
- run brenda_download_reaction.py to retrieve reaction information
- run brenda_reaction_data_filter.py to parse reaction text files and generate unique reaction entries
- run brenda_match_ec_reaction.py to assign a unique reaction to each EC number
- run brenda_match_pH_temperature.py to extract pH and temperature information
- run brenda_data_break_uniprotID.py to resolve enzyme complexes and assign representative UniProt identifiers
- run brenda_get_sub_pro_clean_ion.py to clean substrates and products and remove ions
- run brenda_get_sub_pro_clean_ion_version2.py to re-parse and correct reaction stoichiometric coefficients

### Data merging
- run merge_brenda_sabio.py to merge the curated BRENDA and SABIO-RK datasets
- run expand_seq.py to retrieve protein sequences from UniProt based on UniProt identifiers
- run match_mutant_find_chain.py to match enzyme chains and annotate mutation information
- run cut_duplication.py to remove duplicate entries between databases
- run expand_smile.py to retrieve canonical SMILES for substrates and products
- run atomic_conservation.py to perform atomic mass conservation checks on reactions

## Reaction representation, feature construction, and model training
For reaction representation and feature construction:
- run reaction_mapper.py to perform local atom mapping for reaction equations
- run reaction_indice.py to extract atom-mapping indices for downstream processing
- run reaction_feature.py to generate DeepChem-based features for substrates and products in each reaction
- run combined_reactions.py to merge reaction features into unified reaction representations

For protein and pocket feature preparation:
- run prepare_input_esm1b_pH_T.py to encode protein sequences using the ESM-1b model
- run combine_pocket_esm1b_to_data.py to combine protein pocket features with ESM-1b sequence embeddings

For model training:
- run train.py to train the deep learning model using GPU acceleration

