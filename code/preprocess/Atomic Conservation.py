import json
from rdkit import Chem
from collections import Counter
from code.tools import write_utf_8_sig


with open('../../data/preprocess/expand_smile2_clear_version.json', 'r', encoding='utf-8-sig') as f:
    data = json.load(f)


def remove_semicolons(smiles):
    return smiles.replace(";", "")


def get_atomic_counts(smiles):
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return None
    atoms = mol.GetAtoms()
    atom_counts = Counter([atom.GetSymbol() for atom in atoms])
    return atom_counts



def check_atom_conservation(substrate_smiles, product_smiles):
    # substrate_smiles = remove_semicolons(substrate_smiles)
    # product_smiles = remove_semicolons(product_smiles)
    substrate_counts = Counter()
    product_counts = Counter()
    for sub in substrate_smiles:
        # if sub == '[H+]':
        #     continue
        # print(sub)
        _count = get_atomic_counts(sub)
        if _count is None:
            return False, Counter(), Counter()
        substrate_counts += _count
    for pro in product_smiles:
        # if pro == '[H+]':
        #     continue
        _count = get_atomic_counts(pro)
        if _count is None:
            return False, Counter(), Counter()
        product_counts += _count

    return substrate_counts == product_counts, substrate_counts, product_counts

match = 0
unMatch = 0
res = {}

with open('../../data/preprocess/atomic_conservation.txt', 'w') as f:
    for key, value in data.items():
        substrate_smiles = value["SubstrateSmiles"].split(';')
        product_smiles = value["ProductSmiles"].split(';')

        is_conserved, substrate_counts, product_counts = check_atom_conservation(substrate_smiles, product_smiles)

        if not is_conserved:
            f.write(f"unmatch Key: {key}\t")
            f.write(f"Substrate: {substrate_counts}\t")
            f.write(f"Product: {product_counts}\t")
            f.write("Atom conservation is NOT maintained.\n")
            unMatch += 1
        else:

            match += 1
            res[key] = value

write_utf_8_sig('../../data/preprocess/atomic_match_results.json', res)
