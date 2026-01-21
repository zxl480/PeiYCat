import json
import hashlib
import numpy as np
from code.tools import read_utf_8_sig, write_utf_8_sig


sabio_brenda_data = read_utf_8_sig('../../data/preprocess/match_mutant_find_chain.json')


def normalize_order(data):
    """Normalize the order of semi-colon-separated strings."""
    return ';'.join(sorted(data.split(';')))


def hash_entry(entry):
    """Generate a hash for a given entry, ignoring EntryID and normalizing substrate order."""
    entry_copy = entry.copy()

    # Normalize Substrate and Substrate_Smiles order
    if "Substrate" in entry_copy:
        entry_copy["Substrate"] = normalize_order(entry_copy["Substrate"])
    # if "Substrate_Smiles" in entry_copy:
    #     entry_copy["Substrate_Smiles"] = normalize_order(entry_copy["Substrate_Smiles"])

    if "Product" in entry_copy:
        entry_copy["Product"] = normalize_order(entry_copy["Product"])
    # if "Product_Smiles" in entry_copy:
    #     entry_copy["Product_Smiles"] = normalize_order(entry_copy["Product_Smiles"])

    # Convert to JSON string and generate a hash
    entry_json = json.dumps(entry_copy, sort_keys=True)
    return hashlib.md5(entry_json.encode()).hexdigest()

sabio_brenda_res = {}
for key, _entry in sabio_brenda_data.items():
    if "mutant" in _entry:
        tmp = {
            'ECNumber': _entry['ECNumber'],
            'Substrate': _entry['Substrate'],
            'Organism': _entry['Organism'],
            'Product': _entry['Product'],
            'temperature': str(_entry['pH']),
            'pH': str(_entry['pH']),
            'UniprotID': _entry['UniprotID'],
            'sequence': _entry['sequence'],
            'chain': _entry['chain'],
            'mutant': _entry['mutant'],
        }
        entry_hash = hash_entry(tmp)
    elif "mutant" not in _entry:
        tmp = {
            'ECNumber': _entry['ECNumber'],
            'Substrate': _entry['Substrate'],
            'Organism': _entry['Organism'],
            'Product': _entry['Product'],
            'temperature': str(_entry['pH']),
            'pH': str(_entry['pH']),
            'UniprotID': _entry['UniprotID'],
            'sequence': _entry['sequence'],
            'chain': _entry['chain'],

        }
        entry_hash = hash_entry(tmp)


    if entry_hash in sabio_brenda_res.keys():
        _entry['kcat'] = str(np.sqrt(float(_entry['kcat']) * float(sabio_brenda_res[entry_hash]['kcat'])))

    sabio_brenda_res[entry_hash] = _entry

res = {}
for key, entry in sabio_brenda_res.items():
    res[key] = entry

write_utf_8_sig('../../data/preprocess/cut_duplication.json', res)


