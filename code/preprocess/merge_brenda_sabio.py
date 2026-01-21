import json
import hashlib
import numpy as np
from code.tools import read_utf_8_sig, write_utf_8_sig


sabio_data = read_utf_8_sig('../../data/preprocess/sabio_data_get_sub_pro_clean_ion.json')
brenda_data = read_utf_8_sig('../../data/preprocess/brenda_data_get_sub_pro_clean_ion_2.json')


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

brenda_res = {}
for key, _entry in brenda_data.items():
    if _entry['kcat_max'] != '-':
        _entry['kcat'] = str(np.sqrt(float(_entry['kcat']) * float(_entry['kcat_max'])))
    _entry.pop('kcat_max')

    tmp = {
        'ECNumber': _entry['ECNumber'],
        'Substrate': _entry['Substrate'],
        'commentary': _entry['commentary'],
        'Organism': _entry['Organism'],
        'Product': _entry['Product'],
        'temperature': str(_entry['pH']),
        'pH': str(_entry['pH']),
    }
    for _t in _entry.keys():
        if _t.startswith('UniprotID'):
            tmp[_t] = _entry[_t]

    entry_hash = hash_entry(tmp)

    if entry_hash in brenda_res.keys():
        _entry['kcat'] = str(np.sqrt(float(_entry['kcat']) * float(brenda_res[entry_hash]['kcat'])))
    _entry['unit'] = 's^(-1)'
    if 'mutant' in _entry['commentary']:
        _entry['Type'] = 'mutant'
    else:
        _entry['Type'] = 'wildtype'
    brenda_res[entry_hash] = _entry

sabio_res = {}
for key, _entry in sabio_data.items():
    if 'kcat' not in _entry.keys():
        continue

    if 'unit' not in _entry:
        _entry['unit'] = 's^(-1)'
    _entry.pop('EntryID')
    tmp = {
        'ECNumber': _entry['ECNumber'],
        'Substrate': _entry['Substrate'],
        'EnzymeType': _entry['EnzymeType'],
        'Organism': _entry['Organism'],
        'Product': _entry['Product'],
        'temperature': str(_entry['pH']),
        'pH': str(_entry['pH']),
    }
    for _t in _entry.keys():
        if _t.startswith('UniprotID'):
            tmp[_t] = _entry[_t]
    entry_hash = hash_entry(tmp)
    # print(key)
    if entry_hash in sabio_res.keys():
        _entry['kcat'] = str(np.sqrt(float(_entry['kcat']) * float(sabio_res[entry_hash]['kcat'])))
    if 'mutant' in _entry['EnzymeType']:
        _entry['Type'] = 'mutant'
    else:
        _entry['Type'] = 'wildtype'

    sabio_res[entry_hash] = _entry

res = {}
for key, entry in brenda_res.items():
    res[key] = entry


for key, entry in sabio_res.items():
    res[key] = entry

write_utf_8_sig('../../data/preprocess/brenda_sabio_process.json', res)



