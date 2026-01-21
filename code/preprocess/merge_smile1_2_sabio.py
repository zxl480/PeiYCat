import os
from code.tools import read_utf_8_sig, write_utf_8_sig, request_uniprot
import requests
import time
from rdkit import Chem

smilesCachePath1 = '../../data/preprocess/smiles_cache_1.json'
smilesCachePath2 = '../../data/preprocess/smiles_cache_2.json'
smilesCachePath3 = '../../data/preprocess/sabio_smile_Cache.json'
resPath = '../../data/preprocess/merge_smile1_2_sabio_clear_version.json'
smilesCache1 = read_utf_8_sig(smilesCachePath1)
smilesCache2 = read_utf_8_sig(smilesCachePath2)
smilesCache3 = read_utf_8_sig(smilesCachePath3)

res = smilesCache1

for key, val in smilesCache2.items():
    if key in res:
        if res[key] == '':
            res[key] = val
        else:
            continue
    else:
        res[key] = val


for key, val in smilesCache3.items():
    if key in res:
        if res[key] == '':
            res[key] = val
        else:
            continue
    else:
        res[key] = val
count = 0
result = {}
for key, val in res.items():
    mol = Chem.MolFromSmiles(val)
    if mol is not None:
        val = Chem.MolToSmiles(mol, isomericSmiles=False)
    else:
        val = ''
    if val != "":
        count += 1

    result[key] = val

write_utf_8_sig(resPath, result)
