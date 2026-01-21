import requests
import json
import os
from code.tools import read_utf_8_sig, write_utf_8_sig, request_uniprot
import requests
import time

smilesCachePath = '../../data/preprocess/merge_smile1_2_sabio_clear_version.json'

resPath = '../../data/preprocess/expand_smile2_clear_version.json'
expandSeqDataPath = '../../data/preprocess/cut_duplication.json'
smilesCache = read_utf_8_sig(smilesCachePath)
seqData = read_utf_8_sig(expandSeqDataPath)
if os.path.exists(resPath):
    res = read_utf_8_sig(resPath)
else:
    res = {}

noFind = 0
count = 0
for key, val in seqData.items():
    if key in res.keys():
        print('skip {}'.format(key))
        continue
    print(key)
    SubstrateSmiles = []
    for sub in val['Substrate'].split(';'):
        if sub in smilesCache.keys():
            SubstrateSmiles.append(smilesCache[sub])
        else:
            SubstrateSmiles.append("")
    if '' in SubstrateSmiles or len(SubstrateSmiles) == 0:
        noFind += 1
        continue
    ProductSmiles = []
    for pro in val['Product'].split(';'):
        if pro in smilesCache.keys():
            ProductSmiles.append(smilesCache[pro])
        else:
            ProductSmiles.append("")

    if '' in ProductSmiles or len(ProductSmiles) == 0:
        noFind += 1
        continue

    val['SubstrateSmiles'] = ';'.join(SubstrateSmiles)
    val['ProductSmiles'] = ';'.join(ProductSmiles)
    res[key] = val
    count += 1

write_utf_8_sig(resPath, res)


