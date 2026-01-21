import json
import hashlib
import numpy as np
import os
from code.tools import read_utf_8_sig, write_utf_8_sig, request_uniprot
import requests

expandSeqPath = '../../data/preprocess/expand_seq.json'
uniprotCachePath = '../../data/preprocess/uniprot_entry_cache.json'
processPath = '../../data/preprocess/brenda_sabio_process.json'
uniprot_entry_cache = read_utf_8_sig(uniprotCachePath)
process = read_utf_8_sig(processPath)
if os.path.exists(expandSeqPath):
    res = read_utf_8_sig(expandSeqPath)
else:
    res = {}

skip = 0
count = 0
uniprot_entry_cache_len = len(uniprot_entry_cache)

for key, value in process.items():
    if key in res:
        print("skip: " + key)
        continue

    uni_id = value['UniprotID']
    seq_id = 'sequence'
    if uni_id not in uniprot_entry_cache.keys():
        uniprot_entry_cache[uni_id] = request_uniprot(uni_id)

    if uniprot_entry_cache[uni_id] is None or 'sequence' not in uniprot_entry_cache[uni_id].keys():
        skip += 1
        value[seq_id] = ''
        continue
    value[seq_id] = uniprot_entry_cache[uni_id]['sequence']['value']
    res[key] = value
    print(key + '---' + uni_id)
    count += 1
    if count % 2000 == 0:
        write_utf_8_sig(expandSeqPath, res)
    if len(uniprot_entry_cache) - uniprot_entry_cache_len >= 1000:

        write_utf_8_sig(uniprotCachePath, uniprot_entry_cache)
        uniprot_entry_cache_len = len(uniprot_entry_cache)


if uniprot_entry_cache_len != len(uniprot_entry_cache):
    write_utf_8_sig(uniprotCachePath, uniprot_entry_cache)
write_utf_8_sig(expandSeqPath, res)

