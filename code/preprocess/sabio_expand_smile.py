import requests
import os
from code.tools import read_utf_8_sig, write_utf_8_sig, request_uniprot
import json


def get_smile(rId):
    QUERY_URL = 'https://sabiork.h-its.org/sabioRestWebServices/searchReactionParticipants'

    # input: SabioReactionID
    # valid output fields: "fields[]":["Name","Role","SabioCompoundID","ChebiID","PubChemID","KeggCompoundID", "InChI","Smiles"]
    query = {"SabioReactionID": rId, "fields[]": ["Name", "Smiles"]}

    response = requests.get(QUERY_URL, params=query)
    if response.status_code == 200:
        return response.text
    else:
        print(f"Error: Received status code {response.status_code}")
        return None

# d = get_smile('1954')
# print(d)


EcPath = '../../data/preprocess/sabio_Ec_data'
files = os.listdir(EcPath)
smileCaCheFile = '../../data/preprocess/sabio_smile_Cache.json'
reactionIdCacheFile = '../../data/preprocess/reactionIdCache.json'
if os.path.exists(smileCaCheFile):
    sabio_smile_cache = read_utf_8_sig(smileCaCheFile)
else:
    sabio_smile_cache = {}
if os.path.exists(reactionIdCacheFile):
    reactionIdCache = read_utf_8_sig(reactionIdCacheFile)
else:
    reactionIdCache = []
count = 0

for file in files:
    _path = os.path.join(EcPath, file)
    print(_path)
    with open(_path, 'r', encoding='utf-8-sig') as f:
        fileData = f.readlines()
    for line in fileData[1:]:
        _line = line.replace('\n', '').split('\t')
        reactionId = _line[-1]
        if reactionId in reactionIdCache:
            continue
        smiles = get_smile(reactionId)
        reactionIdCache.append(reactionId)
        if smiles is not None:
            smiles = smiles.split('\n')
            for smile in smiles[1:]:
                _smile = smile.split('\t')
                if len(_smile) != 2 or _smile[1] == 'null':
                    continue
                name, _smile = _smile
                sabio_smile_cache[name] = _smile.strip()
                # print("{}-----{}".format(name, sabio_smile_cache[name]))
                count += 1
                if count % 1000 == 0:
                    write_utf_8_sig(smileCaCheFile, sabio_smile_cache)
                    write_utf_8_sig(reactionIdCacheFile, reactionIdCache)

write_utf_8_sig(smileCaCheFile, sabio_smile_cache)
write_utf_8_sig(reactionIdCacheFile, reactionIdCache)



