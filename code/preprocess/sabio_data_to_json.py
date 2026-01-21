import os
import json
import requests
from bs4 import BeautifulSoup


def get_reaction(reactionId, ecNumber):
    QUERY_URL = 'https://sabiork.h-its.org/sabioRestWebServices/searchReactionDetails'

    # input: SabioReactionID
    # valid output fields: "fields[]":["Name","Role","SabioCompoundID","ChebiID","PubChemID","KeggCompoundID", "InChI","Smiles"]
    # example
    query = {"SabioReactionID": str(reactionId), "fields[]": ["SabioReactionID", "ReactionEquation", "ECNumber",]}


    # note that you can also retrieve the data for ALL reaction participants by performing a wildcard search
    # for this wildcard query, the output field 'SabioReactionID' is added by default
    # query = {"SabioReactionID":"*", "fields[]":["SabioCompoundID","Name","Role"]}
    request = requests.get(QUERY_URL, params=query)
    if request.status_code == 200:
        _d = request.text.split('\n')
        for __d in _d:
            __d = __d.split('\t')
            if reactionId == __d[0] and ecNumber == __d[2]:
                return __d[1]
    else:
        return None

def get_reaction_from_sabio(reaction_id):
    url = "http://sabio.h-its.org/reacdetails.jsp?reactid={}".format(reaction_id)

    headers = {
        "User-Agent": "Mozilla/5.0 (Windows NT 10.0; Win64; x64) AppleWebKit/537.36 (KHTML, like Gecko) Chrome/91.0.4472.124 Safari/537.36"
    }


    response = requests.get(url, headers=headers)

    if response.status_code == 200:

        soup = BeautifulSoup(response.content, 'html.parser')


        reaction = soup.find('td', id='StochiometricEquation')

        if reaction:
            return reaction.get_text(strip=True).replace('<->', '=')
        else:
            return "Reaction information not found."
    else:
        return f"Failed to fetch the page, status code: {response.status_code}"

EcPath = '../../data/preprocess/sabio_Ec_data'
files = os.listdir(EcPath)
skip_len = 0
skip_no_id = 0
skip_no_value = 0
skip_no_reaction = 0
reactionCaCheFile = '../../data/preprocess/reactionIdCache.json'
if os.path.exists(reactionCaCheFile):
    with open(reactionCaCheFile, 'r') as f:
        reaction_cache = json.load(f)
else:
    reaction_cache = {}
count = 0
resPath = '../../data/preprocess/sabio_data_to_json-filter_unit.json'
if os.path.exists(resPath):
    with open(resPath, 'r') as f:
        res = json.load(f)
else:
    res = {}

# reaction = get_reaction('12133', '2.5.1.31')

for file in files:
    kcat_path = os.path.join(EcPath, file)
    print(kcat_path)

    with open(kcat_path, 'r', encoding='utf-8-sig') as f:
        data = f.readlines()

    title = data[0].split('\t')
    for i in range(1, len(data)):

        line = data[i].replace('\n', '').split('\t')
        if len(line) != 18:
            skip_len += 1
            continue

        if line[0] == '':
            skip_no_id += 1
            continue
        if line[11] == '':
            skip_no_value += 1
            continue
        # columns_to_keep = [
        #     "UniprotID", "EntryID", "ECNumber", "EnzymeType", "Organism",
        #     "Reaction", "Substrate", "Product", "parameter.type",
        #     "parameter.startValue", "parameter.unit", "temperature", "pH"
        # ]
        EntryID = line[1]

        if line[14] not in ["M", "s^(-1)", "M^(-1)*s^(-1)", "mol*s^(-1)*mol^(-1)"]:
            continue

        if line[14] == 'mol*s^(-1)*mol^(-1)':
            line[14] = "s^(-1)"

        if line[15] == '-':
            line[15] = 'RT'
        if line[16] == '-':
            line[16] = 'NT'

        if EntryID not in res.keys():
            reactionId = line[17]
            if reactionId in reaction_cache.keys():
                reaction = reaction_cache[reactionId]
            else:
                reaction = get_reaction_from_sabio(line[17])
                reaction_cache[reactionId] = reaction
            if reaction is None:
                skip_no_reaction += 1
                continue
            res[EntryID] = {
                "UniprotID": line[0],
                "EntryID": line[1],
                "ECNumber": line[2],
                "EnzymeType": line[3],
                "Organism": line[4],
                "Reaction": reaction,
                "Substrate": line[6],
                "Product": line[7],
                "temperature": line[15],
                "pH": line[16],
                "Parameter": [
                    {
                        "type": line[9],
                        "value": line[11],
                        "unit": line[14],
                    }
                ]

            }
        else:
            res[EntryID]['Parameter'].append({
                "type": line[9],
                "value": line[11],
                "unit": line[14],
            })

        count += 1
        if count % 500 == 0:
            with open(reactionCaCheFile, 'w') as f:
                f.write(json.dumps(reaction_cache, indent=2))
            with open(resPath, 'w') as f:
                f.write(json.dumps(res, indent=2))

with open(resPath, 'w') as f:
    f.write(json.dumps(res, indent=2))

with open(reactionCaCheFile, 'w') as f:
    f.write(json.dumps(reaction_cache, indent=2))

