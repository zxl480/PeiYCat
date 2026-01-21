import json
import hashlib
import os

def hash_entry(entry):
    entry_copy = entry.copy()
    # Convert to JSON string and generate a hash
    entry_json = json.dumps(entry_copy, sort_keys=True)
    return hashlib.md5(entry_json.encode()).hexdigest()

res = {}
files = os.listdir('../../data/preprocess/brenda_Ec_data')
for file in files:
    with open(os.path.join('../../data/preprocess/brenda_Ec_data', file), 'r', encoding='utf-8-sig') as f:
        data = f.readlines()

    for line in data:
        line = line.replace('\n', '').replace('�', '°').split('\t')

        # print(line)
        if len(line) < 10:
            continue
        # EC, turnovernumber,turnovernumberMax, substrate, commentary, organism, primaryAccession
        ec = line[0]
        turnover_number = line[1]
        turnover_number_max = line[2]
        substrate = line[3]
        commentary = line[4]
        organism = line[5]
        uniprotId = line[6]

        if uniprotId == '' or turnover_number == 'additional information' or uniprotId == '-':
            continue
        hashKey = hash_entry({
            "ec": ec,
            "uniprotId": uniprotId,
            'Organism': organism
        })
        if hashKey not in res:
            res[hashKey] = [{
                "ECNumber": ec,
                "kcat": turnover_number,
                "kcat_max": turnover_number_max,
                "Substrate": substrate,
                "commentary": commentary,
                "Organism": organism,
                "UniprotID": uniprotId
            }]
        else:
            res[hashKey].append({
                "ECNumber": ec,
                "kcat": turnover_number,
                "kcat_max": turnover_number_max,
                "Substrate": substrate,
                "commentary": commentary,
                "Organism": organism,
                "UniprotID": uniprotId
            })


with open('../../data/preprocess/brenda_Ec_data_filtered.json', 'w', encoding='utf-8-sig') as f:
    f.write(json.dumps(res, indent=2))


