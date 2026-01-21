import json
import hashlib
import os

def hash_entry(entry):
    """Generate a hash for a given entry, ignoring EntryID and normalizing substrate order."""
    entry_copy = entry.copy()
    # Convert to JSON string and generate a hash
    entry_json = json.dumps(entry_copy, sort_keys=True)
    return hashlib.md5(entry_json.encode()).hexdigest()

res = {}
files = os.listdir('../../data/preprocess/brenda_Reaction_data')
for file in files:
    with open(os.path.join('../../data/preprocess/brenda_Reaction_data', file), 'r', encoding='utf-8-sig') as f:
        data = f.readlines()

    for line in data:
        line = line.replace('\n', '').replace('�', '°').split('\t')

        # print(line)
        if len(line) < 11:
            continue
        # ECNumber，Substrates, CommentarySubstrates, Organism,
        # PrimaryAccessionNo., Products, Commentary(
        #     Products), Reversibility
        ec = line[0]
        proteinName = line[1]
        reaction = "{} = {}".format(line[2], line[6])
        if "?" in reaction:
            continue
        Substrate = ';'.join(line[2].split(' + '))
        Product = ';'.join(line[6].split(' + '))
        organism = line[4]
        uniprotId = line[5]

        if uniprotId == '' or uniprotId == '-':
            continue
        hashKey = hash_entry({
            "ec": ec,
            "uniprotId": uniprotId,
            'Organism': organism
        })
        if hashKey not in res:
            res[hashKey] = [{
                "ECNumber": ec,
                "Substrate": Substrate,
                "Product": Product,
                "Organism": organism,
                "UniprotID": uniprotId,
                "proteinName": proteinName,
                "Reaction": reaction,
            }]
        else:
            _d = {
                "ECNumber": ec,
                "Substrate": Substrate,
                "Product": Product,
                "Organism": organism,
                "UniprotID": uniprotId,
                "proteinName": proteinName,
                "Reaction": reaction,
            }
            if _d not in res[hashKey]:
                res[hashKey].append(_d)
            else:
                pass


with open('../../data/preprocess/brenda_reaction_data_filtered.json', 'w', encoding='utf-8-sig') as f:
    f.write(json.dumps(res, indent=2))

