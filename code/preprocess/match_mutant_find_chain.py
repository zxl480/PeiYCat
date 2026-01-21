import re
from code.tools import read_utf_8_sig, write_utf_8_sig, request_uniprot, match_mutant


expandSeqPath = '../../data/preprocess/expand_seq.json'
expandSeq = read_utf_8_sig(expandSeqPath)

resPath = '../../data/preprocess/match_mutant_find_chain.json'
uniprotCachePath = '../../data/preprocess/uniprot_entry_cache.json'
uniprotCache = read_utf_8_sig(uniprotCachePath)

res = {}

for key, item in expandSeq.items():
    print(key)

    uniprotId = item['UniprotID']
    sequenceId = 'sequence'
    chainId = 'chain'
    uni_data = uniprotCache[uniprotId]
    if uni_data is None or 'features' not in uni_data.keys():
        item[chainId] = item[sequenceId]
        continue
    item[chainId] = item[sequenceId]
    for feature in uni_data['features']:
        if feature['type'] == 'Chain':
            start = feature['location']['start']['value']
            end = feature['location']['end']['value']
            if start == 'null' or start is None:
                start = 0
            else:
                start = int(start) - 1
            if end == 'null' or end is None:
                end = len(item[sequenceId]) + 3
            else:
                end = int(end)
            item[chainId] = item[sequenceId][start:end]
            break
    res[key] = item


_res = {}
for key, item in res.items():
    print(key)

    if item['Type'] == 'wildtype':
        _res[key] = item
    else:
        uniprotId = item['UniprotID']
        sequenceId = 'sequence'
        chainId = 'chain'
        if 'commentary' in item.keys():
            EnzymeType = item['commentary']
        else:
            EnzymeType = item['EnzymeType']
        seqUnmatched, seqIdMatch = match_mutant(item[sequenceId], EnzymeType)
        if not seqUnmatched:
            item['mutant'] = seqIdMatch
        chainUnmatched, chainIdMatch = match_mutant(item[chainId], EnzymeType)
        mutantId = 'mutant'
        if not chainUnmatched:
            item[mutantId] = chainIdMatch
        if not seqUnmatched and chainUnmatched:
            uni_data = uniprotCache[uniprotId]
            item[mutantId] = seqIdMatch
            if uni_data is None or 'features' not in uni_data.keys():
                continue
            for feature in uni_data['features']:
                if feature['type'] == 'Chain':
                    start = feature['location']['start']['value']
                    end = feature['location']['end']['value']
                    if start == 'null' or start is None:
                        start = 0
                    else:
                        start = int(start) - 1
                    if end == 'null' or end is None:
                        end = len(item[sequenceId]) + 3
                    else:
                        end = int(end)
                    item[mutantId] = seqIdMatch[start:end]
                    break
    _res[key] = item


filtered_res = {}
for key, item in _res.items():
    if "mutant" in item:
        if len(item["mutant"]) > 1000:
            continue
    elif "chain" in item:
        if len(item["chain"]) > 1000:
            continue
    filtered_res[key] = item


__res = {}
for key, item in filtered_res.items():
    if item['Type'] == 'wildtype':
        __res[key] = item
    else:
        co = 0
        for _k in item.keys():
            if _k.startswith('mutant'):
                co += 1
        if co != 0:
            __res[key] = item

short_sequences = 0
for key, value in filtered_res.items():
    if "mutant" in value:
        if len(value["mutant"]) < 50:
            print(f"Entry {key} has a mutant sequence shorter than 50: {value['mutant']}")
            short_sequences += 1
    elif "chain" in value:
        if len(value["chain"]) < 50:
            print(f"Entry {key} has a chain sequence shorter than 50: {value['chain']}")
            short_sequences += 1

print(f"Total entries with sequences shorter than 50: {short_sequences}")

write_utf_8_sig(resPath, __res)












