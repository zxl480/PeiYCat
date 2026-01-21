import json
import re

with open('../../data/preprocess/brenda_merged_ec_reaction.json', 'r', encoding='utf-8-sig') as f:
    merged = json.load(f)
res = {}
for key, value in merged.items():
    commentary = value['commentary']
    ph_pattern = r"pH\s(\d+(\.\d+)?)"
    ph_match = re.search(ph_pattern, value['commentary'])
    if ph_match:
        value['pH'] = float(ph_match.group(1))
    else:
        value['pH'] = 'NT'

    t_pattern = r"(\d+\s*)°C"
    t_match = re.search(t_pattern, value['commentary'])
    if t_match:
        value['temperature'] = float(t_match.group(1))
    else:
        value['temperature'] = 'RT'
    res[key] = value

with open('../../data/preprocess/brenda_match_pH_temperature.json', 'w', encoding='utf-8-sig') as f:
    json.dump(res, f, ensure_ascii=False, indent=2)




