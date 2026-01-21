import json

with open('../../data/preprocess/sabio_data_to_json-filter_unit.json', 'r') as f:
    data = json.load(f)

res = {}

for key, value in data.items():
    _r = {}
    for k, v in value.items():
        if k == 'Parameter':
            continue
        _r[k] = v
    p_keys = []
    for p in value['Parameter']:
        # if p['type'] == 'kcat':
        #     _r['kcat'] = p['value']
        #     _r['unit'] = p['unit']
        #     break
        p_keys.append(p['type'].lower())

    if 'kcat' in p_keys:
        for p in value['Parameter']:
            if p['type'].lower() == 'kcat' and p['unit'] == 's^(-1)':
                _r['kcat'] = p['value']
                _r['unit'] = p['unit']
                break
    elif 'kcat/km' in p_keys and 'km' in p_keys:
        kcat_km = None
        km = None
        for p in value['Parameter']:
            if p['type'].lower() == 'kcat/km':
                kcat_km = p
            if p['type'].lower() == 'km':
                km = p
        _r['kcat'] = str(float(kcat_km['value']) * float(km['value']))
        _r['unit'] = "s^(-1)"
    else:
        continue
    res[key] = _r

with open('../../data/preprocess/sabio_filter_Kcat_unify_unit.json', 'w') as f:
    json.dump(res, f, indent=2)

