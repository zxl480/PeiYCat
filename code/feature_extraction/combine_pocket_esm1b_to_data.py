import json
import pickle

with open('./prepared_input_5_esm1b_amber_2_pH_T.json', 'r') as f:
    esm1b_data = json.load(f)

with open('./prepared_pocket_min_distance.json', 'r') as f:
    pocket_data = json.load(f)

# print(1)

new_data = []
for esm1b in esm1b_data:
    key = esm1b['reaction_id']
    print(key)
    if key not in pocket_data:
        print('not exist')
        continue
    for _key, _value in pocket_data[key].items():
        esm1b[_key] = _value

    new_data.append(esm1b)

print(len(new_data))
with open('./prepared_combine_pocket_esm1b_4A_min_distance_pH_T.json', 'w') as f:
    # pickle.dump(new_data, f
    f.write(json.dumps(new_data))
