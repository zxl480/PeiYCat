import json


def read_json(file_path):
    with open(file_path, 'r', encoding='utf-8-sig') as f:
        return json.load(f)


def write_json(data, file_path):
    with open(file_path, 'w', encoding='utf-8-sig') as f:
        json.dump(data, f, indent=4, ensure_ascii=False)



def merge_jsons(json1, json2):
    result = {}

    for hash_key, entries1 in json1.items():

        if hash_key in json2:
            entries2 = json2[hash_key]

            for entry1 in entries1:
                substrate1 = entry1.get("Substrate", "")

                for entry2 in entries2:
                    substrate2 = entry2.get("Substrate", "")


                    if substrate1 and substrate1 in substrate2:

                        new_entry = {
                            "ECNumber": entry1["ECNumber"],
                            "UniprotID": entry1["UniprotID"],
                            "Organism": entry1["Organism"],
                            "kcat": entry1["kcat"],
                            "unit": "s^(-1)",
                            "commentary": entry1.get("commentary", ""),
                            "Reaction": entry2.get("Reaction", entry2.get("reaction", "")),
                            "Substrate": entry2["Substrate"],
                            "Product": entry2["Product"]
                        }


                        result.setdefault(hash_key, []).append(new_entry)

    return result


json1_path = "../../data/preprocess/brenda_Ec_data_filtered.json"
json2_path = "../../data/preprocess/brenda_reaction_data_filtered.json"
output_path = "../../data/preprocess/brenda_merged_ec_reaction.json"


if __name__ == "__main__":

    json1_data = read_json(json1_path)
    json2_data = read_json(json2_path)


    merged_data = merge_jsons(json1_data, json2_data)


    write_json(merged_data, output_path)

    print(f"Merged data has been written to {output_path}")
