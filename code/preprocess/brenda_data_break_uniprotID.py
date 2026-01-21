import json
import re

uniprot_pattern = re.compile(r"(?i)[,\s;]+|(?<!\w)(?i)and(?!\w)")

if __name__ == "__main__":
    input_file = "../../data/preprocess/brenda_match_pH_temperature.json"
    cache_file = "../../data/preprocess/uniprot_entry_cache.json"
    output_file = "../../data/preprocess/brenda_data_only_one_uniprotId.json"

    with open(input_file, 'r', encoding='utf-8-sig') as f:
        _data = json.load(f)

    with open(cache_file, 'r', encoding='utf-8-sig') as f:
        cache = json.load(f)

    # 过滤 kcat ≤ 0 的数据
    data = {k: v for k, v in _data.items() if float(v.get("kcat", 0)) > 0}

    result = {}
    for key, entry in data.items():

        uniprot_ids = list(set(uniprot_pattern.split(entry["UniprotID"])))
        enzyme_ids = []

        for uid in uniprot_ids:
            entry_info = cache.get(uid, {})
            if entry_info is None:
                continue

            protein_desc = entry_info.get("proteinDescription", {})


            name_info = protein_desc.get("recommendedName", {}).get("fullName", {})
            protein_name = name_info.get("value", "")


            if not protein_name:
                submission_names = protein_desc.get("submissionNames", [])
                if isinstance(submission_names, list) and submission_names:

                    name_info = submission_names[0].get("fullName", {})
                    protein_name = name_info.get("value", "")

            if "ase" in protein_name.lower():
                enzyme_ids.append(uid)

        if enzyme_ids:
            del entry["UniprotID"]
            final_uniprot = min(enzyme_ids)
            result[key] = {"UniprotID": final_uniprot, **entry}
            continue
        if len(uniprot_ids) == 1:
            result[key] = entry
            continue
        print(key)


    with open(output_file, 'w', encoding='utf-8') as f:
        json.dump(result, f, indent=4, ensure_ascii=False)

