
import json
import re
from code.tools import write_utf_8_sig, read_utf_8_sig

def update_substrate_product(json_data):
    coefficient_updates = 0
    for entry_id, entry in json_data.items():
        if entry_id =="5926":
            print(111)
        reaction = entry.get("Reaction", "")
        substrate = entry.get("Substrate", "")
        product = entry.get("Product", "")


        updated_substrate = []
        for sub in substrate.split(";"):
            # match = re.search(rf'(?:^|\s)(\d*)\s*{re.escape(sub)}(?:\s|$)', sub)
            match = re.search(rf'^\s*(?:(\d+)\s+)?(.+?)\s*$', sub.strip())
            if match:

                count = int(match.group(1)) if match.group(1) else 1
                if count > 1:  # 如果系数大于1，表示有修改
                    coefficient_updates += 1
                updated_substrate.extend([match.group(2)] * count)
            else:

                updated_substrate.append(sub)
        entry["Substrate"] = ";".join(updated_substrate)


        updated_product = []
        for prod in product.split(";"):
            # match = re.search(rf'(?:^|\s)(\d*)\s*{re.escape(prod)}(?:\s|$)', prod.strip())
            match = re.search(rf'^\s*(?:(\d+)\s+)?(.+?)\s*$', prod.strip())
            if match:
                # print("ersuror  -", prod)

                count = int(match.group(1)) if match.group(1) else 1
                if count > 1:  # 如果系数大于1，表示有修改
                    coefficient_updates += 1
                updated_product.extend([match.group(2)] * count)
            else:

                updated_product.append(prod)
        entry["Product"] = ";".join(updated_product)

    return json_data, coefficient_updates


def clean_and_filter_data(json_data):

    ions_to_remove = [
        "H+", "H2", "Cu+[side 1]", "Cu+[side 2]", "Ca2+", "Ca2+/cis", "Ca2+/trans",
        "Mg2+[side 1]", "Mg2+[side 2]", "Zn2+/out", "Zn2+/in", "Li+/in", "Li+/out", "n H+/in",
        "n H+/out", "H+[side1]", "H+[side2]", "n Li+/in", "n Li+/out", "Na+/in", "Na+/out", "Ca2+"
    ]


    remove_reaction_items = [
        "Cl-", "HCl", "reduced [2Fe-2S] Ferredoxin", "Fe2+", "Fe3+", "Fe4+", "Fe5+", "Fe-coproporphyrin III", "Cu+",
        "Na+", "Na2S", "I-", "Mn2+", "MnO2+", "NH3", "SH-",
        "Hg2+", "Hg", "H2Se", "Mg2+", "C", "N", "K+", "Li+", "Zn2+", "Al3+", "Co3+", "Co", "Cr3+", "Cr6+",

    ]

    filtered_data = {}
    unknown_count = 0  #
    ion_updates = 0
    removed_reaction_count = 0
    removed_large_sub_pro_count = 0

    for entry_id, entry in json_data.items():
        # 删除包含 'Unknown' 的数据
        if "Unknown" in entry.get("Substrate", "") or "Unknown" in entry.get("Product", ""):
            unknown_count += 1
            continue

        substrate = entry.get("Substrate", "")
        original_substrate = substrate.split(";")
        cleaned_substrate = [sub for sub in original_substrate if sub not in ions_to_remove]
        ion_updates += len(original_substrate) - len(cleaned_substrate)
        entry["Substrate"] = ";".join(cleaned_substrate)


        product = entry.get("Product", "")
        original_product = product.split(";")
        cleaned_product = [prod for prod in original_product if prod not in ions_to_remove]
        ion_updates += len(original_product) - len(cleaned_product)
        entry["Product"] = ";".join(cleaned_product)


        reaction = entry.get("Reaction", "")
        substrate = entry.get("Substrate", "")
        product = entry.get("Product", "")

        if any(item in substrate.split(";") for item in remove_reaction_items) or \
                any(item in product.split(";") for item in remove_reaction_items):
            removed_reaction_count += 1
            print(entry_id)
            continue


        if len(cleaned_substrate) > 5 or len(cleaned_product) > 5:
            removed_large_sub_pro_count += 1
            print("--------------------------------")
            print(entry_id)
            print(reaction)
            continue

        filtered_data[entry_id] = entry

    return filtered_data, unknown_count, ion_updates, removed_reaction_count, removed_large_sub_pro_count



with open('../../data/preprocess/brenda_data_only_one_uniprotId.json', 'r', encoding='utf-8') as f:
    json_data = json.load(f)


updated_data, coefficient_updates = update_substrate_product(json_data)


cleaned_data, unknown_count, ion_updates, removed_reaction_count, removed_large_sub_pro_count = clean_and_filter_data(
    updated_data)


with open("../../data/preprocess/brenda_data_get_sub_pro_clean_ion_2.json", 'w') as f:
    json.dump(cleaned_data, f, indent=4)

