import json
import requests


def write_utf_8_sig(file_name, data):
    print("------saving{}".format(file_name))
    with open(file_name, 'w', encoding='utf-8-sig') as file:
        json.dump(data, file, indent=2, ensure_ascii=False)
    print("------finish")


def read_utf_8_sig(file_name):
    with open(file_name, 'r', encoding='utf-8-sig') as file:
        data = json.load(file)
    return data


def request_uniprot(uniprot_entry_id):
    try:
        request_url = 'https://rest.uniprot.org/uniprotkb/{}.json'.format(uniprot_entry_id)
        response = requests.get(request_url)
        if response.status_code == 200:
            r = json.loads(response.content)
        else:
            r = None
    except Exception as e:
        print("uniprot_Id:{}failed".format(uniprot_entry_id))
        print(e.args)

        r = None
    return r