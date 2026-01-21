import requests
import os

# Extract EC number list from ExPASy, which is a repository of information relative to the nomenclature of enzymes.
def eclist():
    with open('./enzyme.dat', 'r') as outfile:
        lines = outfile.readlines()
    exists = []
    for p in os.listdir('./data/brenda_Reaction_data'):
        p = p.replace('.txt', '')
        exists.append(p)
    ec_list = list()

    for line in lines:
        if line.startswith('ID'):
            ec = line.strip().split('  ')[1]
            if ec not in exists:
                ec_list.append(ec)

    print(len(ec_list))
    return ec_list


# print(eclist())
def reaction_info(allEC):


    i = 0
    for EC in allEC:
        i += 1
        print('This is %d ----------------------------' % i)
        print(EC)
        QUERY_URL = ('https://www.brenda-enzymes.org/result_download.php?a=37&RN=&RNV=1&os=&pt=&FNV=&tt=&SYN'
                     '=&Textmining=&W[1]={}&T[1]=1&W[2]=&W[3]=&T[3]=2&V[3]=1&W[4]=&T[4]=2&W[5]=Homo+sapiens&T['
                     '5]=2&W[6]=&T[6]=2&V[6]=1&W[7]=&T[7]=2&V[7]=1&W[10]=&T[10]=2&V[10]=1&nolimit=1').format(EC)

        response = requests.get(QUERY_URL)
        results = response.text
        # print(results)
        print('---------------------------------------------')

        if results:
            with open('../../data/preprocess/brenda_Reaction_data/%s.txt' % EC, 'w', encoding='utf-8-sig') as ECfile:
                ECfile.write(results)


if __name__ == '__main__':
    allEC = eclist()
    reaction_info(allEC)



