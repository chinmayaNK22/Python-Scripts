## The script requires Proteome Database of interest in FASTA format and Peptides list from Proteome Discoverer

import re
from itertools import islice

def Fetch_Alternative_Accession(fasta_file, Peptides):
    dicts = {}
    data = open(fasta_file).readlines()
    for n, line in enumerate(data):
        if line.startswith('>'):
            data[n] = line
            a = line
            dicts[a] = []
        else:
            data[n] = line.rstrip()
            dicts[a].append(data[n])

    wirte_file = open('Test.txt', 'w')
    seq_single_line = ""
    for k,v in dicts.items():
            a = ''.join(map(str.strip,v))
            seq_single_line = seq_single_line + a.replace('+', '').rstrip()
       
    with open(Peptides) as i:
        for s in islice(i, 1, None):
            split_s = s.split('\t')
            if split_s[2] != seq_single_line:
                if split_s[2].split('.')[1] not in seq_single_line:
                    print (s.rstrip())

                    return
