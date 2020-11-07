write1 = open('Human_EGF_Proteins.fasta', 'w')

dicts = {}
data = open('HumanRefSeq_GCF_000001405.39_GRCh38.p13_protein_refseq109_formated.fasta').readlines()
for n, line in enumerate(data):
    if line.startswith('>'):
        data[n] = line
        a = line
        dicts[a] = []
    else:
        data[n] = line.rstrip()
        dicts[a].append(data[n])

for k,v in dicts.items():
    try:
        split_k = k.split('|')
        v1 = ''.join(map(str.strip,v))
        if "epidermal growth factor" in k:
            write1.write(split_k[0]+'|'+split_k[1]+'|'+split_k[2]+'|'+ split_k[3] + '|' + split_k[4]+ v1 + '\n')
    except:
        print (k)
write1.close()
