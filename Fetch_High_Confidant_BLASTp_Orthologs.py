#### extract high confidant identifications from BLASTp ####

from itertools import islice

dicts ={}
with open('Coconut_predicted_proteome_BLASTp.txt') as f:
    for i in f:
        split_i = i.split('\t')
        head = split_i[0]
        if head.startswith('# Query'):
            query = head.rstrip()
            dicts[query] = []
        elif head.startswith('maker_'):
            dicts[query].append(split_i)

write1 = open('Coconut_predicted_proteome_BLASTp_formatted.txt', 'w')
write2 = open('Coconut_predicted_proteome_BLASTp_formatted_low_query_coverage.txt', 'w')
for k, v in dicts.items():
    for j in islice(v,1):
        if j[2] >= '80' and j[-1] >= '60':
            blast = "\t".join(str(x) for x in j)
            write1.write(blast)
        else:
            if j[-1] <= '60' or j[2] <= '80':
                blast1 = "\t".join(str(x) for x in j)
                write2.write(blast1)

write1.close()
write2.close()
