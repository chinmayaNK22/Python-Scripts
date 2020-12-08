import read_fasta_file

fasta = read_fasta_file.read_fasta('Input_fasta_file.fasta')

dicts = {}
with open('Standard_Codon.txt') as f:
    for i in f:
        split_i = i.split('\t')
        codon = split_i[0]
        aa = split_i[1]
        if codon not in dicts:
            dicts[codon] = aa

def translate(instring):
    translated = ''
    j = ''
    for j in range(len(instring) - 2):
        #if instring[j+2] != len(instring) - 1:
            if (j % 3 ==0):
                frame = instring[j] + instring[j+1] + instring[j+2]
                if frame in dicts:
                    translated = translated  + dicts[frame]
                    j = j
    return translated

ATGC = {'A':'T', 'T':'A', 'C':'G', 'G':'C'} 
write1 = open('Output_fasta_file.fasta', 'w')
uniq_id = 1
for rows in fasta:
#### Forward translation ####
    a = 0
    b = len(rows[1].rstrip())
    for k in range(3):
        seq = rows[1].rstrip()[a:b]
        translated_seq = translate(seq)
        seq_final = translated_seq.replace('*', 'X')
        seq_id = 1
        for j in seq_final.split('X'):
            if len(j) >= 10:
                if 'K' in j or 'R' in j:
                    position2 = len(j) + translated_seq.index(j)
                    accession = 'CA' + '|' + str(uniq_id) + str(seq_id) + '_' + 'f' + str(a)
                    write1.write('>' + accession + ' ' + rows[0] + '#' + 'f' + str(a) + '#' + str((((translated_seq.index(j)+ 1)*3)-2)+a) + ':' +  str((position2 * 3)+a) + '\n' + j + '\n')
                    #print('>' + accession + ' ' + rows[0] + '#' + 'f' + str(a) + '#' + str(((translated_seq.index(j)+ 1)*3)-2) + ':' +  str(position2 * 3) + '\n' + j + '\n')
                    seq_id+= 1 
                    uniq_id+= 1
        a = a + 1

#### Reverse translation ####
    strand2 = ''.join(idx if idx not in ATGC else ATGC[idx] for idx in rows[1])
    reverse = strand2[::-1]
    c = 0
    d = len(reverse.rstrip())
    for m in range(3):
        reverse_seq = reverse.rstrip()[c:d]
        translated_reverse_seq = translate(reverse_seq)
        reverse_seq_final = translated_reverse_seq.replace('*', 'X')
        reverse_seq_id = 1
        seq_len = len(reverse.rstrip()) + 1
        for pep in reverse_seq_final.split('X'):
            if len(pep) >= 10:
                if 'K' in pep or 'R' in pep:
                    position2 = len(translated_reverse_seq) - (len(pep) + translated_reverse_seq.index(pep))
                    accession = 'CA' + '|' + str(uniq_id) + str(reverse_seq_id) + '_' + 'r' + str(c)
                    if (seq_len % 3 == 0):
                        write1.write('>' + accession + ' ' + rows[0] + '#' + 'r' + str(c) + '#' + str(((position2 * 3)-c)+3) + ':' +  str((((len(translated_reverse_seq)-(translated_reverse_seq.index(pep)))*3)+2)-c) + '\n' + pep + '\n')
                        #print ('>' + accession + ' ' + rows[0] + '#' + 'r' + str(c) + '#' + str(((position2 * 3)-c)+3) + ':' +  str((((len(translated_reverse_seq)-(translated_reverse_seq.index(pep)))*3)+2)-c) + '\n' + pep + '\n')
                    else:
                        write1.write('>' + accession + ' ' + rows[0] + '#' + 'r' + str(c) + '#' + str(((position2 * 3)-c)+2) + ':' +  str((((len(translated_reverse_seq)-(translated_reverse_seq.index(pep)))*3)+1)-c) + '\n' + pep + '\n')
                        #print ('>' + accession + ' ' + rows[0] + '#' + 'r' + str(c) + '#' + str(((position2 * 3)-c)+2) + ':' +  str((((len(translated_reverse_seq)-(translated_reverse_seq.index(pep)))*3)+1)-c) + '\n' + pep + '\n')
                    reverse_seq_id+= 1
        c = c + 1

write1.close()
