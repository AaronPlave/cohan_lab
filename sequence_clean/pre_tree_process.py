def import_fasta(filename):
    f = open(filename, 'r')
    f2= f.readlines()
    return f2

def read_fasta(fasta_list):
    dct = {}
    index = 0
    
    for line in fasta_list:
        if line[0] == '>':
            #add line as a string as key #1**
            dct[line[1:].strip()] = fasta_list[index + 1].strip()
        index += 1
    return dct

def replace_irregular(fasta_dict):
    """Replaces any nuc that is not A,T,G or C with X"""
    for strain in fasta_dict:
        pos = 0
        for nuc in fasta_dict[strain]:
            if nuc not in ['A','T','G','C','-']:
                new = list(fasta_dict[strain])
                new[pos] = "X"
                new2 = "".join(new)
                fasta_dict[strain] = new2
            pos += 1
    return fasta_dict

def to_file(filename,dct):
    new_fasta = ''
    for k in dct:
        new_fasta += '>' + k + '\n'
        new_fasta += dct[k] + '\n'
    f_out = open(filename, 'w')
    f_out.write(new_fasta)
    f_out.close()

fasta_string = import_fasta('Edited and shortened 3Sept_Aligned 31Aug.fas')
fasta_dict = read_fasta(fasta_string)
new_fasta_dict = replace_irregular(fasta_dict)
to_file('(MOD)Edited and shortened 3Sept_Aligned 31Aug.fas',new_fasta_dict)

