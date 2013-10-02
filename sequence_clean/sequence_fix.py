import closest_rel
import random
import time


##T5c10 has underscore, have to specify some params in closest_rel so things
##like that are accepted..

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

def percent_gap(dct,index):
    """calculates percentage of strains that also have a gap at the same index
    DO I COUNT THE STRAIN WE'RE LOOKING AT?"""
    # start = time.time()
    num_strains = len(dct)
    gaps = 0
    for s in dct:
        #try statement to catch the cases at the ends of seqs where 
        #some seqs will have gaps and other won't have anything
        try:
            if dct[s][index] == "-":
                gaps += 1

        #COUNTING HAVING NO NUCLEOTIDE AS A GAP
        except IndexError: 
            gaps += 1

    # print "%/ gap took ",time.time()-start
    return float(gaps)/num_strains

def remove_column(dct,index):
    """Removes the nucleotide at the given index for
        all sequences"""
    # start = time.time()
    ## COULD OPTIMIZE BY MAYBE NOT REJOINING? KEEP LISTS UNTIL END? 
    # print "DOING NEW"
    for s in dct:
        seq = list(dct[s])
        try:
            seq.pop(index)
        except:
            print "PASSSSSSSSSSSSSSSSSSSSSSSSSSSSSSING"
            pass
        dct[s] = "".join(seq)
    # print "Remove column took", time.time()-start
    return dct

def weighted_seq_choose(strain_dct,strain_node,rels,index,level=0):
    ##LEVEL is how many upward requests have been requested. Idea is that you
    ##store each level of ancestors and only get more if levels > len(rels)
    ##This is to avoid the problem of using all the levels for each replacement
    ##since that sacrifices closeness and makes some of the replacements less acceptable

    start = time.time()
    nucs = {}
    if not rels[0]:
        print "NO RELATIVES!!!!!!!!!!!!!!!!!!!!!!!!!"
    # print rels,"RELATIVES"
    print rels,"RELS"
    print level,"LEVEL",rels[level]
    for r in rels[level]:
        # if r[0] == True:
        try:
            nucleotide = strain_dct[r[0]][index]
            if nucleotide in nucs:
                nucs[nucleotide] += 1
            else:
                nucs[nucleotide] = 1
        except:
            print r[0],"IS A PROBLEM STRAIN, KEY ERROR!"

        # else:
        #     print "FALSE",r
    # print nucs,"NUCS DCT"

    pop = []
    for z in nucs:
        if z in ["A","T","G","C"]:
            for i in range(nucs[z]):
                pop.append(z)
    if not pop:
        print "Nothing was an A,T,G,C, have to go up the tree#############################3a"  
        ##Have to get the parent of the parent of the strain
        ##HOW MANY TIMES SHOULD THIS BE REPEATED??
        ##MAKE SURE THIS IS ACTUALLY WORKING/GOING UP THE TREE MORE THAN ONCE
        # try:
        try:
            print "already have"
            rels_next = rels[level+1]
            print "got next, need parent now"
            parent2 = closest_rel.get_parent(strain_node)
            print "got parent2"
        except:
            print "need new"
            rels_next,parent2 = closest_rel.get_parents_parent_relatives(strain_node)  
            rels.append(rels_next)
        print parent2,strain_node 
        print rels,"I AM NEW RELS AT LEVEL",level
        return weighted_seq_choose(strain_dct,parent2,rels,index,level+1)
        # except:
        #     ######WHAT SHOULD HAPPEN IN THIS CASE? Second pass?##########
        #     print "No parent"
        #     return None
    # print pop,"POP"
    choice = random.choice(pop)
    print rels
    print "weighted_seq_choose took",time.time()-start
    return choice,rels

def replacer(newick_file,dct):
    tree = closest_rel.tree_from_string(newick_file)
    orig_dct = dct.copy()
    bad_list = []
    s_num = 0
    for strain in dct.keys():
        print "STRAIN NUMBER",s_num
        index = 0
        print "STRAIN ("+strain+")",dct[strain]
        rels = [closest_rel.closest_relative(tree,strain)]
        if type(rels[0]) == tuple:
            #for case of just one relative
            a = [[]]
            a[0].append(rels[0])
            rels = a
        sequence_copy = str(dct[strain])
        strain_node = ""
        for nuc in sequence_copy:
            if nuc == '-':
                p_gaps = percent_gap(dct,index)
                if p_gaps <= .15:
                    if strain_node == "":
                        strain_node = closest_rel.strain_to_node(tree,strain)
                    print "Strain_node!!",strain_node

                    closest,rels = weighted_seq_choose(orig_dct,strain_node,rels,index)
                    ###Could optimize this bit here, save closest relatives and just do op once, index later
                    ###BUT now have to create a hierarchy of rel levels so that each replacement will
                    ###only jump up IF necessary, since otherwise you sacrifice closeness.
                    if closest:
                        old_seq = list(dct[strain])
                        old_seq[index] = closest
                        dct[strain] = "".join(old_seq)
                        # print "DCTCLS: ",dct[strain]
                    else:
                        print "FOUND NO CLOSEST RELATIVE EVER...?"
                    index += 1
                else:
                    # print "REMOVING COLUMN"
                    # print "a:",dct[strain]
                    dct = remove_column(dct,index)
                    orig_dct = remove_column(orig_dct,index)
                    # print dct == orig_dct
                    # print "b:",dct[strain]
            elif nuc not in ["A","T","G","C"]:
                # print "DOING REPLACEMENT"
                strain_node = closest_rel.strain_to_node(tree,strain)
                print "Strain_node!!",strain_node
                closest,rels = weighted_seq_choose(orig_dct,strain_node,rels,index)
                if closest:
                    old_seq = list(dct[strain])
                    old_seq[index] = closest
                    dct[strain] = "".join(old_seq)
                index += 1  
            else:
                index += 1  
        s_num += 1

    return dct

def to_file(dct):
    new_fasta = ''
    for k in dct:
        new_fasta += '>' + k + '\n'
        new_fasta += dct[k] + '\n'
    f_out = open('SEQ_R_FASTA.fas', 'w')
    f_out.write(new_fasta)
    f_out.close()

#main calls
x = time.time()
newick_file = 'edit_test.nwk'
fasta_string = import_fasta('1Edited and shortened 3Sept_Aligned 31Aug.fas')
fasta_dict = read_fasta(fasta_string)
fasta_dict_copy = fasta_dict.copy()
new_fasta_dict = replacer(newick_file,fasta_dict_copy)
z = time.time()
print "RUN TIME:",z-x
to_file(new_fasta_dict)
### When doing replacement, get all the non '-' sequences from the relatives
### FROM the VANILLA sequence list and do a weighted probability to determine
### which nucleotide to use. IF all gaps, do a search from the node above the primary node 
### 
###

