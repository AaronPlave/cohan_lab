import os
import dendropy
from operator import itemgetter
import time


###OPTIMIZATION-- don't remake the tree every time!! Just pass it.
#####QUESTION -- should i recalc parents every time? If i don't right now
##it'll affect the weighted average sequence, which might not be the best thing.
#could just do a list of parent levels and go through each level as needed?



def tree_from_string(newick_file):
    return dendropy.Tree(stream=open(newick_file),
        schema="newick")

def down_search(node,dist):
    if node.is_leaf():
        return node.taxon.label,dist
    else:
        dist += 1
        return sorted(to_tuples(flatten([down_search(node.child_nodes()[0],dist),
            down_search(node.child_nodes()[1],dist)])),key=itemgetter(1))


def flatten(l, ltypes=(list, tuple)):
    """TAKEN FROM http://rightfootin.blogspot.com/2006/09/more-on-python-flatten.html
        Removes nesting from down_search"""
    ltype = type(l)
    l = list(l)
    i = 0
    while i < len(l):
        while isinstance(l[i], ltypes):
            if not l[i]:
                l.pop(i)
                i -= 1
                break
            else:
                l[i:i + 1] = l[i]
        i += 1
    return ltype(l)
    
def to_tuples(lst):
    """Takes every other n, n+1 and puts them into a tuple. I.e. ['a',2,'b',3]
    -> [('a',2),('b',3)]"""
    # print lst,len(lst)
    lst2 = []
    for i in range(0,len(lst),2):
        lst2.append((lst[i],lst[i+1]))
    return lst2

def get_parent(strain_node):
    """Gets parent of strain_node"""
    return strain_node.parent_node

def get_parents_parent_relatives(strain_node):
    """
    Returns the next closest relatives of the given strain
    for cases where there are no acceptable closest relatives
    returned by closest_relative, i.e. the sister node + its
    children. This function goes up to the parent of the 
    parent of the given strain and then searches down for 
    relatives
    """
    # tree = tree_from_string(newick_file)
    # primary_node = tree.find_node_with_taxon_label(strain)
    # print primary_node.label
    parent2 = strain_node.parent_node.parent_node
    return down_search(parent2,0),parent2


def strain_to_node(tree,strain):
    """Returns the strain object in the tree"""
    start = time.time()
    # print "strain to node took",time.time()-start
    return tree.find_node_with_taxon_label(strain)

# ASHDBASKDNAOUSD
# Make sure to handle cases where - is less than %15 but then
def closest_relative(tree,strain):
    """Returns closest relatives given a tree and a strain"""
    primary_node = tree.find_node_with_taxon_label(strain)
    # print "is PRIMARY node leaf?:",primary_node.is_leaf()
    if not primary_node:
        return (False,"PRIMARY node "+strain+" could not be found in the tree")

    if primary_node.is_leaf():
        sister_node = primary_node.sister_nodes()[0]
        return down_search(sister_node,0)
    else:
        return (False,"PRIMARY node "+primary_node+" is NOT a leaf!")

# newick_file = '31Aug aligned noninterleaved.phy_phyml_tree.nwk'
newick_file = "test.nwk"
# # # strain = "T17B4" #FAILS HERE, think this may be the outgroup and therefore 
# strain = "M11C10"


# print get_parents_parent_relatives(newick_file,'Out')

# tree = tree_from_string(newick_file)
# closest = closest_relative(newick_file,strain)
# print closest
# print closest
# if closest[0] == True:
#   print "Closest relative to",strain, "was found to be",closest
# else:
#   print "ERROR:",closest[1]

# ls = []
# ls2 = []
# print 'len tree = ',len(tree.nodes())
# print "LEAF NODES = " ,len(tree.leaf_nodes())
# for n in tree.leaf_nodes():
#   rel = closest_relative(tree,n.taxon.label)
#   if rel[0] == False:
#       ls.append(rel)
#   else:
#       ls2.append(rel)

# print len(ls2)


# print len(tree.leaf_nodes())


###COULD DO: if sister not a leaf, take children of sister if one of them is a leaf
#, else keep travelling down, maybe go up..?
# a = tree.find_node_with_taxon_label(strain)
# print a.sister_nodes()[0].child_nodes()[1].taxon.label
