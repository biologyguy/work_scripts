#!/usr/bin/python3
from Bio import Phylo
import sys
import statistics 

summ=[]

def getlimit():
    global summ
    stdevv = statistics.stdev(summ)
    meann = statistics.mean(summ)
    return stdevv+meann

    
    
if __name__ == '__main__':
    if len(sys.argv) < 4 :
        print ("#######################################################################################################")
        print ("# error: not enough arguments in the command")
        print ("#######")
        print ("# goal: trims off branches that are above 1 standard deviation")
        print ("#######")
        print ("#example command:")
        print ("# python3 tree_parse.py test_tree.tre new_tree.txt chopped_branches.txt")
        print ("#######")
        print ("# 1st argument is the tree that you want to parse")
        print ("# 2nd argument is the output file with newly parsed tree")
        print ("# 3rd argument is the output file with list of all the branches that were trimmed off")
        print ("######################################################################################################")
    else :
        args =sys.argv #renamed argument array
        tree = Phylo.read(args[1], "newick")
        f = open(args[3], 'w')
        for clade in tree.find_clades():
            summ.append(clade.branch_length)
        for clade in tree.find_clades():
            if(clade.branch_length>getlimit()):
                try :
                    stringy="%s, %s \n" % (clade.name,clade.branch_length)
                    f.write(stringy)
                    tree.collapse(target=clade)
                except ValueError:
                    print ("WARNING! Value Error!")
                    print (clade.name,clade.branch_length)            
        Phylo.write(tree, args[2], 'newick')
        f.close()
