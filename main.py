from Bio.Alphabet import generic_dna 
from Bio.Seq import Seq 
from Bio.SeqRecord import SeqRecord 
from Bio.Align import MultipleSeqAlignment 
from Bio.Align.Applications import MafftCommandline                                                                                 
from Bio import AlignIO
from Bio import Phylo 
from Bio.Phylo.TreeConstruction import DistanceCalculator 
from Bio.Phylo.TreeConstruction import DistanceTreeConstructor 

import pylab

from Bio import Entrez  
from Bio import pairwise2  
from Bio import SeqIO  
from Bio.SubsMat import MatrixInfo as matlist  
import time  
from tqdm import tqdm 
 
start_time = time.time()  
Entrez.email = "Your.Name.Here@example.org"  
handle = Entrez.esearch(db="nucleotide", retmax=1000, term="txid2697049[Organism] and complete")  
record = Entrez.read(handle)  
handle.close()  
my_alignments = [] #align1, align2, align3]
print("Number of record ", len(record["IdList"]))  
for IdSeq in record["IdList"]:  
    handle = Entrez.efetch(db="nucleotide", id=IdSeq, rettype="gb", retmode="text")  
    record = SeqIO.read(handle, "genbank");
    align1 = MultipleSeqAlignment([ SeqRecord(record.seq, id=record.id) ])
                                                                                                                                       
    my_alignments.append(align1)

    print("Size of Sequence ", IdSeq," ", len(record.seq))  



 
#align1 = MultipleSeqAlignment( 
#    [ SeqRecord(SeqAll[i].seq, id=SeqAll[i].id)  for i in range(2,10) ] 
        #)                                        
#my_alignments = [align1, align2, align3]

AlignIO.write(my_alignments, "my_example.fasta", "fasta")

in_file   = "my_example.fasta" 


print('\nStarting MAFFT\n===================')  
mafft_cline = MafftCommandline(input=in_file, thread=25) 
stdout, stderr = mafft_cline() 
with open("aligned.fasta", "w") as handle: 
    handle.write(stdout) 
from Bio import AlignIO 
 
print('\nMAFFT Ended !!! \n===================')  

aln = AlignIO.read("aligned.fasta", "fasta") 

calculator = DistanceCalculator('identity') 

dm = calculator.get_distance(aln) 
  
# Print the distance Matrix 
print('\nDistance Matrix\n===================') 
print(dm)                                                                                                                                

# Construct the phylogenetic tree using UPGMA algorithm 
constructor = DistanceTreeConstructor() 
tree = constructor.upgma(dm) 
 
# Draw the phylogenetic tree 
Phylo.draw(tree) 
 
# Print the phylogenetic tree in the terminal 
print('\nPhylogenetic Tree\n===================') 
Phylo.draw_ascii(tree)
Phylo.draw(tree)
pylab.savefig('apaf.png')
print("--- %s seconds ---" % (time.time() - start_time))
# 348
#--- 1466.3146817684174 seconds ---



>>> import sys
>>> Phylo.write(tree, sys.stdout, "phyloxml")

tree = Phylo.read('tree1.xml', 'phyloxml') 
   ...: print(tree)

Phylo.write([tree], 'example-both.tree', 'newick')  
