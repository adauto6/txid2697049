from Bio import Entrez  
from Bio import pairwise2  
from Bio import SeqIO  
from Bio.SubsMat import MatrixInfo as matlist  
import time
from tqdm import tqdm 
 
Entrez.email = "Your.Name.Here@example.org"  
handle = Entrez.esearch(db="nucleotide", retmax=3, term="txid2697049[Organism] and complete")  
record = Entrez.read(handle)  
handle.close()  
SeqAll = []  
print("Number of record ", len(record["IdList"]))  
for IdSeq in record["IdList"]:  
    handle = Entrez.efetch(db="nucleotide", id=IdSeq, rettype="gb", retmode="text")  
    record = SeqIO.read(handle, "genbank");  
    print("Size of Sequence ", IdSeq," ", len(record.seq))  
    SeqAll.append(record)  
  
start_time = time.time()  
matrix = matlist.blosum62  
#alignments = pairwise2.align.globalxx(SeqAll[0].seq, SeqAll[2].seq)  
for i in tqdm(range(1)): 
    alignments = pairwise2.align.globaldx(SeqAll[0].seq, SeqAll[2].seq,matrix)  
print("--- %s seconds ---" % (time.time() - start_time))                                                                                 

########### RESULTS ###############
  
#Number of record  3
#Size of Sequence  1827919722   29870
#Size of Sequence  1827919707   29851
#Size of Sequence  1827910768   29882
#  0                                                                                                           | 0/1 [00:00<?, ?it/s]

#100%|███████████████████████████████████████████████████████████████████████████████████████████████████| 1/1 [09:03<00:00, 543.90s/it]
#--- 543.9143180847168 seconds ---



#start_time = time.time()  
#   ...: matrix = matlist.blosum62  
#   ...: #alignments = pairwise2.align.globalxx(SeqAll[0].seq, SeqAll[2].seq)  
#   ...: for i in tqdm(range(1)): 
#   ...:     alignments02 = pairwise2.align.globaldx(SeqAll[0].seq, SeqAll[2].seq,matrix)  
#   ...: print("--- %s seconds ---" % (time.time() - start_time))        
#
#
#
#align1 = MultipleSeqAlignment( 
#    ...:     [ 
#    ...:         SeqRecord(Seq(alignments01[0][0],generic_dna), id="Alpha"), 
#    ...:         SeqRecord(Seq(alignments01[1][0], generic_dna), id="Beta"), 
#    ...:         SeqRecord(Seq(alignments02[0][0][0:len(alignments01[0][0])], generic_dna), id="Gamma"), 
#    ...:     ] 
#    ...: ) 
#
#
#calculator = DistanceCalculator('identity') 
#    ...:                                                                                                                                           
#
#In [66]: dm = calculator.get_distance(align1) 
#    ...:  
#    ...: # Print the distance Matrix 
#    ...: print('\nDistance Matrix\n===================') 
#    ...: print(dm)                                                                                                                                 
#
#Distance Matrix
#===================
#Alpha   0
#Beta    0.0     0
#Gamma   0.649745683308794       0.649745683308794       0
#        Alpha   Beta    Gamma
