import os

#create log file
log = open('miniProject.log', 'w')
#path to where data is 
mydir = '/homes/aahmad9/MiniProject/'

#1
###put data from file into list
f = open(mydir + 'sample.txt', 'r')
data = list(f.read().strip().split('\n'))
f.close()

#get data from NCBI
#command for getting data
get = 'wget'
#run on command line
for i in data:
    os.system(get + ' ' + str(i))

#get SRR numbers
SRR = []
for i in data:
    SRR.append(i[-10:])
#make paired end fastq data
paired = 'fastq-dump -I --split-files'
#run on command line
for i in SRR:
    os.system(paired + ' ' + i)


#2
#extract CDS from genbank and put into fasta file
from Bio import Entrez
from Bio import SeqIO
Entrez.email = 'aahmad9@luc.edu'
handle = Entrez.efetch(db='nucleotide',id='EF999921.1',rettype='gb', retmode='text')
sequence = SeqIO.read(handle,'gb')
count = 0
for f in sequence.features:
    if f.type == "CDS":
        count = count + 1
        SeqIO.write(sequence, "HCMV_cds.fasta", "fasta")

#build index with kallisto
kallisto = 'kallisto index -i HCMV_cds.idx HCMV_cds.fasta'
os.system(kallisto)
#write to log file 
log.write('The HCMV genome (EF99921) has ' + str(count) + ' ' + 'CDS.' + '\n')


