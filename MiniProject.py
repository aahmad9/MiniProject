import os


#create log file
os.system('touch miniProject.log')
#path to where data is 
mydir = '/homes/aahmad9/MiniProject/'

#1
###put data from file into list
f = open(mydir + 'test.txt', 'r')
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
handle = Entrez.efetch(db='nucleotide',id='EF999921.1',rettype='gbwithparts', retmode='text')
record = SeqIO.read(handle,'gb')
output = open('HCMV_cds.fasta', 'w')
count = 0
for f in record.features:
    if f.type == "CDS":
        count = count + 1
        f_name = f.qualifiers["protein_id"]
        f_seq = f.extract(record.seq)
        #FASTA Output
        output.write(">" + str(f_name) + "\n" + str(f_seq) + "\n" )
output.close()

#build index with kallisto
kallisto = 'kallisto index -i HCMV_cds.idx HCMV_cds.fasta'
os.system(kallisto)
#write to log file
line = 'The HCMV genome (EF99921) has '  + str(count) + ' ' + 'CDS.' + '\n'
os.system('echo  ' + " '" + line + "' >> miniProject.log")

#3
#quantify each TPM of each CDS in each transcriptome using kallisto
kallisto = 'time kallisto quant -i HCMV_cds.idx -o'
one = '_1.fastq'
two = '_2.fastq'
os.system('mkdir results')
for i in SRR:
    os.system(kallisto + ' results/' + i + ' -b 30 ' + '-t 4 ' + i + one + ' '
              + i + two)
#Run R script  
os.system('Rscript mp.R')

    



