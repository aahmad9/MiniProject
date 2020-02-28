import os


#path to where output file will be stored
mydir = '/homes/aahmad9/MiniProject/miniProject_Ayesha_Ahmad/'
#create file for outputs
os.system('mkdir miniProject_Ayesha_Ahmad')
#copy files to new directory
os.system('cp test.txt ' + mydir)
os.system('cp sampleinfo.txt ' + mydir)
os.system('cp mp.R ' + mydir)
#change directory to new output file
os.chdir(mydir)
#create log file
os.system('touch miniProject.log')

#1
#put data from file into list
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
#make sample data of only the first 50000 lines of fastq files
one = '_1.fastq'
two = '_2.fastq'
for i in SRR:
    os.system('head -n 50000 ' + i + one + ' > ' + i + '_s_1.fastq')
    os.system('head -n 50000 ' + i + two + ' > ' + i + '_s_2.fastq')


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
os.system('mkdir ' + mydir + 'results')
o = '_s_1.fastq'
t = '_s_2.fastq'
for i in SRR:
    os.system(kallisto + mydir + 'results/' + i + ' -b 30 ' + '-t 4 ' + i + o + ' '
              + i + t)
#Run R script  
os.system('Rscript mp.R')


#4
#create index to map to
os.system('bowtie2-build HCMV_cds.fasta HCMV_CDS' )
#Map reads
bt = 'bowtie2 --quiet -x HCMV_CDS -1 '
for i in SRR:
    os.system(bt + i + o + ' -2 ' + i + t + ' -S ' + i + 'map.sam --al-conc-gz ' + i + '_mapped_%.fq.gz')
#write reads to file
m1 = '_mapped_1.fq.gz'
for i in SRR:
    if i == SRR[0]:
        os.system("cat " + i + o + " | echo ' Donor 1 (2dpi) had' $((`wc -l`/4)) ' read pairs before Bowtie filtering' >> miniProject.log")
        os.system("zcat " + i + m1 + " | echo ' and ' $((`wc -l`/4)) 'read pairs after' >> miniProject.log")
    elif i == SRR[1]:
        os.system("cat " + i + o + " | echo ' Donor 1 (6dpi) had' $((`wc -l`/4)) ' read pairs before Bowtie filtering' >> miniProject.log")
        os.system("zcat " + i + m1 + " | echo ' and ' $((`wc -l`/4)) 'read pairs after' >> miniProject.log")
    elif i == SRR[2]:
        os.system("cat " + i + o + " | echo ' Donor 3 (2dpi) had' $((`wc -l`/4)) ' read pairs before Bowtie filtering' >> miniProject.log")
        os.system("zcat " + i + m1 + " | echo ' and ' $((`wc -l`/4)) 'read pairs after' >> miniProject.log")
    elif i == SRR[3]:
        os.system("cat " + i + o + " | echo ' Donor 3 (6dpi) had' $((`wc -l`/4)) ' read pairs before Bowtie filtering' >> miniProject.log")
        os.system("zcat " + i + m1 + " | echo ' and ' $((`wc -l`/4)) 'read pairs after' >> miniProject.log")
    else:
        break

#5
#Run spades
spades = 'spades -k 55,77,99,127 -t 4 --only-assembler --pe1-1 SRR5660030_mapped_1.fq.gz --pe1-2 SRR5660030_mapped_2.fq.gz --pe2-1 SRR5660033_mapped_1.fq.gz --pe2-2 SRR5660033_mapped_2.fq.gz --pe3-1 SRR5660044_mapped_1.fq.gz --pe3-2 SRR5660044_mapped_2.fq.gz --pe4-1 SRR5660045_mapped_1.fq.gz --pe4-2 SRR5660045_mapped_2.fq.gz -o Assembly1'
os.system( "'" + spades + "'" )
os.system('echo ' + "'" + spades + "' >> miniProject.log")

#6 and 7
#calculate number of contigs with length > 1000
#read contigs.fasta file and put in list
c =[]
h = open(mydir + 'Assembly1/contigs.fasta')
for record in SeqIO.parse(h, 'fasta'):
    c.append(record)
h.close()

#count number of contigs with length greater than 1000
c1 =0
bp = 0
contigs = []
for contig in c:
    if len(contig.seq) > 1000:
        if len(contig.seq) > 1000:
            c1 = c1 + 1
            bp= bp+ len(contig.seq)
            contigs.append(contig)
    else:
        continue

#write to log file 
os.system(' echo There are ' + "'"  + str(c1) + "'" + ' contigs > 1000 in the assembly >> miniProject.log')
os.system('echo There are ' + "'" + str(bp) + "'" + ' bp in the assembly >> miniProject.log')

    
#8
#rewrite the contents of contigs.fasta to only include the contigs with length > 1000
h2 = open('contigs2.fasta', 'w')
N = 'NNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNN'
for contig in contigs:
    h2.write(str(contig.seq) + N)
h2.close()

#9
#Blast results
from Bio.Blast import NCBIWWW
fasta_string = open('contigs2.fasta').read()
result_handle = NCBIWWW.qblast('blastn', 'nt', fasta_string, entrez_query = "txid10292[Organism:exp]")
with open("myblast.xml", "w") as out_handle:
    out_handle.write(result_handle.read())
result_handle.close()
#read results
from Bio.Blast import NCBIXML
rh = open('myblast.xml', 'r')
blast_records = NCBIXML.parse(rh)

p = open('tophits.txt', 'w')
p.write('seq_title align_len  number_HSPs  topHSP_ident topHSP_gaps  topHSP_bits topHSP_expect \n')
for i, blast_record in enumerate(blast_records):
    if i == 10: break
    for alignment in blast_record.alignments:
        for hsp in alignment.hsps:
            p.write(str(alignment.title) + '\t' + str(alignment.length) + '\t' + str(hsp.num_alignments) + '\t' + str(hsp.gaps) + '\t' + str(hsp.bits) + '\t' + str(hsp.expect) + '\n')
p.close()

os.system("sed -n '1,11p' tophits.txt > th.txt")
os.system('cat th.txt >> miniProject.log')

print('You are done!')











    

        



    

        
    

