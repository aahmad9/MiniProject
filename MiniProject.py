import os

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
##run on command line
for i in SRR:
    os.system(paired + ' ' + i)
