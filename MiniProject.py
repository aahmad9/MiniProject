import os

#path to where data is 
mydir = '/Users/ayeshaahmad/Desktop/COMP_388/MiniProject/'

#1
#put data from file into list
f = open(mydir + 'sample.txt', 'r')
data = list(f.read().strip().split('\n'))
f.close()

#get data from NCBI
#command for getting data
get = 'wget'
#run on command line
for i in data:
    os.system(get + ' ' + str(i))



