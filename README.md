# MiniProject

## Description 
This is a python wrapper to automate the execution of Bioinformatics tools on transcriptomes. 

## Software Needed 
1. SRA Toolkit: https://trace.ncbi.nlm.nih.gov/Traces/sra/sra.cgi?view=software
2. Kallisto: https://pachterlab.github.io/kallisto/download
3. Sleuth: https://pachterlab.github.io/sleuth/download
4. Bowtie2: http://bowtie-bio.sourceforge.net/bowtie2/manual.shtml#obtaining-bowtie-2
5. SPAdes: http://cab.spbu.ru/files/release3.14.0/manual.html

## Files Needed 
```
MiniProject.py
mp.R
sampleinfo.txt
test.txt
```
## Sample Data
Sample data is from four HCMV transcriptomes 2- and 6-days post infection. The ```test.txt``` data provides the urls from which to retrieve the SRA files for these four transcriptomes. The wrapper will not go through the entire fastq file that is generated from these SRA files, rather it will look at the first 50,000 lines, which is a step that is done in ```MiniProject.py```

## Instructions 
You will need the wget command to retrieve the SRA files from NCBI. Change the mydir variable in both the ```MiniProject.py``` and ```mp.R``` script to the directory in which all the downloaded software is located. All files that are generated from the wrapper are output into a folder called miniProject_Ayesha_Ahmad, so add the directory path to where that file will be located in both mydir variables. The wrapper will not execute properly if all the software is not downloaded to the same directory. Additionally, all the files must be in the same directory as well. Run the python wrapper directly on command line with the command ```python3 MiniProject.py``` and it will call the R code as well. 


