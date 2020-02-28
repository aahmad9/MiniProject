#load sleuth 
library(sleuth)

#directory 
mydir <- "/home/aahmad9/MiniProject/miniProject_Ayesha_Ahmad"
#specify where kallisto results are stored
sample_id <- dir(file.path(mydir,"results"))
#list of paths to kallisto results
kal_dirs <- file.path(mydir, "results", sample_id)
#load table describing experimental design
s2c <- read.table(file.path(mydir, "sampleinfo.txt"), header = TRUE, stringsAsFactors=FALSE)
s2c <- dplyr::select(s2c, sample = run_accession, condition)
s2c <- dplyr::mutate(s2c, path = kal_dirs)

#initalize sleuth object 
so <- sleuth_prep(s2c, extra_bootstrap_summary = TRUE)

#fit full model 
so <- sleuth_fit(so, ~condition, 'full')

#fit reduced model
so <- sleuth_fit(so, ~1, 'reduced')

#perform test
so <- sleuth_lrt(so, 'reduced', 'full')

#examine results
sleuth_table <- sleuth_results(so, 'reduced:full', 'lrt', show_all = FALSE)
sleuth_significant <- dplyr::filter(sleuth_table, qval <= 0.05)

#write to log file
t <- sleuth_significant[,c(1,4,2,3)]
write.table(t, 'miniProject.log', append=TRUE, sep='\t', dec = ".", 
row.names=TRUE, col.names=TRUE)
