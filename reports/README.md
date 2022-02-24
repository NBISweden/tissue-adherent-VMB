# Report templates

# NB!
  # In order to preserve the exact results as submitted to paper
  # you have to use the objects found in the data folder.
  # namly the two phyloseq objects from Gabriellas run that have been cleaned up
  # the number of samples changes slightly the normalization so that when the samples are removed after 
  # normalization they still reflect a larger dataset
  # im also not sure if the dada2 script is repoducible
  # might be differences also due to running all and a subset of the samples a s well as 
  # different package versions
  
  paths_p <- c('../../Gabriella_repo/data/phyloseq_boston_r1.RDS',   # CVL V3 (2 replicates of one tissue + 27 CVLv2)
             '../../Gabriella_repo/data/phyloseq_boston_r2.RDS')   # Tissue V3 & CVL V2')  
             
The result are reporiducible from the 01_taxonomy_improved.Rmd script, however you have to use the phyloseq objects which contains all the 111 samples from CLVL v3 and 96 samples from Tissue V3. 
