
library(MSnID)
library(PNNL.DMS.utils)
library(data.table)
library(dplyr)
library(MSnSet.utils)
library(tidyverse)
library(PlexedPiper)

msnid <- read_msgf_data_from_DMS(6520, param_file = "MSGFPlus_Tryp_MetOx_20ppmParTol.txt")


msnid <- apply_filter(msnid, "!grepl('Contaminant', accession)")
msnid <- apply_filter(msnid, "!isDecoy")

#msgf as data frame 
psm.df <- psms(msnid) 

#single top protein based on all spectral matches 
#Sobic.010G160700 is most abundant for all runs (by a lot)
#https://www.greenphyl.org/cgi-bin/sequence.cgi?p=id&sequence_accession=Sobic.010G160700.1&mmode=full1585862730349142975404041#tab-gene-sequence
psm.df %>%
   group_by(Dataset, Protein) %>%
   tally() %>%
   group_by(Dataset) %>%
   slice_max(n) %>%
   View()