

#start at .mzML and run through MSGF on PC 
#loops over files and performs msgf_plus search on each file, convert from .mzid to tsv in separate loop 
#uses the default mods.txt from msgf download, need to optimize further 

# Define the MSGFPlus parameters
memory_allocation <- "8000M"  # Adjust this according to your system's memory, 4000M means 4 GB

tolerance <- "10ppm"

frag_mode <- "0"
# [-m FragmentationMethodID] (Fragmentation Method; Default: 0)
# 0: As written in the spectrum or CID if no info
# 1: CID
# 2: ETD
# 3: HCD
# 4: UVPD

instrument <- "1" 
# [-inst InstrumentID] (Instrument ID; Default: 0)
# 0: Low-res LCQ/LTQ
# 1: Orbitrap/FTICR/Lumos
# 2: TOF
# 3: Q-Exactive
enzyme <- "1"
# [-e EnzymeID] (Enzyme ID; Default: 1)
# 0: Unspecific cleavage
# 1: Trypsin
# 2: Chymotrypsin
# 3: Lys-C
# 4: Lys-N
# 5: glutamyl endopeptidase
# 6: Arg-C
# 7: Asp-N
# 8: alphaLP
# 9: no cleavage

isotope_error <- "-1,2"
# [-ti IsotopeErrorRange] (Range of allowed isotope peak errors; Default: 0,1)
# Takes into account the error introduced by choosing a non-monoisotopic peak for fragmentation.
# The combination of -t and -ti determines the precursor mass tolerance.
# E.g. "-t 20ppm -ti -1,2" tests abs(ObservedPepMass - TheoreticalPepMass - n * 1.00335Da) < 20ppm for n = -1, 0, 1, 2.

num_tolerable_termini <- "1"
# [-ntt 0/1/2] (Number of Tolerable Termini; Default: 2)
# When EnzymeID is 1 (trypsin),
# 2: Only search for fully-tryptic peptides
# 1: Search for semi-tryptic and fully-tryptic peptides
# 0: Non-tryptic search
# Number of tolerable termini

decoy_database <- "1"
# [-tda 0/1] (Target decoy strategy; Default: 0)
# 0: Don't use a decoy database
# 1: Search with a decoy database (forward + reverse proteins)

min_length <- "6"
max_length <- "50" #can also specify a min and max charge, see syntax document 
num_results <- "1"
#  [-n NumMatchesPerSpec] (Number of matches per spectrum to be reported; Default: 1)

num_threads <- "7"
# [-thread NumThreads] (Number of concurrent threads to be executed; Default: Number of available cores)
# This is best set to the number of physical cores in a single NUMA node.
# Generally a single NUMA node is 1 physical processor.
# The default will try to use hyperthreading cores, which can increase the amount of time this process will take.
# This is because the part of Scoring param generation that is multithreaded is also I/O intensive.


# Define the paths to your files and the Java executable
java_exe <- "C:/Java/jre1.8.0_441/bin/java.exe"
msgfplus_jar <- "C:/MSGFPlus_v20240326/MSGFPlus.jar"
database_file <- "C:/FASTA/ID_008654_8377E8C8.fasta" #Sorghum_bicolor_JGI_v5.1_2024-05-06; Contaminants_2006-08-01; Tryp_Pig_Bov
mods_file <- "C:/MSGFPlus_v20240326/MSGFPlus_FungalEffectors.txt" #N term acetyl, Met Ox, samples are not alkylated 

# Get all .mzML files
mzML_files <- unique(list.files(path = "C:/BRAVE", pattern = "\\.mzML$", full.names = TRUE, recursive = TRUE))

for (x in mzML_files) {
   print(x)
   
   input_file <- paste0("C:/BRAVE/", basename(x))
   output_file <- paste0("C:/BRAVE/output/", sub("\\.mzML$", "", basename(x)), ".mzid")
   
   # Construct the command
   cmd <- sprintf('"%s" -Xmx%s -jar %s -s %s -o %s -d %s -t %s -m %s -inst %s -e %s -ti %s -ntt %s -tda %s -minLength %s -maxLength %s -n %s -thread %s -mod %s',
                  java_exe, memory_allocation, msgfplus_jar, input_file, output_file, database_file, tolerance, frag_mode, instrument, enzyme, isotope_error, num_tolerable_termini, decoy_database, min_length, max_length, num_results, num_threads, mods_file)
   cmd
   
   # Run the command using system2 with error handling
   tryCatch({
      system2(command = "cmd.exe", args = c("/c", cmd), stdout = TRUE, stderr = TRUE)
   }, error = function(e) {
      cat("Error while processing", x, "\n")
      cat("Error message:", e$message, "\n")
      next
   })
}

#convert .mzid output from msgf to tsvs
#####################################################################

mzidlist <- unique(list.files(path = "C:/BRAVE/output", pattern = "\\.mzid$", full.names = TRUE, recursive = TRUE))

for(x in unique(mzidlist)){
   print(x)
   
   exe_path <- "C:\\Users\\ives435\\Documents\\net48\\MzidToTsvConverter.exe"
   mzid_path <- paste0("C:\\BRAVE\\output\\", basename(x))
   
   
   args <- c(mzid_path, "-showDecoy", "0", "-unroll", "0")
   
   system2(exe_path, args)
   
}


#read tsvs back in to inspect results 

tsvlist <- unique(list.files(path = "C:/BRAVE/output", pattern = "\\.tsv$", full.names = TRUE, recursive = TRUE))


data <- data.frame()

for(x in unique(tsvlist)){
   print(x)
   
   current_data <- read.delim(x)
   
   data <- rbind(data, current_data)
   
}

library(tidyverse)

# Function to read and parse the FASTA file
# parse_fasta <- function(fasta_file) {
#   # Read the FASTA file
#   fasta_lines <- readLines(fasta_file)
# 
#   # Initialize an empty data frame to store parsed data
#   parsed_df <- data.frame(
#     identifier = character(),
#     sp = character(),
#     id = character(),
#     name = character(),
#     description = character(),
#     stringsAsFactors = FALSE
#   )
# 
#   # Parse each line
#   for (line in fasta_lines) {
#     if (startsWith(line, ">")) {
#       # Split the identifier line at each instance of "|"
#       split_line <- strsplit(sub("^>", "", line), split = "\\|")[[1]]
# 
#       # Extract columns
#       sp <- split_line[1]
#       id <- split_line[2]
#       description <- split_line[3]
# 
#       # Extract additional information from the description
#       name <- strsplit(description, " ", fixed = TRUE)[[1]][1]
# 
#       # Add to the data frame
#       parsed_df <- rbind(parsed_df, data.frame(identifier = line, sp = sp, id = id, name = name, description = description))
#     }
#   }
# 
#   return(parsed_df)
# }
# 
# parsed_fasta <- parse_fasta("C:/FASTA/HomosapiensPlusAAP9.fasta") 

#there both CID and ETD runs so it correctly parsed the scan types 
data2 <- data %>%
   filter(!grepl("rev_sp", Protein)) %>%
   filter(EValue < 0.1) #checked scan at EValue ~0.03 and there are 3 y-ions 

write.csv(data2, file="C:/Nakai/gluc_hits.csv")



data3 <- data %>%
   filter(!grepl("rev_sp", Protein)) %>%
   filter(grepl("1302", Peptide))