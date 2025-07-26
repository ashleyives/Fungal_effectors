
library(MSnID)
library(PNNL.DMS.utils)
library(data.table)
library(dplyr)
library(MSnSet.utils)
library(tidyverse)
library(PlexedPiper)

#NOTE PHYTOZYME ACCESSIONS HAVE THIS FORMAT Accession.splice variant.protein, E.G. Sobic.010G160700.1.p

msnid <- read_msgf_data_from_DMS(6810, param_file = "MSGFPlus_Tryp_MetOx_20ppmParTol.txt")
msnid <- filter_msgf_data(msnid, "peptide", fdr.max = 0.01)

msnid <- apply_filter(msnid, "!grepl('Contaminant', accession)")
msnid <- apply_filter(msnid, "!isDecoy")

#msgf as data frame 
psm.df <- psms(msnid) 

#many of the IDs are from TrEMBL and overlap with known swiss-prot IDs, fill filter and remove TrEMBL IDs
psm.df <- psm.df %>%
   filter(grepl("ingel", Dataset)) %>%
   filter(grepl("SORBI", Protein)) %>%
   separate(Protein, sep="\\|", into=c("db", "Accession", "uniprotkb")) 
   # filter(db == "sp")

#single top protein based on all spectral matches 
#Sobic.010G160700 is most abundant for all runs (by a lot)
#https://www.greenphyl.org/cgi-bin/sequence.cgi?p=id&sequence_accession=Sobic.010G160700.1&mmode=full1585862730349142975404041#tab-gene-sequence
psm.df %>%
   group_by(Dataset, Accession) %>%
   tally() %>%
   group_by(Dataset) %>%
   slice_max(n) %>%
   View()

percent <- psm.df %>%
   group_by(Dataset, Accession) %>%
   # Count occurrences per group
   tally(name = "n") %>% 
   # Ungroup to calculate Dataset-level totals
   ungroup() %>%
   # Calculate total `n` per Dataset
   group_by(Dataset) %>%
   mutate(total = sum(n)) %>% 
   # Calculate percentage of `n` relative to total
   mutate(percentage = (n / total) * 100) %>%
   mutate(genotype = case_when(grepl("btx623", Dataset )~ "btx623", 
                               grepl("sc112114", Dataset )~ "sc112114"))%>%
   mutate(Elution = case_when(grepl("250mM", Dataset )~ "250mM", 
                              grepl("500mM", Dataset )~ "500mM", 
                               grepl("750mM", Dataset )~ "750mM")) 
  
psm <- percent %>%
   # Filter only the top 10 Proteins based on percentage
   arrange(desc(percentage)) %>%
   slice(1:5) %>% # Select the top 10 rows
   ggplot(aes(x = reorder(Accession, -percentage), y = percentage, fill = Elution)) +
   # geom_bar(stat = "identity", position = "dodge") + # Unstacked bars
   geom_col(position = position_dodge2(preserve = "single"))+
   facet_wrap(~ genotype, nrow = 1, scales = "free") +
   scale_fill_manual(values = c("#000000", "#E69F00", "#56B4E9")) + # Manually specify colors
   labs(x = "Accession", y = "Percentage of PSMs", title = "In-gel digestion of ~130kDa bands") +
   scale_y_continuous(expand = c(0, 0), limits = c(0,40))+
   theme_bw(base_size = 24)+
   theme(axis.text.x = element_text(angle = 45, hjust = 1))
psm
   
   ggsave(plot = psm, filename= "C:/BRAVE/figures/percentpsm_uniprot.png", 
          scale=2, 
          width = 12,
          height = 6, #maxheight with caption is 9.167 in 
          dpi = 600,
          units = c("in"))
   
#analysis of total 750 mM elutions 
#########################################################################################################
   
msnid <- read_msgf_data_from_DMS(6810, param_file = "MSGFPlus_Tryp_MetOx_20ppmParTol.txt")
msnid <- filter_msgf_data(msnid, "peptide", fdr.max = 0.01)

msnid <- apply_filter(msnid, "!grepl('Contaminant', accession)")
msnid <- apply_filter(msnid, "!isDecoy")

#msgf as data frame 
psm.df <- psms(msnid) 

psm.df <- psm.df %>%
   filter(grepl("elution", Dataset))  %>%
   filter(grepl("SORBI", Protein)) %>%
   separate(Protein, sep="\\|", into=c("db", "Accession", "uniprotkb")) %>%
   group_by(Dataset, Accession) %>%
   tally() %>%
   filter(n > 3) #must have more than 3 PSMs 
   
# List of unique proteins from each dataset
proteins_bt <- psm.df %>%
   filter(Dataset == "bt750_kdn64000_elution_01Apr25_Ned_BEHCoA-25-02-14") %>%
   pull(Accession)

proteins_sc <- psm.df %>%
   filter(Dataset == "sc750_kdn64000_elution_01Apr25_Ned_BEHCoA-25-02-14") %>%
   pull(Accession)

# Find common and unique proteins
common_proteins <- intersect(proteins_bt, proteins_sc) # Proteins shared between the two datasets
unique_to_bt <- setdiff(proteins_bt, proteins_sc)      # Proteins unique to the "bt750" dataset
unique_to_sc <- setdiff(proteins_sc, proteins_bt)      # Proteins unique to the "sc750" dataset

export <- bind_rows(
   data.frame(Accession = common_proteins, condition = "Both"), # Common proteins
   data.frame(Accession = unique_to_bt, condition = "btx623"),  # Unique to bt750
   data.frame(Accession = unique_to_sc, condition = "sc112114") # Unique to sc750
)

write.csv(export, file = "C:/BRAVE/750mmelution_overlap.csv", row.names = FALSE)

# Print results
cat("Common Proteins:\n", common_proteins, "\n")
cat("Proteins Unique to btx623:\n", unique_to_bt, "\n")
cat("Proteins Unique to sc112114:\n", unique_to_sc, "\n")

# Load library
library(VennDiagram)

# Generate a Venn Diagram using proteins from both datasets
venn.plot <- venn.diagram(
   x = list(
      btx623 = proteins_bt,       # Proteins from "bt750"
      sc112114 = proteins_sc        # Proteins from "sc750"
   ),
   filename = NULL,            # Prevent saving as a file
   col = c("blue", "red"),     # Colors for the outlines
   fill = c("lightblue", "pink"), # Fill colors
   alpha = 0.5,                # Transparency
   cex = 1.5,                  # Text size for labels
   cat.cex = 1.5,              # Text size for dataset names
   cat.col = c("blue", "red")  # Colors for dataset names
)

# Plot the Venn Diagram
grid.draw(venn.plot)
   
#GO analysis of proteins in sc versus bt 
#######################################################################

library(biomaRt)
library(gprofiler2)

# marts <- listMarts(host="https://plants.ensembl.org")
# listDatasets(useMart(biomart="plants_mart",host="https://plants.ensembl.org")) %>% View()
# 
# # Use the Ensembl Plants BioMart
# mart <- useMart(biomart = "plants_mart",
#                 host="https://plants.ensembl.org",
#                 dataset = "sbicolor_eg_gene")
# 
# listFilters(mart) %>% View() #find what you want
# 
# # Query BioMart
# results <- getBM(attributes = c("ensembl_gene_id", "uniprotsptrembl", "uniprotswissprot"),
#                  filters = "uniprotswissprot",
#                  values = "Sorghum bicolor",
#                  mart = mart)

list <- export %>%
   pull(Accession)

#to copy and paste into gprofiler website 
writeLines(list)
writeLines(common_proteins)
writeLines(c(unique_to_bt, common_proteins))
writeLines(c(unique_to_sc, common_proteins))

gostres <- gprofiler2::gost(query = list,
                            organism = "sbicolor",
                            ordered_query = F, 
                            multi_query = FALSE, 
                            significant = FALSE, #will show all data not just p < .05
                            exclude_iea = F, 
                            measure_underrepresentation = FALSE, evcodes = TRUE, 
                            user_threshold = 0.05, correction_method = "g_SCS", 
                            domain_scope = "annotated", 
                            # custom_bg = bg, 
                            numeric_ns = "", sources = c("GO"), as_short_link = FALSE)

hist(gostres$result$p_value)

#gost has a multi-test correction, see correction_method, but I can graphically see in distribution of p values that it would benefit from additional BH
gostresdf <- gostres$result %>% 
   as.data.frame() %>%
   mutate(adj_p = p.adjust(p_value, method  = "BH")) %>%
   filter(adj_p < 0.05)

#upset plot of overlap 
####################################################################
   
# library(UpSetR)
# library(tidyverse)
#    
# protein_conditions <- percent %>%
#    # Create the condition column ("genotype_elution")
#    mutate(condition = paste(genotype, Elution, sep = "_")) %>%
#    # Reduce to relevant columns
#    select(Protein, condition) %>%
#    # Mark presence of Protein in condition with 1
#    mutate(presence = 1) %>%
#    # Pivot to create the binary presence/absence matrix
#    pivot_wider(names_from = condition, 
#                values_from = presence, 
#                values_fill = 0) # Fill missing values with 0
#    
#    # Convert to binary matrix (presence/absence of Protein in each condition)
#    binary_matrix <- protein_conditions[3:8] # Remove the column with Protein names
#    binary_matrix <- binary_matrix %>%
#       mutate(across(everything(), as.numeric)) # Ensure all values are numeric (0/1)
#    
#    # Convert the row names to make proteins as identifiers
#    binary_matrix <- as.data.frame(binary_matrix)
#    # rownames(binary_matrix) <- protein_conditions$Protein
#    
#    # Generate the upset plot
#    upset(binary_matrix,
#          sets = colnames(binary_matrix),
#          order.by = "freq",
#          main.bar.color = "dodgerblue")
#    




