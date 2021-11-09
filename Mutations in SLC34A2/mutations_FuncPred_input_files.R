library(readr)
library(dplyr)
library(stringr)
library(tidyr)
library(tibble)
library(readxl)

# IMPORT
# cbioportal data
mutations_cbio <- read.table("cBioPortal/mutations.txt", sep = "\t")

# genie data
mutations_genie <- read_tsv("Genie/mutations.txt")

# ICGC data
icgc_files <- list.files("ICGC/")[2:12]
for (i in 1:length(icgc_files)) {
  
  name <- str_sub(icgc_files[i], 30, 36)
  name <- str_replace(name,"\\-", "\\_")
  mutations_icgc <- read_tsv(paste0("ICGC/", icgc_files[i])) %>%
    filter(gene_affected == "ENSG00000157765" & consequence_type == "missense_variant") %>%
    select(project_code, icgc_donor_id, aa_mutation)
  names(mutations_icgc) <- c("STUDY_ID", "SAMPLE_ID", "SLC34A2")
  if (nrow(mutations_icgc) > 0) {
    assign(name, mutations_icgc)
  }
  print(i)
  
}

# COMBINE DATA
mutations2 <-
  rbind(
    mutations_cbio,
    mutations_genie,
    BRCA_EU,
    COCA_CN,
    GACA_JP,
    LICA_CN,
    LUSC_KR,
    LICA_FR,
    MALY_DE
  )

# FILTER MUTATIONS
mutations_only <- mutations2 %>%
  filter(!is.na(SLC34A2))  %>%
  mutate(codons = SLC34A2)
mutations3 <- mutations_only

# SEPARATE MUTATIONS IF PATIENT HAVE MORE THAN ONE
for (i in 1:nrow(mutations_only)) {
  
  if (str_detect(mutations_only[i, 3], "\\s")) {
    
    string <- unlist(mutations_only[i, 3])
    string_mod <- unlist(strsplit(string, " "))
    mutations3[i, 4] <- string_mod[1]
    mutations3 <- rbind(mutations3, mutations3[i, ])
    mutations3[nrow(mutations3), 4] <- string_mod[2]
    
    if (length(string_mod) == 3) {
      mutations3 <- rbind(mutations3, mutations3[i, ])
      mutations3[nrow(mutations3), 4] <- string_mod[3]
    }
    
  }
  
}

# DELETE DUPLICATES
mutations3$codons <- unlist(mutations3$codons)
mutations3 <- distinct(mutations3)
mutations3 <- mutations3 %>%
  filter(!str_detect(codons, "SLC34A2|MUTATED"))

# RENAME AND SELECT COLUMNS
mutations_only_2 <- mutations3 %>%
  mutate(SLC34A2 = "SLC34A2")
names(mutations_only_2) <-
  c("STUDY_ID", "Sample_ID", "Hugo_Symbol", "Protein_Change")
mutations_only_2 <- mutations_only_2 %>%
  select(Hugo_Symbol, Sample_ID, Protein_Change)
# SAVE LIST FOR LOLLIPOP PLOT (CBIOPORTAL TOOL)
write_tsv(mutations_only_2, "Mutations in SLC34A2/Result_tables/3studies_lollipop_plot.txt")

# INPUT FILES FOR FUNCTIONAL PREDICTION INSTRUMENTS

# PROVEAN / SIFT
mutations_for_provean <- mutations_only_2 %>%
  mutate(for_provean = Protein_Change) %>%
  filter(!str_detect(for_provean, "\\*|\\=|SLC34A2|MUTATED|splice|del|dup|fusion|\\?")) %>%
  mutate(ref_aa = str_sub(for_provean, 1, 1)) %>%
  mutate(alt_aa = str_sub(for_provean, -1)) %>%
  mutate(codon = str_extract_all(for_provean, '\\d+', simplify = T)) %>%
  mutate(ens_gene = "ENSP00000371483") #!!!
for_provean <- unite(mutations_for_provean, for_provean, ens_gene, codon, ref_aa, alt_aa, sep = ', ')
for_provean <- for_provean[4]
write_tsv(for_provean, "Mutations in SLC34A2/Result_tables/for_provean.tsv", col_names = F)

# Polyphen-2
for_pp2 <- mutations_for_provean %>%
  mutate(Hugo_Symbol = "NPT2B_HUMAN")
for_pp2 <- unite(for_pp2, for_pp2, Hugo_Symbol, codon, ref_aa, alt_aa, sep = ' ')
for_pp2 <- for_pp2[1]
write_tsv(for_pp2, "Mutations in SLC34A2/Result_tables/for_pp2.tsv", col_names = F)

# PANTHER-PSEP
for_panther <- as.data.frame(mutations_for_provean$Protein_Change)
write_tsv(for_panther, "Mutations in SLC34A2/Result_tables/for_panther.tsv", col_names = F)

# FATHMM
for_fathmm <- mutations_for_provean %>%
  select(Protein_Change, ens_gene)
for_fathmm <- unite(for_fathmm, for_fathmm, ens_gene, Protein_Change, sep = ' ')
write_tsv(for_fathmm, "Mutations in SLC34A2/Result_tables/for_fathmm.tsv", col_names = F)

# Mutation Assessor
for_assessor <- mutations_for_provean %>%
  mutate(Hugo_Symbol = "NPT2B_HUMAN")
for_assessor <- unite(for_assessor, for_assessor, Hugo_Symbol, for_provean, sep = ' ')
for_assessor <- for_assessor[1]
write_tsv(for_assessor, "Mutations in SLC34A2/Result_tables/for_assessor.tsv", col_names = F)
