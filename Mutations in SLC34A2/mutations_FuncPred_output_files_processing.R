library(readr)
library(dplyr)
library(stringr)
library(tidyr)
library(ggplot2)
library(tibble)
library(readxl)

# provean
provean <- read_tsv("provean.result.tsv")

# pp2
pp2 <- read_tsv("pph2-full.result.txt")

# PANTHER-PSEP
panther <- read_delim("panther_results.txt", delim = "\t")

# FATHMM
fathmm <- read_delim("fathmm_results.txt", delim = "\t")

# Mutation Assessor
ma <- read_csv("MutationAssessor.output.csv")

#REVEL
#Mutation Taster

# joining table
pp2_pro_del_table <- provean %>%
  full_join(pp2,
            by = c(
              "POSITION" = "o_pos",
              "RESIDUE_REF" = "o_aa1",
              "RESIDUE_ALT" = "o_aa2"
            )) %>%
  select(-1) %>%
  distinct()

pp2_pro_del_table <- pp2_pro_del_table %>%
  unite('Protein_Change', RESIDUE_REF, POSITION, RESIDUE_ALT, sep = "")

full_del_table <- pp2_pro_del_table %>%
  full_join(panther, by = c('Protein_Change' = "substitution")) %>%
  full_join(fathmm, by = c('Protein_Change' = "Substitution"))  %>%
  full_join(ma, by = c('Protein_Change' = "AA variant"))

full_del_table_score <- full_del_table %>%
  mutate(score_prov = ifelse(
    `PREDICTION (cutoff=-2.5)` == 'Deleterious' &
      !is.na(`PREDICTION (cutoff=-2.5)`),
    1,
    0
  )) %>%
  mutate(score_sift = ifelse(
    `PREDICTION (cutoff=0.05)` == 'Damaging' &
      !is.na(`PREDICTION (cutoff=0.05)`),
    1,
    0
  )) %>%
  mutate(score_pp2 = ifelse(prediction == 'probably damaging' &
                              !is.na(prediction), 1, 0)) %>%
  mutate(score_panther = ifelse(message == 'probably damaging' &
                                  !is.na(message), 1, 0)) %>%
  mutate(score_fathmm = ifelse(Prediction == 'CANCER' &
                                 !is.na(Prediction), 1, 0)) %>%
  mutate(score_ma = ifelse(`Func. Impact` == 'high' &
                             !is.na(`Func. Impact`), 1, 0)) %>%
  mutate(score_sum = score_prov + score_sift + score_pp2 + score_panther + score_fathmm + score_ma)

filter_del_table <- full_del_table_score %>%
  filter(score_sum >= 4) %>%
  distinct(
    Protein_Change,
    `PREDICTION (cutoff=-2.5)`,
    `PREDICTION (cutoff=0.05)`,
    prediction,
    message,
    Prediction,
    `Func. Impact`,
    .keep_all = TRUE
  )

filter_del_table_for_paper <- filter_del_table %>%
  select(2:5, 8:9, 20, 24, 67, 68, 72, 73, 89, 90, 15:17)

# full_del_table_for_lollipop <- full_del_table %>%
#   select(2:3)
# 
# del_table_for_lollipop <- mutations_only_2 %>%
#   semi_join(full_del_table_for_lollipop)
# write_tsv(del_table_for_lollipop,
#           "3studies_lollipop__del_only_plot.txt")

gnomad_all <- read_csv("SLC34A2_gnomad.csv")
gnomad_protein <- gnomad_all %>%
  filter(!is.na(`Protein Consequence`)) %>%
  filter(Annotation == 'missense_variant') %>%
  mutate(codon = str_extract_all(`Protein Consequence`, '\\d+', simplify = T)) %>%
  select(3, 10, 16)

aa_acr <- read_excel("amino acids acr.xlsx") %>%
  select(2:3)

# ?????????? ? Gnomad
del_table_join <- filter_del_table_for_paper %>%
  left_join(aa_acr, by = c("aa1" = "Letter code 2")) %>%
  rename(ref_acr = `Letter code 1`) %>%
  left_join(aa_acr, by = c("aa2" = "Letter code 2")) %>%
  rename(alt_acr = `Letter code 1`) %>%
  unite('acr_protein', ref_acr, pos, alt_acr, sep = "") %>%
  mutate(acr_protein = paste0('p.', acr_protein))

gnomad_del_table <- del_table_join %>%
  left_join(gnomad_protein, by = c("acr_protein" = "Protein Consequence"))
write.csv(gnomad_del_table, 'gnomad_del_table.csv')
