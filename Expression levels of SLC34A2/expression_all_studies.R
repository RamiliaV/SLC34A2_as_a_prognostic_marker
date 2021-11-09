library(readr)
library(dplyr)
library(stringr)
library(tidyr)
library(ggplot2)
library(tibble)
library(readxl)
library(rstatix)
library(survival)
library(ggfortify)
library(ggpubr)
library(survminer)

# IMPORT
# cbioportal data
clinical_data_cbio <-
  read.table("cBioPortal/combined_study_clinical_data.tsv", sep = "\t")
names(clinical_data_cbio) <- clinical_data_cbio[1,]
clinical_data_cbio <- clinical_data_cbio[-1,]
expression_cbio <- read.table("cBioPortal/mRNA expression z-scores relative to diploid samples (RNA Seq V2 RSEM).txt", sep = "\t")
names(expression_cbio) <- expression_cbio[1,]
expression_cbio <- expression_cbio[-1,]

# AE data
expression_ae <- read_csv2("ArrayExpress/slc34a2_AE_2.csv")

# survival plots - cbioportal
expression_data_new <- expression_cbio %>%
  mutate(SLC34A2_group = ifelse(SLC34A2 > 1, 1, 0))

clinical_data_join <- clinical_data_cbio %>%
  select(`Sample ID`, `Overall Survival Status`, `Overall Survival (Months)`)
full_data <- expression_data_new %>%
  left_join(clinical_data_join, by = c("SAMPLE_ID"="Sample ID"))
full_data <- full_data[complete.cases(full_data),]
full_data$Overall_Free <- as.numeric(str_sub(full_data$`Overall Survival Status`, end=1))
full_data$Overall_Months <- full_data$`Overall Survival (Months)`

# km_trt_fit <- survfit(Surv(Disease_Months, Disease_Free) ~ SLC34A2_group, data=full_data)
# surv_plot <- ggsurvplot(
#   km_trt_fit,
#   data = full_data,
#   size = 1,                 # change line size
#   palette =
#     c("#E7B800", "#2E9FDF"),# custom color palettes
#   conf.int = TRUE,          # Add confidence interval
#   pval = TRUE,              # Add p-value
#   risk.table = TRUE,        # Add risk table
#   risk.table.col = "strata",# Risk table color by groups
#   legend.labs =
#     c("Low Expression", "High Expression"),    # Change legend labels
#   risk.table.height = 0.25, # Useful to change when you have multiple groups
#   ggtheme = theme_bw(),      # Change ggplot2 theme
#   xlab = "Time in months",
#   title = "SLC34A2"
# )
# 
# surv_plot

# Boxplot - cbio - all samples
# data <- expression_cbio %>%
#   select(-2) %>%
#   gather(Gene, expression, -STUDY_ID) %>%
#   filter(!is.na(expression))
# data$STUDY_ID <- sub("\\_.*", "", data$STUDY_ID)
# acronyms <- read_excel("SLC43A2_all_tables.xlsx", sheet = "for_plots") %>%
#   select(2,3)
# acronyms$For_table <- tolower(acronyms$For_table)
# data <- data %>%
#   left_join(acronyms, by = c("STUDY_ID" = "For_table"))
# 
# stat_test <- data %>%
#   na.omit() %>%
#   dunn_test(expression ~ Study) %>%
#   filter(p.adj < 0.05)
# 
# groups <- data.frame(groups = c(stat_test$group1, stat_test$group2))
# groups_1 <- groups %>%
#   group_by(groups) %>%
#   summarise(n = n()) %>%
#   arrange(desc(n))
# 
# data %>% 
#   ggplot(aes(x = Study, y = expression)) + 
#   geom_boxplot() +
#   theme(axis.text.x = element_text(angle = 90, vjust = 0.5, size = 11), 
#         axis.text = element_text(size = 13)) +
#   labs(x = 'TCGA Study', y = 'Normalized expression') +
#   scale_x_discrete(labels = function(x) str_wrap(x, width = 20)) +
#   coord_cartesian(ylim = c(-1.2, 2.5))

# Boxplot - AE
data_samples <- read_delim("ArrayExpress/E-MTAB-3732.sdrf.txt", delim = "\t")
samples_healthy_new <- data_samples %>%
  filter(`Characteristics[disease]` == "normal") %>%
  select(1) %>%
  unlist()
samples_cancer_new <- data_samples %>%
  filter(str_detect(`Characteristics[disease]`, "cancer|carcinoma|lymphoma|blastoma|leiomyoma|myeloma|meningioma|sarcoma|leukaemia|melanoma|CCRCC|astrocytoma|glioma")) %>%
  select(1) %>%
  unlist()
data_samples_part <- data_samples %>% select(1, 7)

data_expr_updated <- expression_ae %>%
  select(-1) %>%
  gather("sample", "Expression", -`HGNC symbol`) %>%
  mutate(`HGNC symbol` = "SLC34A2")

results_type <- data_expr_updated %>%
  mutate(tissue_type = ifelse(sample %in% samples_healthy_new, "normal", ifelse(sample %in% samples_cancer_new, "cancer", "-"))) %>%
  filter(tissue_type != '-') %>%
  left_join(data_samples_part, by = c('sample' = 'Source Name')) %>%
  filter(`Characteristics[organism part]` != '-')

summary <- results_type %>%
  group_by(`Characteristics[organism part]`, tissue_type) %>%
  count()

all_tissue_type <- unique(summary$`Characteristics[organism part]`)

tissue_type_plot <- c()
for (k in 1:length(all_tissue_type)) {
  
  type_t <- all_tissue_type[k]
  found <- summary %>%
    filter(`Characteristics[organism part]` == type_t) %>%
    filter(n > 5)
  n_found <-  nrow(found)
  if (n_found > 1) {
    
    tissue_type_plot <- c(tissue_type_plot, type_t)
    
  }
  
}

summary_plot <- summary %>%
  filter(`Characteristics[organism part]` %in% tissue_type_plot)

results_plot <- results_type %>%
  filter(`Characteristics[organism part]` %in% tissue_type_plot)
results_plot <- results_plot %>%
  mutate(`Characteristics[organism part]` = ifelse(`Characteristics[organism part]` == "mammary gland", "breast", `Characteristics[organism part]`))

list_results <- data.frame(matrix(data = NA, nrow = nrow(summary_plot)/2, ncol = 3))
names(list_results) <- c("location", "pvalue_l", "pvalue_g")
list_results$location <- unique(summary_plot$`Characteristics[organism part]`)

for (i in 1:length(unique(summary_plot$`Characteristics[organism part]`))) {
  
  loc <- unique(summary_plot$`Characteristics[organism part]`)[i]
  loc_table <- results_plot %>%
    filter(`Characteristics[organism part]` == loc)
  
  result1 <- wilcox.test(Expression ~ tissue_type, data = loc_table, alternative = "l", paired = F)$p.value
  result1 <- p.adjust(result1)
  result2 <- wilcox.test(Expression ~ tissue_type, data = loc_table, alternative = "g", paired = F)$p.value
  result2 <- p.adjust(result2)
  list_results[i, 2] <- result1
  list_results[i, 3] <- result2
  
}

list_results_stat_l <- list_results %>% filter(pvalue_l <= 0.01)
list_results_stat_g <- list_results %>% filter(pvalue_g <= 0.01)

results_plot %>% 
  ggplot(aes(x = `Characteristics[organism part]`, y = Expression, fill = tissue_type)) + 
  geom_boxplot() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.3)) +
  labs(x = 'Site', y = 'Normalized expression', title = 'mRNA expression of SLC34A2, ArrayExpress') +
  scale_fill_discrete(name = 'Tissue type', labels = c('Tumor tissue', 'Normal tissue'))

site_info <- read_csv2("Expression levels of SLC34A2/site_info.csv")

expression_cbio_div <- expression_cbio %>%
  mutate(STUDY_ID = if_else(str_detect(STUDY_ID, 'tcga'), toupper(sub("\\_.*", "", STUDY_ID)), STUDY_ID)) 

expression_location <- expression_cbio_div %>%
  left_join(site_info, by = c('STUDY_ID' = 'TCGA'))

expression_loc_list <- expression_location %>%
  group_split(Site)

expression_loc_names <- expression_location %>%
  group_keys(Site) %>%
  unlist()

names(expression_loc_list) <- expression_loc_names

full_data <- full_data %>%
  mutate(STUDY_ID = if_else(str_detect(STUDY_ID, 'tcga'), toupper(sub("\\_.*", "", STUDY_ID)), STUDY_ID)) 

for (j in 1:length(expression_loc_list)) {
  
  expression_study <- expression_loc_list[[j]]
  
  study <- unique(expression_study$STUDY_ID)
  full_data_study <- full_data %>%
      filter(STUDY_ID %in% study)  
  
  full_data_study$SLC34A2_group <- factor(full_data_study$SLC34A2_group, labels = c("Low expression", "High expression"))
  ggplot(full_data_study, aes(x = SLC34A2_group, y = SLC34A2)) +
    geom_boxplot(fill = 'gray') +
    scale_y_log10() +
    labs(x = "") +
    theme_classic() +
    theme(text = element_text(size=20))
    
  if (nrow(full_data_study) > 0) {
    
    km_trt_fit_study <- survfit(Surv(Overall_Months, Overall_Free) ~ SLC34A2_group, data=full_data_study)
    pvalue <- surv_pvalue(km_trt_fit_study)$pval
    
    if (pvalue < 0.05) {
      
      surv_plot_study <- ggsurvplot(
        km_trt_fit_study,
        data = full_data_study,
        size = 1,                 # change line size
        palette =
          c("#E7B800", "#2E9FDF"),# custom color palettes
        conf.int = TRUE,          # Add confidence interval
        pval = F,              # Add p-value
        risk.table = TRUE,        # Add risk table
        risk.table.col = "strata",# Risk table color by groups
        legend.labs =
          c("Low Expression", "High Expression"),    # Change legend labels
        risk.table.height = 0.25, # Useful to change when you have multiple groups
        ggtheme = theme_bw(base_size = 35),      # Change ggplot2 theme
        fontsize = 8,
        xlab = "Overall Survival (Months)",
        ylab = "Overall Survival",
        title = paste0("SLC34A2 - ", names(expression_loc_list[j]), " study")
      )
      
      ggsave(filename = paste0(names(expression_loc_list[j]), "_SLC34A2_surv_expr2_paper.png"), plot = print(surv_plot_study, newpage = FALSE), device = "png",
               width = 8,
               height = 9)
      
    }
    
    
  }
  
}

