# **Characterization of SLC34A2 as a potential prognostic marker of tumor diseases**

The main goal of this study is to consider SLC34A2 as a potential prognostic marker of tumor diseases using mutational, expression, and survival data of cancer studies publicly available online. We have collected data from 4 databases (cBioPortal, TCGA; cBioPortal, Genie; ICGC; ArrayExpress). 111283 samples were categorized by 27 tumor locations. 

## Data Collection and Preparation
Data were obtained from 4 databases: cBioPortal, AACR Project Genie, The Insertion-al Cancer Genome Consortium (ICGC), ArrayExpress (data downloaded on 12/21/2020). Mutational and expression data were collected for SLC34A2 gene and mRNA. cBioPortal was used to access data from 32 most recent TCGA studies. The TCGA data consisted of mutational, expression and clinical data. From the AACR Project Genie database 33 studies were utilized; only mutational data were collected. The ICGC data portal was used to obtain 6 studies, also only mutational data were collected. From the ArrayExpress database the E-MTAB-3732 study was used to obtain expression data from tumor and relatively healthy samples divided by 24 groups. All collected studies and ArrayEx-press samples were categorized by 27 tumor locations: adrenal gland, biliary tract, blad-der, bone, bowel, breast, CNS and brain, cervix, esophagus and stomach, eye, head and neck, kidney, liver, lung, lymphoid, myeloid, ovary, pancreas, pleura, prostate, skin, soft tissue, testis, thymus, thyroid, uterus, and various tumors. Sample counts and study grouping is summarized in Supplementary Table S1. Different IDs other than the HGNC symbol (SLC34A2) were used: AFFY HG U133A 2 probe ID − 204124_at, and Ensembl Protein ID − ENSP00000371483, TIGRFAMs ID - TIGR01013.
## Prediction of functional impact of mutation
Analysis was conducted on SLC34A2 mutational data from 3 databases: cBioPortal, AACR Project Genie, ICGC. Only missense and indel (insertion or deletion) SLC34A2 missense and indel mutations were analyzed with tools for predicting the functional sig-nificance of mutations. Following tools were used to analyze missense mutations: PROVEAN, SIFT, PolyPhen-2, Panther-PSEP, FATHMM, Mutation Assessor. To consider a missense mutation as pathogenic it needed to be reported as such by at least five of the variant effect prediction tools. Indel mutations were analyzed with PROVEAN tool. Determination of highly conserved regions of SLC34A2 protein were performed with Conserved Domains and Protein Classification resource.
## Analysis of expression levels in various tumors
Comparison of expression levels of SLC34A2 in relatively healthy and tumor sam-ples was carried out on ArrayExpress SLC34A2 expression data. Wilcoxon test (p <0.05) was performed to compare healthy and tumor samples.
## Survival analysis
Survival analysis was performed using the Kaplan-Meyer estimate (p <0.05). For this analysis only cBioPortal TCGA studies were used. Tumor samples were divided into groups by level of SLC34A2 mRNA expression (upregulation were considered as 2 stand-ard deviation above mean) and by the presence or absence of SLC34A2 mutation.
