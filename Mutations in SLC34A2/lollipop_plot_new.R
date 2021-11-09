library(g3viz)
library(readr)

mutation.csv <- file.path("Mutations in SLC34A2/Result_tables/SLC43A2_lollipop_final.csv")

mutation.dat <- readMAF(mutation.csv,
                        gene.symbol.col	= "Hugo_Symbol",
                        variant.class.col = "Mutation_Type",
                        protein.change.col = "Protein_Change",
                        sep = ";")  # column-separator of csv file

plot.options <- g3Lollipop.theme(theme.name = "default",
                                 title.text = "SLC34A2",
                                 y.axis.label = "# of SLC34A2 Mutations")

g3Lollipop(mutation.dat,
           gene.symbol = "SLC34A2",
           plot.options = plot.options,
           output.filename = "cbioportal_theme")
