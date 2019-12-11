The folder contains the data used in Driver_IRW method.
There are four subfolders BRCA, HNSC, KIRC and THCA which represent four different types of cancer.
In each of subfolder, four datasets are listed.
Take the BRCA as the example, "mutation" means the processed mutation data integrated by single nucleotide variants and copy number variations.
"tumor_exp_no0" represent the expression data of tumor samples. Here, the genes with 0 in all samples are removed.
"BRCA_seeds" represents the known BRCA-related driver genes (seeds) of BRCA extracted from DisGeNet [1] and CGC [2] database.
"BRCA_diff_Dawn_net_TCGAhub" represents the processed network which is integrated by differential co-expression network and PPI network from DawnRank method [3].

References: [1] Piñero, J., Bravo, À., Queralt-Rosinach, N., Gutiérrez-Sacristán, A., Deu-Pons, J., Centeno, E., et al. (2016).
DisGeNET: a comprehensive platform integrating information on human disease-associated genes and variants. Nucleic acids research, gkw943.
[2] Sondka, Z., Bamford, S., Cole, C. G., Ward, S. A., Dunham, I., & Forbes, S. A. (2018).
The COSMIC Cancer Gene Census: describing genetic dysfunction across all human cancers. Nature Reviews Cancer, 1.
[3] Hou, J. P., & Ma, J. (2014). DawnRank: discovering personalized driver genes in cancer. Genome Medicine, 6(7), 56.
