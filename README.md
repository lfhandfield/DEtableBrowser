# DEtableBrowser
R Shiny interface for browsing results from Differential expression results

===============

  This repository is defining a shiny server that allows to browse differential expression results produce by a pipeline that analysed single cell expression analysis in a specific framework relating to organoids with CripsR introduced mutation associated with Familial case of Alzheimer disease. The tool requires the results to be stored in a specific nested list structure within .rds files, which were produced by a snakemake pipeline that computed differential expression using both deseq2 and a modified wilcoxon test. The scripts used to produce such file is found in the Rscripts subfolder.
  
  In this experiment, mutation were introduce in either neurons or microglia, or both, and then cells were cocultured for 90 days; hence, on might look for DE gene spefific to mutation in neurons, while the background is either wild type microglia, or microglia with similar mutation. Alternatively, the mutation queried could be in microglia instead, and a total of 3 different mutation were considered. As such, many comparisons were possible, each yeilding DE genes for each celltype defined, justifying the need of a browsing tool to navigate throught the results.
  
  Though one might use a different mean to define DE genes, the approach used here is to find gene that are differentially expressed in specific celltypes which were annotated in the input dataset. Genes relevent to a given celltype are first defined as gene whose expression is higher than what is exprected from ambiant RNA, which is captured within "empty droplets", which are capture with cell barcode with little RNA. Pseudo mini-bulks are made from empty droplets and cell of a given cell type, using replicate experiements to allow DEseq2 to quantify inner/outer class variance. With such lists automatically defined, one can then evaluate the significance of DE genes bycorrect for multiple test hypothesis only on this gene subset, allowing only detect genes relevant to the celltype queried.
  
  In order to populate such results a snakemake pipeline of R scripts was used, and populated every combination of DE tests considered for this experiment, and the shine server menus were built according to the selected comparisons.
  

  

  

  
  




