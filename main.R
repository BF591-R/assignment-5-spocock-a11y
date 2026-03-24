library('tidyverse')
library('SummarizedExperiment')
library('DESeq2')
library('biomaRt')
library('testthat')
library('fgsea')

#' Function to generate a SummarizedExperiment object with counts and coldata
#' to use in DESeq2
#'
#' @param csv_path (str): path to the file verse_counts.tsv
#' @param metafile (str): path to the metadata sample_metadata.csv
#' @param selected_times (list): list of sample timepoints to use
#' 
#'   
#' @return SummarizedExperiment object with subsetted counts matrix
#'   and sample data. Ensure that the timepoints column used as input 
#'   to the model design has 'vP0' set as the reference factor level. Your 
#'   colData dataframe should have columns named samplename and timepoint.
#' @export
#'
#' @examples se <- make_se('verse_counts.tsv', 'sample_metadata.csv', c('vP0', 'vAd'))
make_se <- function(counts_csv, metafile_csv, selected_times) {
  #counts section
  counts <- read_delim(counts_csv, delim = "\t") %>%
    pivot_longer(cols = -gene,
                 names_to = "sample",
                 values_to = "count") %>%
    mutate(timepoint2 = str_sub(sample, 1, 3)) %>%
    filter(timepoint2 %in% selected_times) %>%
    dplyr::select(-timepoint2) %>%
    pivot_wider(names_from = sample,
                values_from = count) %>%
    column_to_rownames("gene")
  
  #sample metadata
  metafile <- read_delim(metafile_csv, delim = ",") %>%
    filter(timepoint %in% selected_times) %>%
    dplyr::select(samplename, timepoint) %>%
    mutate(timepoint = relevel(factor(timepoint), ref = "vP0")) %>%       #used claude for this
    column_to_rownames("samplename")
  
  #se object
  se <- SummarizedExperiment(
    assays  = list(counts = as.matrix(counts)),
    colData = metafile
  )                                          
  
  return(se)
}
se <- make_se('data/verse_counts.tsv', 'data/sample_metadata.csv', c('vP0', 'vAd'))
se

#' Function that runs DESeq2 and returns a named list containing the DESeq2
#' results as a dataframe and the dds object returned by DESeq2
#'
#' @param se (obj): SummarizedExperiment object containing counts matrix and
#' coldata
#' @param design: the design formula to be used in DESeq2
#'
#' @return list with DESeqDataSet object after running DESeq2 and results from
#'   DESeq2 as a dataframe
#' @export
#'
#' @examples results <- return_deseq_res(se, ~ timepoint)
return_deseq_res <- function(se, design) {
  dds <- DESeqDataSet(se, design = design)
  dds <- DESeq(dds)
  res <- results(dds)
  res_df <- as.data.frame(res)
  return(list(dds = dds, results = res_df))
}

results <- return_deseq_res(se, ~ timepoint) 
head(results$results)

#' Function that takes the DESeq2 results dataframe, converts it to a tibble and
#' adds a column to denote plotting status in volcano plot. Column should denote
#' whether gene is either 1. Significant at padj < .10 and has a positive log
#' fold change, 2. Significant at padj < .10 and has a negative log fold change,
#' 3. Not significant at padj < .10. Have the values for these labels be UP,
#' DOWN, NS, respectively. The column should be named `volc_plot_status`. Ensure
#' that the column name for your rownames is called "genes". 
#'
#' @param deseq2_res (df): results from DESeq2 
#' @param padj_threshold (float): threshold for considering significance (padj)
#'
#' @return Tibble with all columns from DESeq2 results and one additional column
#'   labeling genes by significant and up-regulated, significant and
#'   downregulated, and not significant at padj < .10.
#'   
#' @export
#'
#' @examples labeled_results <- label_res(res, .10)
label_res <- function(deseq2_res, padj_threshold) {
  deseq2_res <- as_tibble(deseq2_res, rownames = "genes") %>%
    filter(!is.na(padj)) %>%
    mutate(volc_plot_status = case_when(
      padj < padj_threshold & log2FoldChange > 0 ~ "UP",
      padj < padj_threshold & log2FoldChange < 0 ~ "DOWN",
      padj >= padj_threshold ~ "NS"))
  return(deseq2_res)
}
labeled_results <- label_res(results$results, .10)
labeled_results
#' Function to plot the unadjusted p-values as a histogram
#'
#' @param labeled_results (tibble): Tibble with DESeq2 results and one additional
#' column denoting status in volcano plot
#'
#' @return ggplot: a histogram of the raw p-values from the DESeq2 results
#' @export
#'
#' @examples pval_plot <- plot_pvals(labeled_results)
plot_pvals <- function(labeled_results) {
  plot <- ggplot(labeled_results, aes(x = pvalue)) +
    geom_histogram()
  return(plot)
}
pval_plot <- plot_pvals(labeled_results)
pval_plot
#' Function to plot the log2foldchange from DESeq2 results in a histogram
#'
#' @param labeled_results (tibble): Tibble with DESeq2 results and one additional
#' column denoting status in volcano plot
#' @param padj_threshold (float): threshold for considering significance (padj)
#'
#' @return ggplot: a histogram of log2FC values from genes significant at padj 
#' threshold of 0.1
#' @export
#'
#' @examples log2fc_plot <- plot_log2fc(labeled_results, .10)
plot_log2fc <- function(labeled_results, padj_threshold) {
  plot <- ggplot(labeled_results, aes(x=log2FoldChange)) +
    geom_histogram()
  return(plot)
}
log2fc_plot <- plot_log2fc(labeled_results, .10)
log2fc_plot
#' Function to make scatter plot of normalized counts for top ten genes ranked
#' by ascending padj
#'
#' @param labeled_results (tibble): Tibble with DESeq2 results and one
#'   additional column denoting status in volcano plot
#' @param dds_obj (obj): The object returned by running DESeq (dds) containing
#' the updated DESeqDataSet object with test results
#' @param num_genes (int): Number of genes to plot
#'
#' @return ggplot: a scatter plot with the normalized counts for each sample for
#' each of the top ten genes ranked by ascending padj
#' @export
#'
#' @examples norm_counts_plot <- scatter_norm_counts(labeled_results, dds, 10)
scatter_norm_counts <- function(labeled_results, dds_obj, num_genes){
  top_n_genes <- labeled_results %>%
    arrange(padj) %>%
    dplyr::select(genes) %>%
    head(num_genes) %>%
    as.list()
  
  norm_counts <- counts(dds_obj, normalized = TRUE) %>%
    as_tibble(rownames = "genes") 
  
  results_filtered <- labeled_results %>%
    filter(genes %in% top_n_genes$genes) %>%
    left_join(norm_counts, by = c("genes" = "genes")) %>%
    pivot_longer(
      cols = c(vP0_1, vP0_2, vAd_1, vAd_2),
      names_to = "sample",
      values_to = "normalized_count"
    ) %>%
    mutate(timepoint = str_sub(sample, 1, 3))
  
  plot <- ggplot(results_filtered, aes(x = genes, y = normalized_count, color = timepoint)) +
    geom_point() 
  return(plot)
}
norm_counts_plot <- scatter_norm_counts(labeled_results, results$dds, 10)
norm_counts_plot

#' Function to generate volcano plot from DESeq2 results
#'
#' @param labeled_results (tibble): Tibble with DESeq2 results and one
#'   additional column denoting status in volcano plot
#'
#' @return ggplot: a scatterplot (volcano plot) that displays log2foldchange vs
#'   -log10(padj) and labeled by status
#' @export
#'
#' @examples volcano_plot <- plot_volcano(labeled_results)
#' 
plot_volcano <- function(labeled_results) {
  plot <- labeled_results %>%
    mutate(sig = -log10(padj)) %>%
    ggplot(aes(x = log2FoldChange, y = sig, color = volc_plot_status)) +
    geom_point() 
  return(plot)
}
volcano_plot <- plot_volcano(labeled_results)
volcano_plot



#' Function to generate a named vector ranked by log2FC descending
#'
#' @param labeled_results (tibble): Tibble with DESeq2 results and one
#'   additional column denoting status in volcano plot
#' @param id2gene_path (str): Path to the file containing the mapping of
#' ensembl IDs to MGI symbols
#'
#' @return Named vector with gene symbols as names, and log2FoldChange as values
#' ranked in descending order
#' @export
#'
#' @examples rnk_list <- make_ranked_log2fc(labeled_results, 'data/id2gene.txt')

make_ranked_log2fc <- function(labeled_results, id2gene_path) {
  gene_table <- read_delim(id2gene_path, delim = "\t", col_names = c("ENSEMBL_ID","symbols")) %>%
    right_join(labeled_results, by = c("ENSEMBL_ID" = "genes")) %>%
    as.data.frame() %>%
    dplyr::select(symbols, log2fc = log2FoldChange) %>%
    filter(!is.na(log2fc)) %>%
    arrange(desc(log2fc)) %>% 
    deframe()
  return(gene_table)
}

rnk_list <- make_ranked_log2fc(labeled_results, 'data/id2gene.txt')
rnk_list


#' Function to run fgsea with arguments for min and max gene set size
#'
#' @param gmt_file_path (str): Path to the gene sets of interest in GMT format
#' @param rnk_list (named vector): Named vector generated previously with gene 
#' symbols and log2Fold Change values in descending order
#' @param min_size (int): Minimum number of genes in gene sets to be allowed
#' @param max_size (int): Maximum number of genes in gene sets to be allowed
#'
#' @return Tibble of results from running fgsea
#' @export
#'
#' @examples fgsea_results <- run_fgsea('data/m2.cp.v2023.1.Mm.symbols.gmt', rnk_list, 15, 500)
run_fgsea <- function(gmt_file_path, rnk_list, min_size, max_size) {
  pathways <- gmtPathways(gmt_file_path)
  
  results <- fgsea(
    pathways = pathways,
    stats    = rnk_list,
    minSize  = min_size,  
    maxSize  = max_size) %>%
    as_tibble()
  
  return(results)
}

fgsea_results <- run_fgsea('data/m2.cp.v2023.1.Mm.symbols.gmt', rnk_list, 15, 500)
fgsea_results

#' Function to plot top ten positive NES and top ten negative NES pathways
#' in a barchart
#'
#' @param fgsea_results (tibble): the fgsea results in tibble format returned by
#'   the previous function
#' @param num_paths (int): the number of pathways for each direction (top or
#'   down) to include in the plot. Set this at 10.
#'
#' @return ggplot with a barchart showing the top twenty pathways ranked by positive
#' and negative NES
#' @export
#'
#' @examples fgsea_plot <- top_pathways(fgsea_results, 10)
top_pathways <- function(fgsea_results, num_paths){
  top_n <- fgsea_results %>%
    arrange(desc(NES)) %>%
    head(num_paths)
  bottom_n <- fgsea_results %>%
    arrange(NES) %>%
    head(num_paths)
  all_pathways_plot <- bind_rows(top_n, bottom_n) %>%
    arrange(desc(NES)) %>% 
    mutate(sig = case_when(padj <= 0.05 ~ "Significant", TRUE ~ "Not Significant")) %>%
    ggplot(aes(x = pathway, y = NES)) +
      geom_col() +
      coord_flip() +
    theme(axis.text.y = element_text(size = 4))
   
  return(all_pathways_plot)
} 

fgsea_plot <- top_pathways(fgsea_results, 10)
fgsea_plot




