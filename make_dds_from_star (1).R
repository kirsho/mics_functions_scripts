library(DESeq2)

make_dds_from_star <- function(mapped_dir,
                               design_csv,
                               condition_col = "Condition",
                               strand_mode = c("unstranded", "stranded_forward", "stranded_reverse")
                               ) {
  
  strand_mode <- match.arg(strand_mode)
  
  # 1. Lister les fichiers de comptes
  count_files <- list.files(mapped_dir, pattern = "ReadsPerGene.out.tab$", full.names = TRUE)
  if (length(count_files) == 0) {
    stop("Aucun fichier *ReadsPerGene.out.tab* trouvé dans ", mapped_dir)
  }
  
  message("Fichiers de comptes trouvés :")
  print(basename(count_files))
  
  # 2. Choix de la colonne selon le mode
  col_index <- switch(strand_mode,
                      unstranded       = 2,
                      stranded_forward = 3,
                      stranded_reverse = 4)
  
  read_star_counts <- function(f) {
    tab <- read.table(f, header = FALSE, stringsAsFactors = FALSE)
    srr <- sub("_.*", "", basename(f))  # SRR6825339_ReadsPerGene.out.tab -> SRR6825339
    counts <- tab[, c(1, col_index)]
    colnames(counts) <- c("gene_id", srr)
    return(counts)
  }
  
  # 3. Lire tous les fichiers et FUSIONNER avec Reduce (base R)
  count_list <- lapply(count_files, read_star_counts)
  
  # merge
  counts_merged <- Reduce(function(x, y) merge(x, y, by = "gene_id", all = TRUE),
                          count_list)
  
  # 4. Nettoyer : enlever lignes spéciales STAR
  star_special <- c("N_unmapped", "N_multimapping", "N_noFeature", "N_ambiguous")
  counts_merged <- counts_merged[!counts_merged$gene_id %in% star_special, ]
  
  # 5. Matrice de comptes
  rownames(counts_merged) <- counts_merged$gene_id
  counts_mat <- as.matrix(counts_merged[, -1])
  
  # 6. Importer design
  design <- read.csv(design_csv, stringsAsFactors = FALSE)
  rownames(design) <- design$Run
  
  # 7. Échantillons communs seulement
  common_samples <- intersect(colnames(counts_mat), rownames(design))
  if (length(common_samples) == 0) {
    stop("Aucun échantillon commun trouvé!")
  }
  
  counts_mat <- counts_mat[, common_samples, drop = FALSE]
  design     <- design[common_samples, , drop = FALSE]
  
  # 8. Vérifications
  if (!condition_col %in% colnames(design)) {
    stop("Colonne '", condition_col, "' absente du design.")
  }
  
  design[[condition_col]] <- factor(design[[condition_col]])
  
  # 9. Créer DESeqDataSet
  dds <- DESeqDataSetFromMatrix(countData = counts_mat,
                                colData   = design,
                                design    = as.formula(paste("~", condition_col)))
  
  # 10. Filtrage
  #keep <- rowSums(counts(dds)) >= min_total_count
  #dds <- dds[keep, ]
  
  message("DESeqDataSet créé: ", nrow(dds), " gènes, ", ncol(dds), " échantillons")
  return(dds)
}
