
##Functions scores from https://github.com/olapuentesantana/easier_manuscript

#' Compute Expanded Immune signature
#'
#' \code{compute_ayersEI} computes Expanded Immune signature score as the arithmetic mean of genes included
#' in the Expanded Immune signature (Ayers et al., JCI, 2017)
#'
#' @importFrom stats na.omit 
#'
#' @param RNA.tpm numeric matrix with rows=genes and columns=samples Log2 transformed
#'
#' @return numeric matrix with rows=samples and columns=Expanded Immune signature score
#'
#' @export
#-------------------------------------------------------------------------------------------------------------

compute.Ayers_expIS <- function(RNA.tpm){
  # Literature genes
  Ayers_expIS.read <- c("GZMB", "GZMK", "CXCR6", "CXCL10", "CXCL13", "CCL5", "STAT1","CD3D", "CD3E",
                        "CD2", "IL2RG" , "NKG7", "HLA-E", "CIITA","HLA-DRA", "LAG3", "IDO1", "TAGAP")
  match_Ayers_expIS.genes <- match(Ayers_expIS.read, rownames(RNA.tpm))
  
  if (anyNA(match_Ayers_expIS.genes)){
    warning(c("differenty named or missing signature genes : \n", paste(Ayers_expIS.read[!Ayers_expIS.read %in% rownames(RNA.tpm)], collapse = "\n")))
    match_Ayers_expIS.genes <- stats::na.omit(match_Ayers_expIS.genes)
  }
  
  # Log2 transformation:
  log2.RNA.tpm <- RNA.tpm
  
  # Subset log2.RNA.tpm
  sub_log2.RNA.tpm  <- log2.RNA.tpm[match_Ayers_expIS.genes, ]
  
  # Calculation: average of the included genes for Expanded Immune signature
  score <- apply(sub_log2.RNA.tpm, 2, mean)
  
  message("Ayers_expIS score computed")
  return(data.frame(Ayers_expIS = score, check.names = FALSE))
}


#' Compute cytolytic activity score
#'
#' \code{compute_CYT} computes cytolytic activity score as the geometric mean of immune cytolytic genes
#' (Rooney et al., 2015).
#'
#' @importFrom stats na.omit
#'
#' @param RNA.tpm numeric matrix with rows=genes and columns=samples
#'
#' @return numeric matrix with rows=samples and columns=cytolytic activity score
#'
#' @export
#-------------------------------------------------------------------------------------------------------------

compute.CYT <- function(RNA.tpm){
  
  # Literature genes
  CYT.read <- c("GZMA", "PRF1")
  match_CYT.genes <- match(CYT.read, rownames(RNA.tpm))
  
  if (anyNA(match_CYT.genes)){
    warning(paste0("differenty named or missing signature genes : \n", paste(CYT.read[!CYT.read %in% rownames(RNA.tpm)], collapse = "\n")))
    match_CYT.genes <- stats::na.omit(match_CYT.genes)
  }
  
  # Subset RNA.tpm
  subset_RNA.tpm <- RNA.tpm[match_CYT.genes, ]
  
  # Calculation: geometric mean (so-called log-average) [TPM, 0.01 offset]
  score <- as.matrix(apply(subset_RNA.tpm + 0.01, 2, function(X) exp(mean(log(X)))))
  
  message("CYT score computed")
  return(data.frame(CYT = score, check.names = FALSE))
}

#' Compute Davoli immune signature
#'
#' \code{compute_davoliIS} computes Davoli immune signature as the arithmetic mean of cytotoxic
#' immune infiltrate signature genes, after rank normalization (Davoli et al., 2017).
#'
#' @importFrom stats na.omit
#'
#' @param RNA.tpm numeric matrix with rows=genes and columns=samples
#'
#' @return numeric matrix with rows=samples and columns=Davoli immune signature
#'
#' @export
#-------------------------------------------------------------------------------------------------------------

compute.Davoli_IS <- function(RNA.tpm){
  
  # Literature genes
  Davoli_IS.read <- c("CD247", "CD2", "CD3E", "GZMH", "NKG7", "PRF1", "GZMK")
  match_Davoli_IS.genes <- match(Davoli_IS.read, rownames(RNA.tpm))
  
  if (anyNA(match_Davoli_IS.genes)){
    warning(c("differenty named or missing signature genes : \n", paste(Davoli_IS.read[!Davoli_IS.read %in% rownames(RNA.tpm)], collapse = "\n")))
    match_Davoli_IS.genes <- stats::na.omit(match_Davoli_IS.genes)
  }
  
  # Log2 transformation:
  log2.RNA.tpm <- RNA.tpm 
  
  # Subset log2.RNA.tpm
  sub_log2.RNA.tpm <- log2.RNA.tpm[match_Davoli_IS.genes, ]
  
  # Calculate rank position for each gene across samples
  ranks_sub_log2.RNA.tpm <- apply(sub_log2.RNA.tpm, 1, rank)
  
  # Get normalized rank by divided
  ranks_sub_log2.RNA.tpm.norm <- (ranks_sub_log2.RNA.tpm - 1)/(nrow(ranks_sub_log2.RNA.tpm) - 1)
  
  # Calculation: average of the expression value of all the genes within-sample
  score <- apply(ranks_sub_log2.RNA.tpm.norm, 1, mean)
  
  message("Davoli_IS score computed")
  return(data.frame(Davoli_IS = score, check.names = FALSE))
}

#' Compute IFNy signature score
#'
#' \code{compute_ayersIFNy} computes IFNy signature score as the arithmetic mean of genes included
#' in the IFN-γ signature (Ayers et al., JCI, 2017)
#'
#' @importFrom stats na.omit
#'
#' @param RNA.tpm numeric matrix with rows=genes and columns=samples
#'
#' @return numeric matrix with rows=samples and columns=IFNy signature score
#'
#' @export
#-------------------------------------------------------------------------------------------------------------

compute.IFNy <- function(RNA.tpm){
  
  # Literature genes
  IFNy.read <- c("IFNG", "STAT1", "CXCL9", "CXCL10", "IDO1", "HLA-DRA")
  match_IFNy.genes <- match(IFNy.read, rownames(RNA.tpm))
  
  if (anyNA(match_IFNy.genes)){
    warning(paste0("differenty named or missing signature genes : \n", IFNy.read[!IFNy.read %in% rownames(RNA.tpm)]))
    match_IFNy.genes <- stats::na.omit(match_IFNy.genes)
  }
  
  # Log2 transformation:
  log2.RNA.tpm <- RNA.tpm 
  
  # Subset log2.RNA.tpm
  sub_log2.RNA.tpm <- log2.RNA.tpm[match_IFNy.genes, ]
  
  # Calculation: average of the included genes for the IFN-y signature
  score <- apply(sub_log2.RNA.tpm, 2, mean)
  
  message("IFNy score computed")
  return(data.frame(IFNy = score, check.names = FALSE))
}


#' Compute MSI score
#'
#' \code{compute_MSI} computes MSI score by applying logical comparison of MSI-related gene pairs
#' (Fu et al., 2019).
#'
#' @importFrom stats na.omit
#'
#' @param RNA.tpm numeric matrix with rows=genes and columns=samples
#'
#' @return numeric matrix with rows=samples and columns=MSI score
#'
#' @export
#-------------------------------------------------------------------------------------------------------------

compute.MSI <- function(RNA.tpm){
  
  # Literature genes: * (CCRN4L in tcga, NOCT approved symbol)
  MSI.basis <- data.frame(Gene_1 = c("HNRNPL","MTA2","CALR","RASL11A","LYG1", "STRN3", "HPSE",
                                     "PRPF39","NOCT","AMFR"),
                          Gene_2 = c("CDC16","VGF","SEC22B","CAB39L","DHRS12", "TMEM192", "BCAS3",
                                     "ATF6","GRM8","DUSP18"))
  MSI.read <- unique(as.vector(as.matrix(MSI.basis))) # 20 genes
  
  # Some genes might have other name: case for "CCRN4L", it's called "NOCT", be carefull
  if (any(rownames(RNA.tpm) %in% "CCRN4L")){
    cat("Gene name changed: NOCT is approved symbol, not CCRN4L","\n")
    rownames(RNA.tpm)[rownames(RNA.tpm) %in% "CCRN4L"] <- "NOCT"
  }
  
  # Subset RNA.tpm
  match_F_1 <- match(as.character(MSI.basis[,1]), rownames(RNA.tpm))
  match_F_2 <- match(as.character(MSI.basis[,2]), rownames(RNA.tpm))
  
  if (anyNA(c(match_F_1,match_F_2))){
    warning(c("differenty named or missing signature genes : \n", paste(MSI.read[!MSI.read %in% rownames(RNA.tpm)], collapse = "\n")))
  }
  
  # Initialize variables
  F_pair_expr_A <- matrix(0, nrow(MSI.basis), ncol(RNA.tpm))
  F_pair_expr_B <- matrix(0, nrow(MSI.basis), ncol(RNA.tpm))
  MSI.matrix <- matrix(0, nrow(MSI.basis), ncol(RNA.tpm)) ; colnames(MSI.matrix) <- colnames(RNA.tpm)
  remove_pairs <- vector("list", length = ncol(RNA.tpm)) ; names(remove_pairs) <- colnames(RNA.tpm)
  score <- vector("numeric", length = ncol(RNA.tpm)) ; names(score) <- colnames(RNA.tpm)
  
  # Log2 transformation:
  log2.RNA.tpm <- as.data.frame(RNA.tpm)
  
  # Calculation:
  F_pair_expr_A <- log2.RNA.tpm[match_F_1, ]
  F_pair_expr_B <- log2.RNA.tpm[match_F_2, ]
  
  if(anyNA(F_pair_expr_A + F_pair_expr_B)) {
    remove_pairs <- as.vector(which(is.na(rowSums(F_pair_expr_A + F_pair_expr_B) == TRUE)))
  }
  
  MSI.matrix <- F_pair_expr_A > F_pair_expr_B
  if(anyNA(MSI.matrix)){
    score <- colSums(MSI.matrix, na.rm = TRUE)
    score <- (score * nrow(MSI.matrix)) / (nrow(MSI.matrix) - length(remove_pairs))
  }else{
    score <- colSums(MSI.matrix)
  }
  
  message("MSI score computed")
  return(data.frame(MSI = score, check.names = FALSE))
}

#' Compute Roh immune score
#'
#' \code{compute_rohIS} computes Roh immune score as the geometric-mean of immune score genes
#' (Roh et al., 2017).
#'
#' @importFrom stats na.omit
#'
#' @param RNA.tpm numeric matrix with rows=genes and columns=samples
#'
#' @return numeric matrix with rows=samples and columns=Roh immune score
#'
#' @export
#-------------------------------------------------------------------------------------------------------------

compute.Roh_IS <- function(RNA.tpm){
  
  # Literature genes
  Roh_IS.read <- c("GZMA", "GZMB", "PRF1", "GNLY", "HLA-A", "HLA-B", "HLA-C", "HLA-E", "HLA-F",
                   "HLA-G", "HLA-H", "HLA-DMA", "HLA-DMB", "HLA-DOA", "HLA-DOB", "HLA-DPA1",
                   "HLA-DPB1", "HLA-DQA1", "HLA-DQA2", "HLA-DQB1", "HLA-DRA", "HLA-DRB1",
                   "IFNG", "IFNGR1", "IFNGR2", "IRF1", "STAT1", "PSMB9", "CCR5", "CCL3", "CCL4",
                   "CCL5", "CXCL9", "CXCL10", "CXCL11", "ICAM1", "ICAM2", "ICAM3", "ICAM4", "ICAM5", "VCAM1")
  match_Roh_IS.genes <- match(Roh_IS.read, rownames(RNA.tpm))
  
  if (anyNA(match_Roh_IS.genes)){
    warning(c("differenty named or missing signature genes : \n", paste(Roh_IS.read[!Roh_IS.read %in% rownames(RNA.tpm)], collapse = "\n")))
    match_Roh_IS.genes <- stats::na.omit(match_Roh_IS.genes)
  }
  
  # Subset RNA.tpm
  sub_gene.tpm <- RNA.tpm[match_Roh_IS.genes, ]
  
  # Pseudocount of 0.01 for all genes
  sub_gene.tpm <- sub_gene.tpm + 0.01
  
  # Pseudocount of 1 for genes with 0 expr
  if(any(sub_gene.tpm == 0)) sub_gene.tpm[sub_gene.tpm == 0] <- sub_gene.tpm[sub_gene.tpm == 0] + 1
  
  # Calculation: geometric mean (so-called log-average) [TPM, 0.01 offset]
  score <- apply(sub_gene.tpm, 2, function(X) exp(mean(log(X))))
  
  message("Roh_IS computed score")
  return(data.frame(Roh_IS = score, check.names = FALSE))
}

#' Compute tertiary lymphoid structures signature
#'
#' \code{compute_TLS} computes TLS signature as the geometric-mean of TLS signature genes
#' (Cabrita et al., 2020).
#'
#' @importFrom stats na.omit
#'
#' @param RNA.tpm numeric matrix with rows=genes and columns=samples
#'
#' @return numeric matrix with rows=samples and columns=TLS signature
#'
#' @export
#-------------------------------------------------------------------------------------------------------------

compute.TLS <- function(RNA.tpm){
  
  # Literature genes
  TLS.read <- c("CD79B", "CD1D", "CCR6", "LAT", "SKAP1", "CETP", "EIF1AY", "RBP5", "PTGDS")
  match_TLS.read <- match(TLS.read, rownames(RNA.tpm))
  
  if (anyNA(match_TLS.read)){
    warning(c("differenty named or missing signature genes : \n", paste(TLS.read[!TLS.read %in% rownames(RNA.tpm)], collapse = "\n")))
    match_TLS.read <- stats::na.omit(match_TLS.read)
  }
  
  # Subset RNA.tpm
  sub_gene.tpm <- RNA.tpm[match_TLS.read, ]
  
  # Calculation: geometric mean (so-called log-average) [TPM, 1 offset]
  geom_mean <- apply(sub_gene.tpm, 2, function(X) exp(mean(X)))
  
  message("TLS score computed")
  return(data.frame(TLS = geom_mean, check.names = FALSE))
}

#' Compute T cell-inflamed signature score
#'
#' \code{compute_ayersTcellInfl} computes T cell-inflamed signature score by taking a weighted sum of
#'  the housekeeping normalized values of the T cell-inflamed signature genes
#'
#' @importFrom stats na.omit
#'
#' @param RNA.tpm numeric matrix with rows=genes and columns=samples
#'
#' @return numeric matrix with rows=samples and columns=T cell-inflamed signature score
#'
#' @export
#-------------------------------------------------------------------------------------------------------------

compute.Tcell_inflamed <- function(RNA.tpm){
  
  # Literature genes
  Tcell_inflamed.read <- c("CCL5", "CD27", "CD274", "CD276", "CD8A", "CMKLR1", "CXCL9", "CXCR6", "HLA-DQA1",
                           "HLA-DRB1", "HLA-E", "IDO1", "LAG3", "NKG7", "PDCD1LG2", "PSMB10", "STAT1", "TIGIT")
  Housekeeping.read <- c("STK11IP", "ZBTB34", "TBC1D10B", "OAZ1", "POLR2A", "G6PD", "ABCF1", "NRDE2", "UBB", "TBP", "SDHA") # C14orf102 = NRDE2
  weights <- data.frame(CCL5=0.008346, CD27=0.072293, CD274=0.042853, CD276=-0.0239, CD8A=0.031021 ,CMKLR1=0.151253, CXCL9=0.074135,
                        CXCR6=0.004313, `HLA-DQA1`=0.020091, `HLA-DRB1`=0.058806, `HLA-E`=0.07175, IDO1=0.060679, LAG3=0.123895, NKG7=0.075524, PDCD1LG2=0.003734,
                        PSMB10=0.032999, STAT1=0.250229, TIGIT=0.084767, check.names = FALSE)
  
  # Some genes might have other name: case for "C14orf102", it's called "NRDE2", be careful
  if (any(rownames(RNA.tpm) %in% "C14orf102")){
    cat("Gene name changed: NRDE2 is approved symbol, not C14orf102","\n")
    rownames(RNA.tpm)[rownames(RNA.tpm) %in% "C14orf102"] <- "NRDE2"
  }
  
  match_genes.housekeeping <- match(Housekeeping.read, rownames(RNA.tpm))
  match_genes.predictors <- match(Tcell_inflamed.read, rownames(RNA.tpm))
  
  if (anyNA(c(match_genes.housekeeping, match_genes.predictors))){
    tmp <- c(Tcell_inflamed.read, Housekeeping.read)
    warning(c("differenty named or missing signature genes : \n", paste(tmp[!tmp %in% rownames(RNA.tpm)], collapse = "\n")))
    match_genes.housekeeping <- stats::na.omit(match_genes.housekeeping)
    match_genes.predictors <- stats::na.omit(match_genes.predictors)
  }
  
  # Log2 transformation:
  log2.RNA.tpm <- RNA.tpm 
  
  # Subset log2.RNA.tpm
  ## housekeeping
  log2.RNA.tpm.housekeeping <- log2.RNA.tpm[match_genes.housekeeping, ]
  ## predictors
  log2.RNA.tpm.predictors <- log2.RNA.tpm[match_genes.predictors, ]
  weights <- weights[,rownames(log2.RNA.tpm.predictors)]
  
  # Housekeeping normalization
  average.log2.RNA.tpm.housekeeping <- apply(log2.RNA.tpm.housekeeping, 2, mean)
  log2.RNA.tpm.predictors.norm <- sweep(log2.RNA.tpm.predictors, 2, average.log2.RNA.tpm.housekeeping, FUN = "-")
  
  # Calculation: weighted sum of the normalized predictor gene values
  tidy <- match(rownames(log2.RNA.tpm.predictors.norm), colnames(weights))
  score <- t(log2.RNA.tpm.predictors.norm[tidy,]) %*% t(weights)
  
  message("Tcell_inflamed score computed")
  return(data.frame( Tcell_inflamed = score, check.names = FALSE))
}

#' Compute chemokine score
#'
#' \code{compute_chemokine} computes chemoine score as the PC1 score that results from applying PCA
#' to z-score expression of 12 chemokine genes (Messina et al., 2012).
#'
#' @importFrom stats na.omit prcomp
#'
#' @param RNA.tpm numeric matrix with rows=genes and columns=samples
#'
#' @return numeric matrix with rows=samples and columns=chemokine score
#'
#' @export
#-------------------------------------------------------------------------------------------------------------

compute.chemokines <- function(RNA.tpm){
  # Literature genes
  chemokines.read <- c("CCL2", "CCL3", "CCL4", "CCL5", "CCL8", "CCL18", "CCL19", "CCL21",
                       "CXCL9", "CXCL10", "CXCL11", "CXCL13")
  match_chemokines.genes <- match(chemokines.read, rownames(RNA.tpm))
  
  if (anyNA(match_chemokines.genes)){
    warning(c("differenty named or missing signature genes : \n", paste(chemokines.read[!chemokines.read %in% rownames(RNA.tpm)], collapse = "\n")))
    match_chemokines.genes <- stats::na.omit(match_chemokines.genes)
  }
  
  # Log2 transformation:
  log2.RNA.tpm <- RNA.tpm 
  
  # Subset gene_expr
  sub_log2.RNA.tpm <- log2.RNA.tpm[match_chemokines.genes, ]
  
  sub_log2.RNA.tpm = sub_log2.RNA.tpm[apply(sub_log2.RNA.tpm,1,sum)>0,]
  # calculation: using PCA (Z-score calculated within prcomp)
  chemokine.pca <- stats::prcomp(t(sub_log2.RNA.tpm), center = TRUE, scale = TRUE)
  score <- chemokine.pca$x[, 1]
  
  message("Chemokines score computed")
  return(data.frame(chemokines = score, check.names = FALSE))
}


#-------------------------------------------------------------------------------------------------------------


#' Compute the expression of the immune checkpoints genes
#'
#' \code{computation_ICB_genes} computes the scores for the immune checkpoint genes.
#'
#' @export
#'
#' @param RNA.tpm numeric matrix with data
#'
#' @return List with the expression of the immune checkpoint genes
#'
#-------------------------------------------------------------------------------------------------------

compute_ICB_genes <- function(RNA.tpm){
  
  # Extract position genes for GZMA and PRF1
  tmp <- match(c("CD274","CTLA4","PDCD1"), rownames(RNA.tpm))
  
  # PDL-1 calculation
  PDL1_expr = RNA.tpm[tmp[1],] ; rownames(PDL1_expr) <- "PDL1"
  
  # CTLA-4 calculation
  CTLA4_expr = RNA.tpm[tmp[2],] ; rownames(CTLA4_expr) <- "CTLA4"
  
  # PD-1 calculation
  PD1_expr = RNA.tpm[tmp[3],] ; rownames(PD1_expr) <- "PD1"
  
  ICB_genes_expr <- list(PDL1 = PDL1_expr , CTLA4 = CTLA4_expr , PD1 =PD1_expr)
  message("ICB genes expression computed")
  return(ICB_genes_expr)
}



#' \code{computation_gold_standards} computes the scores for the gold standards required by the user
#'
#' @export
#'
#' @param RNA.tpm numeric matrix with data
#' @param list_gold_standards string with gold standards names
#' @param cancertype string character
#'
#' @return List with the scores of all the gold standards specified.
#'
#-------------------------------------------------------------------------------------------------------
# Input: Transcriptomics data as tpm values
# A list with the names of the scores to be computed has to be provided
# Output: Gold standards scores
#-------------------------------------------------------------------------------------------------------
compute_gold_standards <- function(RNA.tpm){
  list_gold_standards = list_gold_standards=c("CYT","Roh_IS","Davoli_IS", "IFNy","Ayers_expIS","Tcell_inflamed")
  RNA.tpm = as.data.frame(RNA.tpm)
  gold.standards <- sapply(list_gold_standards, function(X){
    
    if ("CYT" == X) {
      
      # calculate Cytolytic activity #
      CYT <- t(compute.CYT(RNA.tpm))
      return(list(CYT))
      
    }else if("Roh_IS" == X) {
      
      # calculate roh immune signature #
      Roh_IS <- t(compute.Roh_IS(RNA.tpm))
      return(list(Roh_IS))
      
    }else if("chemokines" == X) {
      
      # calculate chemokine signature #
      chemokines <- t(compute.chemokines(RNA.tpm))
      return(list(chemokines))
      
    }else if("Davoli_IS" == X) {
      
      # calculate davoli cytotoxic immune signature #
      Davoli_IS <- t(compute.Davoli_IS(RNA.tpm))
      return(list(Davoli_IS))
      
    }else if("IFNy" == X) {
      
      # calculate ayers IFNy #
      IFNy <- t(compute.IFNy(RNA.tpm))
      return(list(IFNy))
      
    }else if("Ayers_expIS" == X) {
      
      # calculate ayers expanded immune signature #
      Ayers_expIS <- t(compute.Ayers_expIS(RNA.tpm))
      return(list(Ayers_expIS))
      
    }else if("Tcell_inflamed" == X) {
      
      # calculate ayers T cell inflamed signature #
      Tcell_inflamed <- t(compute.Tcell_inflamed(RNA.tpm))
      return(list(Tcell_inflamed))
      
    }else if("MSI" == X) {
      
      # calculate MSI signature #
      MSI <- t(compute.MSI(RNA.tpm))
      return(list(MSI))
      
    }else if("RIR" == X) {
      
      # calculate MSI signature #
      RIR <- t(compute.RIR(RNA.tpm))
      return(list(RIR))
      
    }else if("TLS" == X) {
      
      # calculate MSI signature #
      TLS <- t(compute.TLS(RNA.tpm))
      return(list(TLS))
      
    }
    
  })
  
  gold.standards = lapply(gold.standards,function(x){t(x) %>% as.data.frame() %>% rownames_to_column("condition")})
  gold.standards = Reduce(function(x,y) {left_join(x,y ,by="condition")}, gold.standards)
  rownames(gold.standards) = gold.standards$condition 
  gold.standards$condition = NULL
  
  gold.standards = data.frame(minMax(gold.standards))
  
  return(gold.standards)
  
}


