# Functions for processing deconvolution features for LungPredict1 paper

# author: Marcelo Hurtado
# email: marcelo.hurtado@inserm.fr
# organization: INSERM CRCT - Pancaldi team 21
# place: Toulouse, France

#' Processing of deconvolution features
#' 
#' \code{compute.deconvolution.analysis} Deconvolution analysis to reduce data dimensionality by applying iterative correlations techniques  
#' 
#' @param deconvolution.mat Deconvolution matrix from GEMDeCan output (samples X features).
#' @return list of Deconvolution matrix with subgroupped features (samples X features), groups formed by Pearson corr, groups formed by proportionality and subgroups by cell type
#' 
#' -------------------------------------------------------------------------------------------------------------
#' 
compute.deconvolution.analysis <- function(deconvolution, corr, zero = 0.9, high_corr = 0.9, seed = NULL){
  deconvolution.mat = deconvolution
  
  #####Unsupervised filtering 
  
  #Remove high zero number features
  cat(paste0("Removing features with high zero number ", round(zero*100,2), "%...............................................................\n\n"))
  deconvolution.mat = deconvolution.mat[, colSums(deconvolution.mat == 0, na.rm=TRUE) < round(zero*nrow(deconvolution.mat)) , drop=FALSE]
  diff_colnames <- setdiff(colnames(deconvolution), colnames(deconvolution.mat))
  zero_features <- deconvolution[, diff_colnames]
  
  #Remove low_variance features
  variance = remove_low_variance(deconvolution.mat)
  deconvolution.mat = variance[[1]]
  low_variance_features = variance[[2]]
  
  #Scale deconvolution features by columns for making them comparable between cell types (0-1). 
  cat("Scaling deconvolution features for comparison between cell types...............................................................\n\n")
  for (i in 1:ncol(deconvolution.mat)) {
    deconvolution.mat[,i] = deconvolution.mat[,i]/max(deconvolution.mat[,i])
  } 
  
  #####Cell types split 
  cat("Splitting deconvolution features per cell type...............................................................\n\n")
  cells_types = compute.cell.types(deconvolution.mat)
  cells = cells_types[[1]]
  cells_discarded = cells_types[[2]]
  
  ######Pairwise correlation filtering (Highly correlated variables >0.9) within cell types
  cat("Finding group of features with high correlation between each other...............................................................\n\n")
  features_high_corr = list()
  j = 1
  for (i in 1:length(cells)) {
    data = cells[[i]]
    if(is.null(ncol(data))==T){
      cells[[i]] = data
    }else if(ncol(data)>1){
      data = removeCorrelatedFeatures(data, high_corr, names(cells)[i], seed)
      cells[[i]] = data[[1]]
      if(length(data[[2]])>0 && is.null(data[[3]])==F){
        features_high_corr[[j]] = data[[2]]
        names(features_high_corr)[j] = data[[3]]
        j = j+1
      }
    }
  }
  
  #####Subgrouping of deconvolution features
  res = list()
  groups = list()
  groups_similarity = list()
  groups_discard = list()
  for (i in 1:length(cells)) {
    x = compute_subgroups(cells[[i]], file_name = names(cells)[i], thres_corr = corr)
    res = c(res, x[1])
    groups = c(groups, x[2])
    groups_similarity = c(groups_similarity, x[3])
    groups_discard = c(groups_discard, x[4])
  }
  names_cells = c("B.cells", "B.naive", "B.memory", "Macrophages.cells", "Macrophages.M0", "Macrophages.M1", "Macrophages.M2", "Monocytes", "Neutrophils", "NK.cells", "NK.activated", "NK.resting", "NKT.cells", "CD4.cells", "CD4.memory.activated",
                  "CD4.memory.resting", "CD4.naive", "CD8.cells", "T.cells.regulatory", "T.cells.non.regulatory","T.cells.helper", "T.cells.gamma.delta", "Dendritic.cells", "Dendritic.activated", "Dendritic.resting", "Cancer", "Endothelial",
                  "Eosinophils", "Plasma.cells", "Myocytes", "Fibroblasts", "Mast.cells", "Mast.activated", "Mast.resting", "CAF")
  
  names(res) = names_cells
  names(groups) = names_cells
  names(groups_similarity) = names_cells
  names(groups_discard) = names_cells
  
  #####Preparing output
  dt = c()
  for (i in 1:length(res)) {
    dt = c(dt, res[[i]])
  }
  dt = data.frame(dt)
  rownames(dt) = rownames(deconvolution.mat)
  
  #####Create and export table with subgroups
  
  #Count number of subgroups - Linear-based
  idx = c()
  for (i in 1:length(groups)){
    if(length(groups[[i]])>0){
      for (j in 1:length(groups[[i]])){
        idx = c(idx, names(groups[[i]])[[j]])
      } 
    }
  }  
  data.groups = data.frame(matrix(nrow = length(idx), ncol = 2)) #Create table
  colnames(data.groups) = c("Cell_subgroups", "Methods-signatures")
  data.groups$Cell_subgroups = idx #Assign subgroups 
  
  #Save methods corresponding to each subgroup
  contador = 1
  for (i in 1:length(groups)){
    if(length(groups[[i]])>0){
      for (j in 1:length(groups[[i]])){
        data.groups[contador,2] = paste(groups[[i]][[j]], collapse ="\n")
        contador = contador + 1
      } 
    }
  }
  
  
  #Count number of subgroups - Proportionality-based
  idy = c()
  for (i in 1:length(groups_similarity)){
    if(length(groups_similarity[[i]])>0){
      for (j in 1:length(groups_similarity[[i]])){
        idy = c(idy, names(groups_similarity[[i]])[[j]])
      } 
    }
  }
  data.groups.similarity = data.frame(matrix(nrow = length(idy), ncol = 2)) #Create table
  colnames(data.groups.similarity) = c("Cell_subgroups", "Methods-signatures")
  data.groups.similarity$Cell_subgroups = idy #Assign subgroups 
  
  #Save methods corresponding to each subgroup
  contador = 1
  for (i in 1:length(groups_similarity)){
    if(length(groups_similarity[[i]])>0){
      for (j in 1:length(groups_similarity[[i]])){
        data.groups.similarity[contador,2] = paste(groups_similarity[[i]][[j]], collapse ="\n")
        contador = contador + 1
      } 
    }
  }
  
  #Save data to export
  data.output = rbind(data.groups.similarity, data.groups)
  write.csv(dt, 'Results/Deconvolution_after_subgrouping.csv')
  write.csv(data.output, 'Results/Cell_subgroups.csv', row.names = F)
  
  message("Deconvolution features subgroupped")
  
  results = list(dt, res, groups, groups_similarity, groups_discard, zero_features, low_variance_features, cells_discarded, features_high_corr)
  names(results) = c("Deconvolution matrix", "Deconvolution groups per cell types", "Deconvolution groups - Linear-based correlation", "Deconvolution groups - Proportionality-based correlation",
                     "Discarded groups with equal method", "Discarded features with high number of zeros", "Discarded features with low variance", "Discarded cell types",
                     "High correlated deconvolution groups (>0.9) per cell type")
  return(results)
  
}

compute.deconvolution.preprocessing = function(deconv){
  cat("Preprocessing deconvolution features...............................................................\n\n")

  #Convert mcp and xcell features to proportions by row-scaling 
  for (i in 1:nrow(deconv)) {
    deconv[,grep("MCP", colnames(deconv))][i,] = deconv[,grep("MCP", colnames(deconv))][i,]/sum(deconv[,grep("MCP", colnames(deconv))][i,])
    deconv[,grep("XCell", colnames(deconv))][i,] = deconv[,grep("XCell", colnames(deconv))][i,]/sum(deconv[,grep("XCell", colnames(deconv))][i,])
  } 
  
  ##### Edit cell names for consistency across features
  ##### Macrophages (M0, M1, M2)
  Macrophages = deconv[,grep("acrophage", colnames(deconv))]
  M0 = deconv[,grep("M0", colnames(deconv))]
  M1 = deconv[,grep("M1", colnames(deconv))]
  M2 <- deconv[,grep("M2", colnames(deconv))]
  if(length(grep("LM22", colnames(M2)))>0){M2 <- M2[,-grep("LM22", colnames(M2))]}else{M2 <- M2} 
  test = deconv[,grep("LM22", colnames(deconv))]
  test = test[,grep("Macrophages.M2", colnames(test))]
  M2 = cbind(M2, test)
  Macrophages = Macrophages[,-which(colnames(Macrophages)%in%c(colnames(M0), colnames(M1), colnames(M2)))]
  deconv = deconv[,-which(colnames(deconv)%in%c(colnames(Macrophages), colnames(M0), colnames(M1), colnames(M2)))]
  
  colnames(Macrophages) = stringr::str_replace(colnames(Macrophages), "Macrophages", "Macrophages.cells") 
  colnames(Macrophages) = stringr::str_replace(colnames(Macrophages), "Macrophage(?!.)", "Macrophages.cells") 
  colnames(M0) = stringr::str_replace(colnames(M0), "Macrophages_M0", "Macrophages.M0") 
  colnames(M0) = stringr::str_replace(colnames(M0), "_M0", "_Macrophages.M0")
  colnames(M1) = stringr::str_replace(colnames(M1), "Macrophages_M1", "Macrophages.M1") 
  colnames(M1) = stringr::str_replace(colnames(M1), "Macrophages_M1", "Macrophages.M1")
  colnames(M1) = stringr::str_replace(colnames(M1), "Macrophage_M1", "Macrophages.M1")
  colnames(M1) = stringr::str_replace(colnames(M1), "_M1", "_Macrophages.M1")
  colnames(M2) = stringr::str_replace(colnames(M2), "Macrophage_M2", "Macrophages.M2") 
  colnames(M2) = stringr::str_replace(colnames(M2), "Macrophages_M2", "Macrophages.M2")
  colnames(M2) = stringr::str_replace(colnames(M2), "_M2", "_Macrophages.M2")
  
  ##### Monocytes
  Monocytes = deconv[,grep("Mono|mono", colnames(deconv))]
  deconv = deconv[,-which(colnames(deconv)%in%colnames(Monocytes))]
  colnames(Monocytes) = stringr::str_replace(colnames(Monocytes), "Monocytic_lineage", "Monocytes") 
  colnames(Monocytes) = stringr::str_replace(colnames(Monocytes), "Monocyte(?!s)", "Monocytes") 
  colnames(Monocytes) = stringr::str_replace(colnames(Monocytes), "Mono(?!cytes)", "Monocytes")
  colnames(Monocytes) = stringr::str_replace(colnames(Monocytes), "Mono(?!cytes)", "Monocytes")

  ##### Neutrophils
  Neutrophils <- deconv[,grep("Neu", colnames(deconv))]
  deconv = deconv[,-which(colnames(deconv)%in%colnames(Neutrophils))]
  
  colnames(Neutrophils) = stringr::str_replace(colnames(Neutrophils), "Neutrophil(?!s)", "Neutrophils") 
  colnames(Neutrophils) = stringr::str_replace(colnames(Neutrophils), "Neu(?!trophils)", "Neutrophils") 
  
  ###NK cells
  NK = deconv[,grep("NK", colnames(deconv))]
  NKT = NK[,grep("NKT", colnames(NK))]
  NK.activated <- NK[,grep("activated", colnames(NK), value = TRUE)]
  NK.resting <- NK[,grep("resting", colnames(NK), value = TRUE)]
  NK = NK[,-which(colnames(NK)%in%c(colnames(NK.activated), colnames(NK.resting), colnames(NKT)))]  
  deconv = deconv[,-which(colnames(deconv)%in%c(colnames(NK), colnames(NK.activated), colnames(NK.resting), colnames(NKT)))]
  
  colnames(NK) = stringr::str_replace(colnames(NK), "NK(?!.)", "NK.cells")
  colnames(NK) = stringr::str_replace(colnames(NK), "NK_cells", "NK.cells")
  colnames(NK) = stringr::str_replace(colnames(NK), "NK_cell", "NK.cells")
  colnames(NKT) = stringr::str_replace(colnames(NKT), "NKT_", "NKT.")
  colnames(NK.activated) = stringr::str_replace(colnames(NK.activated), "NK.cells.activated", "NK.activated")
  colnames(NK.activated) = stringr::str_replace(colnames(NK.activated), "NK.cells_activated", "NK.activated")
  colnames(NK.resting) = stringr::str_replace(colnames(NK.resting), "NK.cells.resting", "NK.resting")
  colnames(NK.resting) = stringr::str_replace(colnames(NK.resting), "NK.cells_resting", "NK.resting")
  
  ###CD4 cells
  CD4 <- deconv[,grep("CD4", colnames(deconv))]
  CD4.memory.activated = CD4[,grep("activated", colnames(CD4))]
  CD4.memory.resting = CD4[,grep("resting", colnames(CD4))]
  CD4.naive = CD4[,grep("naive", colnames(CD4))]
  CD4.non.regulatory = CD4[,grep("regulatory", colnames(CD4))]
  CD4 = CD4[,-which(colnames(CD4)%in%c(colnames(CD4.memory.activated), colnames(CD4.memory.resting), colnames(CD4.naive), colnames(CD4.non.regulatory)))]  
  deconv = deconv[,-which(colnames(deconv)%in%c(colnames(CD4), colnames(CD4.memory.activated), colnames(CD4.memory.resting), colnames(CD4.naive), colnames(CD4.non.regulatory)))]
  
  colnames(CD4) = stringr::str_replace(colnames(CD4), "CD4", "CD4.cells")
  colnames(CD4) = stringr::str_replace(colnames(CD4), "T.cells.CD4.cells", "CD4.cells")
  colnames(CD4.memory.activated) = stringr::str_replace(colnames(CD4.memory.activated), "CD4_memory_activated", "CD4.memory.activated")
  colnames(CD4.memory.resting) = stringr::str_replace(colnames(CD4.memory.resting), "CD4_memory_resting", "CD4.memory.resting")
  colnames(CD4.naive) = stringr::str_replace(colnames(CD4.naive), "CD4_naive", "CD4.naive")
  colnames(CD4.naive) = stringr::str_replace(colnames(CD4.naive), "CD4._naive", "CD4.naive")
  colnames(CD4.naive) = stringr::str_replace(colnames(CD4.naive), "T.cells.CD4.naive", "CD4.naive")
  colnames(CD4.naive) = stringr::str_replace(colnames(CD4.naive), "T_cells_CD4.naive", "CD4.naive")
  colnames(CD4.non.regulatory) = stringr::str_replace(colnames(CD4.non.regulatory), "T_cell_CD4._.non.regulatory.", "T.cells.non.regulatory")

  ####CD8
  CD8 <- deconv[,grep("CD8", colnames(deconv))]

  deconv = deconv[,-which(colnames(deconv)%in%colnames(CD8))]
  
  colnames(CD8) = stringr::str_replace(colnames(CD8), "T_cells_CD8", "CD8.cells")
  colnames(CD8) = stringr::str_replace(colnames(CD8), "T_cell_CD8", "CD8.cells")
  colnames(CD8) = stringr::str_replace(colnames(CD8), "CD8_T_cells", "CD8.cells")
  colnames(CD8) = stringr::str_replace(colnames(CD8), "T.cells.CD8", "CD8.cells")
  colnames(CD8) = stringr::str_replace(colnames(CD8), "CD8(?!.)", "CD8.cells")
  colnames(CD8) = stringr::str_replace(colnames(CD8), "CD8.cells.", "CD8.cells")

  ##### Regulatory T cells 
  Tregs = deconv[,grep("regs", colnames(deconv))]
  deconv = deconv[,-which(colnames(deconv)%in%colnames(Tregs))]
  
  colnames(Tregs) = stringr::str_replace(colnames(Tregs), "T_cell_regulatory_.Tregs.", "T.cells.regulatory")
  colnames(Tregs) = stringr::str_replace(colnames(Tregs), "T.cells.regulatory..Tregs.", "T.cells.regulatory")
  colnames(Tregs) = stringr::str_replace(colnames(Tregs), "Tregs", "T.cells.regulatory")
  
  ##### Helper T cells 
  Thelper = deconv[,grep("helper", colnames(deconv))]
  deconv = deconv[,-which(colnames(deconv)%in%colnames(Thelper))]
  
  colnames(Thelper) = stringr::str_replace(colnames(Thelper), "T.cells.follicular.helper", "T.cells.helper")
  colnames(Thelper) = stringr::str_replace(colnames(Thelper), "T_cells_follicular_helper", "T.cells.helper")
  
  ##### Gamma delta T cells 
  Tgamma = deconv[,grep("gamma", colnames(deconv))]
  deconv = deconv[,-which(colnames(deconv)%in%colnames(Tgamma))]
  
  colnames(Tgamma) = stringr::str_replace(colnames(Tgamma), "T_cells_gamma_delta", "T.cells.gamma.delta")
  colnames(Tgamma) = stringr::str_replace(colnames(Tgamma), "T_cell_gamma_delta", "T.cells.gamma.delta")
  
  ##### Dendritic cells (activated, resting)
  Dendritic = deconv[,grep("endritic", colnames(deconv))]
  Dendritic.activated = Dendritic[,grep("activated", colnames(Dendritic))]
  Dendritic.resting = Dendritic[,grep("resting", colnames(Dendritic))]
  Dendritic = Dendritic[,-which(colnames(Dendritic)%in%c(colnames(Dendritic.activated), colnames(Dendritic.resting)))]
  deconv = deconv[,-which(colnames(deconv)%in%c(colnames(Dendritic), colnames(Dendritic.activated), colnames(Dendritic.resting)))]
  
  colnames(Dendritic) = stringr::str_replace(colnames(Dendritic), "Myeloid_dendritic_cells", "Dendritic.cells")
  colnames(Dendritic) = stringr::str_replace(colnames(Dendritic), "Myeloid_dendritic_cell", "Dendritic.cells")
  colnames(Dendritic) = stringr::str_replace(colnames(Dendritic), "Dendritic_cells", "Dendritic.cells")
  
  colnames(Dendritic.activated) = stringr::str_replace(colnames(Dendritic.activated), "dendritic_cell_activated", "Dendritic.activated.cells")
  colnames(Dendritic.activated) = stringr::str_replace(colnames(Dendritic.activated), "Dendritic.cells.activated", "Dendritic.activated.cells")
  colnames(Dendritic.activated) = stringr::str_replace(colnames(Dendritic.activated), "Dendritic_cells_activated", "Dendritic.activated.cells")
  colnames(Dendritic.resting) = stringr::str_replace(colnames(Dendritic.resting), "Dendritic.cells.resting", "Dendritic.resting.cells")
  colnames(Dendritic.resting) = stringr::str_replace(colnames(Dendritic.resting), "Dendritic_cells_resting", "Dendritic.resting.cells")
  
  ##### CAF cells 
  CAF = deconv[,grep("CAF|Cancer_associated_fibroblast", colnames(deconv))]
  deconv = deconv[,-which(colnames(deconv)%in%colnames(CAF))]
  
  colnames(CAF) = stringr::str_replace(colnames(CAF), "Cancer_associated_fibroblast", "CAF")
  colnames(CAF) = stringr::str_replace(colnames(CAF), "CAFs", "CAF")
  
  ##### Cancer cells
  Cancer = deconv[,grep("ancer", colnames(deconv))]
  deconv = deconv[,-which(colnames(deconv)%in%colnames(Cancer))]
  
  colnames(Cancer) = stringr::str_replace(colnames(Cancer), "cancer", "Cancer")
  colnames(Cancer) = stringr::str_replace(colnames(Cancer), "Cancer.cells", "Cancer")
  
  malignant = deconv[,grep("alignant", colnames(deconv))]
  deconv = deconv[,-which(colnames(deconv)%in%colnames(malignant))]
  
  colnames(malignant) = stringr::str_replace(colnames(malignant), "Malignant", "Cancer")
  colnames(malignant) = stringr::str_replace(colnames(malignant), "Cancer_cells", "Cancer")
  colnames(malignant) = stringr::str_replace(colnames(malignant), "Cancer.cells", "Cancer")
  Cancer = cbind(Cancer, malignant)
  
  ##### Endothelial cells
  Endothelial = deconv[,grep("dothelial", colnames(deconv))]
  deconv = deconv[,-which(colnames(deconv)%in%colnames(Endothelial))]
  
  colnames(Endothelial) = stringr::str_replace(colnames(Endothelial), "Endothelial_cells", "Endothelial")
  colnames(Endothelial) = stringr::str_replace(colnames(Endothelial), "Endothelial.cells", "Endothelial")
  colnames(Endothelial) = stringr::str_replace(colnames(Endothelial), "Endothelial_cell", "Endothelial")
  
  ##### Eosinophils cells
  Eosinophils = deconv[,grep("osinophil", colnames(deconv))]
  deconv = deconv[,-which(colnames(deconv)%in%colnames(Eosinophils))]
  colnames(Eosinophils) = stringr::str_replace(colnames(Eosinophils), "Eosinophil(?!.)", "Eosinophils")
  
  ##### Plasma cells
  Plasma = deconv[,grep("lasma", colnames(deconv))]
  deconv = deconv[,-which(colnames(deconv)%in%colnames(Plasma))]
  
  colnames(Plasma) = stringr::str_replace(colnames(Plasma), "plasma(?!.)", "Plasma.cells")
  colnames(Plasma) = stringr::str_replace(colnames(Plasma), "Plasma_cells", "Plasma.cells")

  ##### Myocytes cells
  Myocytes = deconv[,grep("yocytes", colnames(deconv))]
  deconv = deconv[,-which(colnames(deconv)%in%colnames(Myocytes))]
  #colnames(Myocytes) = stringr::str_replace(colnames(Myocytes), "Myocytes", "Myocytes")
  
  ##### Fibroblasts cells :1 column
  Fibroblasts = data.frame(deconv[,grep("ibroblast", colnames(deconv))])
  colnames(Fibroblasts) = colnames(deconv)[grep("ibroblast", colnames(deconv))]
  deconv = deconv[,-which(colnames(deconv)%in%colnames(Fibroblasts))]
  #colnames(Fibroblasts) = stringr::str_replace(colnames(Fibroblasts), "Fibroblasts", "Fibroblasts")
  
  ##### Mast cells
  Mast = deconv[,grep("Mast", colnames(deconv))]
  Mast.activated = Mast[,grep("activated", colnames(Mast))]
  Mast.resting = Mast[,grep("resting", colnames(Mast))]
  Mast = Mast[,-which(colnames(Mast)%in%c(colnames(Mast.activated), colnames(Mast.resting)))]  
  deconv = deconv[,-which(colnames(deconv)%in%c(colnames(Mast), colnames(Mast.activated), colnames(Mast.resting)))]
  
  colnames(Mast) = stringr::str_replace(colnames(Mast), "Mast_cell(?!.)", "Mast.cells")
  colnames(Mast) = stringr::str_replace(colnames(Mast), "Mast_cells", "Mast.cells")
  colnames(Mast.activated) = stringr::str_replace(colnames(Mast.activated), "Mast.cells.activated", "Mast.activated.cells")
  colnames(Mast.activated) = stringr::str_replace(colnames(Mast.activated), "Mast_cells_activated", "Mast.activated.cells")
  colnames(Mast.resting) = stringr::str_replace(colnames(Mast.resting), "Mast.cells.resting", "Mast.resting.cells")
  colnames(Mast.resting) = stringr::str_replace(colnames(Mast.resting), "Mast_cells_resting", "Mast.resting.cells")
  
  ##### B cells
  B.naive = deconv[,grep("naive", colnames(deconv))]
  B.memory = deconv[,grep("memory", colnames(deconv))]
  B = deconv[,-which(colnames(deconv)%in%c(colnames(B.naive), colnames(B.memory)))] #"B" will include also the cells haven't been re-named, as it is the last cell
  colnames(B) = stringr::str_replace(colnames(B), "B_cells", "B.cells")
  colnames(B) = stringr::str_replace(colnames(B), "B_cell", "B.cells")
  colnames(B) = stringr::str_replace(colnames(B), "B_lineage", "B.cells")
  colnames(B) = stringr::str_replace(colnames(B), "_B(?!.)", "_B.cells")
  colnames(B.naive) = stringr::str_replace(colnames(B.naive), "B.cells.naive", "B.naive.cells")
  colnames(B.naive) = stringr::str_replace(colnames(B.naive), "B_cells_naive", "B.naive.cells")
  colnames(B.naive) = stringr::str_replace(colnames(B.naive), "B_cell_naive", "B.naive.cells")
  colnames(B.memory) = stringr::str_replace(colnames(B.memory), "B.cells.memory", "B.memory.cells")
  colnames(B.memory) = stringr::str_replace(colnames(B.memory), "B_cells_memory", "B.memory.cells")
  colnames(B.memory) = stringr::str_replace(colnames(B.memory), "B_cell_memory", "B.memory.cells")
  
  
  cell_types = cbind(B, B.naive, B.memory, Macrophages, M0, M1, M2, Monocytes, Neutrophils, NK, NK.activated, NK.resting, NKT, CD4, CD4.memory.activated,
                     CD4.memory.resting, CD4.naive, CD4.non.regulatory, CD8, Tregs, Thelper, Tgamma, Dendritic, Dendritic.activated, Dendritic.resting, Cancer, 
                     Endothelial, Eosinophils, Plasma, Myocytes, Fibroblasts, Mast, Mast.activated, Mast.resting, CAF)
  
  cat("Checking consistency in deconvolution cell fractions across patients...............................................................\n\n")
  
  combinations = c("Quantiseq", "MCP", "XCell", "Epidish_BPRNACan_",  "Epidish_BPRNACanProMet", "Epidish_BPRNACan3DProMet", "Epidish_CBSX.HNSCC.scRNAseq", "Epidish_CBSX.Melanoma.scRNAseq",
                   "Epidish_CBSX.NSCLC.PBMCs.scRNAseq", "Epidish_CCLE.TIL10", "Epidish_TIL10", "Epidish_LM22", "DeconRNASeq_BPRNACan_", "DeconRNASeq_BPRNACanProMet",
                   "DeconRNASeq_BPRNACan3DProMet", "DeconRNASeq_CBSX.HNSCC.scRNAseq", "DeconRNASeq_CBSX.Melanoma.scRNAseq", "DeconRNASeq_CBSX.NSCLC.PBMCs.scRNAseq", "DeconRNASeq_CCLE.TIL10",
                   "DeconRNASeq_TIL10", "DeconRNASeq_LM22", "CBSX_BPRNACan_", "CBSX_BPRNACanProMet", "CBSX_BPRNACan3DProMet", "CBSX_CBSX.HNSCC.scRNAseq",
                   "CBSX_CBSX.Melanoma.scRNAseq", "CBSX_CBSX.NSCLC.PBMCs.scRNAseq", "CBSX_CCLE.TIL10", "CBSX_TIL10", "CBSX_LM22")

  for(i in combinations){
    mat = cell_types[,grep(i, colnames(cell_types))]
    print(paste("Total sum across samples of combination", i, "is", round(sum(mat[1, ]), 2)))
  }
  
  return(cell_types)
}

#' Splitting deconvolution features by cell types
#' 
#' \code{compute.cell.types} Split deconvolution features into matrices for each cell type (features X samples)
#' 
#' @param data Deconvolution matrix from GEMDeCan output (samples X features).
#' @return list of cell types matrices (samples X features).
#' 
#' @details Cell types included: B, Macrophages, M0, M1, M2, Monocytes, Neutrophils, NK, NK.activated, NK.resting, NKT, CD4, CD4.activated, CD4.resting, CD8, Tregs, Dendritic, 
#' Dendritic.activated, Dendritic.resting, Cancer, Endothelial, CAF
#' 
#' -------------------------------------------------------------------------------------------------------------
#' 
compute.cell.types = function(data){
  ##### B cells
  B = grep("B.cells", colnames(data))
  B = data[, B, drop = FALSE]
  ##### B naive
  B.naive = grep("B.naive", colnames(data))
  B.naive = data[, B.naive, drop = FALSE]
  ##### B memory
  B.memory = grep("B.memory", colnames(data))
  B.memory = data[, B.memory, drop = FALSE]  
  ##### Macrophages (M0, M1, M2)
  Macrophages = grep("Macrophages.cells", colnames(data))
  Macrophages = data[, Macrophages, drop = FALSE]  
  M0 = grep("Macrophages.M0", colnames(data))
  M0 = data[, M0, drop = FALSE]  
  M1 = grep("Macrophages.M1", colnames(data))
  M1 = data[, M1, drop = FALSE]  
  M2 = grep("Macrophages.M2", colnames(data))
  M2 = data[, M2, drop = FALSE]  
  ##### Monocytes
  Monocytes = grep("Monocytes", colnames(data))
  Monocytes = data[, Monocytes, drop = FALSE]   
  ##### Neutrophils
  Neutrophils = grep("Neutrophils", colnames(data))
  Neutrophils = data[, Neutrophils, drop = FALSE]     
  ##### NK cells (activated, resting)
  NK = grep("NK.cells", colnames(data))
  NK = data[, NK, drop = FALSE]   
  NK.activated = grep("NK.activated", colnames(data))
  NK.activated = data[, NK.activated, drop = FALSE]   
  NK.resting = grep("NK.resting", colnames(data))
  NK.resting = data[, NK.resting, drop = FALSE]   
  ##### NKT cells
  NKT = grep("NKT.cells", colnames(data))
  NKT = data[, NKT, drop = FALSE]   
  ##### CD4 cells (activated, resting)
  CD4 = grep("CD4.cells", colnames(data))
  CD4 = data[, CD4, drop = FALSE]   
  #memory = CD4[,grep("memory", colnames(CD4))]
  #helper = CD4[,grep("Th", colnames(CD4))]
  #CD4 = CD4[,-which(colnames(CD4)%in%c(colnames(memory), colnames(helper)))]
  CD4.memory.activated = grep("CD4.memory.activated", colnames(data))
  CD4.memory.activated = data[, CD4.memory.activated, drop = FALSE]   
  CD4.memory.resting = grep("CD4.memory.resting", colnames(data))
  CD4.memory.resting = data[, CD4.memory.resting, drop = FALSE]   
  CD4.naive = grep("CD4.naive", colnames(data))
  CD4.naive = data[, CD4.naive, drop = FALSE]
  ##### CD8 cells
  CD8 = grep("CD8.cells", colnames(data))
  CD8 = data[, CD8, drop = FALSE]
  #naive = grep("naive", colnames(CD8))
  #naive = CD8[, naive, drop = FALSE]
  #memory = grep("memory", colnames(CD8))
  #memory = CD8[, memory, drop = FALSE]
  #CD8 = CD8[,-which(colnames(CD8)%in%c(colnames(memory), colnames(naive)))]
  ##### Regulatory T cells 
  Tregs = grep("T.cells.regulatory", colnames(data))
  Tregs = data[, Tregs, drop = FALSE]
  ##### Non regulatory T cells 
  T.non.regs = grep("T.cells.non.regulatory", colnames(data))
  T.non.regs = data[, T.non.regs, drop = FALSE]
  ##### Helper T cells 
  Thelper = grep("T.cells.helper", colnames(data))
  Thelper = data[, Thelper, drop = FALSE]
  ##### Gamma delta T cells 
  Tgamma = grep("T.cells.gamma.delta", colnames(data))
  Tgamma = data[, Tgamma, drop = FALSE]
  ##### Dendritic cells (activated, resting)
  Dendritic = grep("Dendritic.cells", colnames(data))
  Dendritic = data[, Dendritic, drop = FALSE]
  Dendritic.activated = grep("Dendritic.activated", colnames(data))
  Dendritic.activated = data[, Dendritic.activated, drop = FALSE]
  Dendritic.resting = grep("Dendritic.resting", colnames(data))
  Dendritic.resting = data[, Dendritic.resting, drop = FALSE]
  ##### Cancer cells
  Cancer = grep("Cancer", colnames(data))
  Cancer = data[, Cancer, drop = FALSE]
  ##### Endothelial cells
  Endothelial = grep("Endothelial", colnames(data))
  Endothelial = data[, Endothelial, drop = FALSE]
  ##### Eosinophils cells
  Eosinophils = grep("Eosinophils", colnames(data))
  Eosinophils = data[, Eosinophils, drop = FALSE]
  ##### Plasma cells
  Plasma = grep("Plasma.cells", colnames(data))
  Plasma = data[, Plasma, drop = FALSE]
  ##### Myocytes cells
  Myocytes = grep("Myocytes", colnames(data))
  Myocytes = data[, Myocytes, drop = FALSE]
  ##### Fibroblasts cells
  Fibroblasts = grep("Fibroblasts", colnames(data))
  Fibroblasts = data[, Fibroblasts, drop = FALSE]
  ##### Mast cells
  Mast = grep("Mast.cells", colnames(data))
  Mast = data[, Mast, drop = FALSE]
  Mast.activated = grep("Mast.activated", colnames(data))
  Mast.activated = data[, Mast.activated, drop = FALSE]
  Mast.resting = grep("Mast.resting", colnames(data))
  Mast.resting = data[, Mast.resting, drop = FALSE]
  ##### CAF cells 
  CAF = grep("CAF", colnames(data))
  CAF = data[, CAF, drop = FALSE]
  
  #####Output list
  cell_types = list(B, B.naive, B.memory, Macrophages, M0, M1, M2, Monocytes, Neutrophils, NK, NK.activated, NK.resting, NKT, CD4, CD4.memory.activated, CD4.memory.resting, CD4.naive,
                    CD8, Tregs, T.non.regs, Thelper, Tgamma, Dendritic, Dendritic.activated, Dendritic.resting, Cancer, Endothelial, Eosinophils, Plasma, Myocytes, Fibroblasts, Mast, Mast.activated,
                    Mast.resting, CAF)
  
  names(cell_types) = c("B.cells", "B.naive", "B.memory", "Macrophages.cells", "Macrophages.M0", "Macrophages.M1", "Macrophages.M2", "Monocytes", "Neutrophils", "NK.cells", "NK.activated", "NK.resting", "NKT.cells", "CD4.cells", "CD4.memory.activated",
                        "CD4.memory.resting", "CD4.naive", "CD8.cells", "T.cells.regulatory", "T.cells.non.regulatory","T.cells.helper", "T.cells.gamma.delta", "Dendritic.cells", "Dendritic.activated", "Dendritic.resting", "Cancer", "Endothelial",
                        "Eosinophils", "Plasma.cells", "Myocytes", "Fibroblasts", "Mast.cells", "Mast.activated", "Mast.resting", "CAF")
  
  ####Discarded cell types
  cell_types_matrix = cbind(B, B.naive, B.memory, Macrophages, M0, M1, M2, Monocytes, Neutrophils, NK, NK.activated, NK.resting, NKT, CD4, CD4.memory.activated, CD4.memory.resting, CD4.naive,
                      CD8, Tregs, T.non.regs, Thelper, Tgamma, Dendritic, Dendritic.activated, Dendritic.resting, Cancer, Endothelial, Eosinophils, Plasma, Myocytes, Fibroblasts, Mast, Mast.activated,
                      Mast.resting, CAF)
  
  cell_types_discarded = data[,!(colnames(data)%in%colnames(cell_types_matrix)), drop = F]
  
  return(list(cell_types, cell_types_discarded))
  
}

#' Remove highly correlated features
#' 
#' \code{removeCorrelatedFeatures} Remove one of two highly correlated features above a threshold 
#' 
#' @param data Matrix (samples X features).
#' @param threshold Cutoff to define above which corr number you want to consider highly correlated (default = 0.9).
#' @return Matrix with only one feature from a pair of highly correlated variables.
#' 
#' -------------------------------------------------------------------------------------------------------------
#' 
removeCorrelatedFeatures <- function(data, threshold, name, n_seed) {
  set.seed(n_seed)
  features_high_corr = c()
  cell_name = c()
  # Compute correlation matrix
  corr_matrix <- cor(data)
  # Find highly correlated features
  contador = 1
  while(nrow(corr_matrix)>0){
    feature = data.frame(corr_matrix[1, , drop = FALSE]) #Extract first row feature
    feature = feature %>%                                #Take only high corr above threshold
      mutate_all(~ifelse(. > threshold, ., NA)) %>%
      select_if(~all(!is.na(.)))
    
    corr_matrix = corr_matrix[-which(rownames(corr_matrix)%in%colnames(feature)),-which(colnames(corr_matrix)%in%colnames(feature)), drop = F] #Remove already joined features
    
    if(ncol(feature)>1){
      keep = colnames(feature)[sample(ncol(feature), size = 1)] #From high corr group, keep only one feature
      
      print(paste0("Highly correlated features (r>", threshold,"): ", paste(colnames(feature), collapse = ', ')))
      cat(paste0("Keeping only feature: ", keep, "\n\n")) 
      
      if(length(features_high_corr)>0){
        features_high_corr = c(features_high_corr, colnames(feature))
      }else{
        features_high_corr = colnames(feature)
      }
      
      feature = feature[,-which(colnames(feature)%in%keep), drop = F] 

      if(contador==1){
        new_data <- data[, -which(colnames(data)%in%colnames(feature)), drop = F] #Remove rest of the features from original data
      }else{
        new_data <- new_data[, -which(colnames(new_data)%in%colnames(feature)), drop = F] #Remove rest of the features from original data
      }
      contador = contador + 1
      cell_name = name
    }else{
      if(contador == 1){ 
        new_data = data
      }else{ #If it already started the loop
        new_data = new_data
      }
    }
  }
  
  if(length(cell_name)==0){
      cell_name = NULL
  }

  return(list(new_data, features_high_corr, cell_name))
}

#' Post-filtering after subgroups computation
#' 
#' \code{remove_subgroups} Remove subgroups that have the same method across different signatures 
#' 
#' @param groups List of Groups of features within cell types.
#' @return List of position of groups which have features of same method.
#' 
#' -------------------------------------------------------------------------------------------------------------
#' 
remove_subgroups = function(groups){
  lis = c()
  for (pos in 1:length(groups)){
    x = c()
    if(length(groups[[pos]])!=0){
      for (i in 1:length(groups[[pos]])) {
        x =  c(x,str_split(groups[[pos]][[i]], "_")[[1]][[1]])
      }
      if(length(unique(x)) == 1){
        lis = c(lis, pos)
      }
    }
  }
  
  return(lis)
} 

#' Post-filtering after subgroups computation
#' 
#' \code{compute_subgroups} Remove subgroups that have the same method across different signatures 
#' 
#' @param data Cell type deconvolution matrix from GEMDeCan output after splitting (samples X features).
#' @param thres_similarity Threshold for grouping features by proportionality.
#' @param thres_corr Threshold for grouping features by Pearson correlation.
#' @param thres_change Accepted variation between two subgroups to decide if keep subgrouping or finish the iteration.
#' @param file_name Name of the file.
#' 
#' @return List of subgroups of features within cell types. (data frame with subgroups, subgroups of corr, subgroups by similarity)
#' 
#' -------------------------------------------------------------------------------------------------------------
#' 
compute_subgroups = function(deconvolution, thres_corr, file_name){
  data = data.frame(deconvolution)
  cell_subgroups = list()
  cell_groups_similarity = list()
  cell_groups_discard = list()
  if (ncol(data) < 2) {
    warning("Deconvolution features with less than two columns for subgrouping (skipping)\n")
    return(list(data, cell_subgroups, cell_groups_similarity, cell_groups_discard))
  }else{
    #################### Proportionality-based correlation   
    is_similar <- function(value1, value2, threshold) {return(abs(value1 - value2) <= threshold)}
    similarity_matrix <- matrix(FALSE, nrow = ncol(data), ncol = ncol(data), dimnames = list(names(data), names(data)))
    for (col1 in names(data)) {
      for (col2 in names(data)) {
        similarity <- all(mapply(is_similar, data[[col1]], data[[col2]], MoreArgs = list(0.05))) #similarity threshold = 0.05
        similarity_matrix[col1, col2] <- similarity
      }
    }
    get_upper_tri <- function(cormat){
      cormat[lower.tri(cormat, diag = T)]<- NA
      return(cormat)
    }
    upper_tri <- get_upper_tri(similarity_matrix)
    x <- melt(upper_tri) %>%
      na.omit() %>%
      mutate_all(as.character)
    indice = 1
    subgroup = list()
    vec = unique(x$Var1)
    while(length(vec)>0){
      sub = x[which(x$Var1%in%vec[1]),] 
      sub = sub[which(sub$value==T),]
      if(nrow(sub)!=0){
        subgroup[[indice]] = c(vec[1], sub$Var2)
        x = x[-which(x$Var1%in%subgroup[[indice]]),] #Variable 1
        x = x[-which(x$Var2%in%subgroup[[indice]]),] #Variable 2
        vec = vec[-which(vec%in%subgroup[[indice]])]
        indice = indice + 1
      }else{
        indice = indice
        vec = vec[-1]
      }
    }
    
    if(length(subgroup)!=0){
      for (i in 1:length(subgroup)){ 
        names(subgroup)[i] = paste0(file_name, "_Subgroup.Similarity.", i) #Name subgroups
      }
      lis = remove_subgroups(subgroup) #Map subgroups with same method
      if(length(lis)>0){
        cell_groups_discard = subgroup[lis]
        subgroup = subgroup[-lis] #Remove subgroups if all subgroupped features belong to the same method
      }
      
      if(length(subgroup)!=0){  #check if after removal of subgroups with equal method, you still have subgroups
        cell_groups_similarity = subgroup
        data_sub = c()
        for(i in 1:length(cell_groups_similarity)){ #Create data frame with features subgroupped
          sub = data.frame(data[,colnames(data)%in%cell_groups_similarity[[i]]]) #Map features that are inside each subgroup from input (deconvolution)
          sub$median = rowMedians(as.matrix(sub), useNames = FALSE) #Compute median of subgroups across patients 
          data_sub = data.frame(cbind(data_sub, sub$median)) #Save median in a new data frame
          colnames(data_sub)[i] = names(cell_groups_similarity)[i]
          name = colnames(data)[which(!(colnames(data)%in%cell_groups_similarity[[i]]))]
          data = data[,-which(colnames(data)%in%cell_groups_similarity[[i]])] #Remove from deconvolution features that are subgrouped
          if(ncol(data.frame(data))==1){
            data = as.data.frame(data)
            colnames(data)[1] = name
          }
        }
        
        rownames(data_sub) = rownames(data) #List of patients
        data_sub = data.frame(data_sub[,colnames(data_sub)%in%names(cell_groups_similarity)])
        colnames(data_sub) = names(cell_groups_similarity)
        
        data = cbind(data, data_sub) #Join subgroups in deconvolution file
      }else{
        cell_groups_similarity = list()
      }
      k = 2
    }else{
      k = 3
    }
    if(ncol(data) == 1){ #everything is already subgroupped
      return(list(data, cell_subgroups, cell_groups_similarity, cell_groups_discard))  
    }
    
    #################### Linear-based correlation
    if(k==2 | k==3){
      terminate = FALSE
      iteration = 1
      while (terminate == FALSE) {
        corr_df <- correlation(data.matrix(data))
        vec = colnames(data)
        indice = 1
        subgroup = list()
        data_sub = c() 
        while(length(vec)>0){ #Keep running until no features are left
          if(vec[1] %in% corr_df$measure1){ #Check if feature still no-grouped
            tab = corr_df[corr_df$measure1 == vec[1],] #Take one feature against the others
            tab = tab[tab$r>thres_corr,] #Select features corr above the threshold 
            if(nrow(tab)!=0){ #If algorithm found features above corr
              subgroup[[indice]] = c(vec[1], tab$measure2) #Save features as subgroup
              idx = which(corr_df$measure1 %in% subgroup[[indice]])
              if(length(idx)>0){corr_df = corr_df[-idx,]} #Remove features already subgroupped 
              idy = which(corr_df$measure2 %in% subgroup[[indice]])
              if(length(idy)>0){corr_df = corr_df[-idy,]} #Remove features already subgroupped
              vec = vec[-which(vec%in%subgroup[[indice]])] #Remove feature already subgroupped from vector  
              indice = indice + 1
            }else{ #Condition when there is no correlation above the threshold (features no subgroupped)
              corr_df = corr_df[-which(corr_df$measure1 == vec[1]),] #Remove variable from corr matrix to keep subgrouping the others
              if(length(which(corr_df$measure2==vec[1]))>0){corr_df = corr_df[-which(corr_df$measure2 == vec[1]),]}
              vec = vec[-1] #Remove variable from vector to keep analyzing the others 
              indice = indice #Not increase index cause no subgroup appeared
            }
          }else{ #If feature is not in corr matrix it means that there is no any significant correlation against it and other features 
            vec = vec[-1] #Remove variable from vector to keep analyzing the others
            indice = indice  #Not increase index cause no subgroup appeared
          }
        }
        
        if(length(subgroup)!=0){
          for (i in 1:length(subgroup)){ #Name subgroups
            names(subgroup)[i] = paste0(file_name, "_Subgroup.", i, ".Iteration.", iteration)
          }
          ###Check whenever some subgroups belong to the same method
          if(iteration == 1){
            idx = remove_subgroups(subgroup) #Map subgroups with same method
            if(length(idx)>0){
              if(length(cell_groups_discard)>0){
                cell_groups_discard = c(cell_groups_discard, subgroup[idx])
                duplica = which(duplicated(cell_groups_discard)) #check if there are subgroups duplicated discarded
                if(length(duplica)>0){
                  cell_groups_discard = cell_groups_discard[-duplica]
                }
              }
              else{
                cell_groups_discard = subgroup[idx]
              }
              subgroup = subgroup[-idx] #Remove subgroups if all subgroupped features belong to the same method
            }
          }
          
          if(length(subgroup)!=0){ #check if after removal of subgroups with equal method, you still have subgroups (when iteration == 1)
            #Take median expression of subgroups
            for(i in 1:length(subgroup)){ #Create data frame with features subgroupped
              sub = data.frame(data[,colnames(data)%in%subgroup[[i]]]) #Map features that are inside each subgroup from input (deconvolution)
              sub$median = rowMedians(as.matrix(sub), useNames = FALSE) #Compute median of subgroup across patients 
              data_sub = data.frame(cbind(data_sub, sub$median)) #Save median in a new data frame
              colnames(data_sub)[i] = names(subgroup)[i]
              name = colnames(data)[which(!(colnames(data)%in%subgroup[[i]]))]
              data = data.frame(data[,-which(colnames(data)%in%subgroup[[i]])]) #Remove from deconvolution features that are subgrouped
              if(ncol(data.frame(data))==1){
                data = as.data.frame(data)
                colnames(data)[1] = name
              }
            }
            
            rownames(data_sub) = rownames(data) #List of patients
            
            if(iteration == 1){ #Save what is inside the first subgroups
              cell_subgroups = subgroup 
              data_sub = data.frame(data_sub[,colnames(data_sub)%in%names(cell_subgroups)])
              colnames(data_sub) = names(cell_subgroups)
            }else{
              for (i in 1:length(subgroup)) {
                cell_subgroups[[length(cell_subgroups)+1]] = subgroup[[i]]
                names(cell_subgroups)[length(cell_subgroups)] = names(subgroup)[i]
              }
            }
            
            if(ncol(data)!=0){
              data = cbind(data, data_sub)
            }else{
              data = data_sub #data will be 0 if all deconvolution features have been subgroupped
              terminate = TRUE
            }
            iteration = iteration + 1
          }else{
            terminate = TRUE #when the only subgroup that keep grouping is composed from the same method
          }
          
        }else{
          terminate = TRUE
        }
      }
      
      if(is.null(data_sub)==FALSE){
        data = cbind(data, data_sub)
      }
      
      idx = which(duplicated(t(data)))
      if(length(idx)>0){
        names = colnames(data)[idx]
        data = data[,-idx, drop = F]
        colnames(data) = names 
      }
    }
    
    return(list(data, cell_subgroups, cell_groups_similarity, cell_groups_discard))
  }
  
}

#' Pearson correlation for features
#' 
#' \code{correlation} Perform pairwise correlations between two matrices and return organized table of only significant corr
#' 
#' @param data Matrix (samples X features).
#' @return Matrix with corr, p-value, abs_corr and significant correlations.
#' 
#' -------------------------------------------------------------------------------------------------------------
#' 
correlation <- function(data) {
  
  M <- Hmisc::rcorr(as.matrix(data), type = "pearson")
  Mdf <- map(M, ~data.frame(.x))
  
  corr_df = Mdf %>%
    map(~rownames_to_column(.x, var="measure1")) %>%
    map(~pivot_longer(.x, -measure1, names_to = "measure2")) %>%
    bind_rows(.id = "id") %>%
    pivot_wider(names_from = id, values_from = value) %>%
    dplyr::rename(p = P) %>%
    mutate(sig_p = ifelse(p < .05, T, F),
           p_if_sig = ifelse(sig_p, p, NA),
           r_if_sig = ifelse(sig_p, r, NA)) 
  
  corr_df = na.omit(corr_df)  #remove the ones that are the same TFs (pval = NA)
  corr_df <- corr_df[which(corr_df$sig_p==T),]  #remove not significant
  corr_df <- corr_df[order(corr_df$r, decreasing = T),]
  corr_df$AbsR =  abs(corr_df$r)
  
  return(corr_df)
  
}

remove_low_variance <- function(data) {
  vars <- apply(data, 2, var)
  threshold = summary(vars)[[2]]
  cat("Checking distribution of variances...............................................................\n\n")
  cat("Chosen threshold is:", threshold, "\n\n")
  cat("Saving variance distribution plot in Results/ folder\n\n")
  low_variance <- which(vars < threshold)
  cat("Removing", length(low_variance), "features with variance across samples below this threshold...............................................................\n\n")
  
  pdf("Results/Distribution_variances_deconvolution.pdf")
  hist(vars, main = "Distribution of deconvolution variances across samples\nRemoving features below threshold (low variance)", xlab = "Variance", col = "skyblue", border = "white", xlim = range(vars))
  lines(density(vars), col = "red", lwd = 2)
  legend("topright", legend = c("Density", paste("Threshold =", round(threshold, 5))), col = c("red", "orange"), lty = c(1, 2), lwd = c(2, 2))
  
  # Shade region below threshold
  abline(v = threshold, col = "orange", lwd = 2, lty = 2)
  x <- density(vars)$x
  y <- density(vars)$y
  polygon(c(min(x[vars < threshold]), x[vars < threshold], max(x[vars < threshold])), 
          c(0, y[vars < threshold], 0), col = adjustcolor("orange", alpha.f = 0.3), border = NA)
  dev.off()
  
  data_filt = data[,-low_variance, drop = F]
  low_var_features = data[,low_variance, drop = F]
  
  res = list(data_filt, low_var_features)
  return(res)
}

computeQuantiseq <- function(TPM_matrix) {
  require(immunedeconv)
  TPM_matrix = TPM_matrix[rownames(TPM_matrix)%in%rownames(immunedeconv::dataset_racle$expr_mat),] #To avoid problems regarding gene names (quantiseq error)
  
  quantiseq = immunedeconv::deconvolute(TPM_matrix, "quantiseq", tumor = T) %>% 
    column_to_rownames("cell_type") %>%
    t()
  
  colnames(quantiseq) = paste0("Quantiseq_", colnames(quantiseq))
  colnames(quantiseq) <- colnames(quantiseq) %>%
    str_replace_all(., " ", "_")
  
  return(quantiseq)
}

computeMCP <- function(TPM_matrix, genes_path) {
  require(MCPcounter)
  genes <- read.table(paste0(genes_path, "/MCPcounter/MCPcounter-genes.txt"), sep = "\t", stringsAsFactors = FALSE, header = TRUE, colClasses = "character", check.names = FALSE)
  mcp <- MCPcounter.estimate(TPM_matrix, genes = genes, featuresType = "HUGO_symbols", probesets = NULL) %>%
    t()
  
  colnames(mcp) = paste0("MCP_", colnames(mcp))
  colnames(mcp) <- colnames(mcp) %>%
    str_replace_all(., " ", "_")
  
  return(mcp)
}

computeXCell <- function(TPM_matrix) {
  require(immunedeconv)
  
  xcell = immunedeconv::deconvolute(TPM_matrix, "xcell") %>% 
    column_to_rownames("cell_type") %>%
    t()
  
  colnames(xcell) = paste0("XCell_", colnames(xcell))
  colnames(xcell) <- colnames(xcell) %>%
    str_replace_all(., " ", "_")
  
  return(xcell)
}

computeCBSX = function(TPM_matrix, signature_file, name, password, name_signature){
  require(omnideconv)
  set_cibersortx_credentials(name, password)
  cbsx = omnideconv::deconvolute_cibersortx(TPM_matrix, signature_file)
  
  colnames(cbsx) = paste0("CBSX_", name_signature, "_", colnames(cbsx))
  colnames(cbsx) <- colnames(cbsx) %>%
    str_replace_all(., " ", "_")
  
  return(cbsx)
}

computeEpiDISH = function(TPM_matrix, signature_file, name_signature){
  require(EpiDISH)
  epi <- epidish(TPM_matrix, as.matrix(signature_file), method = "RPC", maxit = 200)
  epidish = epi$estF
  
  colnames(epidish) = paste0("Epidish_", name_signature, "_", colnames(epidish))
  colnames(epidish) <- colnames(epidish) %>%
    str_replace_all(., " ", "_")
  
  return(epidish)
}

computeDeconRNASeq = function(TPM_matrix, signature_file, name_signature){
  require(DeconRNASeq)
  decon <- DeconRNASeq(TPM_matrix, data.frame(signature_file))
  deconRNAseq = decon$out.all
  rownames(deconRNAseq) = colnames(TPM_matrix)
  
  colnames(deconRNAseq) = paste0("DeconRNASeq_", name_signature, "_", colnames(deconRNAseq))
  colnames(deconRNAseq) <- colnames(deconRNAseq) %>%
    str_replace_all(., " ", "_")
  
  return(deconRNAseq)
}


compute_methods_variable_signature = function(TPM_matrix, signatures, methods = c("CBSX", "Epidish", "DeconRNASeq"), exclude = NULL, cbsx.name, cbsx.token){
  
  db=list.files(signatures, full.names = T, pattern = "\\.txt$")
  name_exclude = c()
  if(is.null(methods)==F){
    cat("\nThe following method-signature combinations are going to be calculated...............................................................\n")
    
    cat("\nMethods\n")
    for (method in methods) {
      cat("* ", method, "\n", sep = "")
    }
    cat("\nSignatures\n")
    for (i in 1:length(db)) {
      name = str_split(basename(db[[i]]), "\\.")[[1]][1]
      cat("* ", name, "\n", sep = "")
      if(is.null(exclude)==F && name %in% exclude){
        name_exclude = c(name_exclude, name)
      }
    }
    
    if(length(name_exclude)>0){
      cat("\nExcluding signatures: ", paste0(name_exclude, collapse = ", "), "\n")
    }
    
    deconvolution = list()
    
    if("CBSX" %in% methods){
      if(is.null(cbsx.name)==T || is.null(cbsx.token)==T){
        cat("\nYou select to run CBSX but no credentials were found")
        cat("\nPlease set your credentials in the function for running CibersortX")
        stop()
      }
    }
    
    for (i in 1:length(db)) {
      signature <- read.delim(db[[i]], row.names=1)
      signature_name = str_split(basename(db[[i]]), "\\.")[[1]][1]
      if(!is.null(exclude) && signature_name %in% exclude) {
        next
      }else{
        if("DeconRNASeq"%in%methods){
          cat("\nRunning DeconRNASeq...............................................................\n\n")
          deconrnaseq <- computeDeconRNASeq(TPM_matrix, signature, signature_name)}
        if("Epidish"%in%methods){
          cat("\nRunning Epidish...............................................................\n\n")
          epidish <- computeEpiDISH(TPM_matrix, signature, signature_name)}
        if("CBSX"%in%methods){
          cat("\nRunning CBSX...............................................................\n\n")
          cbsx <- computeCBSX(TPM_matrix, signature, cbsx.name, cbsx.token, signature_name)}
        
        combined_data <- NULL
        if (exists("deconrnaseq")) {
          combined_data <- deconrnaseq
        }
        if (exists("epidish")) {
          if (is.null(combined_data)) {
            combined_data <- epidish
          } else {
            combined_data <- cbind(combined_data, epidish)
          }
        }
        if (exists("cbsx")) {
          if (is.null(combined_data)) {
            combined_data <- cbsx
          } else {
            combined_data <- cbind(combined_data, cbsx)
          }
        }
        
        deconvolution[[i]] <- combined_data
      }
    }
    
    deconv = do.call(cbind, deconvolution)
    
    return(deconv)
  }else{
    cat("\nNo methods to be calculated using variable signatures.")
    return(NULL)
  }
  
}


compute.deconvolution <- function(counts.matrix, methods = c("Quantiseq", "MCP", "xCell", "CBSX", "Epidish", "DeconRNASeq"), signatures_exclude = NULL, normalized = F, 
                                  credentials.mail = NULL, credentials.token = NULL){
  
  path_signatures = 'src/signatures'
  
  if(normalized == T){
    cat("Performing TPM normalization log transformed...............................................................\n\n")
    TPM_matrix = data.frame(ADImpute::NormalizeTPM(counts.matrix, log=T)) 
  }else{
    TPM_matrix = data.frame(counts.matrix)
  }

  cat("Running deconvolution using the following methods...............................................................\n\n")
  for (method in methods) {
    cat("* ", method, "\n", sep = "")
  }
  
  if("Quantiseq" %in% methods){
    cat("\nRunning Quantiseq...............................................................\n")
    quantiseq = computeQuantiseq(TPM_matrix)}
  if("MCP" %in% methods){
    cat("\nRunning MCPCounter...............................................................\n")
    mcp = computeMCP(TPM_matrix, path_signatures)}
  if("xCell" %in% methods){
    xcell = computeXCell(TPM_matrix)
    cat("\nRunning XCell...............................................................\n")}
  
  default_sig = c("Quantiseq", "MCP", "xCell")
  methods = methods[!(methods %in% default_sig)]
  if(length(methods) == 0){
    methods = NULL
  }
  deconv_sig = compute_methods_variable_signature(TPM_matrix, signatures = path_signatures, method = methods, exclude = signatures_exclude, cbsx.name = credentials.mail, cbsx.token = credentials.token)
  
  deconv_default <- NULL
  if (exists("quantiseq")) {
    deconv_default <- quantiseq
  }
  if (exists("mcp")) {
    if (is.null(deconv_default)) {
      deconv_default <- mcp
    } else {
      deconv_default <- cbind(deconv_default, mcp)
    }
  }
  if (exists("xcell")) {
    if (is.null(deconv_default)) {
      deconv_default <- xcell
    } else {
      deconv_default <- cbind(deconv_default, xcell)
    }
  }
  
  if(is.null(deconv_sig)){
    all_deconvolution_table = deconv_default
  }else{
    all_deconvolution_table = cbind(deconv_default, deconv_sig)
  }
  
  deconvolution = compute.deconvolution.preprocessing(data.frame(all_deconvolution_table))
  
  return(deconvolution)
  
}

