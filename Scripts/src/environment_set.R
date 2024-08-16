
##Functions for analyze data for LungPredict1 paper

# author: Marcelo Hurtado
# email: marcelo.hurtado@inserm.fr
# organization: INSERM CRCT - Pancaldi team 21
# place: Toulouse, France

##Load libraries

libraries_set <- function(){
  suppressMessages(library("BiocManager"))
  suppressMessages(library("devtools"))
  suppressMessages(library("remotes"))
  suppressMessages(library("decoupleR"))
  suppressMessages(library("OmnipathR"))
  suppressMessages(library("tidyr"))
  suppressMessages(library("dplyr"))
  suppressMessages(library("matrixStats"))
  suppressMessages(library("org.Hs.eg.db"))
  suppressMessages(library("ReactomePA"))
  suppressMessages(library("WGCNA"))
  suppressMessages(library("reshape2"))
  suppressMessages(library("purrr"))
  suppressMessages(library("tidygraph"))
  suppressMessages(library("stringr"))
  suppressMessages(library("tibble"))
  suppressMessages(library("gplots"))
  suppressMessages(library("ggplot2"))
  suppressMessages(library("AnnotationDbi"))
  suppressMessages(library("DESeq2"))
  suppressMessages(library("RColorBrewer"))
  suppressMessages(library("pheatmap"))
  suppressMessages(library("ggfortify"))
  suppressMessages(library("viper"))
  suppressMessages(library("msigdbr"))
  suppressMessages(library("GSVA"))
  suppressMessages(library("clusterProfiler"))
  suppressMessages(library("Hmisc"))
  suppressMessages(library("DOSE"))
  suppressMessages(library("fgsea"))
  suppressMessages(library("ggpubr"))
  suppressMessages(library("ComplexHeatmap"))
  suppressMessages(library("ggstatsplot"))
  suppressMessages(library("dendextend"))
  suppressMessages(library("PCAtools"))
  suppressMessages(library("stats"))
  suppressMessages(library("Boruta"))
  suppressMessages(library("survival"))
  suppressMessages(library("survminer"))
  suppressMessages(library("rms"))
  suppressMessages(library("igraph"))
  suppressMessages(library("enrichplot"))
  suppressMessages(library("uuid"))
}

libraries_set()

source("src/cell_deconvolution.R") #Call deconvolution functions

dir.create(file.path(getwd(), "Results"))

### Functions

vsd_normalization = function(raw.counts, coldata){
  dds <- DESeqDataSetFromMatrix(countData = round(raw.counts),
                                colData = coldata,
                                design = ~1)
  dds <- estimateSizeFactors(dds)
  dds <- estimateDispersions(dds)
  normalized.counts <- counts(dds, normalized=TRUE)
  
  return(normalized.counts)
}

remove.outliers = function(counts){
  sampleDists <- dist(t(counts))
  sampleDistMatrix <- as.matrix(sampleDists)
  colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)
  pdf("Results/SampleDistanceMatrix.pdf", width = 16, height = 14)
  print(pheatmap(sampleDistMatrix,
           clustering_distance_rows=sampleDists,
           clustering_distance_cols=sampleDists,
           col=colors))
  dev.off()
  sampleTree = hclust(dist(t(counts)), method = "average");
  pdf("Results/SamplesClustering.pdf", width = 16)
  plot(sampleTree, main = "Sample clustering to detect outliers", sub="", xlab="", cex.lab = 1.5, cex.axis = 1.5, cex.main = 2)
  dev.off()
  response = readline(prompt = "Do you have outliers? (Y|N): ")
  if(response == 'Y'){
    cutheight = as.numeric(readline(prompt = "Cut height parameter: "))
    pdf("Results/SamplesClusteringOutliers.pdf", width = 16)
    plot(sampleTree, main = "Sample clustering to detect outliers", sub="", xlab="", cex.lab = 1.5, cex.axis = 1.5, cex.main = 2)
    abline(h = cutheight, col = "red")
    dev.off()
    clust = cutreeStatic(sampleTree, cutHeight = cutheight, minSize = 10) # Determine cluster under the line
    table(clust) # clust 1 contains the samples we want to keep.
    keepSamples = (clust==1)
    counts = counts[,keepSamples]
  }
  
  return(counts)
}

compute_pca_analysis <- function(data, coldata, trait, trait2 = NULL, ncomp = 5, cortype = "pearson"){
  #data: genes as rows and samples as columns
  pca_res <- prcomp(t(data))
  coldata[,colnames(coldata)%in%trait] = factor(coldata[,colnames(coldata)%in%trait])
  print(autoplot(pca_res, data = coldata, colour = trait, shape = trait2, size  = 3))
  
  p <- PCAtools::pca(data, metadata = coldata, removeVar = 0.1)
  
  peigencor  <- eigencorplot(p,
                             components = PCAtools::getComponents(p, 1:ncomp),
                             metavars = colnames(coldata),
                             col = c('white', 'cornsilk1', 'gold', 'forestgreen', 'darkgreen'),
                             cexCorval = 1.2,
                             fontCorval = 2,
                             posLab = 'all',
                             rotLabX = 45,
                             scale = TRUE,
                             main = paste("Principal component", cortype, "r^2 clinical correlates"),
                             plotRsquared = TRUE,
                             corFUN = cortype,
                             corUSE = 'pairwise.complete.obs',
                             corMultipleTestCorrection = 'BH',
                             signifSymbols = c('****', '***', '**', '*', ''),
                             signifCutpoints = c(0, 0.0001, 0.001, 0.01, 0.05, 1))
  
  print(peigencor)
}

compute.TFs.activity <- function(RNA.counts, TF.collection = "CollecTRI", min_targets_size = 5, tfs.pruned = FALSE){
  
  tfs2viper_regulons <- function(df){
    regulon_list <- split(df, df$source)
    regulons <- lapply(regulon_list, function(regulon) {
      tfmode <- stats::setNames(regulon$mor, regulon$target)
      list(tfmode = tfmode, likelihood = rep(1, length(tfmode)))
    })
    return(regulons)}
  
  if(tfs.pruned==T){
    cat("Pruned TFs is set to TRUE. Please specify the maximun size of targets allowed/n")
    max_size_targets = as.numeric(readline(prompt = "Maximun size of TFs-targets: "))
  }
  
  if(TF.collection == "CollecTRI"){
    net = decoupleR::get_collectri(organism = 'human', split_complexes = F)
    net_regulons = tfs2viper_regulons(net)
  } else if(TF.collection == "Dorothea"){
    net = decoupleR::get_dorothea(organism = 'human', levels = c("A", "B", "C", "D"))
    net_regulons = tfs2viper_regulons(net)
  } 
  
  if(TF.collection == "ARACNE"){
    cat("For ARACNE analysis you need to specify the path of your network file. Remember this file should be a 3 columns text file, with regulator in the first column, target in the second and mutual information in the third column")
    network_file = readline(prompt = "Path for network file from aracne (no quotes): ")
    net_regulons <- aracne2regulon(network_file, as.matrix(RNA.counts), format = "3col")
  }
  
  if(tfs.pruned == TRUE){
    net_regulons = pruneRegulon(net_regulons, cutoff = max_size_targets)
  }
  
  sample_acts <- viper(as.matrix(RNA.counts), net_regulons, minsize = min_targets_size, verbose=F, method = "scale")
  message("TFs scores computed")
  
  return(data.frame(t(sample_acts)))
  
}

compute.pathway.activity <- function(RNA.tpm, gene_sets = NULL){

  RNA.tpm = as.matrix(RNA.tpm)
  #Get universe
  paths <- get_progeny(organism = 'human', top = 500)
  # Run mlm
  progeny <- run_mlm(mat=RNA.tpm, net=paths, .source='source', .target='target', .mor='weight', minsize = 5)
  # Transform to wide matrix
  sample_acts_progeny <- progeny %>%
    pivot_wider(id_cols = 'condition', names_from = 'source',
                values_from = 'score') %>%
    column_to_rownames('condition') %>%
    as.matrix()
  
  if(is.null(gene_sets)==F){
    
    cat("Computing GSVA analysis using provided gene sets.....................................................\n")
    
    gsva_results <- gsva(
      RNA.tpm,
      gene_sets,
      method = "gsva",
      kcdf = "Gaussian",
      min.sz = 1,
      mx.diff = TRUE,
      verbose = TRUE
    )
    sample_acts_hallmarks <- data.frame(scale(t(gsva_results)))
    return(list(sample_acts_progeny, sample_acts_hallmarks))
  }

  # Scale per feature
  sample_acts_progeny <- data.frame(scale(sample_acts_progeny))
  
  message("Pathways scores computed")
  return(sample_acts_progeny)
  
}

compute.modules.enrichment <- function(RNA.tpm, hub_tfs){
  net = get_collectri(organism = 'human', split_complexes = F) #Get universe  
  res = list()
  pathways = list()
  
  #Pathway enrichment using target genes
  for (i in 1:length(hub_tfs[[1]])) {
    color = names(hub_tfs[[1]])[i]
    res[[i]] = module_enrich(as.matrix(RNA.tpm), color, hub_tfs, net)
    names(res)[i] = color
    pathways[[i]] = res[[i]]@result[["Description"]]
    names(pathways)[i] = color
  }

  ItemsList <- venn(pathways, show.plot = FALSE)
  x = attributes(ItemsList)
  
  #Keep only unique pathways 
  for (i in 1:length(hub_tfs[[1]])) {
    color = names(hub_tfs[[1]])[i]
    res[[i]]@result = res[[i]]@result[which(res[[i]]@result[["Description"]]%in%x$intersections[[color]]),] #take unique pathways per module
    if(nrow(res[[i]]@result)==0){
      print(paste0("No enrichment for module ", color))
    }else{
    pdf(paste0("Results/Enrichment by Reactome\nModule ", color))
    print(dotplot(res[[i]],  title=paste0("Enrichment by Reactome\nModule ", color)))
    dev.off()
    }
  }
  
}

module_enrich = function(tpm.counts, module_color, hub_genes, tfs_universe){

  targets = tfs_universe[tfs_universe$source %in% hub_genes[[1]][[module_color]],] #Extract targets from TFs
  targets = unique(targets$target) #Keep only unique targets from TFs
  
  targets_genes = tpm.counts[rownames(tpm.counts)%in%targets,] #Extract gene expression from targets
  targets_genes = targets_genes[order(rowVars(targets_genes), decreasing = T),][1:round(0.2*nrow(targets_genes)),] #Keep only highly variable targets (20%)
  
  entrz <- AnnotationDbi::select(org.Hs.eg.db, keys = rownames(targets_genes), columns = "ENTREZID", keytype = "SYMBOL") #Change to EntrezID
  universe = AnnotationDbi::select(org.Hs.eg.db, keys = rownames(tpm.counts), columns = "ENTREZID", keytype = "SYMBOL") #Change to EntrezID
  
  reac <- enrichPathway(gene    = entrz$ENTREZID,
                        organism     = 'human',
                        universe = universe$ENTREZID,
                        pvalueCutoff = 0.05)
  
  reac@result = reac@result[reac@result$p.adjust<0.05,]
  
  return(reac)
  
}

identify_hub_TFs <- function(datExpr, moduleEigengenes, moduleColors, MM_thresh = 0.8, degree_thresh = 0.9) {
  # Calculate Module Membership (MM)
  moduleMemberships <- sapply(unique(moduleColors), function(module) {
    genesInModule <- which(moduleColors == module)
    eigengene <- moduleEigengenes[, paste0("ME", module)]
    cor(t(datExpr[genesInModule, ]), eigengene)
  })
  
  # Calculate adjacency matrices and degrees
  adjacencyList <- lapply(unique(moduleColors), function(module) {
    genesInModule <- which(moduleColors == module)
    moduleData <- datExpr[genesInModule, ]
    adjacency <- cor(t(moduleData))
    adjacency[lower.tri(adjacency)] <- 0
    adjacency
  })
  
  names(adjacencyList) <- unique(moduleColors)
  
  moduleDegrees <- lapply(unique(moduleColors), function(module) {
    adjacency <- adjacencyList[[module]]
    degree <- rowSums(adjacency)
    names(degree) <- rownames(datExpr)[which(moduleColors == module)]
    degree
  })
  names(moduleDegrees) <- unique(moduleColors)
  
  # Identify hub genes
  hubGenesList <- list()
  hubGenesData <- data.frame()  # Initialize an empty data frame
  
  allDegrees <- numeric()  
  allMemberships <- numeric() 
  
  for (module in unique(moduleColors)) {
    genesInModule <- which(moduleColors == module)
    degrees <- moduleDegrees[[module]]
    memberships <- moduleMemberships[[module]]
    
    # Update the global lists of degrees and memberships
    allDegrees <- c(allDegrees, degrees)
    allMemberships <- c(allMemberships, memberships)
    
    # Calculate the cutoff for the top 10% by degree
    degreeCutoff <- quantile(degrees, degree_thresh)
    
    # Plot distribution of Degrees for the current module
    degree_plot_data <- data.frame(Degree = degrees)

    degree_plot <- ggplot(degree_plot_data, aes(x = Degree)) +
      geom_histogram(binwidth = 5, fill = module, color = "black", alpha = 0.6) +
      geom_vline(xintercept = degreeCutoff, linetype = "dashed", color = "red") +
      labs(
        title = paste("Distribution of Gene Degrees in Module", module),
        x = "Degree",
        y = "Frequency"
      ) +
      theme_minimal()
    
    print(degree_plot)

    # Plot distribution of Module Membership for the current module
    membership_plot_data <- data.frame(ModuleMembership = memberships)
    
    membership_plot <- ggplot(membership_plot_data, aes(x = ModuleMembership)) +
      geom_histogram(binwidth = 0.05, fill = module, color = "black", alpha = 0.6) +
      geom_vline(xintercept = MM_thresh, linetype = "dashed", color = "red") +
      labs(
        title = paste("Distribution of Module Membership in Module", module),
        x = "Module Membership",
        y = "Frequency"
      ) +
      theme_minimal()
    
    print(membership_plot)
    # Identify hub genes (top 10% by degree)
    hubGenes <- names(degrees[degrees >= degreeCutoff])
    if (length(hubGenes) == 0) next
    
    # Identify hub genes (>0.8 MM)
    highMMGenes <- which(memberships > MM_thresh)
    if (length(highMMGenes) == 0) next
    
    # Filter hub genes (degree) to only include those with high MM 
    finalHubGenes <- intersect(hubGenes, rownames(memberships)[highMMGenes])
    if (length(finalHubGenes) == 0) next
    
    hubGenesList[[module]] = finalHubGenes
    
    # Create a dataframe with detailed information for hub genes only
    moduleData <- data.frame(
        Module = module,
        ModuleMembership = unname(memberships[finalHubGenes,]),
        Degree = degrees[finalHubGenes],
        stringsAsFactors = FALSE)
    
    # Append to the combined dataframe
    hubGenesData <- rbind(hubGenesData, moduleData)
  }
  
  return(list(hubGenes = hubGenesList, detailedData = hubGenesData))
}

compute.WTCNA <- function(TFs.matrix, network.type = "signed", clustering.method = "ward.D2", minMod = 15, corr_mod = 0.9, cor_type = "p"){
  
  #####Choose parameter for scale-free network topology
  powers = c(c(1:10), seq(from = 12, to=20, by=2))
  sft = pickSoftThreshold(TFs.matrix, powerVector = powers, verbose = 5, networkType = network.type)
  pdf("Results/Soft Threshold")
  plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
             xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit,signed R^2",type="n",
             main = paste("Scale independence"))
  text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2], labels=powers,cex=0.9,col="red");
  abline(h=0.90,col="red")
  dev.off()
  softPower = as.numeric(readline(prompt = "Enter soft thresholding parameter: "))
  if(is.na(softPower)==T){
    softPower = as.numeric(readline(prompt = "Enter soft thresholding parameter: "))
  }
  #####Co-expression matrix using nodes adjacency and topological overlapping nodes
  adjacency = adjacency(TFs.matrix, power =softPower, type=network.type, corFnc = "cor", corOptions = list(use = cor_type))
  TOM = TOMsimilarity(adjacency, TOMType = network.type)
  dissTOM = 1-TOM
  
  #####Unsupervised hierarchical clustering using dissimilarity matrix
  geneTree = hclust(dist(dissTOM), method = clustering.method)
  dynamicMods = cutreeDynamic(dendro = geneTree, distM = dissTOM,
                              deepSplit = 2, pamRespectsDendro = FALSE,
                              minClusterSize = minMod);
  dynamicColors = labels2colors(dynamicMods)
  pdf("Results/Gene dendrogram and module colors")
  plotDendroAndColors(geneTree, dynamicColors, "Dynamic Tree Cut",
                            dendroLabels = FALSE, hang = 0.03,
                            addGuide = TRUE, guideHang = 0.05,
                            main = "Gene dendrogram and module colors")
  dev.off()
  
  #####Calculate eigenvectors from modules
  MEList = moduleEigengenes(TFs.matrix, colors = dynamicColors, scale = F) #Data already scale
  MEs = MEList$eigengenes
  MEs = orderMEs(MEs)
  
  #####Merge modules that are significantly correlated
  mergeModules = function(data, colors, corr){
    df = correlation(data)
    print(df)
    idx = which(round(df$r,2) > corr)
    if(length(idx)>0){
      for(i in seq(1, length(idx), by=2)){ 
        module1 = df$measure1[idx[i]]
        module2 = df$measure2[idx[i]]
        if((module1 %in% colnames(data)) && (module2 %in% colnames(data))){
          colors[which(colors%in%c(substring(module1, 3), substring(module2, 3)))] = substring(module1, 3)
          data <- data %>%
            mutate(new_column = rowMeans(dplyr::select(., module1, module2))) %>%
            dplyr::rename(module1 = new_column) %>%
            dplyr::select(., -module1, -module2) 
        }
      }
    }
    
    return(list(data, colors))
  }
  
  print(paste0("Merging modules significantly correlated with ", corr_mod, "........"))
  merge = mergeModules(MEs, dynamicColors, corr_mod)
  MEs = merge[[1]]
  dynamicColors = merge[[2]]
  modtfs = list()
  for(i in 1:ncol(MEs)){
    tfs = colnames(TFs.matrix)
    modules = c(substring(names(MEs), 3))
    inModule = is.finite(match(dynamicColors,modules[i]))
    modtfs[[i]] = tfs[inModule]
  }
  names(modtfs) = modules
  
  pdf("Results/Gene dendrogram and module colors after merging")
  plotDendroAndColors(geneTree, dynamicColors, "Dynamic Tree Cut",
                      dendroLabels = FALSE, hang = 0.03,
                      addGuide = TRUE, guideHang = 0.05,
                      main = "Gene dendrogram and module colors")
  dev.off()
     
  TFspropVar = propVarExplained(TFs.matrix, dynamicColors, MEs, corFnc = "cor", corOptions = "use = 'p'") 

  output = list(MEs, dynamicColors, modtfs, TFspropVar)
  names(output) = c("TFs module matrix", "TFs colors", "TFs per module", "Proportion of variance")
  
  contador = 1
  tfs_modules = data.frame(matrix(nrow = length(modtfs), ncol = 2))
  colnames(tfs_modules) = c("TFs module", "Composition")
  for (i in 1:length(modtfs)) {
    tfs_modules[contador,1] = names(modtfs)[i]
    tfs_modules[contador,2] = paste(modtfs[[i]], collapse = ",")
    contador = contador + 1
  }
  write.csv(tfs_modules, 'Results/TFs_modules.csv', row.names = F)
  
  return(output)
  
}

compute.modules.relationship <- function(tfs_network, matB, file_name, width = 8, height = 8, pval=0.05, padj = F, cor_type = "p", return = F, vertical=F){
  
  tfs_network = data.frame(tfs_network)
  matB = data.frame(matB)
  
  ##check if names from both features are the same
  if(all(rownames(tfs_network)==rownames(matB)) == F){
    stop("No equal names, verify the input objects")
  }
  
  
  moduleTraitCor = cor(tfs_network, matB, method = cor_type)
  moduleTraitPvalue = corPvalueStudent(moduleTraitCor, nrow(tfs_network))
  
  rev = which(colSums(moduleTraitPvalue > pval)==nrow(moduleTraitPvalue)) #check if there are features no significant with any module
  
  if(length(rev)>0){
    moduleTraitCor = moduleTraitCor[,-rev]
    moduleTraitPvalue = moduleTraitPvalue[,-rev] 
  }
  
  if(padj == T){
    for (i in 1:ncol(moduleTraitPvalue)) {
      moduleTraitPvalue[,i] = p.adjust(moduleTraitPvalue[,i], method = 'bonferroni')
    }
  }
  
  ##Plot in vertical
  if(vertical == T){
    if(ncol(data.frame(moduleTraitCor))>1){
      #Extract significant features per trait
      sig = list()
      for (i in 1:nrow(moduleTraitPvalue)) {
        sig[[i]] = names(which(signif(moduleTraitPvalue[i,],2)<=pval))
      }
      names(sig) = substring(rownames(moduleTraitPvalue), 3)
      
      if(return == T){
        retu = list(moduleTraitCor, sig)
        return(retu)
      }else{
        d <- dist(t(moduleTraitCor), method = "manhattan")
        hc1 <- hclust(d, method = "ward.D2")
        vec = hc1[["order"]]
        textMatrix = paste(signif(moduleTraitCor, 2), "\n(", signif(moduleTraitPvalue, 2), ")", sep = "")
        dim(textMatrix) = dim(moduleTraitCor)
        idx = which(round(moduleTraitPvalue,2)>pval)
        for (i in idx) {
          textMatrix[i] = NA
        }
        textMatrix = t(textMatrix)
        moduleTraitCor = data.frame(t(moduleTraitCor))
        pdf(paste0("Results/",file_name), width = width, height = height)
        par(mar = c(3, 25, 5, 3))
        labeledHeatmap(Matrix = moduleTraitCor[vec,],
                       xLabels = colnames(moduleTraitCor),
                       yLabels = rownames(moduleTraitCor[vec,]),
                       xLabelsPosition = "top",
                       colors = blueWhiteRed(50),
                       textMatrix = textMatrix[vec,],
                       setStdMargins = F,
                       cex.text = 0.5,
                       zlim = c(-1,1))
        dev.off()
      }}else{
        #Extract significant features per trait
        sig = list()
        for (i in 1:nrow(moduleTraitPvalue)) {
          sig[[i]] = colnames(moduleTraitPvalue)[which(signif(moduleTraitPvalue[i,],2)<=pval)]
        }
        names(sig) = substring(rownames(moduleTraitPvalue), 3)
        if(return == T){
          retu = list(moduleTraitCor, sig)
          return(retu)
        }else{
          textMatrix = paste(signif(moduleTraitCor, 2), "\n(", signif(moduleTraitPvalue, 2), ")", sep = "")
          idx = which(round(moduleTraitPvalue,2)>pval)
          for (i in idx) {
            textMatrix[i] = NA
          }
          textMatrix = t(textMatrix)
          moduleTraitCor = data.frame(t(moduleTraitCor))
          colnames(moduleTraitCor)[1] = colnames(matB)
          pdf(paste0("Results/",file_name), width = width, height = height)
          par(mar = c(25, 15, 3, 3))
          labeledHeatmap(Matrix = moduleTraitCor,
                         xLabels = colnames(moduleTraitCor),
                         yLabels = rownames(moduleTraitCor),
                         xLabelsPosition = "top",
                         colors = blueWhiteRed(50),
                         textMatrix = textMatrix,
                         setStdMargins = F,
                         cex.text = 0.5,
                         zlim = c(-1,1))
          dev.off()
        }}}
  
  
  ###Plot in horizontal
  if(vertical == F){
    if(ncol(data.frame(moduleTraitCor))>1){
      #Extract significant features per trait
      sig = list()
      for (i in 1:nrow(moduleTraitPvalue)) {
        sig[[i]] = names(which(signif(moduleTraitPvalue[i,],2)<=pval))
      }
      names(sig) = substring(rownames(moduleTraitPvalue), 3)
      if(return==T){
        retu = list(moduleTraitCor, sig)
        return(retu)
      }else{
        d <- dist(t(moduleTraitCor), method = "manhattan")
        hc1 <- hclust(d, method = "ward.D2")
        vec = hc1[["order"]]
        textMatrix = paste(signif(moduleTraitCor, 2), "\n(", signif(moduleTraitPvalue, 2), ")", sep = "")
        dim(textMatrix) = dim(moduleTraitCor)
        idx = which(round(moduleTraitPvalue,2)>pval)
        for (i in idx) {
          textMatrix[i] = NA
        }
        moduleTraitCor = data.frame(moduleTraitCor)
        pdf(paste0("Results/",file_name), width = width, height = height)
        par(mar = c(25, 15, 3, 3))
        labeledHeatmap(Matrix = moduleTraitCor[,vec],
                       xLabels = names(moduleTraitCor[,vec]),
                       yLabels = rownames(moduleTraitCor),
                       ySymbols = rownames(moduleTraitCor),
                       colorLabels = FALSE,
                       colors = blueWhiteRed(50),
                       textMatrix = textMatrix[,vec],
                       setStdMargins = FALSE,
                       cex.text = 0.5,
                       zlim = c(-1,1),
                       main = paste("Module-trait relationships"))
        dev.off()
      }}else{
        #Extract significant features per trait
        sig = list()
        for (i in 1:nrow(moduleTraitPvalue)) {
          sig[[i]] = colnames(moduleTraitPvalue)[which(signif(moduleTraitPvalue[i,],2)<=pval)]
        }
        names(sig) = substring(rownames(moduleTraitPvalue), 3)
        
        if(return == T){
          retu = list(moduleTraitCor, sig)
          return(retu)
        }else{
          textMatrix = paste(signif(moduleTraitCor, 2), "\n(", signif(moduleTraitPvalue, 2), ")", sep = "")
          idx = which(round(moduleTraitPvalue,2)>pval)
          for (i in idx) {
            textMatrix[i] = NA
          }
          moduleTraitCor = data.frame(moduleTraitCor)
          colnames(moduleTraitCor)[1] = colnames(matB)
          pdf(paste0("Results/",file_name), width = width, height = height)
          par(mar = c(25, 15, 3, 3))
          labeledHeatmap(Matrix = moduleTraitCor,
                         xLabels = names(moduleTraitCor),
                         yLabels = rownames(moduleTraitCor),
                         ySymbols = rownames(moduleTraitCor),
                         colorLabels = FALSE,
                         colors = blueWhiteRed(50),
                         textMatrix = textMatrix,
                         setStdMargins = FALSE,
                         cex.text = 0.5,
                         zlim = c(-1,1),
                         main = paste("Module-trait relationships"))
          dev.off()
        }}}
}

classify.tfs = function(tfs.modules, col.data, color){
  col.data = col.data %>%
    mutate(Module = tfs.modules[,paste0("ME", color)],
           Module_level = ifelse(Module > summary(Module)[5], 'High',
                                 ifelse(Module < summary(Module)[2], "Low", "na")))
  return(col.data)
}

diff_analysis = function(RNA.counts.normalized, feature, file.name){
  
  ####TFs differential active 
  net = get_collectri(organism = 'human', split_complexes = F) #Get universe  
  
  collectri2viper_regulons <- function(df) {
    regulon_list <- split(df, df$source)
    regulons <- lapply(regulon_list, function(regulon) {
      tfmode <- stats::setNames(regulon$mor, regulon$target)
      list(tfmode = tfmode, likelihood = rep(1, length(tfmode)))
    })
    
    return(regulons)
  }
  
  net_regulons = collectri2viper_regulons(net)
  #net_regulons <- pruneRegulon(net_regulons, 50, adaptive = FALSE, eliminate = TRUE) #Pruned regulons to a maximum of 50 targets to avoid statistic bias
  
  # Generating test and ref data
  test_i <- which(feature == 1)
  ref_i <- which(feature == 2)
  
  mat_test <- as.matrix(RNA.counts.normalized[,test_i])
  mat_ref <- as.matrix(RNA.counts.normalized[,ref_i])
  
  # Generating NULL model (test, reference)
  dnull <- ttestNull(mat_test, mat_ref, per=1000, verbose = F)
  
  # Generating signature
  signature <- rowTtest(mat_test, mat_ref)
  signature <- (qnorm(signature$p.value/2, lower.tail = FALSE) * sign(signature$statistic))
  signature <- na.omit(signature)
  signature <- signature[,1]
  
  # Running msVIPER
  mra <- msviper(signature, net_regulons, dnull, verbose = F, minsize = 5)
  
  # Plot DiffActive TFs
  pdf(paste0("Results/Differential_TFs_", file.name))
  print(plot(mra, mrs=15, cex=1, include = c("expression","activity")))
  dev.off()
  
  tfs_da = mra$es$p.value
  tfs_da = tfs_da[tfs_da < 0.05]
  tfs_names = names(tfs_da)
  
  nes = mra$es$nes
  nes = nes[tfs_names] #only significant
  up = names(nes[nes > 0])
  down = names(nes[nes<0])
  
  #Map TFs targets
  targets_up = net$target[which(net$source%in%up)]
  targets_down = net$target[which(net$source%in%down)]
  common = intersect(targets_up, targets_down)
  if(length(common)!=0){
    targets_up = targets_up[-which(targets_up %in% common)]
    targets_down = targets_down[-which(targets_down %in% common)]    
  }
  
  #Enrichment upregulated TFs
  entrz <- AnnotationDbi::select(org.Hs.eg.db, keys = targets_up, columns = "ENTREZID", keytype = "SYMBOL") #Change to EntrezID
  universe <- AnnotationDbi::select(org.Hs.eg.db, keys = rownames(RNA.counts.normalized), columns = "ENTREZID", keytype = "SYMBOL") #Change to EntrezID
  
  reac <- enrichPathway(gene    = entrz$ENTREZID,
                        organism     = 'human',
                        universe = universe$ENTREZID,
                        pvalueCutoff = 0.05)
  
  kegg <- enrichKEGG(gene    = entrz$ENTREZID,
                     organism     = 'human',
                     universe = universe$ENTREZID,
                     pvalueCutoff = 0.05)
  
  if(nrow(data.frame(kegg))!=0){
    pdf(paste0("Results/ORA_KEGG_UP_TFs_", file.name))
    print(dotplot(kegg))
    dev.off()
  }
  
  if(nrow(data.frame(reac))!=0){
    pdf(paste0("Results/ORA_Reactome_UP_TFs_", file.name))
    print(dotplot(reac))
    dev.off()
  }
  
  return(up)
  
}

anova_test = function(data, trait_name, y_name, pval = 0.05, file_name){
  print("Performing one way ANOVA test...................................")
  data$y = data[,y_name]
  data$trait = as.factor(data[,trait_name])
  model  <- lm(y ~ trait, data = data)
  res.aov <- data %>% rstatix::anova_test(y ~ trait)
  if(round(res.aov$p, 5) <= pval){
    print("Significant test found!")
    pdf(paste0("Results/ANOVA_", file_name))
    print(ggplot(data, aes(x=trait, y=y, fill=trait)) +
            geom_violin(width=0.6) +
            geom_boxplot(width=0.07, color="black", alpha=0.2) +
            scale_fill_brewer() +
            geom_smooth(aes(x=trait, y=y), method = "loess") +
            ylab(paste0("Values for ", y_name)) +
            xlab(paste0("Clinical trait: ", trait_name)) +
            labs(title="One way ANOVA test",
                 subtitle=rstatix::get_test_label(res.aov, detailed = TRUE)) +
            theme(axis.text.x = element_text(angle = 0),
                  axis.title.y = element_text(size = 8, angle = 90)))
    dev.off()
  }
}

minMax <- function(x) {
  #columns: features
  x = data.matrix(x)
  for(i in 1:ncol(x)){
    x[,i] = (x[,i] - min(x[,i], na.rm = T)) / (max(x[,i], na.rm = T) - min(x[,i], na.rm = T))    
  }
  
  return(x)
  
}

check_normal_distribution = function(data_df){
  data_df = data.frame(data_df)
  # Visualize the data
  plot_histograms <- function(column, colname) {
    if(is.numeric(column)) {
      ggplot(data_df, aes_string(x = colname)) +
        geom_histogram(aes(y = ..density..), bins = 30, fill = "blue", alpha = 0.5) +
        geom_density(color = "red") +
        ggtitle(paste("Histogram and Density Plot for", colname))
    }
  }
  
  # Plot histograms 
  for (colname in names(data_df)) {
    print(plot_histograms(data_df[[colname]], colname))
  }
}

is_scaled <- function(x) {
  mean_x <- mean(x)
  sd_x <- sd(x)
  return(abs(mean_x) < 1e-10 && abs(sd_x - 1) < 1e-10)
}

get_best_mfrow <- function(n_plots) {
  rows <- floor(sqrt(n_plots))
  cols <- ceiling(n_plots / rows)
  return(c(rows, cols))
}

identify.cell.groups = function(features, tfs.modules.groups, cor_type = "p", clustering.method = "ward.D2", distance.method = "manhattan", width = 12, height = 18){

  moduleTraitCor = features[[1]]

  names(features[[2]]) = paste0("ME", names(features[[2]])) #To match names of columns from corr matrix
  
  #Gather significant features across TF modules clusters
  features_vec = list()
  for (i in 1:length(tfs.modules.groups)) {
    features_vec[[i]] = unique(unlist(unname(features[[2]][tfs.modules.groups[[i]]])))
    names(features_vec)[i] = names(tfs.modules.groups)[i]
  }
  
  lis.dendrogram = list()
  
  for (i in 1:length(features_vec)){
    TFmoduleTraitcor = moduleTraitCor[,colnames(moduleTraitCor)%in%features_vec[[i]]]
    data_scaled = scale(t(TFmoduleTraitcor))
    ###Dendogram by Module
    d <- dist(data_scaled[,colnames(data_scaled)%in%tfs.modules.groups[[i]]], method = distance.method)
    dendrogram <- hclust(d, method = clustering.method)
    pdf(paste0("Results/Dendogram_cell_types_", names(features_vec)[[i]]), width = width, height = height)
    par(mar = c(5, 2, 4, 35)) #bottom, left, top, right
    plot(as.dendrogram(dendrogram), horiz= T)
    dev.off()
    lis.dendrogram[[i]] = dendrogram
  }
  
  names(lis.dendrogram) = names(features_vec)

  return(lis.dendrogram)
  
}

classify.deconvolution = function(coldata, deconvolution, group){
  deconv = deconvolution[,colnames(deconvolution)%in%group]
  
  #Patients high in group 1 of cells
  vec = c()
  if(is.null(ncol(deconv))==T){
    idx = which(deconv > median(deconv))
    vec = c(vec,idx)
  }else{
    for (i in 1:ncol(deconv)) {
      idx = which(deconv[,i] > median(deconv[,i]))
      vec = c(vec, idx)
    }
  }
  
  #High in all deconv features from group
  pos = which(table(vec) == length(group))
  
  coldata = coldata %>%
    mutate(Cells_level = "Low")
  
  coldata$Cells_level[pos] = "High"
  coldata$Cells_level = factor(coldata$Cells_level)

  return(coldata)
  
}

compute.fisher.test = function(coldata, trait){

  contingency = table(coldata[,"Cells_level"], coldata[,trait])
  test <- fisher.test(contingency)
  
  df = data.frame("Cells_level" = coldata[,"Cells_level"],
                  "Trait" = coldata[,trait])
  
  pdf(paste0("Results/Fisher test for ", trait), width = 12, height = 9)
  ggbarstats(
    df, Cells_level, Trait,
    results.subtitle = FALSE,
    title = paste0("Analysis for ", trait),
    subtitle = paste0(
      "Fisher's exact test", ", p-value = ",
      ifelse(test$p.value < 0.001, "< 0.001", round(test$p.value, 3))
    )
  )
}

remove_equal = function(lis, sort = F){
  #Sorted list to avoid no recognizing vectors with equal composition but different order of cells
  if(sort == T){
    for(i in 1:length(lis)){
      for (j in 1:length(lis[[i]])) {
        lis[[i]][[j]] = sort(lis[[i]][[j]])
      }
    }      
  }
  #Remove cell groups
  for(i in 1:length(lis)){
    exist = c() #Initialize vector of equalities
    rang = seq(1, length(lis))[-i] #Create sequence to iterate all list elements except the one being analyzed
    for (j in rang){
      idx = which(lis[[i]] %in% lis[[j]] == TRUE) #Map all cell groups which already existed
      exist = c(exist, idx) #Save index cluster
    }
    if(length(exist)!=0){
      lis[[i]] = lis[[i]][-unique(exist)] #Remove cell groups that already exist
    }
  }
  
  #Remove dendrograms without elements (length equal 0)
  vec = c()
  for(i in 1:length(lis)){
    if(length(lis[[i]]) == 0){
      vec = c(vec, i)
    }
  }
  
  if(length(vec)>0){
    lis = lis[-vec]
  }
  
  return(lis)
}

compute.TF.network.classification = function(tf.network, pathways.features){
  require(factoextra)
  
  tf.network[[1]] = data.frame(tf.network[[1]])
  pathways.features = data.frame(pathways.features)
  
  moduleTraitCor = cor(tf.network[[1]], pathways.features, method = "p")
  
  ### Find clusters
  silhouette = fviz_nbclust(moduleTraitCor, hcut, method = "silhouette", k.max = nrow(moduleTraitCor)-1)
  k_cluster = as.numeric(silhouette$data$clusters[which.max(silhouette$data$y)])
  pdf(paste0("Results/TFs_modules_Silhouette_scores"))
  print(silhouette)
  dev.off()
  
  hc_modules = hclust(dist(moduleTraitCor), method = "ward.D2")  
  dend_pathways = as.dendrogram(hc_modules)
  
  pdf(paste0("Results/TFs_modules_clusters"))
  plot(dend_pathways, cex = 0.6)
  rect.hclust(hc_modules, k = k_cluster, border = 2:5)
  dev.off()
  
  ### Extract clusters
  sub_grp <- cutree(hc_modules, k = k_cluster)
  
  ### Plot PCA and biplot
  p = fviz_cluster(list(data = moduleTraitCor, cluster = sub_grp))
  pdf(paste0("Results/PCA_TFs_modules_clusters"))
  print(p)
  dev.off()
  
  res.pca <- prcomp(moduleTraitCor,  scale = F)
  p = fviz_pca_biplot(res.pca, label="all", select.var = list(contrib = 6), addEllipses=TRUE, ellipse.level=0.75)
  pdf(paste0("Results/PCA_Biplot_TFs_modules_clusters"))
  print(p)
  dev.off()
  
  # Extract the loadings
  loadings <- res.pca$rotation
  contribution <- (loadings^2)*100
  features = contribution %>%
    data.frame() %>%
    rownames_to_column("Features") %>%
    arrange(desc(PC1))
  
  pdf("Results/Pathways_contribution_pca.pdf", width = 12, height = 8)
  p = ggplot(features, aes(x = reorder(Features, -PC1), y = PC1)) +
    geom_bar(stat = "identity", fill = "skyblue") +
    labs(title = "Contribution of pathways",
         x = "Feature",
         y = "Contribution (%)") +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 90, hjust = 1, size=15))
  print(p)
  dev.off()

  groups = list()
  for (i in 1:length(unique(sub_grp))) {
    groups[[i]] = names(sub_grp)[sub_grp == i]
    names(groups)[i] = paste0(gsub("ME", "", groups[[i]]), collapse = "_")
  }
  
  return(groups)
}

cell.groups.analysis = function(deconvolution, tfs.module.network, coldata = NULL, trait = NULL, cell.dendrograms, cut.height, stat.test, pval = 0.05, prop_var = 0.7, width = 12, height = 14){
  tfs.module.matrix = tfs.module.network[[1]] 
  module_colors = unique(tfs.module.network[[2]])
  if(is.null(trait) == F){
    coldata[,trait] = as.factor(coldata[,trait]) 
  }
  
  for (i in 1:length(cell.dendrograms)) {
    pdf(paste0("Results/Dendrogram_cell_types_", names(cell.dendrograms)[i], "_cut_", cut.height), width = 12, height = 20)
    par(mar = c(5, 2, 4, 35)) #bottom, left, top, right
    plot(as.dendrogram(cell.dendrograms[[i]]), horiz= T)
    dendextend::rect.dendrogram(as.dendrogram(cell.dendrograms[[i]]), h=cut.height, horiz=TRUE)
    dev.off()      
  }
  
  ###Identify cell groups
  cell.groups.dendrograms = list()
  cell.groups = list()
  exists = FALSE
  contador = 1
  for (i in 1:length(cell.dendrograms)) {
    x = list() #Iterator for cell names
    y = list() #Iterator for median of cell groups
    clusters <- dendextend::cutree(cell.dendrograms[[i]], h = cut.height, order_clusters_as_data=FALSE)
    for (j in 1:length(table(clusters))) {
      cells = names(clusters)[clusters==j]
      y[[j]] = cells #Filling vector with cell names
      ###Compute score for each cell group
      if(length(cells)>1){
        pca_group = deconvolution[,colnames(deconvolution) %in% cells, drop = F]
        color = extract_colors(module_colors, names(cell.dendrograms)[i])
        x[[j]] <- compute_composite_score(pca_group, color, tfs.module.matrix)#Give a normal distribution
      }else{
        x[[j]] = deconvolution[,colnames(deconvolution) %in% cells]
      }
      #######################################################
      names(x)[j] = paste0("group_", j)
      names(y)[j] = paste0("group_", j)
    }
    cell.groups[[contador]] = y
    cell.groups.dendrograms[[contador]] = data.frame(x)
    names(cell.groups.dendrograms)[contador] = names(cell.dendrograms)[contador]
    names(cell.groups)[contador] = names(cell.dendrograms)[contador]
    #cell.groups[[contador]] <- sapply(cell.groups[[contador]], function(x) paste0(x, "_" , names(cell.dendrograms)[i])) #Adding name of TF module to each deconvolution feature
    contador = contador + 1
  }
  
  #Remove groups with equal composition across TFs dendrograms
  # cell.groups = remove_equal(cell.groups, sort = T) #cell groups contain information about the composition thus it is only necessary to map it here and not in cell.groups.dendrograms
  # for (i in 1:length(cell.groups.dendrograms)) {
  #   cell.groups.dendrograms[[i]] = cell.groups.dendrograms[[i]][names(cell.groups.dendrograms[[i]]) %in% names(cell.groups[[i]])]
  # }

  #Remove groups composed of 1 feature
  for (i in 1:length(cell.groups)) {
    vec = c()
    for (j in 1:length(cell.groups[[i]])) {
      if(length(cell.groups[[i]][[j]])==1){
        vec = c(vec, j)
      }
    }
    if(length(vec)>0){
      cell.groups[[i]] = cell.groups[[i]][-vec] 
      cell.groups.dendrograms[[i]] = cell.groups.dendrograms[[i]][-vec] 
    }
  }
  
  cell.groups.dendrograms_all = list()
  cell.groups_all = list()
  contador = 1
  ##Naming cell groups 
  for (i in 1:length(cell.groups.dendrograms)) {
    for (j in 1:ncol(cell.groups.dendrograms[[i]])) {
      cell.groups.dendrograms_all[[contador]] = cell.groups.dendrograms[[i]][[j]]
      names(cell.groups.dendrograms_all)[contador] = paste0("Dendrogram_", names(cell.groups.dendrograms)[[i]], ".", names(cell.groups.dendrograms[[i]])[[j]])
      
      cell.groups_all[[contador]] = cell.groups[[i]][[j]]
      names(cell.groups_all)[contador] = paste0("Dendrogram_", names(cell.groups)[[i]], ".", names(cell.groups[[i]])[[j]])
      
      contador = contador + 1
    }
  }
  
  cell.groups.dendrograms_all = data.frame(cell.groups.dendrograms_all)
  rownames(cell.groups.dendrograms_all) = rownames(deconvolution)
  
  #Remove high correlating cell groups between each TF module
  groups = list(cell.groups.dendrograms_all, cell.groups_all)
  output = remove.cell.groups.corr(groups, unique(tfs.module.network[[2]]), threshold = 0.9)
  
  cell.groups.dendrograms_all = output[[1]]
  cell.groups_all = output[[2]]
  
  ##Export clusters in a table and save
  clusters = data.frame(matrix(nrow = length(cell.groups_all), ncol = 2))
  for (i in 1:length(cell.groups_all)) {
    clusters[i,1] = names(cell.groups_all)[[i]]
    clusters[i,2] = paste(cell.groups_all[[i]], collapse ="\n")
  }
  colnames(clusters) = c("Cell groups", "Methods-signatures")
  
  write.csv(clusters, paste0("Results/Cell.groups_all_", cut.height, ".csv"), row.names = F)
  write.csv(cell.groups.dendrograms_all, paste0("Results/Cell.groups.values_all_", cut.height, ".csv"))
  
  ######Statistical test (if asked)
  if(is.null(trait)==F){
    contador = 1
    cell.groups.dendrograms_sig = list()
    cell.groups_sig = list()
    for (i in 1:length(cell.groups.dendrograms)) {
      for (j in 1:ncol(cell.groups.dendrograms[[i]])) {
        if(stat.test == "Fisher"){
          ##Compute Fisher test
          coldata = coldata %>%
            mutate(level = cell.groups.dendrograms[[i]][,j],
                   Cells_level = ifelse(level > summary(level)[3], 'High', 'Low'))
          contingency = table(coldata[,"Cells_level"], coldata[,trait])
          test = fisher.test(contingency)
          ##Extract only significant features
          if(round(test$p.value, 5) <= pval){
            print(paste0("Significant pval after doing Fisher test for dendrogram ", names(cell.groups.dendrograms)[i], " - Cell ", names(cell.groups.dendrograms[[i]])[[j]]))
            df = data.frame("Cells_level" = coldata[,"Cells_level"], "Trait" = coldata[,trait])
            pdf(paste0("Results/Fisher_", trait, "_Dendrogram_", names(cell.groups.dendrograms)[[i]], "_", names(cell.groups.dendrograms[[i]])[[j]]), width = 12, height = 9)
            print(ggbarstats(df, Cells_level, Trait, results.subtitle = F, 
                             title= paste0("Dendrogram_", names(cell.groups.dendrograms)[[i]], " - ", names(cell.groups.dendrograms[[i]])[[j]]), 
                             subtitle = paste0("Fisher's exact test, p-value = ", ifelse(test$p.value < 0.001, "< 0.001", round(test$p.value, 5))))+
                    ggplot2::theme(plot.title = ggplot2::element_text(size=15), axis.text = ggplot2::element_text(size=14), legend.title = ggplot2::element_text(size=14)))
            dev.off()
            
            cell.groups.dendrograms_sig[[contador]] = cell.groups.dendrograms[[i]][[j]]
            names(cell.groups.dendrograms_sig)[contador] = paste0("Dendrogram_", names(cell.groups.dendrograms)[[i]], ".", names(cell.groups.dendrograms[[i]])[[j]])
            
            cell.groups_sig[[contador]] = cell.groups[[i]][[j]]
            names(cell.groups_sig)[contador] = paste0("Dendrogram_", names(cell.groups)[[i]], ".", names(cell.groups[[i]])[[j]])
            contador = contador + 1
          }
        }
        if(stat.test == "Anova"){
          ##Compute Anova test
          data = data.frame("Value" = cell.groups.dendrograms[[i]][,j], "Trait" = coldata[,trait])
          model  <- lm(Value ~ Trait, data = data)
          res.aov <- data %>% rstatix::anova_test(Value ~ Trait)
          ##Extract only significant features
          if(round(res.aov$p, 5) <= pval){
            print(paste0("Significant pval after doing Anova test for dendrogram ", names(cell.groups.dendrograms)[i], " - Cell ", names(cell.groups.dendrograms[[i]])[[j]]))
            pdf(paste0("Results/Anova_", trait, "_Dendrogram_", names(cell.groups.dendrograms)[[i]], "_", names(cell.groups.dendrograms[[i]])[[j]]), width = 12, height = 9)
            print(ggplot(data, aes(x=Trait, y=Value, fill=Trait)) +
                    geom_violin(width=0.6) +
                    geom_boxplot(width=0.07, color="black", alpha=0.2) +
                    scale_fill_brewer() +
                    geom_smooth(aes(x=Trait, y=Value), method = "loess") +
                    xlab(paste0("Clinical trait: ", trait)) +
                    labs(title= paste0("Dendrogram_", names(cell.groups.dendrograms)[[i]], " - ", names(cell.groups.dendrograms[[i]])[[j]]), 
                         subtitle = rstatix::get_test_label(res.aov, detailed = TRUE)) +
                    theme(axis.text.x = element_text(angle = 0),
                          axis.title.y = element_text(size = 8, angle = 90)))
            dev.off()
            
            cell.groups.dendrograms_sig[[contador]] = cell.groups.dendrograms[[i]][[j]]
            names(cell.groups.dendrograms_sig)[contador] = paste0("Dendrogram_", names(cell.groups.dendrograms)[[i]], ".", names(cell.groups.dendrograms[[i]])[[j]])
            
            cell.groups_sig[[contador]] = cell.groups[[i]][[j]]
            names(cell.groups_sig)[contador] = paste0("Dendrogram_", names(cell.groups)[[i]], ".", names(cell.groups[[i]])[[j]])
            contador = contador + 1
          }
        }
      }
    }
    
    if(length(cell.groups_sig)!=0){
      ##Export clusters in a table
      clusters = data.frame(matrix(nrow = length(cell.groups_sig), ncol = 2))
      for (i in 1:length(cell.groups_sig)) {
        clusters[i,1] = names(cell.groups_sig)[[i]]
        clusters[i,2] = paste(cell.groups_sig[[i]], collapse ="\n")
      }
      colnames(clusters) = c("Cell groups", "Methods-signatures")
      
      write.csv(clusters, paste0("Results/Cell.groups_", stat.test, "_" , trait, "_", cut.height, ".csv"), row.names = F)
      
      ###Convert to data frame
      cell.groups.dendrograms_sig = data.frame(cell.groups.dendrograms_sig)
      rownames(cell.groups.dendrograms_sig) = rownames(coldata)
      
      write.csv(cell.groups.dendrograms_sig, paste0("Results/Cell.groups.values_", stat.test, "_" , trait, "_", cut.height, ".csv"))
      
      return(list(cell.groups.dendrograms_sig, cell.groups_sig))
    }else{
      print("No significant cell groups")
    }
  }else{
    return(list(cell.groups.dendrograms_all, cell.groups_all))
  }
  

}

plot.modules.categorical = function(tfs.modules, coldata){
  data = cbind(tfs.modules, coldata)
  for(i in 1:ncol(tfs.modules)){
    for (j in (ncol(tfs.modules)+1):ncol(data)) {
      module <- names(data[i])
      trait <- names(data[j])
      avz <- broom::tidy(aov(data[,i] ~ data[,j], data = data))
      if(avz$p.value[1] < 0.05) {
        pdf(paste0("Results/ANOVA_", module, "-", trait))
        print(ggplot(data, aes(x=data[,j], y=data[,i], fill=data[,j])) +
                geom_violin(width=0.6) +
                geom_boxplot(width=0.07, color="black", alpha=0.2) +
                scale_fill_brewer() +
                geom_smooth(aes(x=data[,j], y=data[,i]), method = "loess") +
                ylab(paste0("Values for ", module)) +
                xlab(paste0("Clinical trait: ", trait)) +
                labs(title="One way ANOVA test",
                     subtitle=paste0("pvalue: ", avz$p.value[1])) +
                theme(axis.text.x = element_text(angle = 0),
                      axis.title.y = element_text(size = 8, angle = 90))+
                scale_fill_discrete(name = trait))
        dev.off()
      }
    }
  }
}

compute.metada.association = function(tfs.modules, coldata, pval = 0.05, width = 20, height = 8){
  ###Association with categorical variables
  coldata_categorical = coldata %>%
    dplyr::select(where(is.character)|where(is.factor))
  
  if(ncol(coldata_categorical)!=0){
    data = cbind(tfs.modules, coldata_categorical)
    pvals = data.frame()
    fvals = data.frame()
    for(i in 1:ncol(tfs.modules)){
      contador = 1
      for (j in (ncol(tfs.modules)+1):ncol(data)) {
        module <- names(data[i])
        trait <- names(data[j])
        avz <- broom::tidy(aov(data[,i] ~ data[,j], data = data))
        pvals[i,contador] = avz$p.value[1]
        fvals[i,contador] = avz$statistic[1]
        contador = contador + 1
        if(avz$p.value[1] < pval) {
          pdf(paste0("Results/ANOVA_", module, "-", trait))
          print(ggplot(data, aes(x=data[,j], y=data[,i], fill=data[,j])) +
                  geom_violin(width=0.6) +
                  geom_boxplot(width=0.07, color="black", alpha=0.2) +
                  scale_fill_brewer() +
                  geom_smooth(aes(x=data[,j], y=data[,i]), method = "loess") +
                  ylab(paste0("Values for ", module)) +
                  xlab(paste0("Clinical trait: ", trait)) +
                  labs(title="One way ANOVA test",
                       subtitle=paste0("F statistic: ", round(avz$statistic[1],3), "\npvalue: ", round(avz$p.value[1], 3))) +
                  theme(axis.text.x = element_text(angle = 0),
                        axis.title.y = element_text(size = 8, angle = 90))+
                  scale_fill_discrete(name = trait))
          dev.off()
        }
      }
    }
    rownames(pvals) = colnames(tfs.modules)
    colnames(pvals) = colnames(coldata_categorical)
    
    rownames(fvals) = colnames(tfs.modules)
    colnames(fvals) = colnames(coldata_categorical)
    
    pvals = as.matrix(pvals)
    textMatrix2 = paste("ANOVA\n(", signif(pvals, 2), ")", sep = "")
    dim(textMatrix2) = dim(pvals)
  }

  ###Association with quantitative variables
  coldata_quantitative = coldata %>%
    dplyr::select(where(is.numeric))
  
  if(ncol(coldata_quantitative)!=0){
    moduleTraitCor = cor(tfs.modules, coldata_quantitative, method = "p");
    moduleTraitPvalue = corPvalueStudent(moduleTraitCor, nrow(tfs.modules))
    
    textMatrix = paste(signif(moduleTraitCor, 2), "\n(", signif(moduleTraitPvalue, 2), ")", sep = "")
    dim(textMatrix) = dim(moduleTraitCor)
    
    if(ncol(coldata_categorical)!=0){
      textMatrix = cbind(textMatrix, textMatrix2)
      moduleTraitPvalue = cbind(moduleTraitPvalue, pvals)
      simulated_corr = matrix(runif(n=nrow(pvals)*ncol(pvals), min=-0.1, max=0.1), nrow = nrow(pvals), ncol = ncol(pvals))
      colnames(simulated_corr) = colnames(pvals)
      moduleTraitCor = cbind(moduleTraitCor, simulated_corr)
    }
    
    idx = which(round(moduleTraitPvalue,2)>pval)
    for (i in idx) {
      textMatrix[i] = NA
    }
    
    pdf("Results/TF.modules_metadata", width = width, height = height)
    par(mar = c(25, 15, 3, 3))
    labeledHeatmap(Matrix = moduleTraitCor,
                   xLabels = colnames(moduleTraitCor),
                   yLabels = rownames(moduleTraitCor),
                   ySymbols = rownames(moduleTraitCor),
                   colorLabels = FALSE,
                   colors = blueWhiteRed(50),
                   textMatrix = textMatrix,
                   setStdMargins = FALSE,
                   cex.text = 0.5,
                   zlim = c(-1,1),
                   main = paste0("Clinical associations with TFs modules\nOnly showing significant associations (pvalue < ", pval, ")"))
    dev.off()
  }

}

feature.selection.boruta <- function(train_data, iterations = 50, fix = TRUE, thres = 0.5, file_name, return = T) {
  library(Boruta)
  library(dplyr)
  results_list <- list()
  imps_list <- list()
  
  cat('\nComputing Boruta algorithm with bootstraps...................................\n\n')
  
  if(fix){
    cat('\nApplying Rough Fix in Boruta...................................\n\n')
  }
  
  for (i in 1:iterations) {
    set.seed(sample.int(100000, 1))
    boruta_output <- Boruta(target ~ ., data = train_data, doTrace = 0)
    
    if (fix) {
      roughFixMod <- TentativeRoughFix(boruta_output)
      boruta_output <- roughFixMod
    }
    
    imps <- attStats(boruta_output)
    results_list[[i]] <- as.character(imps$decision)
    imps_list[[i]] <- imps %>%
      data.frame() %>%
      rownames_to_column("Variable") %>%
      dplyr::select(-decision)
  }
  
  boruta = merge_boruta_results(importance_values = imps_list, decisions = results_list, file_name, iterations = iterations, threshold = thres, return_plot = return)
  
  return(boruta)
}

merge_boruta_results = function(importance_values, decisions, file_name, iterations, threshold = 0.8, return_plot = T){
  
  ### Construct matrix of importance
  combined_importance <- do.call(rbind, importance_values)
  combined_results_long <- combined_importance %>% #Matrix for plotting
    pivot_longer(cols = meanImp, names_to = "Measure", values_to = "Value")
  
  median_df <- combined_importance %>% #Calculate the median for each column, grouped by the variable name
    group_by(Variable) %>%
    dplyr::summarize(across(everything(), median, na.rm = TRUE))
  
  ### Retrieve important and tentatives variables
  combined_results <- do.call(cbind, decisions)
  rownames(combined_results) = median_df$Variable
  decisions_summary <- apply(combined_results, 1, function(x) {
    table(factor(x, levels = c("Confirmed", "Tentative", "Rejected")))
  })
  confirmed_vars <- names(which(decisions_summary["Confirmed",] >= round(threshold*iterations))) 
  tentative_vars <- names(which(decisions_summary["Tentative",] >= round(threshold*iterations))) 
  
  # For plotting
  combined_results_long$Decision = "Rejected"
  combined_results_long$Decision[which(combined_results_long$Variable %in% confirmed_vars)] = "Confirmed"
  combined_results_long$Decision[which(combined_results_long$Variable %in% tentative_vars)] = "Tentative"
  
  mean_order <- median_df %>% #Extract the order of variables for plotting
    arrange(meanImp) %>%
    pull(Variable)
  
  # For result 
  median_df$Decision = "Rejected"
  median_df$Decision[which(median_df$Variable %in% confirmed_vars)] = "Confirmed"
  median_df$Decision[which(median_df$Variable %in% tentative_vars)] = "Tentative"
  
  if(return_plot == T){
    # Plot variable importance boxplots
    pdf(paste0("Results/Boruta_variable_importance_", file_name, ".pdf"), width = 8, height = 12)
    print(ggplot(combined_results_long, aes(x = factor(Variable, levels = mean_order), y = Value, fill = Decision)) +
            geom_bar(stat = "identity", position = "dodge") +
            theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
            coord_flip() +
            labs(x = "Features", y = "Importance", title = paste0("Variable Importance by Boruta after ", iterations, " bootstraps\n", file_name)) +
            scale_fill_manual(values = c("Confirmed" = "green", "Tentative" = "yellow", "Rejected" = "red")) +
            facet_wrap(~ Measure, scales = "free_y"))
    dev.off()
  }
  
  cat('\nCell groups selected across different training sets:\n\n')
  print(median_df[which(median_df$Decision %in% c("Confirmed", "Tentative")),])
  
  return(list(Confirmed = confirmed_vars, Tentative = tentative_vars, Matrix_Importance = median_df))
}

compute.survival.analysis = function(features, survival.data, time_unit, p.value = 0.05, thres = 0.5, max_factors = Inf) {
  n_features <- ncol(features)
  significant_combinations <- list() # To store significant feature combinations
  
  # Generate all possible combinations of the features
  contador = 1
  for (n in 1:min(n_features, max_factors)) {
    combinations <- combn(1:n_features, n, simplify = FALSE)
    
    for (comb in combinations) {
      # Create a formula dynamically based on the combination
      formula <- as.formula(paste("Surv(time, status) ~", paste(colnames(features)[comb], collapse = " + ")))
      
      # Prepare the data for survival analysis
      data_for_model <- data.frame("time" = survival.data$PFS,
                                   "status" = survival.data$DRP_st) 
      
      data_for_model = cbind(data_for_model, features[,comb, drop=F])
      
      # Fit the Cox PH model with the combination of features (cox PH take into account covariates and measure the impact of each variable in the survival time)
      cox <- cph(formula, data = data_for_model)
      data_for_model$CoxPredictors <- cox$linear.predictors #linear predictors is the risk score for each individual in the dataset
      
      # Check that the model is significant as a predictor (maybe not useful?, it gives the same linear.predictos - to be check)
      cphmodel <- coxph(Surv(time, status) ~ CoxPredictors, data = data_for_model)
      data_for_model$CoxPredictors <- cphmodel$linear.predictors
      
      quantiles <- quantile(data_for_model$CoxPredictors, thres) 
      
      # Binarize the Cox model output to draw two KM lines (linear predictors are used to stratify between high-risk and low-risk groups)
      data_for_model$coxHL <- ifelse(cphmodel$linear.predictors >= quantiles, 'High', "Low") 
      
      # Perform Kaplan-Meier based on coxHL
      km_fit <- survfit(Surv(time, status) ~ coxHL, data = data_for_model)
      
      pval <- surv_pvalue(km_fit, data = data_for_model)$pval #Performs log-rank test to see whether both survival curves are significantly different
      
      if (!is.na(pval) && pval < p.value) {
        significant_combinations[[contador]] <- formula
        names(significant_combinations)[contador] = paste0("Formula_", contador)
        
        pdf(paste0("Results/SurvPlot_", names(significant_combinations)[contador]), width = 10, height = 5, onefile = FALSE)
        print(ggsurvplot(km_fit,
                         data = data_for_model,
                         size = 1,
                         palette = c("#E7B800", "#2E9FDF"),
                         conf.int.style = "step",
                         pval = TRUE,
                         risk.table = TRUE,
                         risk.table.col = "strata",
                         legend.labs = c("High", "Low"),
                         risk.table.height = 0.3,
                         ggtheme = theme_grey(),
                         title = paste0("Cox PH for ", names(significant_combinations)[contador]),
                         xlab = paste0("Time to death/recurrence/progression (", time_unit, ")")
        ))
        dev.off()
        contador = contador + 1
      }
    }
  }
  
  if (length(significant_combinations) == 0) {
    print("No significant combinations found.")
  } else {
    return(significant_combinations)
  }
  
}

extract_cells = function(groups){
  names_cells = c("B.cells", "B.naive", "B.memory", "Macrophages.cells", "Macrophages.M0", "Macrophages.M1", "Macrophages.M2", "Monocytes", "Neutrophils", "NK.cells", "NK.activated", 
                  "NK.resting", "NKT.cells", "CD4.cells", "CD4.memory.activated", "CD4.memory.resting", "CD4.naive", "CD8.cells", "T.cells.regulatory", "T.cells.non.regulatory","T.cells.helper", 
                  "T.cells.gamma.delta", "Dendritic.cells", "Dendritic.activated", "Dendritic.resting", "Cancer", "Endothelial", "Eosinophils", "Plasma.cells", "Myocytes", "Fibroblasts",
                  "Mast.cells", "Mast.activated", "Mast.resting", "CAF")
  
  regex_pattern <-  paste0("(", paste(names_cells, collapse = "|"), ")")
  
  extracted_names <- sapply(groups, function(x) {
    match <- regexpr(regex_pattern, x)
    if (match != -1) {
      return(regmatches(x, match))
    } else {
      return(NA)
    }
  })
  
  extracted_names <- unname(extracted_names)
  extracted_names <- unique(na.omit(extracted_names))
  return(extracted_names)
}

extract_colors <- function(module_colors, cell_group_name) {
  matches <- c() # For storing the matches
  for (color in module_colors) {
    match <- regexpr(color, cell_group_name) # Find the position of the match
    if (match != -1) {
      matches <- c(matches, regmatches(cell_group_name, match))
    }
  }
  
  
  if (length(matches) > 0) {
    order <- sapply(matches, function(m) regexpr(m, cell_group_name)) # Sort matches based on their position in the original string to ensure names are the same
    matches <- matches[order(order)]  # Order the matches based on their position
    return(matches)  # Return the ordered matches
  } else {
    return(NA)  # If no matches are found, return NA
  }
  
}

remove.cell.groups.corr <- function(data, colors, threshold = 0.9) {
  
  features_high_corr = c()
  # Compute correlation matrix
  corr_matrix <- cor(data[[1]])
  # Find highly correlated features
  contador = 1
  while(nrow(corr_matrix)>0){
    color_features = c()
    feature = data.frame(corr_matrix[1, , drop = FALSE]) #Extract first row feature
    feature = feature %>%                                #Take only high corr above threshold
      mutate_all(~ifelse(. > threshold, ., NA)) %>%
      select_if(~all(!is.na(.)))
    
    corr_matrix = corr_matrix[-which(rownames(corr_matrix)%in%colnames(feature)),-which(colnames(corr_matrix)%in%colnames(feature)), drop = F] #Remove already joined features
    color_features = list()
    if(ncol(feature)>1){
      for (m in 1:ncol(feature)) {
        color_group = extract_colors(colors, colnames(feature)[m])
        color_features[[m]] = color_group
      }
      
      all_equal <- all(sapply(color_features, function(x) identical(x, color_features[[1]]))) #Check whether highly correlated features belong to the same group of TFs modules
      
      if(all_equal == TRUE){
        new_group_composition = unique(unlist(unname(data[[2]][colnames(feature)])))
        new_group_value = rowMeans(data[[1]][,colnames(data[[1]])%in%colnames(feature)])
        
        print(paste0("Highly correlated features (r>", threshold,"): ", paste(colnames(feature), collapse = ', '), ". Combining."))
        
        if(contador==1){
          #Remove features from original data
          new_data <- data[[1]][, -which(colnames(data[[1]])%in%colnames(feature)), drop = F] 
          new_groups = data[[2]][-which(names(data[[2]]) %in% colnames(feature))]
        }else{
          new_data <- new_data[, -which(colnames(new_data)%in%colnames(feature)), drop = F] 
          new_groups = new_groups[-which(names(new_groups) %in% colnames(feature))]
        }
        
        #Add new combined features 
        new_name = paste0("Dendrogram_",  paste0(color_group, collapse = "_"), ".group_combined_", contador)
        new_data = cbind(new_data, new_group_value)
        colnames(new_data)[length(new_data)] = new_name
        
        new_groups[[length(new_groups)+1]] = new_group_composition
        names(new_groups)[length(new_groups)] = new_name
        
        contador = contador + 1
      }else{
        print(paste0("Highly correlated features between different clusters of TFs modules (r>", threshold,"): ", paste(colnames(feature), collapse = ', '), ". Not joining"))
      }
    }else{
      if(contador == 1){ 
        new_data = data[[1]]
        new_groups = data[[2]]
      }
    }
  }
  
  res = list(new_data, new_groups)
  
  return(res)
}

compute_composite_score = function(cell_group, color_group, tfs.module.matrix, prop_var = 0.7){
  module_group = paste0("ME", color_group) #To match with columns of TFs modules
  pca_group <- cbind(minMax(tfs.module.matrix[,module_group, drop=F]), cell_group) #Combined TF module corresponding to each cell group
  svd_result <- svd(pca_group, nu = min(nrow(pca_group), ncol(pca_group)), nv = min(nrow(pca_group), ncol(pca_group))) #Performs Singular Value Decomposition (SVD)
  singular_values <- svd_result$d #Extract singular values
  nPCs <- sum(cumsum(singular_values^2) / sum(singular_values^2) < prop_var) + 1 #Extract number of PCs necessary to cover the prop_var
  selected_components <- svd_result$u[,1]*-1 #Take eigenvalue and avoiding negatives (refers only to the direction)
  #selected_components = minMax(selected_components) #Make the PCs values positives (between 0-1)
  variance_explained <- (singular_values^2/sum(singular_values^2)) #Extract the % variance explained by each PC
  weights <- variance_explained[1:nPCs]
  #deconv_composite_score = rowSums(selected_components * weights)
  composite_score = selected_components
  #composite_score = minMax(tfs.module.matrix[,module_group]) * deconv_composite_score #Add scaled TFs module i score for considering TF contribution
  #composite_score <- scale(composite_score) #Give a normal distribution
  return(composite_score)
}

find.maximum.iteration = function(cells.groups){
  max_iteration = c()
  for (i in 1:length(cells.groups)){
    if(is.null(names(cells.groups[[i]]))==F){
      iterations <- sapply(names(cells.groups[[i]]), function(x) {
        as.numeric(sub(".*\\.Iteration\\.(\\d+)", "\\1", x))
      })
      local_max = max(iterations)
      max_iteration = c(max_iteration, local_max)
    }
  }
  
  return(max(max_iteration))
}

create_tfs_modules = function(TF.matrix, network_tfs){
  
  tfs.modules = TF.matrix %>%
    t() %>%
    data.frame() %>%
    mutate(Module = "na")
  
  for (i in 1:length(network_tfs[[3]])) {
    tfs.modules$Module[which(rownames(tfs.modules) %in% network_tfs[[3]][[i]])] = names(network_tfs[[3]])[i]
  }
  
  tfs_colors = tfs.modules %>%
    pull(Module) 
  
  MEList = moduleEigengenes(TF.matrix, colors = tfs_colors, scale = F) #Data already scale
  MEs = MEList$eigengenes
  MEs = orderMEs(MEs)
  
  return(MEs)
}

compute_cell_groups_signatures = function(deconv_res, network_res, cell_groups, features, deconv_test, tfs_test){
  
  #Remove colors indicatives from deconvolution features to be able to project them in the raw deconvolution results
  #pattern_colors <- paste0("_(", paste(unique(network_res[[2]]), collapse = "|"), ")$")
  
  # Use the regular expression to remove the suffix colors from the deconvolution features
  # for (i in 1:length(cell_groups[[2]])) {
  #   cell_groups[[2]][[i]] = gsub(pattern_colors, "", cell_groups[[2]][[i]])
  # }

  ################################################################################Simulate TFs module scores
  TF.matrix_simulated = create_tfs_modules(tfs_test, network_res)
  module_colors = unique(network_res[[2]])
  
  ################################################################################Simulate cell subgroups
  #Scale deconvolution features by columns for making them comparable between cell types (0-1). 
  for (i in 1:ncol(deconv_test)) {
    deconv_test[,i] = deconv_test[,i]/max(deconv_test[,i])
  } 
  
  deconv_subgroups <- mapply(c, deconv_res[[3]], deconv_res[[4]], SIMPLIFY = FALSE) #Join cell groups 
  iterations = find.maximum.iteration(deconv_subgroups)
  
  # Create same groups composition
  for (m in 1:iterations) {
    base_groups = list()
    for (i in 1:length(deconv_subgroups)){
      if(length(deconv_subgroups[[i]])!=0){
        idy = grep(paste0("Iteration.",m), names(deconv_subgroups[[i]]))
        if(length(idy)!=0){
          base_groups = append(base_groups, deconv_subgroups[[i]][idy])
        }
      }
    }
    
    deconv_subgroups_values = c()
    for (i in 1:length(base_groups)) {
      deconv_subgroups_values = cbind(deconv_subgroups_values, rowMedians(as.matrix(deconv_test[,base_groups[[i]]]))) #Compute median using base groups
    }
    colnames(deconv_subgroups_values) = names(base_groups)
    deconv_test = cbind(deconv_subgroups_values, deconv_test) # Join cell subgroups and deconv features
    
  }
  
  deconv_test = deconv_test[,colnames(deconv_test)%in%colnames(deconv_res[[1]])]
  
  # Compute composite scores
  idx = which(names(cell_groups[[2]]) %in% features)
  cell_dendrogram = c()
  names = c()
  for (i in 1:length(idx)) {
    pca_cells = deconv_test[,cell_groups[[2]][[idx[i]]]]
    pca_cells <- pca_cells[, apply(pca_cells, 2, function(x) var(x) != 0), drop = F] #Check if we have zero-columns
    name_cell_group = names(cell_groups[[2]][idx[i]])
    color = extract_colors(module_colors, name_cell_group)
    
    if(ncol(pca_cells) > 1){ #Check if there are more than 2 columns
      cell_dendrogram = cbind(cell_dendrogram, compute_composite_score(pca_cells, color, TF.matrix_simulated))
      names = c(names, names(cell_groups[[2]])[idx[i]])
    }
  }
  
  if(is.null(cell_dendrogram)==T){
    print("No composite scores because all features have zero variance.")
  }else{
    colnames(cell_dendrogram) = names
    rownames(cell_dendrogram) = rownames(TF.matrix_simulated)
  }
  
  return(data.frame(cell_dendrogram))
}
