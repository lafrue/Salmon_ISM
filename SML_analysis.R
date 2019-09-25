
# setting the working directory
setwd("/path/to/working/folder/data_sml")

# metadata file to be imported
compo <- "compo.txt"

## for parralel
library('doMC')
registerDoMC(cores = detectCores())
library('phyloseq')
library('metagenomeSeq')
library('matrixStats')

## scripts for SML and plots
source("Rscripts/sml_compo.R")
source("Rscripts/sml_compo_hyperparameter_fitting.R")
source("Rscripts/plot_ml.R")

## the column name of the metadata file for samples ids
s <- "Samplename"

## how many combinations of parameter to try?
hyperparam <- 100

# for each ASVs table (ciliates and bacteria)
for (m in c("ciliates", "bacteria"))
{
  # setting the asv table to be imported
  # table must be tab separated, ASV as rows and samples as columns, and with no taxonomic info in the last columns.
  if (m == "ciliates") tab <- "ASVtable_CIL_noTax.tsv"
  if (m == "bacteria") tab <- "ASVtable_BAC_noTax.tsv"

  # import data 
  comp <- read.table(compo, header=TRUE, sep="\t", dec=",")
  asv <- read.table(tab, header=TRUE, sep="\t", row.names=1)
  
  # for match between ASVs table and metadata
  colnames(asv) == as.character(comp[,s])
  table(colnames(asv) %in% as.character(comp[,s]))
  
  # re-order the matrix
  asv_sort <- asv[,as.character(comp[,s])]
  asv_raw <- t(asv_sort)
  colnames(asv_sort) == as.character(comp[,s])
  
  # removing rare ASVs
  asv <- data.frame(asv_raw[,colSums(asv_raw)>10])
  
  # removing samples with not enough reads 
  asv <- subset(asv, rowSums(asv_raw) > 1000)
  comp <- subset(comp, rowSums(asv_raw) > 1000)
  
  ### CSS NORMALIZATION
  # with metagenomeSeq package
  # building a phyloseq object
  asv_to_norm  <- asv
  comp_to_norm <- comp
  samples_names <- as.character(comp[,s])
  comp_ <- comp_to_norm[,c("Farm", "Station")]
  rownames(comp_) <- samples_names
  rownames(asv_to_norm) <- samples_names
  obj <- phyloseq(otu_table(asv_to_norm, taxa_are_rows = F), sample_data(comp_))
  obj_m <- phyloseq_to_metagenomeSeq(obj)
  p = cumNormStatFast(obj_m)
  dat_mol_norm = cumNorm(obj_m, p = p)
  asv_NORM <- t(MRcounts(dat_mol_norm, norm = TRUE, log = TRUE))
  
  ### full ASV table (asv_NORM) and comp for composition file
  ASV <- as.data.frame(asv_NORM)
  COMP <- cbind(comp, col_plot = as.numeric(comp$Locality))
  
  ############################################################################################################
  ############################################################################################################
  
  ##### SML 
  ### without hyperparameter tuning
  # SVM
  preds_svm <- sml_compo(ASV, COMP, index = "AMBI", algo = "SVM", cross_val = "Locality")
  plot_ml(preds_svm, metadata = COMP, index = "AMBI", aggreg = c("Grab", "Station", "Locality"), title = paste0(m, "_SVM_default"), pdf = T, folder_export = paste0(m, "_SVM_default"))
  # RF
  preds_rf <- sml_compo(ASV, COMP, index = "AMBI", algo = "RF", cross_val = "Locality")
  plot_ml(preds_rf, metadata = COMP, index = "AMBI", aggreg = c("Grab", "Station", "Locality"), title = paste0(m, "_RF_default"), pdf = T, folder_export = paste0(m, "_RF_default"))
  
  ### with hyperparameter tuning
  ### random search for SVM
  grid_svm <- data.frame(array(NA, c(hyperparam,4)))
  dimnames(grid_svm)[[2]] <- c("type", "kernel", "epsilon", "tolerance")
  for (i in 1:hyperparam)
  {
    grid_svm[i,"type"] <- c("eps-regression", "nu-regression")[sample(1:2, 1)]
    grid_svm[i,"kernel"] <- c("linear", "polynomial", "radial", "sigmoid")[sample(1:4, 1)]
    grid_svm[i,"epsilon"] <- runif(1, 0.01, 0.2) # default is 0.1
    grid_svm[i,"tolerance"] <- runif(1, 0.0001, 0.002) # default is 0.001
  }
  
  ### random search for RF
  grid_rf <- data.frame(array(NA, c(hyperparam,3)))
  dimnames(grid_rf)[[2]] <- c("mtry", "splitrule", "min.node.size")
  for (i in 1:hyperparam)
  {
    grid_rf[i,"mtry"] <- as.numeric(c("0.33", "0.5", "0.66")[sample(1:3, 1)])
    grid_rf[i,"splitrule"] <- c("variance", "extratrees", "maxstat")[sample(1:3, 1)]
    grid_rf[i,"min.node.size"] <- as.numeric(c("3", "4", "5", "6", "7")[sample(1:5, 1)])
  }
  
  ## hypeparameter optimization and export of plots and rmse
  # RF
  preds_hyp_rf <- sml_compo_hyper_fit(ASV, COMP, index = "AMBI", algo = "RF", grid = grid_rf)
  plot_ml(preds_hyp_rf$preds, metadata = COMP, index = "AMBI", aggreg = c("Grab", "Station", "Locality"), folder_export = paste0(m,"_RF_100hypOpt"), pdf = T)
  write.table(preds_hyp_rf$RMSE, paste0(m,"_RF_100hypOpt/RMSE.tsv"), sep="\t", row.names = F, quote = F, dec=",")
  # SVM
  preds_hyp_svm <- sml_compo_hyper_fit(ASV, COMP, index = "AMBI", algo = "SVM", grid = grid_svm)
  plot_ml(preds_hyp_svm$preds, metadata = COMP, index = "AMBI", aggreg = c("Grab", "Station", "Locality"), folder_export = paste0(m,"_SVM_100hypOpt"), pdf = T)
  write.table(preds_hyp_svm$RMSE, paste0(m,"_SVM_100hypOpt/RMSE.tsv"), sep="\t", row.names = F, quote = F, dec=",")
  
  ### ASV importance
  mod <- ranger(COMP[,"AMBI"] ~ ., data=ASV, mtry=floor(dim(ASV)[2]/3), num.trees = 300, importance= "impurity", write.forest = T)
  imp <- tail(sort(mod$variable.importance), 100)
  pdf(width = 5, height=10, file = paste0("ASV_importance_",m,".pdf"))
  p <- barplot(imp, horiz = T, xlab="Variable importance", main=paste("ASVs importance for AMBI"), las=2, cex.names = 0.5)
  dev.off()

}








