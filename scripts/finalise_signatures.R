apply_sig_thresholds <- function(sigs=NULL,thresholds=NULL){
  sig_col_names <- colnames(sigs)
  newSigs <- do.call(cbind,lapply(1:ncol(sigs), FUN = function(x){
      sig_col <- sigs[,x]
      sig_col[sig_col <= thresholds[x]] <- 0
      return(sig_col)
    })
  )
  newNormSigs <- t(apply(newSigs,1,FUN = function(x) x / sum(x)))
  colnames(newNormSigs) <- sig_col_names
  return(newNormSigs)
}
			 
cn <- read.table(file = "data/cellLine_segment_ascat_sc_fixed_purity_tCN.tsv",header = T,sep = "\t")
signature_data <- readRDS(file = "data/signatures/4_Exposures_to_TCGA_Signatures_ASCAT_cellLines_penalty70.rds")
sig_thresholds <- read.table("resources/3_Boxplots_per_sig_fullTCGA_1000sims_10pGaussian_10pSamplePoisson.txt",
                             header = T,
                             sep = "\t",
                             row.names = 1)

new_signature_data <- apply_sig_thresholds(sigs = signature_data,thresholds = sig_thresholds$Thresh_ZeroGMM_0.01)
new_signature_data.z <- apply(new_signature_data,MARGIN = 2, FUN = function(x) (x-mean(x))/sd(x))

sig_data <- list(copy_number=cn,
	signatures=new_signature_data,
	signatures.zscore=new_signature_data.z,
	signatures.raw=signature_data)
saveRDS(object = sig_data,file = "data/cellLine_signature_data.rds")
