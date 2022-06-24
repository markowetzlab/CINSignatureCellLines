library(CINSignatureQuantification)

CN.data = readRDS("data/cellLine_segment_ascat_sc_fixed_purity_tCN.rds")

sigs = quantifyCNSignatures(CN.data,experimentName="cell.line",method="drews",cores=1)

write.table(as.matrix(sigs@activities$thresholdAct2),
	file="data/cellLine_signature_data.tsv",
	row.names=T,col.names=T,sep="\t")
saveRDS(sigs,file="data/cellLine_signature_data.rds")
