samples <- gsub(gsub(x = list.dirs(path = "data/ascat_results",full.names=FALSE),pattern = "\\./",replacement = ""),pattern = "_70$",replacement= "")
samples <- samples[-1]

meta <- data.frame()

for(i in samples){
  if(file.exists(paste0("data/ascat_results/",i,"_70/solution.second.Rda"))){
  	load(paste0("data/ascat_results/",i,"_70/solution.second.Rda"))
  	data <- data.frame(sample=i,purity=as.numeric(solution$purity),ploidy=as.numeric(solution$ploidy))
  	write.table(x=data,"data/ascat_results/ascat_sc_profile_statistics.tsv",sep="\t",quote=F,col.names=F,row.names=F,append=T)
   rm(solution)
	}
}
