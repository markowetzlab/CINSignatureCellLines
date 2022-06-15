samples <- gsub(gsub(x = list.dirs(path = "data/ascat_results",full.names=FALSE),pattern = "\\./",replacement = ""),pattern = "_70$",replacement= "")
samples <- samples[-1]

meta <- data.frame()

for(i in samples){
	if(file.exists(paste0("data/ascat_results/",i,"_70/ASCAT.output.second.Rda"))){
    	#print(i)
  		load(paste0("data/ascat_results/",i,"_70/ASCAT.output.second.Rda"))
    	if(length(ascat.output$failedarrays) == 0){
  	  data <- data.frame(sample=i,
			purity=as.numeric(ascat.output$aberrantcellfraction),
			ploidy=as.numeric(ascat.output$ploidy),
			goodnessOffit=as.numeric(ascat.output$goodnessOfFit))
          print(data)
          write.table(x=data,"data/ascat_results/ascat_profile_statistics.tsv",
		sep="\t",
		quote=F,
		col.names=F,
		row.names=F,
		append=T)
    }
   rm(ascat.output)
	}
}

