PEN <- 70
ASCAT_DIR <- "data/ascat_results/"
OUTPUT_SUFFIX <- "ASCAT_ASCAT.sc.logr.fit.txt"

SAMPLE_LIST <- gsub(gsub(x = list.files(ASCAT_DIR,pattern=paste0("*_ASCAT.sc.logr.fit.txt"),recursive=T),pattern=OUTPUT_SUFFIX,replacement=""),pattern="_70\\/",replacement="")

FILE_BOOL <- file.size(paste0(ASCAT_DIR,SAMPLE_LIST,"_",PEN,"/",OUTPUT_SUFFIX)) > 1
SAMPLE_LIST_FILT <- SAMPLE_LIST[FILE_BOOL]

segment_file <- do.call(rbind,
	lapply(SAMPLE_LIST_FILT,FUN = function(x){
		file <- paste0(ASCAT_DIR,x,"_",PEN,"/",OUTPUT_SUFFIX)
		tab <- read.table(file,sep="\t",skip=1,stringsAsFactors=F)
		tab$sample <- rep(x,times=nrow(tab))
		if(ncol(tab) != 8){
			print(x)
			print(head(tab))
			
		}
		return(tab)
			}))
colnames(segment_file) <- c("chr","startpos","endpos","segVal_rounded","segVal","logr","logr.sd","sample")

#segment_file$segVal <- as.numeric(segment_file$nAraw) + as.numeric(segment_file$nBraw)
segment_file <- segment_file[,c("chr","startpos","endpos","segVal","sample")]
colnames(segment_file) <- c("chromosome", "start", "end", "segVal", "sample")

if(file.exists("data/cell_fit_qc_table.txt")){
	filter.qc <- read.table("cell_fit_qc_table.txt",sep="\t",header=T)
	qc.samples <- filter.qc$sample[filter.qc$use == "TRUE"]
} else {
	qc.samples <- segment_file$sample
}
segment_file <- segment_file[segment_file$sample %in% qc.samples,]
segment_file <- segment_file[!is.na(segment_file$segVal),]

write.table(x = segment_file,"data/cellLine_segment_ascat_sc_fixed_purity_tCN.tsv",sep="\t",row.names=F,col.names=T,quote=F)
saveRDS(object = segment_file,file="data/cellLine_segment_ascat_sc_fixed_purity_tCN.rds")
