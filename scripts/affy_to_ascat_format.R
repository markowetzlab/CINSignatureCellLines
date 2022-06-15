# Parse cmdargs
args <- commandArgs(trailingOnly = T)
outDir <- args[1]

# Load joint lrr and baf table
snp.data <- read.table(file = paste0(outDir,"/raw/lrr_baf.txt"),
                       header = TRUE,
                       sep = "\t",
                       check.names = FALSE)

# Assign colnames to each logR and baf file
logRratio.cols <- c("Name","Chr","Position",colnames(snp.data)[grep(pattern = "\\.Log R Ratio",x = colnames(snp.data))]) 
baf.cols <- c("Name","Chr","Position",colnames(snp.data)[grep(pattern = "\\.B Allele Freq",x = colnames(snp.data))])

# Split lrr and baf table into separate logR and Baf files
logRratio <- snp.data[,which(colnames(snp.data) %in% logRratio.cols)]
colnames(logRratio) <- c("name","chrs","pos",gsub(pattern = "\\.CEL\\.Log R Ratio|\\.cel\\.Log R Ratio",replacement = "",x = colnames(logRratio[-c(1:3)])))

baf <- snp.data[,which(colnames(snp.data) %in% baf.cols)]
colnames(baf) <- c("name","chrs","pos",gsub(pattern = "\\.CEL\\.B Allele Freq|\\.cel\\.B Allele Freq",replacement = "",x = colnames(baf[-c(1:3)])))

for(i in colnames(baf)[-c(1:3)]){
  sBAF <- baf[,c("name","chrs","pos",i)]
  sLOGR <- logRratio[,c("name","chrs","pos",i)]
  
  write.table(x = sLOGR,
            file = paste0(outDir,"sample_files/",i,"_tumour_LogR.txt"),
            quote = FALSE,
            sep = "\t",
            row.names = FALSE,
            col.names = TRUE)

  write.table(x = sBAF,
            file = paste0(outDir,"sample_files/",i,"_tumour_BAF.txt"),
            quote = FALSE,
            sep = "\t",
            row.names = FALSE,
            col.names = TRUE)
}

# END
