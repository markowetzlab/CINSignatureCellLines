## ARGUMENTS FROM COMMAND LINE
args <- commandArgs(trailingOnly = TRUE)
ROWN <- as.numeric(args[1])
## Added by Phil
file_list <- read.table("data/cellLine_CEL_file_mapping.tsv",sep="\t",as.is=T,header=T)
SAMPLE <- file_list$fileid[ROWN]
BAFin <- fileBAF <- file_list$baf[ROWN]
LogRin <- fileLOGR <- file_list$logr[ROWN]
GENDER <- "XX"
#GENDER <- ifelse(file_list$Gender[ROWN] == "M","XY","XX")
PENALTY <- as.numeric(file_list$penalty[ROWN])
library(ASCAT.sc)
library(ASCAT)
## /!\ CHANGE THIS /!\ DEFINE WORKDIR
## Set working directory: this is where output dirs
## and files will be written - 1 dir per sample
WORKINGDIR <- "data/ascat_results/"

if(!dir.exists(WORKINGDIR)){
    dir.create(WORKINGDIR)
}

#setwd(WORKINGDIR)

## FUNCTION DEFINITIONS
## Need to overwrite predictGermline with new recipe for AffySNP6
ascat.predictGermlineGenotypes <- function (ASCATobj, platform = "AffySNP6", img.dir = ".", img.prefix = "")
{
    require(ASCAT)
    Homozygous = matrix(nrow = dim(ASCATobj$Tumor_LogR)[1], ncol = dim(ASCATobj$Tumor_LogR)[2])
    colnames(Homozygous) = colnames(ASCATobj$Tumor_LogR)
    rownames(Homozygous) = rownames(ASCATobj$Tumor_LogR)
    if (platform == "Custom10k") {
        maxHomozygous = 0.05
        proportionHetero = 0.59
        proportionHomo = 0.38
        proportionOpen = 0.02
        segmentLength = 20
    }
    else if (platform == "Illumina109k") {
        maxHomozygous = 0.05
        proportionHetero = 0.35
        proportionHomo = 0.6
        proportionOpen = 0.02
        segmentLength = 20
    }
    else if (platform == "IlluminaCytoSNP") {
        maxHomozygous = 0.05
        proportionHetero = 0.28
        proportionHomo = 0.62
        proportionOpen = 0.03
        segmentLength = 100
    }
    else if (platform == "Illumina610k") {
        maxHomozygous = 0.05
        proportionHetero = 0.295
        proportionHomo = 0.67
        proportionOpen = 0.015
        segmentLength = 30
    }
    else if (platform == "Illumina660k") {
        maxHomozygous = 0.05
        proportionHetero = 0.295
        proportionHomo = 0.67
        proportionOpen = 0.015
        segmentLength = 30
    }
    else if (platform == "Illumina700k") {
        maxHomozygous = 0.05
        proportionHetero = 0.295
        proportionHomo = 0.67
        proportionOpen = 0.015
        segmentLength = 30
    }
    else if (platform == "Illumina1M") {
        maxHomozygous = 0.05
        proportionHetero = 0.22
        proportionHomo = 0.74
        proportionOpen = 0.02
        segmentLength = 100
    }
    else if (platform == "Illumina2.5M") {
        maxHomozygous = 0.05
        proportionHetero = 0.21
        proportionHomo = 0.745
        proportionOpen = 0.03
        segmentLength = 100
    }
    else if (platform == "IlluminaOmni5") {
        maxHomozygous = 0.05
        proportionHetero = 0.13
        proportionHomo = 0.855
        proportionOpen = 0.01
        segmentLength = 100
    }
    else if (platform == "Affy10k") {
        maxHomozygous = 0.04
        proportionHetero = 0.355
        proportionHomo = 0.605
        proportionOpen = 0.025
        segmentLength = 20
    }
    else if (platform == "Affy100k") {
        maxHomozygous = 0.05
        proportionHetero = 0.27
        proportionHomo = 0.62
        proportionOpen = 0.09
        segmentLength = 30
    }
    else if (platform == "Affy250k_sty") {
        maxHomozygous = 0.05
        proportionHetero = 0.26
        proportionHomo = 0.66
        proportionOpen = 0.05
        segmentLength = 50
    }
    else if (platform == "Affy250k_nsp") {
        maxHomozygous = 0.05
        proportionHetero = 0.26
        proportionHomo = 0.66
        proportionOpen = 0.05
        segmentLength = 50
    }
    else if (platform == "AffySNP6") { ## CHANGED TO BE MORE STRINGENT
        maxHomozygous = 0.15
        proportionHetero = 0.25
        proportionHomo = 0.67
        proportionOpen = 0.04
        segmentLength = 300
    }
    else if (platform == "AffyOncoScan") {
        maxHomozygous = 0.04
        proportionHetero = 0.355
        proportionHomo = 0.605
        proportionOpen = 0.025
        segmentLength = 30
    }
    else if (platform == "AffyCytoScanHD") {
        maxHomozygous = 0.04
        proportionHetero = 0.32
        proportionHomo = 0.6
        proportionOpen = 0.03
        segmentLength = 100
    }
    else if (platform == "HumanCNV370quad") {
        maxHomozygous = 0.05
        proportionHetero = 0.295
        proportionHomo = 0.67
        proportionOpen = 0.015
        segmentLength = 20
    }
    else if (platform == "HumanCore12") {
        maxHomozygous = 0.05
        proportionHetero = 0.295
        proportionHomo = 0.67
        proportionOpen = 0.015
        segmentLength = 20
    }
    else if (platform == "HumanCoreExome24") {
        maxHomozygous = 0.05
        proportionHetero = 0.175
        proportionHomo = 0.79
        proportionOpen = 0.02
        segmentLength = 100
    }
    else if (platform == "HumanOmniExpress12") {
        maxHomozygous = 0.05
        proportionHetero = 0.295
        proportionHomo = 0.67
        proportionOpen = 0.015
        segmentLength = 100
    }
    else if (platform == "IlluminaOmniExpressExome") {
        maxHomozygous = 0.05
        proportionHetero = 0.35
        proportionHomo = 0.6
        proportionOpen = 0.03
        segmentLength = 100
    }
    else {
        print("Error: platform unknown")
    }
  failedarrays = NULL
    for (i in 1:dim(ASCATobj$Tumor_LogR)[2]) {
        Tumor_BAF_noNA = ASCATobj$Tumor_BAF[!is.na(ASCATobj$Tumor_BAF[,
            i]), i]
        names(Tumor_BAF_noNA) = rownames(ASCATobj$Tumor_BAF)[!is.na(ASCATobj$Tumor_BAF[,
            i])]
        Tumor_LogR_noNA = ASCATobj$Tumor_LogR[names(Tumor_BAF_noNA),
            i]
        names(Tumor_LogR_noNA) = names(Tumor_BAF_noNA)
        chr_noNA = list()
        prev = 0
        for (j in 1:length(ASCATobj$chr)) {
            chrke = ASCATobj$chr[[j]]
            next2 = prev + sum(!is.na(ASCATobj$Tumor_BAF[chrke,
                i]))
            chr_noNA[[j]] = (prev + 1):next2
            prev = next2
        }
        ch_noNA = list()
        prev = 0
        for (j in 1:length(ASCATobj$ch)) {
            chrke = ASCATobj$ch[[j]]
            next2 = prev + sum(!is.na(ASCATobj$Tumor_BAF[chrke,
                i]))
            ch_noNA[[j]] = (prev + 1):next2
            prev = next2
        }
        tbsam = Tumor_BAF_noNA
        bsm = ifelse(tbsam < 0.5, tbsam, 1 - tbsam)
        homoLimit = max(sort(bsm)[round(length(bsm) * proportionHomo)],
            maxHomozygous)
        if (homoLimit > 0.25) {
            failedarrays = c(failedarrays, ASCATobj$samples[i])
        }
        Hom = ifelse(bsm < homoLimit, T, NA)
        Homo = sum(Hom == T, na.rm = T)
        Undecided = sum(is.na(Hom))
        extraHetero = round(min(proportionHetero * length(Tumor_BAF_noNA),
            Undecided - proportionOpen * length(Tumor_BAF_noNA)))
        Hetero = 0
        if (extraHetero > 0) {
            allProbes = 1:length(Tumor_BAF_noNA)
            nonHomoProbes = allProbes[is.na(Hom) | Hom == F]
            lowestDist = NULL
            bsmHNA = bsm
            bsmHNA[!is.na(Hom) & Hom] = NA
            for (chrke in chr_noNA) {
                chrNonHomoProbes = intersect(nonHomoProbes, chrke)
                if (length(chrNonHomoProbes) > 5) {
                  segmentLength2 = min(length(chrNonHomoProbes) -
                    1, segmentLength)
                  chrNonHomoProbesStartWindowLeft = c(rep(NA,
                    segmentLength2), chrNonHomoProbes[1:(length(chrNonHomoProbes) -
                    segmentLength2)])
                  chrNonHomoProbesEndWindowLeft = c(NA, chrNonHomoProbes[1:(length(chrNonHomoProbes) -
                    1)])
                  chrNonHomoProbesStartWindowRight = c(chrNonHomoProbes[2:length(chrNonHomoProbes)],
                    NA)
                  chrNonHomoProbesEndWindowRight = c(chrNonHomoProbes[(segmentLength2 +
                    1):length(chrNonHomoProbes)], rep(NA, segmentLength2))
                  chrNonHomoProbesStartWindowMiddle = c(rep(NA,
                    segmentLength2/2), chrNonHomoProbes[1:(length(chrNonHomoProbes) -
                    segmentLength2/2)])
                  chrNonHomoProbesEndWindowMiddle = c(chrNonHomoProbes[(segmentLength2/2 +
                    1):length(chrNonHomoProbes)], rep(NA, segmentLength2/2))
                  chrLowestDist = NULL
                  for (probeNr in 1:length(chrNonHomoProbes)) {
                    probe = chrNonHomoProbes[probeNr]
                    if (!is.na(chrNonHomoProbesStartWindowLeft[probeNr]) &
                      !is.na(chrNonHomoProbesEndWindowLeft[probeNr])) {
                      medianLeft = median(bsmHNA[chrNonHomoProbesStartWindowLeft[probeNr]:chrNonHomoProbesEndWindowLeft[probeNr]],
                        na.rm = T)
                    }
                    else {
                      medianLeft = NA
                    }
                    if (!is.na(chrNonHomoProbesStartWindowRight[probeNr]) &
                      !is.na(chrNonHomoProbesEndWindowRight[probeNr])) {
                      medianRight = median(bsmHNA[chrNonHomoProbesStartWindowRight[probeNr]:chrNonHomoProbesEndWindowRight[probeNr]],
                        na.rm = T)
                    }
                    else {
                      medianRight = NA
                    }
                    if (!is.na(chrNonHomoProbesStartWindowMiddle[probeNr]) &
                      !is.na(chrNonHomoProbesEndWindowMiddle[probeNr])) {
                      medianMiddle = median(c(bsmHNA[chrNonHomoProbesStartWindowMiddle[probeNr]:chrNonHomoProbesEndWindowLeft[probeNr]],
                        bsmHNA[chrNonHomoProbesStartWindowRight[probeNr]:chrNonHomoProbesEndWindowMiddle[probeNr]]),
                        na.rm = T)
                    }
                    else {
                      medianMiddle = NA
                    }
                    chrLowestDist[probeNr] = min(abs(medianLeft -
                      bsm[probe]), abs(medianRight - bsm[probe]),
                      abs(medianMiddle - bsm[probe]), Inf, na.rm = T)
                  }
                }
                else {
                  chrLowestDist = NULL
                  if (length(chrNonHomoProbes) > 0) {
                    chrLowestDist[1:length(chrNonHomoProbes)] = 1
                  }
                }
                lowestDist = c(lowestDist, chrLowestDist)
            }
            lowestDistUndecided = lowestDist[is.na(Hom[nonHomoProbes])]
            names(lowestDistUndecided) = names(Tumor_LogR_noNA)[nonHomoProbes[is.na(Hom[nonHomoProbes])]]
            sorted = sort(lowestDistUndecided)
            Hom[names(sorted[1:min(length(sorted), extraHetero)])] = F
            Hetero = sum(Hom == F, na.rm = T)
            Homo = sum(Hom == T, na.rm = T)
            Undecided = sum(is.na(Hom))
        }
        png(filename = file.path(img.dir, paste(img.prefix, "tumorSep",
            colnames(ASCATobj$Tumor_LogR)[i], ".png", sep = "")),
            width = 2000, height = 500, res = 200)
        title = paste(paste(colnames(ASCATobj$Tumor_BAF)[i],
            Hetero), Homo)
        ascat.plotGenotypes(ASCATobj, title, Tumor_BAF_noNA,
            Hom, ch_noNA)
        dev.off()
        ## Hom[is.na(Hom)] = T ## COMMENTED OUT
        Homozygous[names(Hom), i] = Hom
    }
    return(list(germlinegenotypes = Homozygous, failedarrays = failedarrays))
}


my.smooth.CNA <- function (x,chrs, ## hack from DNAcopy for smoothing logr
                           smooth.region = 10, outlier.SD.scale = 4, smooth.SD.scale = 2,
                           trim = 0.025)
{
    require(DNAcopy)
    chrom <- chrs
    uchrom <- unique(chrom)
    genomdat <- x
    ina <- which(is.finite(genomdat))
    trimmed.SD <- sqrt(DNAcopy:::trimmed.variance(genomdat[ina], trim))
    outlier.SD <- outlier.SD.scale * trimmed.SD
    smooth.SD <- smooth.SD.scale * trimmed.SD
    k <- smooth.region
    n <- length(genomdat[ina])
    cfrq <- diff(c(which(!duplicated(chrom[ina])), n + 1))
    nchr <- length(cfrq)
    smoothed.data <- .Fortran("smoothLR", as.integer(n),
                              as.double(genomdat[ina]), as.integer(nchr), as.integer(cfrq),
                              sgdat = double(n), as.integer(k), as.double(outlier.SD),
                              as.double(smooth.SD), PACKAGE = "DNAcopy")$sgdat
    x[ina] <- smoothed.data
    x
}


## Function to use ASCAT segments for ASCAT.sc
getTrackForAll.snp6.v2 <- function (logr,
                                    CHRS,
                                    STARTS,
                                    ENDS,
                                    allchr = ALLCHR,
                                    ascat.output)
{
    require(GenomicRanges)
    segraw <- ascat.output$segments_raw
    ## smoothed <- my.smooth.CNA(logr, chrs=CHRS)
    smoothed <- logr
    lT <- lapply(allchr, function(chr) {
        cond <- CHRS == paste0("chr", chr)
        data.frame(start = STARTS[cond], end = ENDS[cond], records = logr[cond],
            fitted = logr[cond], smoothed = smoothed[cond])
    })
    lSe <- lapply(allchr, function(chr) {
        cond <- CHRS == paste0("chr", chr)
        list(starts = STARTS[cond], ends = ENDS[cond])
    })
    getSegs <- function(logr,chr,start,end, segraw)
    {
        if(chr==23) chr <- "X"
        keep1 <- gsub("chr","",segraw[,"chr"])==unique(gsub("chr","",chr))
        gr1 <- GRanges(gsub("chr","",segraw[keep1,"chr"]),
                       IRanges(segraw[keep1,"startpos"],
                               segraw[keep1,"endpos"]))
        gr2 <- GRanges(gsub("chr","",chr),
                       IRanges(start,end))
        ovs <- findOverlaps(gr1,gr2)
        data <- data.frame(chrom=chr,
                           maploc=start,
                           Sample.1=logr)
        output <- data.frame(ID="Sample.1",
                             chrom=chr,
                             loc.start=segraw[keep1,"startpos"],
                             loc.end=segraw[keep1,"endpos"],
                             num.mark=0,
                             seg.mean=NA)
        means <- tapply(1:length(ovs),queryHits(ovs),function(x) mean(logr[subjectHits(ovs)],na.rm=T))
        nmarks <- tapply(1:length(ovs),queryHits(ovs),function(x) length(x))
        output[as.numeric(names(means)),"num.mark"] <- nmarks
        output[as.numeric(names(means)),"seg.mean"] <- means
        list(data=data,
             output=output)
    }
    lSegs <- lapply(1:length(lT), function(x) {
        require(DNAcopy)
        cat(".")
        segments <- getSegs(lT[[x]]$smoothed, chr = paste0(x),
                             lSe[[x]]$starts, lSe[[x]]$ends, segraw)
    })
    names(lSegs) <- paste0(1:length(lT))
    tracks <- list(lCTS = lT, lSegs = lSegs)
    return(tracks)
}

## plotSolution for ASCAT.sc
plotSolution <- function (tracksSingle, purity, ploidy, lwdSeg = 2, ylim = c(0, 8),
                          gamma = 1, ismale = F, isPON = F, colFit = rgb(0.756,0.494, 0.756), ...)
{
    meansSeg <- fitProfile(tracksSingle, purity, ploidy, gamma = gamma,
        ismale = ismale, isPON = isPON)
    tracksSingle <- normaliseByPloidy(tracksSingle)
    breaks <- c(0, cumsum(sapply(tracksSingle$lSegs, function(x) max(x$output$loc.end))/1e+06))
    plot(0, 0, col = rgb(0, 0, 0, 0), xaxt = "n", yaxt = "n",
        xlim = c(0, max(breaks)), xlab = "Genomic Position",
        ylab = "Total copy number", frame = F, ylim = ylim, ...)
    axis(side = 2)
    for (i in 1:length(tracksSingle$lSegs)) {
        segments(tracksSingle$lCTS[[i]]$start/1e+06 + breaks[i],
            transform_bulk2tumour(tracksSingle$lCTS[[i]]$smoothed,
                purity, ploidy, gamma = gamma, ismale = ismale,
                isPON = isPON, isX = i == 23), tracksSingle$lCTS[[i]]$end/1e+06 +
                breaks[i], transform_bulk2tumour(tracksSingle$lCTS[[i]]$smoothed,
                purity, ploidy, gamma = gamma, ismale = ismale,
                isPON = isPON, isX = i == 23), col = rgb(0.7,
                0.7, 0.7, 0.6))
        segments(tracksSingle$lSegs[[i]]$output$loc.start/1e+06 +
            breaks[i], transform_bulk2tumour(sapply(meansSeg[[i]],
            function(x) x$mu), purity, ploidy, gamma = gamma,
            ismale = ismale, isPON = isPON, isX = i == 23), tracksSingle$lSegs[[i]]$output$loc.end/1e+06 +
            breaks[i], transform_bulk2tumour(sapply(meansSeg[[i]],
            function(x) x$mu), purity, ploidy, gamma = gamma,
            ismale = ismale, isPON = isPON, isX = i == 23), lwd = lwdSeg,
            col = rgb(0.4, 0.4, 0.4, 0.4))
        segments(tracksSingle$lSegs[[i]]$output$loc.start/1e+06 +
            breaks[i], round(sapply(meansSeg[[i]], function(x) x$roundmu)),
            tracksSingle$lSegs[[i]]$output$loc.end/1e+06 + breaks[i],
            round(sapply(meansSeg[[i]], function(x) x$roundmu)),
            lwd = 2.5, col = colFit)
    }
    abline(v = breaks, lwd = 1, lty = 2, col = rgb(0.6, 0.6,
        0.6, 0.4))
    abline(h = 1:50, lwd = 1, lty = 2, col = rgb(0.6, 0.6, 0.6,
        0.2))
    text(x = breaks[2:length(breaks)] - 25, y = max(ylim), names(breaks)[2:length(breaks)],
        cex = 0.4)
    mtext(side = 3, paste0("purity=", signif(purity, 2), "; average ploidy=",
        signif(ploidy, 2), "; tumor ploidy=", signif(getTumourPhi(ploidy,
            purity), 2)))
}

## add BAF/LogR and density het SNPs to ASCAT profile
annotateProfiles <- function(profile, BAF, LogR, gg)
{
    library(GenomicRanges)
    grProf <- GRanges(profile[,1],IRanges(as.numeric(profile[,2]),
                                          as.numeric(profile[,3])))
    grBAF <-  GRanges(BAF[,1],IRanges(BAF[,2],BAF[,2]))
    ovs <- findOverlaps(grProf,grBAF)
    sizes <- as.numeric(profile[,"end"])-as.numeric(profile[,"start"])
    profile <- cbind(as.data.frame(profile),
                     sizes=sizes,
                     Nprobes=as.numeric(NA),
                     Nkept=as.numeric(NA),
                     MeanProbes=as.numeric(NA),
                     BAF=as.numeric(NA),
                     LogR=as.numeric(NA))
    BAFs <- tapply(1:length(ovs),queryHits(ovs),function(x)
    {
        keep <- 1:nrow(BAF)%in%subjectHits(ovs)[x]
        keep <- which(keep & gg$germlinegenotypes[,1])
        med <- median(BAF[keep,3],na.rm=T)
        return(med)
    })
    LogRs <- tapply(1:length(ovs),queryHits(ovs),function(x)
    {
        med <- median(LogR[subjectHits(ovs)[x],3],na.rm=T)
        med
    })
    Nprobes <- tapply(1:length(ovs),queryHits(ovs),function(x)
    {
        sum(!is.na(BAF[subjectHits(ovs)[x],3]))
    })
    Nkept <- tapply(1:length(ovs),queryHits(ovs),function(x)
    {
        sum(!is.na(BAF[subjectHits(ovs)[x],3]) & !gg$germlinegenotypes[subjectHits(ovs)[x],1],na.rm=T)
    })
    profile[as.numeric(names(BAFs)),"BAF"] <- BAFs
    profile[as.numeric(names(LogRs)),"LogR"] <- LogRs
    profile[as.numeric(names(Nprobes)),"Nprobes"] <- Nprobes
    profile[as.numeric(names(Nkept)),"Nkept"] <- Nkept
    profile[,"MeanProbes"] <- profile[,"Nprobes"]/sizes
    profile[,"MeanProbesKept"] <- profile[,"Nkept"]/sizes
    profile
}

## useful function for BAF prioritisation
considerBAF <- function(BAF, threshold=.2)
{
    sapply(BAF,function(x) if(is.na(x)) return(NA) else if(abs(x-.5)<threshold) return(x) else return(NA))
}

## main function for rescuing SNPs in the very imbalanced segments
getmoreHetSNPs_where_depleted <- function(gg, nprofile, BAF, BAFPCfed, hetProp=.2)
{
    require(GenomicRanges)
    rownames(BAFPCfed) <- BAFPCfed[,1]
    BAF <- BAF[grepl("SNP",rownames(BAF)),]
    grBAF <- GRanges(BAF[,1],IRanges(BAF[,2],BAF[,2]))
    grProf <- GRanges(nprofile[,1],
                      IRanges(as.numeric(as.character(nprofile[,2])),
                              as.numeric(as.character(nprofile[,3]))))
    ovs <- findOverlaps(grProf,grBAF)
    props <- nprofile$MeanProbesKept/nprofile$MeanProbes
    wwDepleted <- which(props<hetProp)
    okRetrieve <- !is.na(gg$germlinegenotypes[,1])
    for(ww in wwDepleted)
    {
        wwBAFinDepleted <- subjectHits(ovs)[queryHits(ovs)==ww]
        prop <- props[ww]
        medianBAF <- median(BAFPCfed[rownames(BAFPCfed)%in%rownames(BAF)[wwBAFinDepleted],2],na.rm=T)
        try({
            if(!is.null(medianBAF))
                if(!is.na(medianBAF))
                    if(medianBAF>.75)
                    {
                        bbunsorted <- BAF[wwBAFinDepleted,3]
                        bbunsorted <- ifelse(bbunsorted<.5,1-bbunsorted,bbunsorted)
                        orderbb <- order(bbunsorted,decreasing=F)
                        extra <- wwBAFinDepleted[orderbb]
                        inds <- 1:round(length(wwBAFinDepleted)*(hetProp)-prop)
                        keep <- rownames(BAF)[extra]%in%rownames(gg$germlinegenotypes)[okRetrieve]
                        rr <- rownames(BAF)[extra][keep][inds]
                        gg$germlinegenotypes[rownames(gg$germlinegenotypes)%in%rr,1] <- F
                    }
        })
    }
    gg
}
## END FUNCTION DEFINITIONS
## Overwrite fitProfiles - due to NAs introduced by old minimum bin
fitProfile = function (tracksSingle, purity, ploidy, ismale = F, isPON = F,
    gamma = 1)
{
    tracksSingle <- normaliseByPloidy(tracksSingle)
    meansSeg <- lapply(1:length(tracksSingle$lSegs), function(i) {
        out <- tracksSingle$lSegs[[i]]$output
        means <- lapply(1:nrow(out), function(x) {
            isIn <- tracksSingle$lCTS[[i]]$start > out$loc.start[x] &
                tracksSingle$lCTS[[i]]$start <= out$loc.end[x]
            if (sum(isIn) < 1)
                return(list(roundmu = NA, mu = NA, sd = NA, start = out$loc.start[x],
                  end = out$loc.end[x]))
            mu <- median(tracksSingle$lCTS[[i]]$smoothed[isIn],
                na.rm = T)
            sd <- mad(tracksSingle$lCTS[[i]]$smoothed[isIn],
                na.rm = T)
            list(roundmu = transform_bulk2tumour(mu, purity,
                ploidy, gamma = gamma, ismale = ismale, isPON = isPON,
                isX = tracksSingle$lSegs[[i]]$output$chrom[1] %in% c("23", "X")), mu = mu, sd = sd, start = out$loc.start[x],end = out$loc.end[x])
        })
    })
}
## PIPELINE
## STEP1 LOAD RAW DATA & CHANGE WORKDIR TO SAMPLEDIR
BAF <- as.data.frame(data.table::fread(BAFin))
rownames(BAF) <- BAF[,1]; BAF <- BAF[,-c(1)]
LogR <- as.data.frame(data.table::fread(LogRin))
rownames(LogR) <- LogR[,1]; LogR <- LogR[,-c(1)]
OUTDIR <- paste0(WORKINGDIR,SAMPLE,"_",PENALTY,"/")
system(paste0("mkdir ",OUTDIR))
## END STEP1 LOAD RAW DATA & CHANGE WORKDIR TO SAMPLEDIR

## STEP2 FIRST ASCAT RUN
ascat.bc  <-  ascat.loadData(fileLOGR,fileBAF,NULL, NULL, gender=GENDER)
gg <- ascat.predictGermlineGenotypes(ascat.bc, platform = "AffySNP6",img.prefix = paste0(OUTDIR,"tumorSep"))
ascat.bc <- ascat.aspcf(ascat.bc, ascat.gg=gg, penalty = PENALTY, out.prefix = paste0(OUTDIR,"first"))
ascat.plot <- ascat.plotSegmentedData(ascat.bc, img.prefix = paste0(OUTDIR,"first"))
ascat.output <- ascat.runAscat(ascat.bc, img.prefix = paste0(OUTDIR,"first"))
save(ascat.output, file=paste0(OUTDIR,"ASCAT.output.first.Rda"))
save(ascat.bc, file=paste0(OUTDIR,"ASCAT.bc.first.Rda"))
BAFPCfed <- as.data.frame(data.table::fread(paste0(OUTDIR,"first",SAMPLE,".BAF.PCFed.txt")))
## END STEP2 FIRST ASCAT RUN

## STEP2b FIRST ASCAT.sc RUN
duplicates <- unlist(lapply(unique(LogR[,1]),function(x) duplicated(LogR[LogR[,1]==x,2])))
keep <- !is.na(LogR[,3]) & !duplicates
track <- getTrackForAll.snp6.v2(LogR[keep,3],
                                 CHRS=paste0("chr",LogR[keep,1]),
                                 STARTS=LogR[keep,2],
                                 ENDS=LogR[keep,2]+1,
                                 allchr=c(1:22,"X"),
                                 ascat.output=ascat.output)
save(track,file=paste0(OUTDIR,"track.first.Rda"))
GAMMA <- .55
solution <- searchGrid(track,
                      purs=ascat.output$aberrantcellfraction+seq(-0.04, 0.1, 0.01),
                      ploidies=ascat.output$psi+seq(-0.3,0.3,0.01),
                      maxTumourPhi=7,
                      gamma=GAMMA)
save(solution,file=paste0(OUTDIR,"solution.first.Rda"))
profile <- getProfile(fitProfile(track,
                                 purity=solution$purity,
                                 ploidy=solution$ploidy,
                                 gamma=GAMMA),
                      CHRS=c(1:22,"X"))
nprofile <- annotateProfiles(profile, BAF, LogR, gg)
nprofile <- cbind(nprofile,
                  BAFconsider=considerBAF(nprofile[,"BAF"]))
## END STEP2b FIRST ASCAT.sc RUN

## STEP3 RESCUE het SNPS in very imbalanced segments
newgg <- getmoreHetSNPs_where_depleted(gg, nprofile, BAF, BAFPCfed, hetProp=.2)
## END STEP3 RESCUE het SNPS in very imbalanced segments

## STEP4 SECOND ASCAT RUN
ascat.bc <- ascat.aspcf(ascat.bc, ascat.gg=newgg, penalty = PENALTY, out.prefix = paste0(OUTDIR,"second"))
ascat.plot <- ascat.plotSegmentedData(ascat.bc, img.prefix = paste0(OUTDIR,"second"))
ascat.output <- ascat.runAscat(ascat.bc, img.prefix = paste0(OUTDIR,"second"))
save(ascat.output, file=paste0(OUTDIR,"ASCAT.output.second.Rda"))
save(ascat.bc, file=paste0(OUTDIR,"ASCAT.bc.second.Rda"))
## Just added this to write the output from the second ASCAT run to disk
## The run could end here and the remaining code could be skipped
write.table(ascat.output$segments_raw,
            file=paste0(OUTDIR,SAMPLE,"_ASCAT_final_output.tsv"),
            sep="\t",col.names=T,row.names=F,quote=F)
## END STEP4 SECOND ASCAT RUN

## STEP4b SECOND ASCAT.sc RUN
## The following code does not have to be run. It will generate
## ASCAT.sc figures & profiles from the second segmentation run.
duplicates <- unlist(lapply(unique(LogR[,1]),function(x) duplicated(LogR[LogR[,1]==x,2])))
keep <- !is.na(LogR[,3]) & !duplicates
track <- getTrackForAll.snp6.v2(LogR[keep,3],
                                 CHRS=paste0("chr",LogR[keep,1]),
                                 STARTS=LogR[keep,2],
                                 ENDS=LogR[keep,2]+1,
                                 allchr=c(1:22,"X"),
                                 ascat.output=ascat.output)
save(track,file=paste0(OUTDIR,"track.second.Rda"))
GAMMA <- .55
aPur <- signif(ascat.output$aberrantcellfraction,2)
solution <- searchGrid(track,
                      purs=seq(0.95,1.00,0.01),
                      ploidies=ascat.output$psi+seq(-0.3,0.3,0.01),
                      maxTumourPhi=7,
                      gamma=GAMMA)
save(solution,file=paste0(OUTDIR,"solution.second.Rda"))
png(paste0(OUTDIR,"ASCAT_ASCAT.sc.sunrise.second.png"),width=7,height=7, unit="in",res=300)
plotSunrise(solution)
dev.off()
png(paste0(OUTDIR,"ASCATseg_ASCAT.sc.logr.fit.second.png"),width=10,height=4, unit="in",res=300)
plotSolution(track,
             purity=solution$purity,
             ploidy=solution$ploidy,lwdSeg=2,ylim=c(0,15),
             gamma=GAMMA)
dev.off()
load(paste0(OUTDIR,"track.second.Rda"))
load(paste0(OUTDIR,"solution.second.Rda"))
profile <- getProfile(fitProfile(track,
                                 purity=solution$purity,
                                 ploidy=solution$ploidy,
                                 gamma=GAMMA),
                      CHRS=c(1:22,"X"))
                      
write.table(profile,
            file=paste0(OUTDIR,"ASCAT_ASCAT.sc.logr.fit.txt"),
            sep="\t",
            col.names=T,row.names=F,quote=F)
### END STEP4b SECOND ASCAT.sc RUN
#
### END PIPELINE
q(save="no")
