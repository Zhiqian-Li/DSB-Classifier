#### All packages used 
library(ShortRead)
library(plot.matrix)
library(ggplot2)
library(grid)
library(ggpubr)

## Check the size of the RawData
p <- qa("NGS-LZQ-20210303_1_S27_L002_R1_001.fastq.gz")
p[["readCounts"]]

#### Step 1: Importing and reading the FASTQ files
## Importing
dfglst <- list("aa", "bb")
cntt <- 1
forsampstrm <- FastqStreamer("NGS-LZQ-20210303_1_S27_L002_R1_001.fastq.gz", 
                             500000)
while (cntt <= 3) {
  fffqob <- yield(forsampstrm)
  dfglst[[cntt]] <- sread(fffqob)
  cntt <- cntt + 1
}
bgdnstrfor <- c(dfglst[[1]], dfglst[[2]], dfglst[[3]])
totblstfr <- tables(bgdnstrfor, n = 1260)


## Save totblstfr as RDS file
saveRDS(totblstfr, "totblstfr.Rds")
saveRDS(bgdnstrfor, "bgdnstrfor.Rds")

## values for align  =====
dplearef <- read.table("pletop.txt")
drplechar <- dplearef[1, 1]
plegrnchar <- subseq(drplechar, start = 61, end = 83)
ttalnssub <- DNAString(drplechar)
denomval <- length(bgdnstrfor)
ffrnumvec <- unname(totblstfr$top)
pctnumvec <- ffrnumvec/denomval*100

## define tetrmatr maker function  =====
crttetmfun <- function(tmfargval) {
  print(tmfargval)
  ttrklo <- tmfargval
  ttrkhi <- ttrklo + 14
  nnumseqve <- names(totblstfr$top[ttrklo:ttrkhi])
  rnumvec <- as.character(ttrklo:ttrkhi)
  pctchavec <- as.character(round(pctnumvec[ttrklo:ttrkhi], digits = 2))
  ttalsstrpat <- DNAStringSet(c(plegrnchar, drplechar, nnumseqve))
  grpsltaln <- pairwiseAlignment(ttalsstrpat, ttalnssub)
  currmatr <- as.matrix(grpsltaln)
  gralnsc <- round(pid(grpsltaln, type = "PID2"), digits = 0)
  alnscorvec <- as.character(100 - gralnsc[3:17])
  mxrlabscor <- c("*", "$aln", alnscorvec)
  mxrlabrnum <- c("grna", "ref", rnumvec)
  mxrlabrpct <- c("*", "%", pctchavec)
  mxrflblnam <- paste(mxrlabscor, mxrlabrpct, mxrlabrnum, sep = "--")
  rownames(currmatr) <- mxrflblnam
  return(currmatr)
}

## do the big lap with function  =====
bgargind <- seq(1, 1246, 15)
bgmatrlstt <- lapply(bgargind, crttetmfun)


#### Step 2: Making the first Tetra alignment
## define megaplots maker function  =====
megmplofun <- function(pmxargval) {
  print(pmxargval)
  smatrlab <- paste("ple-A-F-F1F", ".matrix", pmxargval, sep = "")
  seltmatr <- bgmatrlstt[[pmxargval]]
  ttflnamlab <- paste("tetraln", pmxargval, ".pdf", sep = "")
  pdf(file = ttflnamlab, width = 8.43, height = 5.13)
  par(mar = c(3.1, 7.6, 4.1, 1.1), xpd = TRUE, adj = 0.9)
  plot(seltmatr,
       main = smatrlab,
       col = c("red2", "turquoise1", "navyblue", "navajowhite2", "white"),
       breaks = c("G", "A", "T", "C", "-"),
       border = NA,
       xlab = "",
       ylab = "",
       axis.row = list(las = 1),
       key = NULL)
  polygon(x = c(64.5, 68.5, 66.5),
          y = c(18.5, 18.5, 17.5),
          col = "black")
  mtext(c("G", "A", "T", "C"),
        side = 3,
        line = 2,
        at = c(15, 25, 35, 45),
        font = 2,
        cex = 1.2)
  pcindlo <- bgargind[pmxargval]
  pcindhi <- pcindlo + 14
  subg1pct <- round(sum(pctnumvec[pcindlo:pcindhi]), digits = 2)
  subg1cha <- as.character(subg1pct)
  rdcornval <- paste(subg1cha, "%", sep = "")
  mtext(rdcornval, side = 3, line = 2, at = -15)
  mtext("of reads", side = 3, line = 1, at = -15)
  polygon(x = c(10, 20, 20, 10), y = c(19.3, 19.3, 18.3, 18.3),
          col = "red2")
  polygon(x = c(20, 30, 30, 20), y = c(19.3, 19.3, 18.3, 18.3),
          col = "turquoise1")
  polygon(x = c(30, 40, 40, 30), y = c(19.3, 19.3, 18.3, 18.3),
          col = "navyblue")
  polygon(x = c(40, 50, 50, 40), y = c(19.3, 19.3, 18.3, 18.3),
          col = "navajowhite2")
  mtext("101",
        side = 3,
        line = 0,
        at = 101,
        font = 2,
        cex = 0.8)
  polygon(x = c(100.5 , 101.5 , 101),
          y = c(17.5, 17.5, 16.5),
          col = "black")
  dev.off()
}


## do the big lap with function  =====
mgmplargin <- c(1:84)
lapply(mgmplargin, megmplofun)

## Output single PDF file with 1246 alleles
qpdf::pdf_combine(input = c("tetraln1.pdf", "tetraln2.pdf", "tetraln3.pdf", "tetraln4.pdf", "tetraln5.pdf", "tetraln6.pdf", "tetraln7.pdf", "tetraln8.pdf", "tetraln9.pdf", "tetraln10.pdf", 
                            "tetraln11.pdf", "tetraln12.pdf", "tetraln13.pdf", "tetraln14.pdf", "tetraln15.pdf", "tetraln16.pdf", "tetraln17.pdf", "tetraln18.pdf", "tetraln19.pdf", "tetraln20.pdf", 
                            "tetraln21.pdf", "tetraln22.pdf", "tetraln23.pdf", "tetraln24.pdf", "tetraln25.pdf", "tetraln26.pdf", "tetraln27.pdf", "tetraln28.pdf", "tetraln29.pdf", "tetraln30.pdf", 
                            "tetraln31.pdf", "tetraln32.pdf", "tetraln33.pdf", "tetraln34.pdf", "tetraln35.pdf", "tetraln36.pdf", "tetraln37.pdf", "tetraln38.pdf", "tetraln39.pdf", "tetraln40.pdf", 
                            "tetraln41.pdf", "tetraln42.pdf", "tetraln43.pdf", "tetraln44.pdf", "tetraln45.pdf", "tetraln46.pdf", "tetraln47.pdf", "tetraln48.pdf", "tetraln49.pdf", "tetraln50.pdf", 
                            "tetraln51.pdf", "tetraln52.pdf", "tetraln53.pdf", "tetraln54.pdf", "tetraln55.pdf", "tetraln56.pdf", "tetraln57.pdf", "tetraln58.pdf", "tetraln59.pdf", "tetraln60.pdf", 
                            "tetraln61.pdf", "tetraln62.pdf", "tetraln63.pdf", "tetraln64.pdf", "tetraln65.pdf", "tetraln66.pdf", "tetraln67.pdf", "tetraln68.pdf", "tetraln69.pdf", "tetraln70.pdf", 
                            "tetraln71.pdf", "tetraln72.pdf", "tetraln73.pdf", "tetraln74.pdf", "tetraln75.pdf", "tetraln76.pdf", "tetraln77.pdf", "tetraln78.pdf", "tetraln79.pdf", "tetraln80.pdf", 
                            "tetraln81.pdf", "tetraln82.pdf", "tetraln83.pdf", "tetraln84.pdf"), output = "TetraAlign.pdf")

#### Step 3: Classification
## Get the WT variants
bgdnstrfor <- readRDS("bgdnstrfor.Rds")
totalreads <- length(bgdnstrfor)
totblstfrpsf1 <- readRDS("totblstfr.Rds")
reads <- unname(totblstfrpsf1$top)
top180seq <- names(totblstfrpsf1$top[1:180])
top180seqstr <- DNAStringSet(top180seq)
saveRDS(top180seqstr, "top180seqstr.Rds")
plewt <- vwhichPDict(wt_l, top180seqstr)
plewtvec <- which(lengths(plewt) > 0)
sort(unique(plewtvec))
plewtreads <- sum(reads[plewtvec])

## Pepper Dictionary search
ple_PEPPER_s_Dic <- readRDS("ple_PEPPER_s_Dic.Rds")
ple_PEPPER_l_Dic <- readRDS("ple_PEPPER_l_Dic.Rds")
pepper <- vwhichPDict(ple_PEPPER_s_Dic, top180seqstr)
plepepper <- which(lengths(pepper) > 0)
plepedicmatvec <- unique(unlist(pepper))
length(plepedicmatvec)
saveRDS(plepepper, "plepepper.Rds")

plepeDic1 <- which(pepper == 1)
plepeDic2 <- which(pepper == 2)
plepeDic3 <- which(pepper == 3)
plepeDic4 <- which(pepper == 4)
plepeDic5 <- which(pepper == 5)
plepeDic6 <- which(pepper == 6)
plepeDic7 <- which(pepper == 7)
plepeDic8 <- which(pepper == 8)
plepeDic9 <- which(pepper == 9)
plepeDic10 <- which(pepper == 10)
plepeDic11 <- which(pepper == 11)
plepeDic12 <- which(pepper == 12)
plepeDic13 <- which(pepper == 13)
plepeDic14 <- which(pepper == 14)
plepeDic15 <- which(pepper == 15)
plepeDic16 <- which(pepper == 16)
plepeDic17 <- which(pepper == 17)
plepeDic18 <- which(pepper == 18)
plepeDic19 <- which(pepper == 19)
plepeDic20 <- which(pepper == 20)
plepeDic21 <- which(pepper == 21)
plepeDic22 <- which(pepper == 22)
plepeDic23 <- which(pepper == 23)
plepeDic24 <- which(pepper == 24)
plepeDic25 <- which(pepper == 25)
plepeDic26 <- which(pepper == 26)
plepeDic27 <- which(pepper == 27)
plepeDic28 <- which(pepper == 28)
plepeDic29 <- which(pepper == 29)
plepeDic30 <- which(pepper == 30)
plepeDic31 <- which(pepper == 31)
plepeDic32 <- which(pepper == 32)
plepeDic33 <- which(pepper == 33)
plepeDic34 <- which(pepper == 34)
plepeDic35 <- which(pepper == 35)
plepeDic36 <- which(pepper == 36)
plepeDic37 <- which(pepper == 37)
plepeDic38 <- which(pepper == 38)
plepeDic39 <- which(pepper == 39)
plepeDic40 <- which(pepper == 40)
plepeDic41 <- which(pepper == 41)
plepeDic42 <- which(pepper == 42)
plepeDic43 <- which(pepper == 43)
plepeDic44 <- which(pepper == 44)
plepeDic45 <- which(pepper == 45)
plepeDic46 <- which(pepper == 46)
plepeDic47 <- which(pepper == 47)
plepeDic48 <- which(pepper == 48)
plepeDic49 <- which(pepper == 49)
plepeDic50 <- which(pepper == 50)
plepeDic51 <- which(pepper == 51)
plepeDic52 <- which(pepper == 52)
plepeDic53 <- which(pepper == 53)
plepeDic54 <- which(pepper == 54)
plepeDic55 <- which(pepper == 55)
plepeDic56 <- which(pepper == 56)
plepeDic57 <- which(pepper == 57)
plepeDic58 <- which(pepper == 58)
plepeDic59 <- which(pepper == 59)
plepeDic60 <- which(pepper == 60)
plepeDic61 <- which(pepper == 61)
plepeDic62 <- which(pepper == 62)
plepeDic63 <- which(pepper == 63)
plepeDic64 <- which(pepper == 64)
plepeDic65 <- which(pepper == 65)
plepeDic66 <- which(pepper == 66)
plepeDic67 <- which(pepper == 67)
plepeDic68 <- which(pepper == 68)
plepeDic69 <- which(pepper == 69)
plepeDic70 <- which(pepper == 70)
plepeDic71 <- which(pepper == 71)
plepeDic72 <- which(pepper == 72)
plepeDic73 <- which(pepper == 73)
plepeppervec <- c(plepeDic1, plepeDic2, plepeDic3, plepeDic4, plepeDic5, plepeDic6, plepeDic7, plepeDic8, plepeDic9, plepeDic10,
                  plepeDic11, plepeDic12, plepeDic13, plepeDic14, plepeDic15, plepeDic16, plepeDic17, plepeDic18, plepeDic19, plepeDic20,
                  plepeDic21, plepeDic22, plepeDic23, plepeDic24, plepeDic25, plepeDic26, plepeDic27, plepeDic28, plepeDic29, plepeDic30,
                  plepeDic31, plepeDic32, plepeDic33, plepeDic34, plepeDic35, plepeDic36, plepeDic37, plepeDic38, plepeDic39, plepeDic40,
                  plepeDic41, plepeDic42, plepeDic43, plepeDic44, plepeDic45, plepeDic46, plepeDic47, plepeDic48, plepeDic49, plepeDic50,
                  plepeDic51, plepeDic52, plepeDic53, plepeDic54, plepeDic55, plepeDic56, plepeDic57, plepeDic58, plepeDic59, plepeDic60,
                  plepeDic61, plepeDic62, plepeDic63, plepeDic64, plepeDic65, plepeDic66, plepeDic67, plepeDic68, plepeDic69, plepeDic70,
                  plepeDic71, plepeDic72, plepeDic73) 

plepe1reads <- sum(reads[plepeDic1])
plepe2reads <- sum(reads[plepeDic2])
plepe3reads <- sum(reads[plepeDic3])
plepe4reads <- sum(reads[plepeDic4])
plepe5reads <- sum(reads[plepeDic5])
plepe6reads <- sum(reads[plepeDic6])
plepe7reads <- sum(reads[plepeDic7])
plepe8reads <- sum(reads[plepeDic8])
plepe9reads <- sum(reads[plepeDic9])
plepe10reads <- sum(reads[plepeDic10])
plepe11reads <- sum(reads[plepeDic11])
plepe12reads <- sum(reads[plepeDic12])
plepe13reads <- sum(reads[plepeDic13])
plepe14reads <- sum(reads[plepeDic14])
plepe15reads <- sum(reads[plepeDic15])
plepe16reads <- sum(reads[plepeDic16])
plepe17reads <- sum(reads[plepeDic17])
plepe18reads <- sum(reads[plepeDic18])
plepe19reads <- sum(reads[plepeDic19])
plepe20reads <- sum(reads[plepeDic20])
plepe21reads <- sum(reads[plepeDic21])
plepe22reads <- sum(reads[plepeDic22])
plepe23reads <- sum(reads[plepeDic23])
plepe24reads <- sum(reads[plepeDic24])
plepe25reads <- sum(reads[plepeDic25])
plepe26reads <- sum(reads[plepeDic26])
plepe27reads <- sum(reads[plepeDic27])
plepe28reads <- sum(reads[plepeDic28])
plepe29reads <- sum(reads[plepeDic29])
plepe30reads <- sum(reads[plepeDic30])
plepe31reads <- sum(reads[plepeDic31])
plepe32reads <- sum(reads[plepeDic32])
plepe33reads <- sum(reads[plepeDic33])
plepe34reads <- sum(reads[plepeDic34])
plepe35reads <- sum(reads[plepeDic35])
plepe36reads <- sum(reads[plepeDic36])
plepe37reads <- sum(reads[plepeDic37])
plepe38reads <- sum(reads[plepeDic38])
plepe39reads <- sum(reads[plepeDic39])
plepe40reads <- sum(reads[plepeDic40])
plepe41reads <- sum(reads[plepeDic41])
plepe42reads <- sum(reads[plepeDic42])
plepe43reads <- sum(reads[plepeDic43])
plepe44reads <- sum(reads[plepeDic44])
plepe45reads <- sum(reads[plepeDic45])
plepe46reads <- sum(reads[plepeDic46])
plepe47reads <- sum(reads[plepeDic47])
plepe48reads <- sum(reads[plepeDic48])
plepe49reads <- sum(reads[plepeDic49])
plepe50reads <- sum(reads[plepeDic50])
plepe51reads <- sum(reads[plepeDic51])
plepe52reads <- sum(reads[plepeDic52])
plepe53reads <- sum(reads[plepeDic53])
plepe54reads <- sum(reads[plepeDic54])
plepe55reads <- sum(reads[plepeDic55])
plepe56reads <- sum(reads[plepeDic56])
plepe57reads <- sum(reads[plepeDic57])
plepe58reads <- sum(reads[plepeDic58])
plepe59reads <- sum(reads[plepeDic59])
plepe60reads <- sum(reads[plepeDic60])
plepe61reads <- sum(reads[plepeDic61])
plepe62reads <- sum(reads[plepeDic62])
plepe63reads <- sum(reads[plepeDic63])
plepe64reads <- sum(reads[plepeDic64])
plepe65reads <- sum(reads[plepeDic65])
plepe66reads <- sum(reads[plepeDic66])
plepe67reads <- sum(reads[plepeDic67])
plepe68reads <- sum(reads[plepeDic68])
plepe69reads <- sum(reads[plepeDic69])
plepe70reads <- sum(reads[plepeDic70])
plepe71reads <- sum(reads[plepeDic71])
plepe72reads <- sum(reads[plepeDic72])
plepe73reads <- sum(reads[plepeDic73])
plepedicmatreads <- c(plepe1reads, plepe2reads, plepe3reads, plepe4reads, plepe5reads, plepe6reads, plepe7reads, plepe8reads, plepe9reads, plepe10reads,
                      plepe11reads, plepe12reads, plepe13reads, plepe14reads, plepe15reads, plepe16reads, plepe17reads, plepe18reads, plepe19reads, plepe20reads,
                      plepe21reads, plepe22reads, plepe23reads, plepe24reads, plepe25reads, plepe26reads, plepe27reads, plepe28reads, plepe29reads, plepe30reads,
                      plepe31reads, plepe32reads, plepe33reads, plepe34reads, plepe35reads, plepe36reads, plepe37reads, plepe38reads, plepe39reads, plepe40reads,
                      plepe41reads, plepe42reads, plepe43reads, plepe44reads, plepe45reads, plepe46reads, plepe47reads, plepe48reads, plepe49reads, plepe50reads,
                      plepe51reads, plepe52reads, plepe53reads, plepe54reads, plepe55reads, plepe56reads, plepe57reads, plepe58reads, plepe59reads, plepe60reads,
                      plepe61reads, plepe62reads, plepe63reads, plepe64reads, plepe65reads, plepe66reads, plepe67reads, plepe68reads, plepe69reads, plepe70reads,
                      plepe71reads, plepe72reads, plepe73reads)

## MMEJ Dictionary search
ple_MMEJ_s_Dic <- readRDS("ple_MMEJ_s_Dic.Rds")
mmej <- vwhichPDict(ple_MMEJ_s_Dic, top180seqstr)
plemmej <- which(lengths(mmej) > 0)
ple_MMEJ_l_Dic <- readRDS("ple_MMEJ_l_Dic.Rds")
plemjadicmatvec <- unique(unlist(mmej))
length(plemjadicmatvec)

plemjaDic1 <- which(vwhichPDict(ple_MMEJ_s_Dic[1], top180seqstr) != 0)
plemjaDic2 <- which(vwhichPDict(ple_MMEJ_s_Dic[2], top180seqstr) != 0)
plemjaDic3 <- which(vwhichPDict(ple_MMEJ_s_Dic[3], top180seqstr) != 0)
plemjaDic4 <- which(vwhichPDict(ple_MMEJ_s_Dic[4], top180seqstr) != 0)
plemjaDic5 <- which(vwhichPDict(ple_MMEJ_s_Dic[5], top180seqstr) != 0)
plemjaDic6 <- which(vwhichPDict(ple_MMEJ_s_Dic[6], top180seqstr) != 0)
plemjaDic7 <- which(vwhichPDict(ple_MMEJ_s_Dic[7], top180seqstr) != 0)
plemjaDic8 <- which(vwhichPDict(ple_MMEJ_s_Dic[8], top180seqstr) != 0)
plemjaDic9 <- which(vwhichPDict(ple_MMEJ_s_Dic[9], top180seqstr) != 0)
plemjaDic10 <- which(vwhichPDict(ple_MMEJ_s_Dic[10], top180seqstr) != 0)
plemjaDic11 <- which(vwhichPDict(ple_MMEJ_s_Dic[11], top180seqstr) != 0)
plemjaDic12 <- which(vwhichPDict(ple_MMEJ_s_Dic[12], top180seqstr) != 0)
plemjaDic13 <- which(vwhichPDict(ple_MMEJ_s_Dic[13], top180seqstr) != 0)
plemjaDic14 <- which(vwhichPDict(ple_MMEJ_s_Dic[14], top180seqstr) != 0)
plemjaDic15 <- which(vwhichPDict(ple_MMEJ_s_Dic[15], top180seqstr) != 0)
plemjaDic16 <- which(vwhichPDict(ple_MMEJ_s_Dic[16], top180seqstr) != 0)
plemjaDic17 <- which(vwhichPDict(ple_MMEJ_s_Dic[17], top180seqstr) != 0)
plemjaDic18 <- which(vwhichPDict(ple_MMEJ_s_Dic[18], top180seqstr) != 0)
plemjaDic19 <- which(vwhichPDict(ple_MMEJ_s_Dic[19], top180seqstr) != 0)
plemjaDic20 <- which(vwhichPDict(ple_MMEJ_s_Dic[20], top180seqstr) != 0)
plemjaDic21 <- which(vwhichPDict(ple_MMEJ_s_Dic[21], top180seqstr) != 0)
plemjaDic22 <- which(vwhichPDict(ple_MMEJ_s_Dic[22], top180seqstr) != 0)
plemjaDic23 <- which(vwhichPDict(ple_MMEJ_s_Dic[23], top180seqstr) != 0)
plemjaDic24 <- which(vwhichPDict(ple_MMEJ_s_Dic[24], top180seqstr) != 0)
plemjaDic25 <- which(vwhichPDict(ple_MMEJ_s_Dic[25], top180seqstr) != 0)
plemjaDic26 <- which(vwhichPDict(ple_MMEJ_s_Dic[26], top180seqstr) != 0)
plemjaDic27 <- which(vwhichPDict(ple_MMEJ_s_Dic[27], top180seqstr) != 0)
plemjaDic28 <- which(vwhichPDict(ple_MMEJ_s_Dic[28], top180seqstr) != 0)
plemjaDic29 <- which(vwhichPDict(ple_MMEJ_s_Dic[29], top180seqstr) != 0)
plemjaDic30 <- which(vwhichPDict(ple_MMEJ_s_Dic[30], top180seqstr) != 0)
plemjaDic31 <- which(vwhichPDict(ple_MMEJ_s_Dic[31], top180seqstr) != 0)
plemjaDic32 <- which(vwhichPDict(ple_MMEJ_s_Dic[32], top180seqstr) != 0)
plemjaDic33 <- which(vwhichPDict(ple_MMEJ_s_Dic[33], top180seqstr) != 0)
plemjaDic34 <- which(vwhichPDict(ple_MMEJ_s_Dic[34], top180seqstr) != 0)
plemjaDic35 <- which(vwhichPDict(ple_MMEJ_s_Dic[35], top180seqstr) != 0)
plemjaDic36 <- which(vwhichPDict(ple_MMEJ_s_Dic[36], top180seqstr) != 0)
plemjaDic37 <- which(vwhichPDict(ple_MMEJ_s_Dic[37], top180seqstr) != 0)
plemjaDic38 <- which(vwhichPDict(ple_MMEJ_s_Dic[38], top180seqstr) != 0)
plemjaDic39 <- which(vwhichPDict(ple_MMEJ_s_Dic[39], top180seqstr) != 0)
plemjaDic40 <- which(vwhichPDict(ple_MMEJ_s_Dic[40], top180seqstr) != 0)
plemjaDic41 <- which(vwhichPDict(ple_MMEJ_s_Dic[41], top180seqstr) != 0)
plemjaDic42 <- which(vwhichPDict(ple_MMEJ_s_Dic[42], top180seqstr) != 0)

plemjavec <- c(plemjaDic1, plemjaDic2, plemjaDic3, plemjaDic4, plemjaDic5, plemjaDic6, plemjaDic7, plemjaDic8, plemjaDic9, plemjaDic10,
              plemjaDic11, plemjaDic12, plemjaDic13, plemjaDic14, plemjaDic15, plemjaDic16, plemjaDic17, plemjaDic18, plemjaDic19, plemjaDic20,
              plemjaDic21, plemjaDic22, plemjaDic23, plemjaDic24, plemjaDic25, plemjaDic26, plemjaDic27, plemjaDic28, plemjaDic29, plemjaDic30,
              plemjaDic31, plemjaDic32, plemjaDic33, plemjaDic34, plemjaDic35, plemjaDic36, plemjaDic37,
              plemjaDic38, plemjaDic39, plemjaDic40, plemjaDic41, plemjaDic42) 

plemja1reads <- sum(reads[plemjaDic1])
plemja2reads <- sum(reads[plemjaDic2])
plemja3reads <- sum(reads[plemjaDic3])
plemja4reads <- sum(reads[plemjaDic4])
plemja5reads <- sum(reads[plemjaDic5])
plemja6reads <- sum(reads[plemjaDic6])
plemja7reads <- sum(reads[plemjaDic7])
plemja8reads <- sum(reads[plemjaDic8])
plemja9reads <- sum(reads[plemjaDic9])
plemja10reads <- sum(reads[plemjaDic10])
plemja11reads <- sum(reads[plemjaDic11])
plemja12reads <- sum(reads[plemjaDic12])
plemja13reads <- sum(reads[plemjaDic13])
plemja14reads <- sum(reads[plemjaDic14])
plemja15reads <- sum(reads[plemjaDic15])
plemja16reads <- sum(reads[plemjaDic16])
plemja17reads <- sum(reads[plemjaDic17])
plemja18reads <- sum(reads[plemjaDic18])
plemja19reads <- sum(reads[plemjaDic19])
plemja20reads <- sum(reads[plemjaDic20])
plemja21reads <- sum(reads[plemjaDic21])
plemja22reads <- sum(reads[plemjaDic22])
plemja23reads <- sum(reads[plemjaDic23])
plemja24reads <- sum(reads[plemjaDic24])
plemja25reads <- sum(reads[plemjaDic25])
plemja26reads <- sum(reads[plemjaDic26])
plemja27reads <- sum(reads[plemjaDic27])
plemja28reads <- sum(reads[plemjaDic28])
plemja29reads <- sum(reads[plemjaDic29])
plemja30reads <- sum(reads[plemjaDic30])
plemja31reads <- sum(reads[plemjaDic31])
plemja32reads <- sum(reads[plemjaDic32])
plemja33reads <- sum(reads[plemjaDic33])
plemja34reads <- sum(reads[plemjaDic34])
plemja35reads <- sum(reads[plemjaDic35])
plemja36reads <- sum(reads[plemjaDic36])
plemja37reads <- sum(reads[plemjaDic37])
plemja38reads <- sum(reads[plemjaDic38])
plemja39reads <- sum(reads[plemjaDic39])
plemja40reads <- sum(reads[plemjaDic40])
plemja41reads <- sum(reads[plemjaDic41])
plemja42reads <- sum(reads[plemjaDic42])

plemjadicmatreads <- c(plemja1reads, plemja2reads, plemja3reads, plemja4reads, plemja5reads, plemja6reads, plemja7reads, plemja8reads, plemja9reads, plemja10reads,
                       plemja11reads, plemja12reads, plemja13reads, plemja14reads, plemja15reads, plemja16reads, plemja17reads, plemja18reads, plemja19reads, plemja20reads,
                       plemja21reads, plemja22reads, plemja23reads, plemja24reads, plemja25reads, plemja26reads, plemja27reads, plemja28reads, plemja29reads, plemja30reads,
                       plemja31reads, plemja32reads, plemja33reads, plemja34reads, plemja35reads, plemja36reads, plemja37reads, plemja38reads, plemja39reads, plemja40reads, 
                       plemja41reads, plemja42reads)

## deletion Wdel+Minidel Dictionary search
ple_DEL_s_Dic <- readRDS("ple_DEL_s_Dic.Rds")
ple_DEL_l_Dic <- readRDS("ple_DEL_l_Dic.Rds")
del <- vwhichPDict(ple_DEL_s_Dic, top180seqstr)
pledel <- which(lengths(del) > 0)
pledeldicmatvec <- unique(unlist(del))
length(pledeldicmatvec)

pledelDic1 <- which(lengths(vwhichPDict(ple_DEL_s_Dic[1], top180seqstr)) >0)
pledelDic2 <- which(lengths(vwhichPDict(ple_DEL_s_Dic[2], top180seqstr)) >0)
pledelDic3 <- which(lengths(vwhichPDict(ple_DEL_s_Dic[3], top180seqstr)) >0)
pledelDic4 <- which(lengths(vwhichPDict(ple_DEL_s_Dic[4], top180seqstr)) >0)
pledelDic5 <- which(lengths(vwhichPDict(ple_DEL_s_Dic[5], top180seqstr)) >0)
pledelDic6 <- which(lengths(vwhichPDict(ple_DEL_s_Dic[6], top180seqstr)) >0)
pledelDic7 <- which(lengths(vwhichPDict(ple_DEL_s_Dic[7], top180seqstr)) >0)
pledelDic8 <- which(lengths(vwhichPDict(ple_DEL_s_Dic[8], top180seqstr)) >0)
pledelDic9 <- which(lengths(vwhichPDict(ple_DEL_s_Dic[9], top180seqstr)) >0)
pledelDic10 <- which(lengths(vwhichPDict(ple_DEL_s_Dic[10], top180seqstr)) >0)
pledelDic11 <- which(lengths(vwhichPDict(ple_DEL_s_Dic[11], top180seqstr)) >0)
pledelDic12 <- which(lengths(vwhichPDict(ple_DEL_s_Dic[12], top180seqstr)) >0)
pledelDic13 <- which(lengths(vwhichPDict(ple_DEL_s_Dic[13], top180seqstr)) >0)
pledelDic14 <- which(lengths(vwhichPDict(ple_DEL_s_Dic[14], top180seqstr)) >0)
pledelDic15 <- which(lengths(vwhichPDict(ple_DEL_s_Dic[15], top180seqstr)) >0)
pledelDic16 <- which(lengths(vwhichPDict(ple_DEL_s_Dic[16], top180seqstr)) >0)
pledelDic17 <- which(lengths(vwhichPDict(ple_DEL_s_Dic[17], top180seqstr)) >0)
pledelDic18 <- which(lengths(vwhichPDict(ple_DEL_s_Dic[18], top180seqstr)) >0)
pledelDic19 <- which(lengths(vwhichPDict(ple_DEL_s_Dic[19], top180seqstr)) >0)
pledelDic20 <- which(lengths(vwhichPDict(ple_DEL_s_Dic[20], top180seqstr)) >0)
pledelDic21 <- which(lengths(vwhichPDict(ple_DEL_s_Dic[21], top180seqstr)) >0)
pledelDic22 <- which(lengths(vwhichPDict(ple_DEL_s_Dic[22], top180seqstr)) >0)
pledelDic23 <- which(lengths(vwhichPDict(ple_DEL_s_Dic[23], top180seqstr)) >0)
pledelDic24 <- which(lengths(vwhichPDict(ple_DEL_s_Dic[24], top180seqstr)) >0)
pledelDic25 <- which(lengths(vwhichPDict(ple_DEL_s_Dic[25], top180seqstr)) >0)
pledelDic26 <- which(lengths(vwhichPDict(ple_DEL_s_Dic[26], top180seqstr)) >0)
pledelDic27 <- which(lengths(vwhichPDict(ple_DEL_s_Dic[27], top180seqstr)) >0)
pledelDic28 <- which(lengths(vwhichPDict(ple_DEL_s_Dic[28], top180seqstr)) >0)
pledelDic29 <- which(lengths(vwhichPDict(ple_DEL_s_Dic[29], top180seqstr)) >0)
pledelDic30 <- which(lengths(vwhichPDict(ple_DEL_s_Dic[30], top180seqstr)) >0)
pledelDic31 <- which(lengths(vwhichPDict(ple_DEL_s_Dic[31], top180seqstr)) >0)
pledelDic32 <- which(lengths(vwhichPDict(ple_DEL_s_Dic[32], top180seqstr)) >0)
pledelDic33 <- which(lengths(vwhichPDict(ple_DEL_s_Dic[33], top180seqstr)) >0)
pledelDic34 <- which(lengths(vwhichPDict(ple_DEL_s_Dic[34], top180seqstr)) >0)
pledelDic35 <- which(lengths(vwhichPDict(ple_DEL_s_Dic[35], top180seqstr)) >0)
pledelDic36 <- which(lengths(vwhichPDict(ple_DEL_s_Dic[36], top180seqstr)) >0)
pledelDic37 <- which(lengths(vwhichPDict(ple_DEL_s_Dic[37], top180seqstr)) >0)
pledelDic38 <- which(lengths(vwhichPDict(ple_DEL_s_Dic[38], top180seqstr)) >0)
pledelDic39 <- which(lengths(vwhichPDict(ple_DEL_s_Dic[39], top180seqstr)) >0)
pledelDic40 <- which(lengths(vwhichPDict(ple_DEL_s_Dic[40], top180seqstr)) >0)
pledelDic41 <- which(lengths(vwhichPDict(ple_DEL_s_Dic[41], top180seqstr)) >0)
pledelDic42 <- which(lengths(vwhichPDict(ple_DEL_s_Dic[42], top180seqstr)) >0)
pledelDic43 <- which(lengths(vwhichPDict(ple_DEL_s_Dic[43], top180seqstr)) >0)
pledelDic44 <- which(lengths(vwhichPDict(ple_DEL_s_Dic[44], top180seqstr)) >0)
pledelDic45 <- which(lengths(vwhichPDict(ple_DEL_s_Dic[45], top180seqstr)) >0)
pledelDic46 <- which(lengths(vwhichPDict(ple_DEL_s_Dic[46], top180seqstr)) >0)
pledelDic47 <- which(lengths(vwhichPDict(ple_DEL_s_Dic[47], top180seqstr)) >0)
pledelDic48 <- which(lengths(vwhichPDict(ple_DEL_s_Dic[48], top180seqstr)) >0)
pledelDic49 <- which(lengths(vwhichPDict(ple_DEL_s_Dic[49], top180seqstr)) >0)
pledelDic50 <- which(lengths(vwhichPDict(ple_DEL_s_Dic[50], top180seqstr)) >0)
pledelDic51 <- which(lengths(vwhichPDict(ple_DEL_s_Dic[51], top180seqstr)) >0)
pledelDic52 <- which(lengths(vwhichPDict(ple_DEL_s_Dic[52], top180seqstr)) >0)
pledelDic53 <- which(lengths(vwhichPDict(ple_DEL_s_Dic[53], top180seqstr)) >0)
pledelDic54 <- which(lengths(vwhichPDict(ple_DEL_s_Dic[54], top180seqstr)) >0)
pledelDic55 <- which(lengths(vwhichPDict(ple_DEL_s_Dic[55], top180seqstr)) >0)
pledelDic56 <- which(lengths(vwhichPDict(ple_DEL_s_Dic[56], top180seqstr)) >0)
pledelDic57 <- which(lengths(vwhichPDict(ple_DEL_s_Dic[57], top180seqstr)) >0)
pledelDic58 <- which(lengths(vwhichPDict(ple_DEL_s_Dic[58], top180seqstr)) >0)
pledelDic59 <- which(lengths(vwhichPDict(ple_DEL_s_Dic[59], top180seqstr)) >0)
pledelDic60 <- which(lengths(vwhichPDict(ple_DEL_s_Dic[60], top180seqstr)) >0)
pledelDic61 <- which(lengths(vwhichPDict(ple_DEL_s_Dic[61], top180seqstr)) >0)
pledelDic62 <- which(lengths(vwhichPDict(ple_DEL_s_Dic[62], top180seqstr)) >0)
pledelDic63 <- which(lengths(vwhichPDict(ple_DEL_s_Dic[63], top180seqstr)) >0)
pledelDic64 <- which(lengths(vwhichPDict(ple_DEL_s_Dic[64], top180seqstr)) >0)
pledelDic65 <- which(lengths(vwhichPDict(ple_DEL_s_Dic[65], top180seqstr)) >0)
pledelDic66 <- which(lengths(vwhichPDict(ple_DEL_s_Dic[66], top180seqstr)) >0)
pledelDic67 <- which(lengths(vwhichPDict(ple_DEL_s_Dic[67], top180seqstr)) >0)
pledelDic68 <- which(lengths(vwhichPDict(ple_DEL_s_Dic[68], top180seqstr)) >0)
pledelDic69 <- which(lengths(vwhichPDict(ple_DEL_s_Dic[69], top180seqstr)) >0)
pledelDic70 <- which(lengths(vwhichPDict(ple_DEL_s_Dic[70], top180seqstr)) >0)
pledelDic71 <- which(lengths(vwhichPDict(ple_DEL_s_Dic[71], top180seqstr)) >0)
pledelDic72 <- which(lengths(vwhichPDict(ple_DEL_s_Dic[72], top180seqstr)) >0)
pledelDic73 <- which(lengths(vwhichPDict(ple_DEL_s_Dic[73], top180seqstr)) >0)
pledelDic74 <- which(lengths(vwhichPDict(ple_DEL_s_Dic[74], top180seqstr)) >0)
pledelDic75 <- which(lengths(vwhichPDict(ple_DEL_s_Dic[75], top180seqstr)) >0)
pledelDic76 <- which(lengths(vwhichPDict(ple_DEL_s_Dic[76], top180seqstr)) >0)
pledelDic77 <- which(lengths(vwhichPDict(ple_DEL_s_Dic[77], top180seqstr)) >0)
pledelDic78 <- which(lengths(vwhichPDict(ple_DEL_s_Dic[78], top180seqstr)) >0)
pledelDic79 <- which(lengths(vwhichPDict(ple_DEL_s_Dic[79], top180seqstr)) >0)
pledelDic80 <- which(lengths(vwhichPDict(ple_DEL_s_Dic[80], top180seqstr)) >0)
pledelDic81 <- which(lengths(vwhichPDict(ple_DEL_s_Dic[81], top180seqstr)) >0)
pledelDic82 <- which(lengths(vwhichPDict(ple_DEL_s_Dic[82], top180seqstr)) >0)
pledelDic83 <- which(lengths(vwhichPDict(ple_DEL_s_Dic[83], top180seqstr)) >0)
pledelDic84 <- which(lengths(vwhichPDict(ple_DEL_s_Dic[84], top180seqstr)) >0)
pledelDic85 <- which(lengths(vwhichPDict(ple_DEL_s_Dic[85], top180seqstr)) >0)
pledelDic86 <- which(lengths(vwhichPDict(ple_DEL_s_Dic[86], top180seqstr)) >0)
pledelDic87 <- which(lengths(vwhichPDict(ple_DEL_s_Dic[87], top180seqstr)) >0)
pledelDic88 <- which(lengths(vwhichPDict(ple_DEL_s_Dic[88], top180seqstr)) >0)
pledelDic89 <- which(lengths(vwhichPDict(ple_DEL_s_Dic[89], top180seqstr)) >0)
pledelDic90 <- which(lengths(vwhichPDict(ple_DEL_s_Dic[90], top180seqstr)) >0)
pledelDic91 <- which(lengths(vwhichPDict(ple_DEL_s_Dic[91], top180seqstr)) >0)
pledelDic92 <- which(lengths(vwhichPDict(ple_DEL_s_Dic[92], top180seqstr)) >0)
pledelDic93 <- which(lengths(vwhichPDict(ple_DEL_s_Dic[93], top180seqstr)) >0)
pledelDic94 <- which(lengths(vwhichPDict(ple_DEL_s_Dic[94], top180seqstr)) >0)
pledelDic95 <- which(lengths(vwhichPDict(ple_DEL_s_Dic[95], top180seqstr)) >0)
pledelDic96 <- which(lengths(vwhichPDict(ple_DEL_s_Dic[96], top180seqstr)) >0)
pledelDic97 <- which(lengths(vwhichPDict(ple_DEL_s_Dic[97], top180seqstr)) >0)
pledelDic98 <- which(lengths(vwhichPDict(ple_DEL_s_Dic[98], top180seqstr)) >0)
pledelDic99 <- which(lengths(vwhichPDict(ple_DEL_s_Dic[99], top180seqstr)) >0)
pledelDic100 <- which(lengths(vwhichPDict(ple_DEL_s_Dic[100], top180seqstr)) >0)
pledelDic101 <- which(lengths(vwhichPDict(ple_DEL_s_Dic[101], top180seqstr)) >0)
pledelDic102 <- which(lengths(vwhichPDict(ple_DEL_s_Dic[102], top180seqstr)) >0)
pledelDic103 <- which(lengths(vwhichPDict(ple_DEL_s_Dic[103], top180seqstr)) >0)
pledelDic104 <- which(lengths(vwhichPDict(ple_DEL_s_Dic[104], top180seqstr)) >0)
pledelDic105 <- which(lengths(vwhichPDict(ple_DEL_s_Dic[105], top180seqstr)) >0)
pledelDic106 <- which(lengths(vwhichPDict(ple_DEL_s_Dic[106], top180seqstr)) >0)
pledelDic107 <- which(lengths(vwhichPDict(ple_DEL_s_Dic[107], top180seqstr)) >0)
pledelDic108 <- which(lengths(vwhichPDict(ple_DEL_s_Dic[108], top180seqstr)) >0)
pledelDic109 <- which(lengths(vwhichPDict(ple_DEL_s_Dic[109], top180seqstr)) >0)
pledelDic110 <- which(lengths(vwhichPDict(ple_DEL_s_Dic[110], top180seqstr)) >0)
pledelDic111 <- which(lengths(vwhichPDict(ple_DEL_s_Dic[111], top180seqstr)) >0)
pledelDic112 <- which(lengths(vwhichPDict(ple_DEL_s_Dic[112], top180seqstr)) >0)
pledelDic113 <- which(lengths(vwhichPDict(ple_DEL_s_Dic[113], top180seqstr)) >0)
pledelDic114 <- which(lengths(vwhichPDict(ple_DEL_s_Dic[114], top180seqstr)) >0)
pledelDic115 <- which(lengths(vwhichPDict(ple_DEL_s_Dic[115], top180seqstr)) >0)
pledelDic116 <- which(lengths(vwhichPDict(ple_DEL_s_Dic[116], top180seqstr)) >0)
pledelDic117 <- which(lengths(vwhichPDict(ple_DEL_s_Dic[117], top180seqstr)) >0)
pledelDic118 <- which(lengths(vwhichPDict(ple_DEL_s_Dic[118], top180seqstr)) >0)
pledelDic119 <- which(lengths(vwhichPDict(ple_DEL_s_Dic[119], top180seqstr)) >0)
pledelDic120 <- which(lengths(vwhichPDict(ple_DEL_s_Dic[120], top180seqstr)) >0)
pledelDic121 <- which(lengths(vwhichPDict(ple_DEL_s_Dic[121], top180seqstr)) >0)
pledelDic122 <- which(lengths(vwhichPDict(ple_DEL_s_Dic[122], top180seqstr)) >0)
pledelDic123 <- which(lengths(vwhichPDict(ple_DEL_s_Dic[123], top180seqstr)) >0)
pledelDic124 <- which(lengths(vwhichPDict(ple_DEL_s_Dic[124], top180seqstr)) >0)
pledelDic125 <- which(lengths(vwhichPDict(ple_DEL_s_Dic[125], top180seqstr)) >0)
pledelDic126 <- which(lengths(vwhichPDict(ple_DEL_s_Dic[126], top180seqstr)) >0)
pledelDic127 <- which(lengths(vwhichPDict(ple_DEL_s_Dic[127], top180seqstr)) >0)
pledelDic128 <- which(lengths(vwhichPDict(ple_DEL_s_Dic[128], top180seqstr)) >0)
pledelDic129 <- which(lengths(vwhichPDict(ple_DEL_s_Dic[129], top180seqstr)) >0)
pledelDic130 <- which(lengths(vwhichPDict(ple_DEL_s_Dic[130], top180seqstr)) >0)

pledelvec <- c(pledelDic1, pledelDic2, pledelDic3, pledelDic4, pledelDic5, pledelDic6, pledelDic7, pledelDic8, pledelDic9, pledelDic10,
              pledelDic11, pledelDic12, pledelDic13, pledelDic14, pledelDic15, pledelDic16, pledelDic17, pledelDic18, pledelDic19, pledelDic20,
              pledelDic21, pledelDic22, pledelDic23, pledelDic24, pledelDic25, pledelDic26, pledelDic27, pledelDic28, pledelDic29, pledelDic30,
              pledelDic31, pledelDic32, pledelDic33, pledelDic34, pledelDic35, pledelDic36, pledelDic37, pledelDic38, pledelDic39, pledelDic40,
              pledelDic41, pledelDic42, pledelDic43, pledelDic44, pledelDic45, pledelDic46, pledelDic47, pledelDic48, pledelDic49, pledelDic50,
              pledelDic51, pledelDic52, pledelDic53, pledelDic54, pledelDic55, pledelDic56, pledelDic57, pledelDic58, pledelDic59, pledelDic60,
              pledelDic61, pledelDic62, pledelDic63, pledelDic64, pledelDic65, pledelDic66, pledelDic67, pledelDic68, pledelDic69, pledelDic70, pledelDic71, pledelDic72, pledelDic73, pledelDic74, pledelDic75,
              pledelDic76, pledelDic77, pledelDic78, pledelDic79, pledelDic80, pledelDic81, pledelDic82, pledelDic83, pledelDic84, pledelDic85,
              pledelDic86, pledelDic87, pledelDic88, pledelDic89, pledelDic90, pledelDic91, pledelDic92, pledelDic93, pledelDic94, pledelDic95, pledelDic96,
              pledelDic97, pledelDic98, pledelDic99, pledelDic100, pledelDic101, pledelDic102, pledelDic103, pledelDic104, pledelDic105, pledelDic106,
              pledelDic107, pledelDic108, pledelDic109, pledelDic110, pledelDic111, pledelDic112, pledelDic113, pledelDic114, pledelDic115, pledelDic116, pledelDic117,
              pledelDic118, pledelDic119, pledelDic120, pledelDic121, pledelDic122, pledelDic123, pledelDic124, pledelDic125, pledelDic126, pledelDic127, pledelDic128, pledelDic129, pledelDic130) 


pledel1reads <- sum(reads[pledelDic1])
pledel2reads <- sum(reads[pledelDic2])
pledel3reads <- sum(reads[pledelDic3])
pledel4reads <- sum(reads[pledelDic4])
pledel5reads <- sum(reads[pledelDic5])
pledel6reads <- sum(reads[pledelDic6])
pledel7reads <- sum(reads[pledelDic7])
pledel8reads <- sum(reads[pledelDic8])
pledel9reads <- sum(reads[pledelDic9])
pledel10reads <- sum(reads[pledelDic10])
pledel11reads <- sum(reads[pledelDic11])
pledel12reads <- sum(reads[pledelDic12])
pledel13reads <- sum(reads[pledelDic13])
pledel14reads <- sum(reads[pledelDic14])
pledel15reads <- sum(reads[pledelDic15])
pledel16reads <- sum(reads[pledelDic16])
pledel17reads <- sum(reads[pledelDic17])
pledel18reads <- sum(reads[pledelDic18])
pledel19reads <- sum(reads[pledelDic19])
pledel20reads <- sum(reads[pledelDic20])
pledel21reads <- sum(reads[pledelDic21])
pledel22reads <- sum(reads[pledelDic22])
pledel23reads <- sum(reads[pledelDic23])
pledel24reads <- sum(reads[pledelDic24])
pledel25reads <- sum(reads[pledelDic25])
pledel26reads <- sum(reads[pledelDic26])
pledel27reads <- sum(reads[pledelDic27])
pledel28reads <- sum(reads[pledelDic28])
pledel29reads <- sum(reads[pledelDic29])
pledel30reads <- sum(reads[pledelDic30])
pledel31reads <- sum(reads[pledelDic31])
pledel32reads <- sum(reads[pledelDic32])
pledel33reads <- sum(reads[pledelDic33])
pledel34reads <- sum(reads[pledelDic34])
pledel35reads <- sum(reads[pledelDic35])
pledel36reads <- sum(reads[pledelDic36])
pledel37reads <- sum(reads[pledelDic37])
pledel38reads <- sum(reads[pledelDic38])
pledel39reads <- sum(reads[pledelDic39])
pledel40reads <- sum(reads[pledelDic40])
pledel41reads <- sum(reads[pledelDic41])
pledel42reads <- sum(reads[pledelDic42])
pledel43reads <- sum(reads[pledelDic43])
pledel44reads <- sum(reads[pledelDic44])
pledel45reads <- sum(reads[pledelDic45])
pledel46reads <- sum(reads[pledelDic46])
pledel47reads <- sum(reads[pledelDic47])
pledel48reads <- sum(reads[pledelDic48])
pledel49reads <- sum(reads[pledelDic49])
pledel50reads <- sum(reads[pledelDic50])
pledel51reads <- sum(reads[pledelDic51])
pledel52reads <- sum(reads[pledelDic52])
pledel53reads <- sum(reads[pledelDic53])
pledel54reads <- sum(reads[pledelDic54])
pledel55reads <- sum(reads[pledelDic55])
pledel56reads <- sum(reads[pledelDic56])
pledel57reads <- sum(reads[pledelDic57])
pledel58reads <- sum(reads[pledelDic58])
pledel59reads <- sum(reads[pledelDic59])
pledel60reads <- sum(reads[pledelDic60])
pledel61reads <- sum(reads[pledelDic61])
pledel62reads <- sum(reads[pledelDic62])
pledel63reads <- sum(reads[pledelDic63])
pledel64reads <- sum(reads[pledelDic64])
pledel65reads <- sum(reads[pledelDic65])
pledel66reads <- sum(reads[pledelDic66])
pledel67reads <- sum(reads[pledelDic67])
pledel68reads <- sum(reads[pledelDic68])
pledel69reads <- sum(reads[pledelDic69])
pledel70reads <- sum(reads[pledelDic70])
pledel71reads <- sum(reads[pledelDic71])
pledel72reads <- sum(reads[pledelDic72])
pledel73reads <- sum(reads[pledelDic73])
pledel74reads <- sum(reads[pledelDic74])
pledel75reads <- sum(reads[pledelDic75])
pledel76reads <- sum(reads[pledelDic76])
pledel77reads <- sum(reads[pledelDic77])
pledel78reads <- sum(reads[pledelDic78])
pledel79reads <- sum(reads[pledelDic79])
pledel80reads <- sum(reads[pledelDic80])
pledel81reads <- sum(reads[pledelDic81])
pledel82reads <- sum(reads[pledelDic82])
pledel83reads <- sum(reads[pledelDic83])
pledel84reads <- sum(reads[pledelDic84])
pledel85reads <- sum(reads[pledelDic85])
pledel86reads <- sum(reads[pledelDic86])
pledel87reads <- sum(reads[pledelDic87])
pledel88reads <- sum(reads[pledelDic88])
pledel89reads <- sum(reads[pledelDic89])
pledel90reads <- sum(reads[pledelDic90])
pledel91reads <- sum(reads[pledelDic91])
pledel92reads <- sum(reads[pledelDic92])
pledel93reads <- sum(reads[pledelDic93])
pledel94reads <- sum(reads[pledelDic94])
pledel95reads <- sum(reads[pledelDic95])
pledel96reads <- sum(reads[pledelDic96])
pledel97reads <- sum(reads[pledelDic97])
pledel98reads <- sum(reads[pledelDic98])
pledel99reads <- sum(reads[pledelDic99])
pledel100reads <- sum(reads[pledelDic100])
pledel101reads <- sum(reads[pledelDic101])
pledel102reads <- sum(reads[pledelDic102])
pledel103reads <- sum(reads[pledelDic103])
pledel104reads <- sum(reads[pledelDic104])
pledel105reads <- sum(reads[pledelDic105])
pledel106reads <- sum(reads[pledelDic106])
pledel107reads <- sum(reads[pledelDic107])
pledel108reads <- sum(reads[pledelDic108])
pledel109reads <- sum(reads[pledelDic109])
pledel110reads <- sum(reads[pledelDic110])
pledel111reads <- sum(reads[pledelDic111])
pledel112reads <- sum(reads[pledelDic112])
pledel113reads <- sum(reads[pledelDic113])
pledel114reads <- sum(reads[pledelDic114])
pledel115reads <- sum(reads[pledelDic115])
pledel116reads <- sum(reads[pledelDic116])
pledel117reads <- sum(reads[pledelDic117])
pledel118reads <- sum(reads[pledelDic118])
pledel119reads <- sum(reads[pledelDic119])
pledel120reads <- sum(reads[pledelDic120])
pledel121reads <- sum(reads[pledelDic121])
pledel122reads <- sum(reads[pledelDic122])
pledel123reads <- sum(reads[pledelDic123])
pledel124reads <- sum(reads[pledelDic124])
pledel125reads <- sum(reads[pledelDic125])
pledel126reads <- sum(reads[pledelDic126])
pledel127reads <- sum(reads[pledelDic127])
pledel128reads <- sum(reads[pledelDic128])
pledel129reads <- sum(reads[pledelDic129])
pledel130reads <- sum(reads[pledelDic130])

pledeldicmatreads <- c(pledel1reads, pledel2reads, pledel3reads, pledel4reads, pledel5reads, pledel6reads, pledel7reads, pledel8reads, pledel9reads, pledel10reads,
                       pledel11reads, pledel12reads, pledel13reads, pledel14reads, pledel15reads, pledel16reads, pledel17reads, pledel18reads, pledel19reads, pledel20reads,
                       pledel21reads, pledel22reads, pledel23reads, pledel24reads, pledel25reads, pledel26reads, pledel27reads, pledel28reads, pledel29reads, pledel30reads,
                       pledel31reads, pledel32reads, pledel33reads, pledel34reads, pledel35reads, pledel36reads, pledel37reads, pledel38reads, pledel39reads, pledel40reads,
                       pledel41reads, pledel42reads, pledel43reads, pledel44reads, pledel45reads, pledel46reads, pledel47reads, pledel48reads, pledel49reads, pledel50reads,
                       pledel51reads, pledel52reads, pledel53reads, pledel54reads, pledel55reads, pledel56reads, pledel57reads, pledel58reads, pledel59reads, pledel60reads,
                       pledel61reads, pledel62reads, pledel63reads, pledel64reads, pledel65reads, pledel66reads, pledel67reads, pledel68reads, pledel69reads, pledel70reads, pledel71reads, pledel72reads, pledel73reads, pledel74reads, pledel75reads,
                       pledel76reads, pledel77reads, pledel78reads, pledel79reads, pledel80reads, pledel81reads, pledel82reads, pledel83reads, pledel84reads, pledel85reads, pledel86reads, pledel87reads, pledel88reads, pledel89reads, pledel90reads, pledel91reads,
                       pledel92reads, pledel93reads, pledel94reads, pledel95reads, pledel96reads, pledel97reads, pledel98reads, pledel99reads, pledel100reads, pledel101reads,
                       pledel102reads, pledel103reads, pledel104reads, pledel105reads,  pledel106reads, pledel107reads, pledel108reads, pledel109reads, pledel110reads, pledel111reads,
                       pledel112reads, pledel113reads, pledel114reads, pledel115reads, pledel116reads, pledel117reads, pledel118reads, pledel119reads, pledel120reads, pledel121reads,
                       pledel122reads, pledel123reads, pledel124reads, pledel125reads, pledel126reads, pledel127reads, pledel128reads, pledel129reads, pledel130reads)

## Getting Insertion class alleles
length(c(plewtvec, plepeppervec, pledelvec, plemjavec))
plevec <- c(plewtvec, plepeppervec, pledelvec, plemjavec)
unique(sort(plevec))
pledna <- totblstfr$top[plevec]
allele <- c(1:180)
plednavec <- intersect(plevec, allele)
nonplednavec <- allele[!allele %in% plednavec]
pleinsvec <- nonplednavec

## Check intersect
intersect(plepeppervec, plemjavec)
intersect(plepeppervec, pledelvec)
intersect(plepeppervec, pleinsvec)
intersect(plemjavec, pledelvec)
intersect(plemjavec, pleinsvec)
intersect(pledelvec, pleinsvec)

## Making a new dataframe
ple_DEL_l_Dic <- readRDS("ple_DEL_l_Dic.Rds")
ple_MMEJ_l_Dic <- readRDS("ple_MMEJ_l_Dic.Rds")
ple_DEL_l_Dic <- readRDS("ple_DEL_l_Dic.Rds")
wt_l <- readRDS("wt_L.Rds")
length(c(plepedicmatvec, plemjadicmatvec, pledeldicmatvec, pleinsvec))
pledna <- c(ple_PEPPER_l_Dic[plepedicmatvec], ple_MMEJ_l_Dic[plemjadicmatvec], ple_DEL_l_Dic[pledeldicmatvec], top180seqstr[pleinsvec])
plepedicmatreads <- plepedicmatreads[!plepedicmatreads %in% (plepedicmatreads == 0)]
plemjadicmatreads <- plemjadicmatreads[!plemjadicmatreads %in% (plemjadicmatreads == 0)]
pledeldicmatreads <- pledeldicmatreads[!pledeldicmatreads %in% (pledeldicmatreads == 0)]
plednareads <- c(plepedicmatreads, plemjadicmatreads, pledeldicmatreads, reads[pleinsvec])
pctnumvecn <-plednareads/totalreads*100
pctchavecn <- as.character(round(pctnumvecn, digits = 2))
length(plepedicmatvec)
length(plemjadicmatvec)
length(pledeldicmatvec)
length(top180seqstr[pleinsvec])
length(plepedicmatreads)
length(plemjadicmatreads)
length(pledeldicmatreads)
length(reads[pleinsvec])

inserttvecn <- rep(NA, 178)
zsumclassvecn <- replace(inserttvecn, c(1:67), "peppr")
zsumclassvecn <- replace(zsumclassvecn, c(68:78), "charmmj")
zsumclassvecn <- replace(zsumclassvecn, c(79:112), "delet")
zsumclassvecn <- replace(zsumclassvecn, c(113:178), "insrt")
frindvecn <- c(1:178)
zclassnewdfrn <- data.frame(frindvecn, zsumclassvecn, pctchavecn, plednareads)
names(zclassnewdfrn) <- c("rank", "classify", "reads ratio", "reads")
zclassnewdfrn$classify <- factor(zclassnewdfrn$classify,
                                levels = c("insrt",
                                           "delet",
                                           "charmmj",
                                           "peppr"))
readspct <- order(zclassnewdfrn$reads, decreasing = TRUE)
zclass <- zclassnewdfrn[readspct, ]
zclass$alleles <- c(1:178)
zclass
saveRDS(zclass, "zclass.Rds")

#### Step 4: Data Visualization
## DSB pattern fingerprints
library(ggplot2)
readbreaks <- c("wt", "peppr", "charmmj", "delet", "insrt")
readsccols <- c("deeppink", "deepskyblue", "deeppink",
                "goldenrod1", "purple")
zclassdfrnplot <- zclass[c(1:50), ]
zclassdfrnplot$allelenew <- c(1:50)
ggplot(data = zclassdfrnplot,
       mapping = aes (x = allelenew, y = classify, fill = classify))+
  geom_tile (size = 6, width = 0.8, height = 1)+
  theme_bw()+
  labs(title = "pattern change between indel classes", 
       x = "Rank of alleles",
       y = "Classify")+
  scale_fill_manual(values = readsccols,
                    breaks = readbreaks,
                    name = "class")+
  theme(aspect.ratio = 0.12,
        panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank())

## Allele landscape plotting
library(scales)
library(ggplot2)
zclass <- readRDS("zclass.Rds")
plesumprofile <- zclass[c(1:50), ]
ggplot(data = plesumprofile) +
  geom_segment(mapping = aes(x = c(1:50),
                          xend = c(1:50),
                          y = 0,
                          yend = reads), 
            color = "skyblue", alpha = 0.8,
            size = 2.0) +
  geom_point(mapping = aes(x = c(1:50),
                           y = reads), 
             color = "red2", 
             size = 0.8) +
  scale_x_continuous(breaks = seq(0, 50, 5)) +
  scale_y_log10(breaks = trans_breaks("log10",
                                      n = 10,
                                      function(x) 10^x)(c(1e1, 1e10)),
                labels = trans_format("log10", math_format(10^.x)))+
  labs(title = "readsrank",
       y = "Number of Reads",
       x = "Rank Index of Specific Sequences") +
  theme(aspect.ratio = 0.35,
        panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank(),
        axis.title.x = element_blank()) 

