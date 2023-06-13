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

#### Step 3: Create deletion skeleton datafr
## get main datafrs for plot  =====
grnumseqvec1 <- names(totblstfr$top[1:180])
ttalsstrpat <- DNAStringSet(grnumseqvec1)
ttalnssub <- DNAString(drplechar)
cbgrpsltaln <- pairwiseAlignment(ttalsstrpat, ttalnssub)
cbmismtdfr <- mismatchTable(cbgrpsltaln)
cbinsrtdfr <- as.data.frame(insertion(cbgrpsltaln))
cbdeletdfr <- as.data.frame(deletion(cbgrpsltaln))
shorskuhldfr <- cbdeletdfr[cbdeletdfr$start < 90, ]
shorskuhldfr$group_name <- NULL

## create deletion skeleton datafr  =====
litdbfvec <- shorskuhldfr[duplicated(shorskuhldfr$group), 1]
wdelsetvec <- unique(litdbfvec)
shorskuhldfr$seldel <- shorskuhldfr$group %in% wdelsetvec
wierdelmdfr <- shorskuhldfr[shorskuhldfr$seldel == TRUE, ]
skelskuhldfr <- shorskuhldfr[shorskuhldfr$seldel == FALSE, ]

## create merged delsdatafr from wierdels  =====
wwidsumval1 <- tapply(wierdelmdfr$width, wierdelmdfr$group, sum)
wstaminval <- tapply(wierdelmdfr$start, wierdelmdfr$group, min)
wstamaxval <- tapply(wierdelmdfr$start, wierdelmdfr$group, max)
wwidval2 <- wstamaxval - wstaminval
wwidtotval <- wwidsumval1 + wwidval2
wcalendval <- wstaminval + wwidtotval - 1
wrddelnewdfr <- data.frame(wdelsetvec, wstaminval, wcalendval, wwidtotval)
names(wrddelnewdfr) <- c("group", "start", "end", "width")

## create zero deletion set datafr  =====
bgcompst <- c(1:180)
zkdelsetvec <- setdiff(bgcompst, shorskuhldfr$group)
cbinsrtdfr$group %in% zkdelsetvec
kaddstavec <- rep(NA, length(zkdelsetvec))
kaddendvec <- rep(NA, length(zkdelsetvec))
kaddwidvec <- rep(0, length(zkdelsetvec))
zkdelnewdfr <- data.frame(zkdelsetvec, kaddstavec, kaddendvec, kaddwidvec)
names(zkdelnewdfr) <- c("group", "start", "end", "width")

## create complete deletion set datafr  =====
skelskuhldfr$seldel <- NULL
comskuhldfr <- rbind(skelskuhldfr, wrddelnewdfr, zkdelnewdfr)
extskuhldfr <- comskuhldfr[order(comskuhldfr$group), ]
sum(extskuhldfr$group == bgcompst)
frindvec <- c(1:180)
frvalvec <- unname(totblstfr$top[1:180])
extskuhldfr$sknewind <- frindvec
extskuhldfr$sknewval <- frvalvec
mixpoisvec <- c(cbinsrtdfr$group, cbmismtdfr[ ,1])
refpoivec <- unique(mixpoisvec)
inserttvec <- rep(NA, 180)
inserttvec <- replace(inserttvec, refpoivec, refpoivec)
extskuhldfr$inserrpois <- inserttvec

## Save Rds and reopen csv
write.csv(extskuhldfr, "TOP180.csv")
saveRDS(extskuhldfr, "TOP180.Rds")
extskuhldfr <- readRDS("TOP180.Rds")


## create plot prep objects  =====
skuhldfr <- extskuhldfr[1:50, ]
rdgrsum <- sum(frvalvec[1:50])
rdgrpct <- round(rdgrsum/denomval*100, digits = 1)
anteegrop <- grid.text(rdgrpct, x = 0.86, y = 0.98,
                       gp = gpar(col = "red3",
                                 fontsize = 12,
                                 fontface = "bold"))
antsgnpop <- grid.text("%", x = 0.94, y = 0.98,
                       gp = gpar(fontsize = 10,
                                 fontface = "bold"))

## create cool ggplot imagge  =====
ggplot(data = skuhldfr) +
  geom_point(mapping = aes(x = 2, y = sknewind, size = sknewval),
             color = "purple2") +
  geom_segment(mapping = aes(x = 8, xend = 155,
                             y = sknewind, yend = sknewind),
               color = "skyblue",
               size = 0.8) +
  geom_segment(mapping = aes(x = start, xend = end,
                             y = sknewind, yend = sknewind),
               na.rm = TRUE,
               color = "black",
               size = 2.0) +
  geom_vline(xintercept = 66.5,
             color = "red2",
             size = 0.5) +
  geom_point(mapping = aes(x = 68, y = inserrpois),
             na.rm = TRUE,
             shape = 23,
             size = 2.4,
             color = "black",
             fill = "chartreuse") +
  geom_label(mapping = aes(x = 155, y = sknewind, label = width),
             size = 2,
             label.padding = unit(0.12, "lines"),
             fontface = "bold") +
  scale_x_continuous(breaks = seq(0, 155, 30)) +
  scale_y_reverse(breaks = seq(0, 50, 5),
                  limits = c(50, 0)) +
  scale_size_area(breaks = seq(1000, 500000, 100000),
                  name = "n Reads",
                  limits = c(1000, 500000)) +
  labs(title = "ple-A-F-F1F",
       subtitle = "summary of lesions: deletions 1-50",
       y = "Index of Sequence Variants",
       x = "Nucleotide Position within Amplicon") +
  theme(aspect.ratio = 1.5,
        panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank(),
        panel.grid.major.y = element_blank(),
        panel.grid.minor.y = element_blank()) +
  annotation_custom(anteegrop) +
  annotation_custom(antsgnpop)

#### Step 4: Classification
## prep selection vecs for peppr class dfrs =====
pepskelvec <- extskuhldfr[extskuhldfr$start == 67 & extskuhldfr$width > 3, 1]
pepoff1vec <- extskuhldfr[extskuhldfr$start == 68 & extskuhldfr$width > 1, 1]
pepoff2vec <- extskuhldfr[extskuhldfr$start == 69 & extskuhldfr$width > 1, 1]
pepoff3vec <- extskuhldfr[extskuhldfr$start == 66 & extskuhldfr$width > 1, 1]
pepoff4vec <- extskuhldfr[extskuhldfr$start == 65 & extskuhldfr$width > 1, 1]

## combine selection vecs for visualize classmuts =====
extskuhldfr[pepskelvec, "inserrpois"]
pepplusvec <- c(pepoff1vec, pepoff2vec)
pepplusvec <- sort(pepplusvec)
pepmnusvec <- c(pepoff3vec, pepoff4vec)
pepmnusvec <- sort(pepmnusvec)

## now tetris matrix up the peppr muts =====
library(plot.matrix)
length(pepplusvec)

## get pepprplus muts for viewing  =====
rnumvec <- as.character(pepplusvec)
pepgrpsltaln <- cbgrpsltaln[pepplusvec]
currpepmatr <- as.matrix(pepgrpsltaln)
rownames(currpepmatr) <- rnumvec
pleref <- unlist(strsplit(drplechar, ""))
blnkrwn <- 30 - length(pepplusvec)
blnkmdat <- rep("-", blnkrwn*155)
blnkmmatrx <- matrix(blnkmdat, nrow = blnkrwn, ncol = 155, byrow = TRUE)
rownames(blnkmmatrx) <- rep("-", blnkrwn)

## get first matrix object and plot  =====
currmatr <- rbind(pleref, currpepmatr, pleref, blnkmmatrx)
pepmatrlab <- paste("Maternal-F", "align.pepplus", sep = "")
par(mar = c(3.1, 4.1, 4.1, 1.6), xpd = TRUE, adj = 0.9)
plot(currmatr,
     main = pepmatrlab,
     col = c("red2", "turquoise1", "navyblue", "navajowhite2", "white"),
     breaks = c("G", "A", "T", "C", "-"),
     border = NA,
     xlab = "",
     ylab = "",
     axis.row = list(las = 1, cex.axis = 0.6, font.axis = 2),
     key = NULL)
polygon(x = c(64.5, 68.5, 66.5),
        y = c(34.5, 34.5, 33),
        col = "black")
mtext(c("G", "A", "T", "C"),
      side = 3,
      line = 2,
      at = c(15, 25, 35, 45),
      font = 2,
      cex = 1.2)
polygon(x = c(10, 20, 20, 10), y = c(36, 36, 34.5, 34.5),
        col = "red2")
polygon(x = c(20, 30, 30, 20), y = c(36, 36, 34.5, 34.5),
        col = "turquoise1")
polygon(x = c(30, 40, 40, 30), y = c(36, 36, 34.5, 34.5),
        col = "navyblue")
polygon(x = c(40, 50, 50, 40), y = c(36, 36, 34.5, 34.5),
        col = "navajowhite2")



### MMEJ
## do selection ops for charmmej dfrs =====
mjaval1 <- extskuhldfr[extskuhldfr$start == 64 & extskuhldfr$end == 69 
                       & extskuhldfr$width > 1, 1]
mjaval2 <- extskuhldfr[extskuhldfr$start == 56 & extskuhldfr$end == 66 
                       & extskuhldfr$width > 1, 1]
mjaval3 <- extskuhldfr[extskuhldfr$start == 62 & extskuhldfr$end == 73 
                       & extskuhldfr$width > 1, 1]
mjaval4 <- extskuhldfr[extskuhldfr$start == 69 & extskuhldfr$end == 98 
                       & extskuhldfr$width > 1, 1]
mjaval5 <- extskuhldfr[extskuhldfr$start == 69 & extskuhldfr$end == 92 
                       & extskuhldfr$width > 1, 1]
mjaval6 <- extskuhldfr[extskuhldfr$start == 67 & extskuhldfr$end == 110 
                       & extskuhldfr$width > 1, 1]
mjaval7 <- extskuhldfr[extskuhldfr$start == 61 & extskuhldfr$end == 66 
                       & extskuhldfr$width > 1, 1]
mjaval8 <- extskuhldfr[extskuhldfr$start == 67 & extskuhldfr$end == 114 
                       & extskuhldfr$width > 0, 1]
mjaval9 <- extskuhldfr[extskuhldfr$start == 68 & extskuhldfr$end == 82 
                       & extskuhldfr$width > 0, 1]
mjaval10 <- extskuhldfr[extskuhldfr$start == 66 & extskuhldfr$end == 98 
                        & extskuhldfr$width > 0, 1]
mjaval11 <- extskuhldfr[extskuhldfr$start == 35 & extskuhldfr$end == 96 
                        & extskuhldfr$width > 0, 1]
mjaval12 <- extskuhldfr[extskuhldfr$start == 58 & extskuhldfr$end == 69 
                        & extskuhldfr$width > 0, 1]
mjaval13 <- extskuhldfr[extskuhldfr$start == 20 & extskuhldfr$end == 71 
                        & extskuhldfr$width > 0, 1]
mjaval14 <- extskuhldfr[extskuhldfr$start == 54 & extskuhldfr$end == 116 
                        & extskuhldfr$width > 0, 1]
mjaval15 <- extskuhldfr[extskuhldfr$start == 63 & extskuhldfr$end == 125 
                        & extskuhldfr$width > 0, 1]
mjaval16 <- extskuhldfr[extskuhldfr$start == 62 & extskuhldfr$end == 103 
                        & extskuhldfr$width > 1, 1]
mjaval17 <- extskuhldfr[extskuhldfr$start == 66 & extskuhldfr$end == 125 
                        & extskuhldfr$width > 1, 1]
mjaval18 <- extskuhldfr[extskuhldfr$start == 63 & extskuhldfr$end == 122 
                        & extskuhldfr$width > 1, 1]
mjaval19 <- extskuhldfr[extskuhldfr$start == 44 & extskuhldfr$end == 67 
                        & extskuhldfr$width > 1, 1]
mjaval20 <- extskuhldfr[extskuhldfr$start == 52 & extskuhldfr$end == 68 
                        & extskuhldfr$width > 0, 1]
mjaval21 <- extskuhldfr[extskuhldfr$start == 37 & extskuhldfr$end == 68 
                        & extskuhldfr$width > 0, 1]
mjaval22 <- extskuhldfr[extskuhldfr$start == 58 & extskuhldfr$end == 67 
                        & extskuhldfr$width > 0, 1]
mjaval23 <- extskuhldfr[extskuhldfr$start == 62 & extskuhldfr$end == 103 
                        & extskuhldfr$width > 1, 1]
mjaval24 <- extskuhldfr[extskuhldfr$start == 64 & extskuhldfr$end == 125 
                        & extskuhldfr$width > 1, 1]
mjaval25 <- extskuhldfr[extskuhldfr$start == 58 & extskuhldfr$end == 92 
                        & extskuhldfr$width > 0, 1]
mjaval26 <- extskuhldfr[extskuhldfr$start == 40 & extskuhldfr$end == 68 
                        & extskuhldfr$width > 0, 1]
mjaval27 <- extskuhldfr[extskuhldfr$start == 47 & extskuhldfr$end == 68 
                        & extskuhldfr$width > 0, 1]
mjaval28 <- extskuhldfr[extskuhldfr$start == 47 & extskuhldfr$end == 87 
                        & extskuhldfr$width > 1, 1]
mjaval29 <- extskuhldfr[extskuhldfr$start == 46 & extskuhldfr$end == 90 
                        & extskuhldfr$width > 1, 1]
mjaval30 <- extskuhldfr[extskuhldfr$start == 41 & extskuhldfr$end == 88 
                        & extskuhldfr$width > 1, 1]
mjaval31 <- extskuhldfr[extskuhldfr$start == 62 & extskuhldfr$end == 116 
                        & extskuhldfr$width > 1, 1]

## modify for charmmej duplicates =====
mjaval1 <- mjaval1[2]
mjaval20 <- mjaval20[1]


## get pepej muts for viewing  =====
mjacharvec <- c(mjaval1, mjaval2, mjaval3, mjaval4,
                mjaval5, mjaval6, mjaval7, mjaval8,
                mjaval9, mjaval10, mjaval11, mjaval12,
                mjaval13, mjaval14, mjaval15, mjaval16,
                mjaval17, mjaval18, mjaval19, mjaval20,
                mjaval21, mjaval22, mjaval23, mjaval24,
                mjaval25, mjaval26, mjaval27, mjaval28,
                mjaval29, mjaval30, mjaval31)
rnumvec <- as.character(mjacharvec)
mjachagrpsltaln <- cbgrpsltaln[mjacharvec]
currmjamatr <- as.matrix(mjachagrpsltaln)
rownames(currmjamatr) <- rnumvec
pleref <- unlist(strsplit(drplechar, ""))

## get first matrix object and plot  =====
currmatr <- rbind(pleref, currmjamatr, pleref)
mjamatrlab <- paste("ple-A-M", "class.mjchar", sep = "")
par(mar = c(3.1, 4.1, 4.1, 1.6), xpd = TRUE, adj = 0.9)
plot(currmatr,
     main = mjamatrlab,
     col = c("red2", "turquoise1", "navyblue", "navajowhite2", "white"),
     breaks = c("G", "A", "T", "C", "-"),
     border = NA,
     xlab = "",
     ylab = "",
     axis.row = list(las = 1, cex.axis = 0.8, font.axis = 2),
     key = NULL)
polygon(x = c(65, 68, 66.5),
        y = c(34.5, 34.5, 33.5),
        col = "black")
mtext(c("G", "A", "T", "C"),
      side = 3,
      line = 2,
      at = c(15, 25, 35, 45),
      font = 2,
      cex = 1.2)
polygon(x = c(10, 20, 20, 10), y = c(36, 36, 35.3, 35.3),
        col = "red2")
polygon(x = c(20, 30, 30, 20), y = c(36, 36, 35.3, 35.3),
        col = "turquoise1")
polygon(x = c(30, 40, 40, 30), y = c(36, 36, 35.3, 35.3),
        col = "navyblue")
polygon(x = c(40, 50, 50, 40), y = c(36, 36, 35.3, 35.3),
        col = "navajowhite2")


### Clean up PEPPER
## fix up pepskel classify list  =====
extskuhldfr[pepskelvec, "inserrpois"]
intersect(pepskelvec, wdelsetvec)
remskellvec <- extskuhldfr[pepskelvec, "inserrpois"]
remskellvec <- remskellvec[! is.na(remskellvec)]
pepskellclvec <- pepskelvec[! pepskelvec %in% remskellvec]

## make clean pep plus and minus classify list  =====
pepkaddvec1 <- pepplusvec
pepkaddvec2 <- pepmnusvec

## make exact pepmut classify list  ===== 
pepktotvec1 <- c(pepskellclvec, pepkaddvec1, pepkaddvec2)
pepktotvec1 <- sort(pepktotvec1)
pepktotvec2 <- intersect(pepktotvec1, mjacharvec)
pepktotvec <- pepktotvec1[!pepktotvec1 %in% pepktotvec2]

## prep second round classify with recycled dfr  =====
rnd1classvec <- c(pepktotvec, mjacharvec)
rnd1classvec <- sort(rnd1classvec)
rnd1clnegvec <- setdiff(bgcompst, rnd1classvec)
secrndkuhldfr <- extskuhldfr[rnd1clnegvec, ]

## prep selection vecs for peppr class dfrs =====
pepvarsecvec <- secrndkuhldfr[(secrndkuhldfr$start == 67 
                               | secrndkuhldfr$start == 68 
                               | secrndkuhldfr$start == 69 
                               | secrndkuhldfr$start == 66 
                               | secrndkuhldfr$start == 65) 
                              & secrndkuhldfr$width > 4, 1]
pepvarsecvec <- sort(pepvarsecvec)

## now tetris matrix up the varpeppr muts =====
library(plot.matrix)
length(pepvarsecvec)

## get varpepej muts for viewing  =====
rnumvec <- as.character(pepvarsecvec)
pepgrpsltaln <- cbgrpsltaln[pepvarsecvec]
currpepmatr <- as.matrix(pepgrpsltaln)
rownames(currpepmatr) <- rnumvec
pleref <- unlist(strsplit(drplechar, ""))
blnkrwn <- 30 - length(pepvarsecvec)
blnkmdat <- rep("-", blnkrwn*155)
blnkmmatrx <- matrix(blnkmdat, nrow = blnkrwn, ncol = 155, byrow = TRUE)
rownames(blnkmmatrx) <- rep("-", blnkrwn)

## get first matrix object and plot  =====
currmatr <- rbind(pleref, currpepmatr, pleref, blnkmmatrx)
pepmatrlab <- paste("ple-A-M", "alngrp1.pepvars", sep = "")
par(mar = c(3.1, 4.1, 4.1, 1.6), xpd = TRUE, adj = 0.9)
plot(currmatr,
     main = pepmatrlab,
     col = c("red2", "turquoise1", "navyblue", "navajowhite2", "white"),
     breaks = c("G", "A", "T", "C", "-"),
     border = NA,
     xlab = "",
     ylab = "",
     axis.row = list(las = 1, cex.axis = 0.6, font.axis = 2),
     key = NULL)
polygon(x = c(64.5, 68.5, 66.5),
        y = c(34.5, 34.5, 33),
        col = "black")
mtext(c("G", "A", "T", "C"),
      side = 3,
      line = 2,
      at = c(15, 25, 35, 45),
      font = 2,
      cex = 1.2)
polygon(x = c(10, 20, 20, 10), y = c(36, 36, 34.5, 34.5),
        col = "red2")
polygon(x = c(20, 30, 30, 20), y = c(36, 36, 34.5, 34.5),
        col = "turquoise1")
polygon(x = c(30, 40, 40, 30), y = c(36, 36, 34.5, 34.5),
        col = "navyblue")
polygon(x = c(40, 50, 50, 40), y = c(36, 36, 34.5, 34.5),
        col = "navajowhite2")

### PEPPER Varients
## prep selection vecs for other dels class dfrs =====
widelvarsvec <- secrndkuhldfr[! secrndkuhldfr$group %in% pepvarsecvec 
                              & (secrndkuhldfr$start > 69 
                                 | secrndkuhldfr$start < 65) 
                              & secrndkuhldfr$width > 8, 1]
widelvarsvec <- sort(widelvarsvec)

## now tetris matrix up the widevar muts =====
library(plot.matrix)
length(widelvarsvec)

## get widevar muts for viewing  =====
rnumvec <- as.character(widelvarsvec)
pepgrpsltaln <- cbgrpsltaln[widelvarsvec]
currpepmatr <- as.matrix(pepgrpsltaln)
rownames(currpepmatr) <- rnumvec
pleref <- unlist(strsplit(drplechar, ""))
blnkrwn <- 50 - length(widelvarsvec)
blnkmdat <- rep("-", blnkrwn*155)
blnkmmatrx <- matrix(blnkmdat, nrow = blnkrwn, ncol = 155, byrow = TRUE)
rownames(blnkmmatrx) <- rep("-", blnkrwn)

## get first matrix object and plot  =====
currmatr <- rbind(pleref, currpepmatr, pleref, blnkmmatrx)
pepmatrlab <- paste("ple-A-M", "alngrp1.diffdels", sep = "")
par(mar = c(3.1, 4.1, 4.1, 1.6), xpd = TRUE, adj = 0.9)
plot(currmatr,
     main = pepmatrlab,
     col = c("red2", "turquoise1", "navyblue", "navajowhite2", "white"),
     breaks = c("G", "A", "T", "C", "-"),
     border = NA,
     xlab = "",
     ylab = "",
     axis.row = list(las = 1, cex.axis = 0.6, font.axis = 2),
     key = NULL)
polygon(x = c(64.5, 68.5, 66.5),
        y = c(54.5, 54.5, 53),
        col = "black")
mtext(c("G", "A", "T", "C"),
      side = 3,
      line = 2,
      at = c(15, 25, 35, 45),
      font = 2,
      cex = 1.2)
polygon(x = c(10, 20, 20, 10), y = c(56, 56, 54.5, 54.5),
        col = "red2")
polygon(x = c(20, 30, 30, 20), y = c(56, 56, 54.5, 54.5),
        col = "turquoise1")
polygon(x = c(30, 40, 40, 30), y = c(56, 56, 54.5, 54.5),
        col = "navyblue")
polygon(x = c(40, 50, 50, 40), y = c(56, 56, 54.5, 54.5),
        col = "navajowhite2")

### Other deletions
## prep selection vecs for mini dels class dfrs =====
mindelvarsvec <- secrndkuhldfr[! secrndkuhldfr$group %in% pepvarsecvec 
                               & (! secrndkuhldfr$group %in% widelvarsvec) 
                               & secrndkuhldfr$width > 0, 1]
mindelvarsvec <- sort(mindelvarsvec)

## now tetris matrix up the minidel muts =====
library(plot.matrix)
length(mindelvarsvec)

## get minidel muts for viewing  =====
rnumvec <- as.character(mindelvarsvec)
pepgrpsltaln <- cbgrpsltaln[mindelvarsvec]
currpepmatr <- as.matrix(pepgrpsltaln)
rownames(currpepmatr) <- rnumvec
pleref <- unlist(strsplit(drplechar, ""))
blnkrwn <- 30 - length(mindelvarsvec)
blnkmdat <- rep("-", blnkrwn*155)
blnkmmatrx <- matrix(blnkmdat, nrow = blnkrwn, ncol = 155, byrow = TRUE)
rownames(blnkmmatrx) <- rep("-", blnkrwn)

## get first matrix object and plot  =====
currmatr <- rbind(pleref, currpepmatr, pleref, blnkmmatrx)
pepmatrlab <- paste("p[le-A-M", "alngrp1.minidels", sep = "")
par(mar = c(3.1, 4.1, 4.1, 1.6), xpd = TRUE, adj = 0.9)
plot(currmatr,
     main = pepmatrlab,
     col = c("red2", "turquoise1", "navyblue", "navajowhite2", "white"),
     breaks = c("G", "A", "T", "C", "-"),
     border = NA,
     xlab = "",
     ylab = "",
     axis.row = list(las = 1, cex.axis = 0.6, font.axis = 2),
     key = NULL)
polygon(x = c(64.5, 68.5, 66.5),
        y = c(34.5, 34.5, 33),
        col = "black")
mtext(c("G", "A", "T", "C"),
      side = 3,
      line = 2,
      at = c(15, 25, 35, 45),
      font = 2,
      cex = 1.2)
polygon(x = c(10, 20, 20, 10), y = c(36, 36, 34.5, 34.5),
        col = "red2")
polygon(x = c(20, 30, 30, 20), y = c(36, 36, 34.5, 34.5),
        col = "turquoise1")
polygon(x = c(30, 40, 40, 30), y = c(36, 36, 34.5, 34.5),
        col = "navyblue")
polygon(x = c(40, 50, 50, 40), y = c(36, 36, 34.5, 34.5),
        col = "navajowhite2")

### Insertions
## prep selection vecs for pure insrts class dfrs =====
purinsvarsvec <- secrndkuhldfr[! secrndkuhldfr$group %in% pepvarsecvec 
                               & (! secrndkuhldfr$group %in% widelvarsvec) 
                               & secrndkuhldfr$width == 0, 1]
purinsvarsvec <- sort(purinsvarsvec)

## now tetris matrix up the pinsrt muts =====
library(plot.matrix)
length(purinsvarsvec)

## get pinsrt muts for viewing  =====
rnumvec <- as.character(purinsvarsvec)
pepgrpsltaln <- cbgrpsltaln[purinsvarsvec]
currpepmatr <- as.matrix(pepgrpsltaln)
rownames(currpepmatr) <- rnumvec
pleref <- unlist(strsplit(drplechar, ""))
blnkrwn <- 30 - length(purinsvarsvec)
blnkmdat <- rep("-", blnkrwn*155)
blnkmmatrx <- matrix(blnkmdat, nrow = blnkrwn, ncol = 155, byrow = TRUE)
rownames(blnkmmatrx) <- rep("-", blnkrwn)

## get first matrix object and plot  =====
currmatr <- rbind(pleref, currpepmatr, pleref, blnkmmatrx)
pepmatrlab <- paste("ple-A-M", "alngrp1.nogapins", sep = "")
par(mar = c(3.1, 4.1, 4.1, 1.6), xpd = TRUE, adj = 0.9)
plot(currmatr,
     main = pepmatrlab,
     col = c("red2", "turquoise1", "navyblue", "navajowhite2", "white"),
     breaks = c("G", "A", "T", "C", "-"),
     border = NA,
     xlab = "",
     ylab = "",
     axis.row = list(las = 1, cex.axis = 0.6, font.axis = 2),
     key = NULL)
polygon(x = c(64.5, 68.5, 66.5),
        y = c(34.5, 34.5, 33),
        col = "black")
mtext(c("G", "A", "T", "C"),
      side = 3,
      line = 2,
      at = c(15, 25, 35, 45),
      font = 2,
      cex = 1.2)
polygon(x = c(10, 20, 20, 10), y = c(35, 35, 34, 34),
        col = "red2")
polygon(x = c(20, 30, 30, 20), y = c(35, 35, 34, 34),
        col = "turquoise1")
polygon(x = c(30, 40, 40, 30), y = c(35, 35, 34, 34),
        col = "navyblue")
polygon(x = c(40, 50, 50, 40), y = c(35, 35, 34, 34),
        col = "navajowhite2")

### Other views for PEPPER
## create plot prep objects  =====
modwidskuhldfr <- extskuhldfr
modwidskuhldfr[ , "sknewval"] == frvalvec
rebreaksvec <- unname(quantile(frvalvec))
septqread <- .bincode(frvalvec, rebreaksvec, include.lowest = TRUE)
modwidskuhldfr$septqtile <- as.factor(septqread)
pepwidskuhldfr <- modwidskuhldfr[pepktotvec, c(1, 4:6, 8)]
rdglabpct <- round(pepwidskuhldfr[ , "sknewval"]/denomval*100, digits = 2)
pepwidskuhldfr$rdpct <- rdglabpct
npepindval <- c(1:length(pepktotvec))
pepwidskuhldfr$newpepind <- npepindval
evevalvec <- npepindval[npepindval %% 2 == 0]
oddvalvec <- npepindval[npepindval %% 2 != 0]
inserttvec <- rep(NA, length(pepktotvec))
eveindvec <- replace(inserttvec, evevalvec, evevalvec)
pepwidskuhldfr$evewidind <- eveindvec
inserttvec <- rep(NA, length(pepktotvec))
oddindvec <- replace(inserttvec, oddvalvec, oddvalvec)
pepwidskuhldfr$oddwidind <- oddindvec
rdgrsum <- sum(pepwidskuhldfr$sknewval)
rdgrpct <- round(rdgrsum/denomval*100, digits = 1)
anteegrop <- grid.text(rdgrpct, x = 0.94, y = 0.98,
                       gp = gpar(col = "grey15",
                                 fontsize = 12,
                                 fontface = "bold"))
antsgnpop <- grid.text("%", x = 0.98, y = 0.98,
                       gp = gpar(fontsize = 10,
                                 fontface = "bold"))
readbreaks <- as.character(c(1:4))
readsccols <- c("blue2", "deepskyblue", "coral1", "darkred")

## create cool ggplot imagge  =====
ggplot(data = pepwidskuhldfr) +
  geom_segment(mapping = aes(x = 77, xend = 82.5,
                             y = newpepind, yend = newpepind,
                             color = septqtile),
               size = 1.6) +
  geom_segment(mapping = aes(x = 0, xend = width,
                             y = newpepind, yend = newpepind),
               na.rm = TRUE,
               color = "black",
               size = 2) +
  geom_vline(xintercept = 0,
             color = "red2",
             size = 0.8) +
  geom_label(mapping = aes(x = 75, y = oddwidind, label = rdpct),
             na.rm = TRUE,
             size = 2,
             label.padding = unit(0.12, "lines"),
             fontface = "bold") +
  geom_label(mapping = aes(x = 78, y = evewidind, label = rdpct),
             na.rm = TRUE,
             size = 2,
             label.padding = unit(0.12, "lines"),
             fontface = "bold") +
  scale_x_continuous(breaks = seq(0, 80, 5),
                     limits = c(0, 83)) +
  scale_y_reverse(breaks = seq(0, 64, 4),
                  limits = c(65, 0)) +
  scale_color_manual(values = readsccols,
                     breaks = readbreaks,
                     guide = guide_legend(reverse = TRUE),
                     name = "Q Read") +
  labs(title = "ple-A-M",
       subtitle = "align group 1.180 : PEPPR deletions 1-65",
       y = "Rank Index of Read Quantity",
       x = "Nucleotide Position Distal to Cut-site") +
  theme(aspect.ratio = 0.6,
        axis.text.y = element_text(size = 8),
        panel.grid.major.x = element_line(size = 0.3),
        panel.grid.minor.x = element_blank(),
        panel.grid.major.y = element_line(color = "grey75", size = 0.5),
        panel.grid.minor.y = element_line(color = "grey85", size = 0.5)) +
  annotation_custom(anteegrop) +
  annotation_custom(antsgnpop)

## create plot prep objects  =====
modwidskuhldfr <- extskuhldfr
modwidskuhldfr[ , "sknewval"] == frvalvec
rebreaksvec <- unname(quantile(frvalvec))
septqread <- .bincode(frvalvec, rebreaksvec, include.lowest = TRUE)
modwidskuhldfr$septqtile <- as.factor(septqread)
pepwidskuhldfr <- modwidskuhldfr[pepktotvec, c(1, 4:6, 8)]
rdglabpct <- round(pepwidskuhldfr[ , "sknewval"]/denomval*100, digits = 2)
pepwidskuhldfr$rdpct <- rdglabpct
pordwidskuhldfr <- pepwidskuhldfr[order(pepwidskuhldfr$width), ]
npepindval <- c(1:length(pepktotvec))
pordwidskuhldfr$newpepind <- rev(npepindval)
evevalvec <- npepindval[npepindval %% 2 == 0]
oddvalvec <- npepindval[npepindval %% 2 != 0]
inserttvec <- rep(NA, length(pepktotvec))
eveindvec <- replace(inserttvec, evevalvec, evevalvec)
pordwidskuhldfr$evewidind <- rev(eveindvec)
inserttvec <- rep(NA, length(pepktotvec))
oddindvec <- replace(inserttvec, oddvalvec, oddvalvec)
pordwidskuhldfr$oddwidind <- rev(oddindvec)
rdgrsum <- sum(pordwidskuhldfr$sknewval)
rdgrpct <- round(rdgrsum/denomval*100, digits = 1)
anteegrop <- grid.text(rdgrpct, x = 0.94, y = 0.98,
                       gp = gpar(col = "grey15",
                                 fontsize = 12,
                                 fontface = "bold"))
antsgnpop <- grid.text("%", x = 0.98, y = 0.98,
                       gp = gpar(fontsize = 10,
                                 fontface = "bold"))
readbreaks <- as.character(c(1:4))
readsccols <- c("blue2", "deepskyblue", "coral1", "darkred")

## create cool ggplot imagge  =====
ggplot(data = pordwidskuhldfr) +
  geom_segment(mapping = aes(x = 77, xend = 82.5,
                             y = newpepind, yend = newpepind,
                             color = septqtile),
               size = 1.6) +
  geom_segment(mapping = aes(x = 0, xend = width,
                             y = newpepind, yend = newpepind),
               na.rm = TRUE,
               color = "black",
               size = 2) +
  geom_vline(xintercept = 0,
             color = "red2",
             size = 0.8) +
  geom_label(mapping = aes(x = 75.5, y = oddwidind, label = sknewind),
             na.rm = TRUE,
             size = 2,
             label.padding = unit(0.12, "lines"),
             fontface = "bold") +
  geom_label(mapping = aes(x = 78, y = evewidind, label = sknewind),
             na.rm = TRUE,
             size = 2,
             label.padding = unit(0.12, "lines"),
             fontface = "bold") +
  scale_x_continuous(breaks = seq(0, 80, 5),
                     limits = c(0, 83)) +
  scale_y_reverse(breaks = seq(0, 64, 4),
                  limits = c(65, 0)) +
  scale_color_manual(values = readsccols,
                     breaks = readbreaks,
                     guide = guide_legend(reverse = TRUE),
                     name = "Q Read") +
  labs(title = "ple-A-M",
       subtitle = "align group 1.180 : PEPPR deletions 1-65",
       y = "Rank Index of Deletion Size",
       x = "Nucleotide Position Distal to Cut-site") +
  theme(aspect.ratio = 0.6,
        axis.text.y = element_text(size = 8),
        panel.grid.major.x = element_line(size = 0.3),
        panel.grid.minor.x = element_blank(),
        panel.grid.major.y = element_line(color = "grey75", size = 0.5),
        panel.grid.minor.y = element_line(color = "grey85", size = 0.5)) +
  annotation_custom(anteegrop) +
  annotation_custom(antsgnpop)

### Step 5: Display data
## first make full classify vec and make dfr  =====
inserttvec <- rep(NA, 180)
zsumclassvec <- replace(inserttvec, mjacharvec, "charmmj")
zsumclassvec <- replace(zsumclassvec, pepktotvec, "peppr")
zsumclassvec <- replace(zsumclassvec, pepvarsecvec, "pepvar")
zsumclassvec <- replace(zsumclassvec, widelvarsvec, "diffdel")
zsumclassvec <- replace(zsumclassvec, mindelvarsvec, "minidel")
zsumclassvec <- replace(zsumclassvec, purinsvarsvec, "insrt")
frindvec <- c(1:180)
zclassnewdfr <- data.frame(frindvec, zsumclassvec)
names(zclassnewdfr) <- c("sknewind", "classify")
zclassnewdfr$classify <- factor(zclassnewdfr$classify,
                                levels = c("insrt",
                                           "minidel",
                                           "diffdel",
                                           "pepvar",
                                           "peppr",
                                           "charmmj"))

## make select label vecs and add to dfr =====
oddmjavec <- c(mjaval5, mjaval3, mjaval9, mjaval10, mjaval12, mjaval8, mjaval1)
evemjavec <- c(mjaval4, mjaval2, mjaval6, mjaval17, mjaval20, mjaval2, mjaval27)
inserttvec <- rep(NA, 180)
inserttvec <- replace(inserttvec, oddmjavec, oddvalvec)
zclassnewdfr$oddmmjpos <- inserttvec
inserttvec <- rep(NA, 180)
inserttvec <- replace(inserttvec, evemjavec, evevalvec)
zclassnewdfr$evemmjpos <- inserttvec
zclassnewdfr$rdpct <- round(extskuhldfr[ , "sknewval"]/denomval*100,
                            digits = 6)

### Plotting 
## ggploting classification
library(ggplot2)
zclassnewdfrplot <- zclassnewdfr[1:100, ]
ggplot(data = zclassnewdfrplot,
       mapping = aes (x = sknewind, y = classify, fill = classify))+
  geom_tile (size = 6, width = 0.8, height = 1)+
  theme_bw()+
  labs(title = "pattern change between indel classes", 
       x = "Rank of alleles",
       y = "Classify")+
  geom_label(mapping = aes(x = sknewind, y = classify, label = oddmmjpos),
             size = 2,
             label.padding = unit(0.12, "lines"),
             fontface = "bold")+
  geom_label(mapping = aes(x = sknewind, y = classify, label = evemmjpos),
             size = 2,
             label.padding = unit(0.12, "lines"),
             fontface = "bold")+
  theme(aspect.ratio = 0.12,
        panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank())

