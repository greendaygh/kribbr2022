## Lecture 09 (0901)

**class 01**

- msa 
- `NC_045512.2` NCBI 다운로드 


```r
library(Biostrings)
library(rentrez)

entrez_dbs()
tmps <- entrez_fetch(db="nuccore", id = "NC_045512.2", rettype = "fasta")
tmps
write.table(tmps, "covid19yu.fasta", quote = F, row.names = F, col.names = F)

download.file(url = "https://raw.githubusercontent.com/greendaygh/kribbr2022/main/covid_table.csv", destfile = "covid_table2.csv")

covid19 <- read.csv("covid_table2.csv")
selid <- covid19$gbacc[1:4]

tmps <- entrez_fetch(db="nuccore", id = selid, rettype = "fasta")
write.table(tmps, "covid19.fasta", quote = F, row.names = F, col.names = F)

```

- msa 설치
- Bioconductor landing page (구글에서 bioconductor msa)
- installation 코드 복사 후 콘솔에서 붙여넣기


```r
## install 
# if (!require("BiocManager", quietly = TRUE))
#     install.packages("BiocManager")
# 
# BiocManager::install("msa")


library(msa)


mycovidseq <- readDNAStringSet("covid19.fasta") 
class(mycovidseq)
mycovidseq
```

- 일부 서열만 가지고 실습
- (팁) 블럭지정 후 ctrl + enter 블럭만 실행 


```r
#subseq
#substr


mycovidsubseq <- subseq(mycovidseq, 1, 2000)

aln <- msa(mycovidsubseq, method="Muscle")
aln
DNAStringSet(aln)

?msa

print(aln, show="complete")
#msaPrettyPrint(aln, output = "pdf")

alphabetFrequency(aln)

```


- IRanges 설치Bioconductor landing page (구글에서 bioconductor IRanges 검색 ) installation 코드 복사 후 콘솔에서 붙여넣기




```r
library(IRanges)
nchar(consensusString(aln))

mym <- IRanges(start=c(1, 73), end=c(36, nchar(consensusString(aln))))
colmask(aln) <- mym
print(aln, show="complete")

alphabetFrequency(aln)
unmasked(aln)

```

- DECIPHER


```r
library(DECIPHER)

aln <- AlignSeqs(mycovidsubseq)
BrowseSeqs(aln)
BrowseSeqs(aln, colWidth = 100)

```

**class 02**

- ggtree (bioconductor), treeio (Bioconductor), ape (CRAN) 설치 
- 1) MSA - 2) Distance - 3) clustering - 4) Tree


```r

aln <- msa(mycovidsubseq, method="Muscle")
class(aln)

?stringDist

mydist <- stringDist(DNAStringSet(aln))
class(mydist)

```


- label 축소
- treeio error 발생시, ggtree 패키지 unload (detach("package:ggtree", unload = TRUE)) 후에 treeio, ggtree 순서로 로딩 


```r
library(treeio)
library(ggtree)


tmps <- strsplit(names(mycovidsubseq), split=" ")
mynames <- sapply(tmps, function(x){x[1]})
names(mycovidsubseq) <- mynames
mycovidsubseq

aln <- msa(mycovidsubseq)
mydist <- stringDist(DNAStringSet(aln))
myclust <- hclust(mydist)
plot(myclust)

myphylo <- as.phylo(myclust)

ggtree(myphylo)  +
  geom_tiplab(color="firebrick")

ggtree(myphylo, layout="circular")  +
  geom_tiplab(color="firebrick")

```

- DECIPHER 버전업 되면서 바뀐부분 체크 


```r
aln <- AlignSeqs(mycovidsubseq)
#mydist <- DistanceMatrix(aln)
# myclust <- IdClusters(aln, cutoff=10, method="NJ", type="dendrogram")
# class(mydist)

```

- as.tibbe


```r
library(tidyverse)

as.tibble(myphylo)
```

- 여러 그림 그리기 


```r

mytree <- rtree(30)
class(mytree)

p <- ggtree(mytree) +
  geom_tiplab()

p

a <- runif(30, 0, 1)
b <- 1 - a
df <- data.frame(mytree$tip.label, a, b)
df2 <- pivot_longer(df, -mytree.tip.label)

p2 <- p + geom_facet(panel="bar", 
                     data=df2, 
                     geom=geom_bar, 
                     mapping = aes(x = value, fill = as.factor(name)),
                     orientation = 'y', 
                     width=0.8,
                     stat = 'identity') +
  xlim_tree(9)

facet_widths(p2, widths = c(1, 2))

```

- Genome 분석  
- `NC_045512.2` NCBI 다운로드 


```r
library(Biostrings)
library(rentrez)


tmps <- entrez_fetch(db="nuccore", id = "NC_045512.2", rettype = "fasta")
tmps
write.table(tmps, "covid19yu.fasta", quote = F, row.names = F, col.names = F)

covid19seq <- readDNAStringSet("covid19yu.fasta")

```


- read genbank


```r
require(genbankr)

tmps <- entrez_fetch(db="nuccore", id = "NC_045512.2", rettype = "gb")
tmps
write.table(tmps, "covid19yu.gb", quote = F, row.names = F, col.names = F)

covidgb <- readGenBank("covid19yu.gb")
covidgb

cds(covidgb)
class(covidgb)
methods(class="GenBankRecord")

covidseq <- getSeq(covidgb)
class(covidgb)

cdsseq <- getSeq(covidseq, cds(covidgb))
class(covidseq)
```

**class 03**



```r
library(IRanges)

ir <- IRanges(start = 1, end = 10)
ir

ir <- IRanges(start = 1:10, width = 10:1)
ir
class(ir)

x <- "GATTGCCCCCCTAG"
y <- unlist(strsplit(x, split=""))
z <- Rle(y)
write.table(z, "tmp.txt")

```


```r
cds(covidgb)
```



```r
library(ggbio)

ir <- ranges(cds(covidgb))
ir

start(ir)
end(ir)
width(ir)

autoplot(ir) +
  theme_bw()

autoplot(ir, aes(fill=width)) +
  theme_bw()




ir <- IRanges(c(1, 10, 100), end=c(20, 50, 110))
autoplot(ir)


disjointBins(ir)
autoplot(disjoin(ir))
autoplot(IRanges::reduce(ir))

```


- ggbio 설치 


```r

IRanges(start = c(1,1,1,1), width = c(100, 15, 30, 4))

```




```r
library(GenomicRanges)
cds(covidgb)

gr <- GRanges(seqnames = "chr1", 
              ranges = IRanges(1:10, 21:30),
              strand = "+")
gr


gr <- GRanges(seqnames = "chr1", 
              ranges = IRanges(1:10, 21:30),
              strand = "+", 
              names = paste0("g", 1:10), 
              score = 1:10)
gr

seqnames(gr)
ranges(gr)
strand(gr)
start(gr)
end(gr)

granges(gr) 
mcols(gr)
seqlengths(gr) <- 50

autoplot(gr)


gr2 <- GRanges(seqnames = "chr2", 
              ranges = IRanges(1:10, 21:30),
              strand = "+", 
              names = paste0("g", 11:20), 
              score = 1:10)
gr3 <- c(gr, gr2)
seqnames(gr3)
seqlengths(gr3)

split(gr3, c(rep("a", 10), rep("b", 10)))

```

- disjoin
- reduce 


**class 04**

- 대장균 지놈 plot 


```r
tmps <- entrez_fetch("nuccore", id="NC_000913.3", rettype="gbwithparts")
write(tmps, "ecoli-mg1655.gb")

ecoligb <- readGenBank("ecoli-mg1655.gb")
ecoli_cds <- cds(ecoligb)
autoplot(ecoli_cds)
seqlengths(ecoli_cds)
ggbio() +
  circle(ecoli_cds, geom="ideo", fill = "gray70") +
  circle(ecoli_cds, geom="scale", size=3) +
  circle(ecoli_cds, geom = "text", aes(label = locus_tag), vjust = 0, size = 3) +
  theme(
    axis.text.x = element_text(angle=90)
  )


gr1 <- granges(ecoli_cds)
gr2 <- granges(ecoli_cds)
mcols(gr2)$test <- rnorm(length(ecoli_cds))

ggplot() + 
  layout_circle(ecoli_cds, geom = "ideo", fill = "gray70", radius = 9, trackWidth = 1) +
  layout_circle(ecoli_cds, geom = "scale", size = 3, trackWidth = 1, scale.n=20) +
  layout_circle(gr1, geom = "rect", color = "steelblue",  radius = 5) +
  layout_circle(gr2, geom = "bar", aes(y=test), size = 3, trackWidth = 1, scale.n=20, radius = 4) 
  

gr1disj <- disjoin(gr1)

ggplot() + 
  layout_circle(ecoli_cds, geom = "ideo", fill = "gray70", radius = 9, trackWidth = 1) +
  layout_circle(ecoli_cds, geom = "scale", size = 3, trackWidth = 1, scale.n=20) +
  layout_circle(gr1disj, geom = "rect", color = "steelblue",  radius = 5) +
  layout_circle(gr2, geom = "bar", aes(y=test), size = 3, trackWidth = 1, scale.n=20, radius = 4) 

```



```r
library(plyranges)

if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install("plyranges")

covidseq <- getSeq(covidgb)
cdsrng <- cds(covidgb)
cdsseq <- getSeq(covidseq, cdsrng)

gcratio <- rowSums(letterFrequency(cdsseq, c("G", "C"))) / nchar(cdsseq)

cdsrng %>% 
  dplyr::select(gene, product) %>% 
  mutate(gcratio)
  

```

- gggene install --> CRAN 
- Rstudio > Packages > install


```r
library(gggenes)

mydf <- data.frame(seqname = seqnames(cdsrng),
           start = start(cdsrng), 
           end = end(cdsrng), 
           strand = case_when(
             as.vector(strand(cdsrng))=="+"~ TRUE,
             as.vector(strand(cdsrng))=="-"~ FALSE
           ),
           gene = mcols(cdsrng)$gene)

ggplot(mydf, aes(xmin = start, 
                 xmax = end, 
                 y = seqname,
                 label = gene, 
                 fill = gene,
                 forward = strand)) +
  geom_gene_arrow(arrowhead_height = grid::unit(12, "mm"),
    arrowhead_width = grid::unit(6, "mm"),
    arrow_body_height = grid::unit(12, "mm")) +
  geom_gene_label() +
  theme(legend.position = "none")
      

```


















