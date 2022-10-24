# DEG Excercise

## Counting with multiple RNAseq datasets

- https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE136101
- https://www.ncbi.nlm.nih.gov/Traces/study/?acc=PRJNA561298&o=acc_s%3Aa 사이트에서 metadata, accession list 파일 다운로드 
- RNAseq 파일 다운로드

~~~
prefetch --option-file SRR_Acc_List.txt
~~~

- fastq 파일 가져오기 

~~~
fastq-dump -X 100000 --split-files SRR10009019/SRR10009019.sralite
fastq-dump -X 100000 --split-files SRR10009020/SRR10009020.sralite
fastq-dump -X 100000 --split-files SRR10009021/SRR10009021.sralite
fastq-dump -X 100000 --split-files SRR10009022/SRR10009022.sralite
fastq-dump -X 100000 --split-files SRR10009023/SRR10009023.sralite
fastq-dump -X 100000 --split-files SRR10009024/SRR10009024.sralite
~~~


```r
dirnames <- dir(file.path("examples", "ecoli"), pattern = "^SRR10")
```

- QC


```r
library(Rfastp)

fastq_report <- rfastp(
  read1 = file.path("examples", "ecoli", "SRR10009019_1.fastq"),
  read2 = file.path("examples", "ecoli", "SRR10009019_2.fastq"), 
  outputFastq = file.path("examples", "ecoli", "filtered_SRR10009019.fastq"))

fastq_report <- rfastp(
  read1 = file.path("examples", "ecoli", "SRR10009020_1.fastq"),
  read2 = file.path("examples", "ecoli", "SRR10009020_2.fastq"), 
  outputFastq = file.path("examples", "ecoli", "filtered_SRR10009020.fastq"))

fastq_report <- rfastp(
  read1 = file.path("examples", "ecoli", "SRR10009021_1.fastq"),
  read2 = file.path("examples", "ecoli", "SRR10009021_2.fastq"), 
  outputFastq = file.path("examples", "ecoli", "filtered_SRR10009021.fastq"))

fastq_report <- rfastp(
  read1 = file.path("examples", "ecoli", "SRR10009022_1.fastq"),
  read2 = file.path("examples", "ecoli", "SRR10009022_2.fastq"), 
  outputFastq = file.path("examples", "ecoli", "filtered_SRR10009022.fastq"))

fastq_report <- rfastp(
  read1 = file.path("examples", "ecoli", "SRR10009023_1.fastq"),
  read2 = file.path("examples", "ecoli", "SRR10009023_2.fastq"), 
  outputFastq = file.path("examples", "ecoli", "filtered_SRR10009023.fastq"))

fastq_report <- rfastp(
  read1 = file.path("examples", "ecoli", "SRR10009024_1.fastq"),
  read2 = file.path("examples", "ecoli", "SRR10009024_2.fastq"), 
  outputFastq = file.path("examples", "ecoli", "filtered_SRR10009024.fastq"))

```


- Reference 서열 준비


```r
library(genbankr)
library(Biostrings)

download.file(url="https://github.com/greendaygh/kribbr2022/raw/main/ecoli-mg1655.gb", destfile="examples/ecoli-mg1655.gb")

ecoli <- readGenBank("examples/ecoli-mg1655.gb")
ecoliseq <- getSeq(ecoli)


maxlen <- 1000000
ecoliseqsub <- subseq(ecoliseq, 1, maxlen)
writeXStringSet(ecoliseqsub, "examples/ecoli/ecolisub.fasta")

```


- index 생성


```r
library(Rsubread)

buildindex(basename = file.path("examples", "ecoli", "ecolimedia"), 
           reference = file.path("examples", "ecoli", "ecolisub.fasta"))
```


- mapping 



```r

alignstat <- align(file.path("examples", "ecoli", "ecolimedia")
                   , readfile1 = file.path("examples", "ecoli", "filtered_SRR10009019.fastq_R1.fastq.gz")
                   , readfile2 = file.path("examples", "ecoli", "filtered_SRR10009019.fastq_R2.fastq.gz")
                   , output_file = file.path("examples", "ecoli", "ecolimedia_19.BAM")
                   , nthreads = 6)



alignstat <- align(file.path("examples", "ecoli", "ecolimedia")
                   , readfile1 = file.path("examples", "ecoli", "filtered_SRR10009020.fastq_R1.fastq.gz")
                   , readfile2 = file.path("examples", "ecoli", "filtered_SRR10009020.fastq_R2.fastq.gz")
                   , output_file = file.path("examples", "ecoli", "ecolimedia_20.BAM")
                   , nthreads = 6)


alignstat <- align(file.path("examples", "ecoli", "ecolimedia")
                   , readfile1 = file.path("examples", "ecoli", "filtered_SRR10009021.fastq_R1.fastq.gz")
                   , readfile2 = file.path("examples", "ecoli", "filtered_SRR10009021.fastq_R2.fastq.gz")
                   , output_file = file.path("examples", "ecoli", "ecolimedia_21.BAM")
                   , nthreads = 6)


alignstat <- align(file.path("examples", "ecoli", "ecolimedia")
                   , readfile1 = file.path("examples", "ecoli", "filtered_SRR10009022.fastq_R1.fastq.gz")
                   , readfile2 = file.path("examples", "ecoli", "filtered_SRR10009022.fastq_R2.fastq.gz")
                   , output_file = file.path("examples", "ecoli", "ecolimedia_22.BAM")
                   , nthreads = 6)


alignstat <- align(file.path("examples", "ecoli", "ecolimedia")
                   , readfile1 = file.path("examples", "ecoli", "filtered_SRR10009023.fastq_R1.fastq.gz")
                   , readfile2 = file.path("examples", "ecoli", "filtered_SRR10009023.fastq_R2.fastq.gz")
                   , output_file = file.path("examples", "ecoli", "ecolimedia_23.BAM")
                   , nthreads = 6)


alignstat <- align(file.path("examples", "ecoli", "ecolimedia")
                   , readfile1 = file.path("examples", "ecoli", "filtered_SRR10009024.fastq_R1.fastq.gz")
                   , readfile2 = file.path("examples", "ecoli", "filtered_SRR10009024.fastq_R2.fastq.gz")
                   , output_file = file.path("examples", "ecoli", "ecolimedia_24.BAM")
                   , nthreads = 6)

#alignstat

```


- Sorting



```r
library(Rsamtools)

sortBam(file = file.path("examples", "ecoli", "ecolimedia_19.BAM")
        , destination = file.path("examples", "ecoli", "sorted_ecolimedia_19.BAM"))

sortBam(file = file.path("examples", "ecoli", "ecolimedia_20.BAM")
        , destination = file.path("examples", "ecoli", "sorted_ecolimedia_20.BAM"))

sortBam(file = file.path("examples", "ecoli", "ecolimedia_21.BAM")
        , destination = file.path("examples", "ecoli", "sorted_ecolimedia_21.BAM"))

sortBam(file = file.path("examples", "ecoli", "ecolimedia_22.BAM")
        , destination = file.path("examples", "ecoli", "sorted_ecolimedia_22.BAM"))

sortBam(file = file.path("examples", "ecoli", "ecolimedia_23.BAM")
        , destination = file.path("examples", "ecoli", "sorted_ecolimedia_23.BAM"))

sortBam(file = file.path("examples", "ecoli", "ecolimedia_24.BAM")
        , destination = file.path("examples", "ecoli", "sorted_ecolimedia_24.BAM"))

indexBam(files = file.path("examples", "ecoli", "sorted_ecolimedia_19.BAM.bam"))
indexBam(files = file.path("examples", "ecoli", "sorted_ecolimedia_20.BAM.bam"))
indexBam(files = file.path("examples", "ecoli", "sorted_ecolimedia_21.BAM.bam"))
indexBam(files = file.path("examples", "ecoli", "sorted_ecolimedia_22.BAM.bam"))
indexBam(files = file.path("examples", "ecoli", "sorted_ecolimedia_23.BAM.bam"))
indexBam(files = file.path("examples", "ecoli", "sorted_ecolimedia_24.BAM.bam"))

```


- Counting



```r
library(GenomicAlignments)
library(plyranges)

ecolicds <- cds(ecoli)
ecolicds_sub <- ecolicds %>% 
  filter(end < maxlen)
seqlengths(ecolicds_sub) <- maxlen


filelist <- c(file.path("examples", "ecoli", "sorted_ecolimedia_19.BAM.bam"),
              file.path("examples", "ecoli", "sorted_ecolimedia_20.BAM.bam"),
              file.path("examples", "ecoli", "sorted_ecolimedia_21.BAM.bam"),
              file.path("examples", "ecoli", "sorted_ecolimedia_22.BAM.bam"),
              file.path("examples", "ecoli", "sorted_ecolimedia_23.BAM.bam"),
              file.path("examples", "ecoli", "sorted_ecolimedia_24.BAM.bam")
              )
mybam <- BamFileList(filelist)
myresult <- summarizeOverlaps(ecolicds_sub, mybam, ignore.strand = T)

class(myresult)
assay(myresult)
rowRanges(myresult)
colData(myresult)
metadata(myresult)

```

- Log transform and normalization
- 가정: 샘플마다 RNA 발현 정도는 유사하며 총 량은 같음. 


```r
library(tidyverse)

mydata <- assay(myresult)
mygenes <- rowRanges(myresult)
boxplot(mydata)

mydatat <- mydata %>% 
  data.frame() %>% 
  rownames_to_column()  %>% 
  pivot_longer(-rowname) %>%
  pivot_wider(names_from = rowname, values_from = value)

mymeanval <- mydatat %>% 
  summarise(across(where(is.numeric), mean))
mysdval <- mydatat %>% 
  summarise(across(where(is.numeric), sd))

mydf <- mymeanval %>% 
  bind_rows(mysdval) %>% 
  as.data.frame
rownames(mydf) <- c("mean", "sd")

mydft <- mydf %>% 
  rownames_to_column() %>% 
  pivot_longer(-rowname) %>% 
  pivot_wider(names_from = rowname, values_from = value)

ggplot(mydft, aes(x=mean, y=sd)) +
  geom_point() +
  scale_x_log10() +
  scale_y_log10()
  
```

- installation DESeq2


```r
if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install("DESeq2")
```

- create DESeq2 object


```r
library(DESeq2)

metaData <- data.frame(Group = c("control", "control", "control", "case", "case", "cases"), row.names = colnames(mydata))
metaData

dds <- DESeqDataSetFromMatrix(countData = mydata,
                      colData = metaData,
                      design = ~Group,
                      rowRanges = mygenes)

dds2 <- DESeq(dds)

counts(dds2)
counts(dds2, normalized =T)

```


- Mean variance


```r

plotDispEsts(dds2)
?plotDispEsts

```


- DESeq results


```r

myres <- results(dds2, contrast = c("Group", "control", "case"))
myres

summary(myres)

plotMA(myres)
```


- Multiple testing
  - Bonferonni = pvalue * total genes tested
  - Benjamini-Hockberg = (Pvalue*total genes)/rank of pvalue





---


<a rel="license" href="http://creativecommons.org/licenses/by-nc-nd/4.0/"><img alt="크리에이티브 커먼즈 라이선스" style="border-width:0" src="https://i.creativecommons.org/l/by-nc-nd/4.0/88x31.png" /></a><br />이 저작물은 <a rel="license" href="http://creativecommons.org/licenses/by-nc-nd/4.0/">크리에이티브 커먼즈 저작자표시-비영리-변경금지 4.0 국제 라이선스</a>에 따라 이용할 수 있습니다.

