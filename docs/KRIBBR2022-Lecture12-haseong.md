## Lecture 12 (1020)

**Class 01**


- 참고 https://greendaygh.github.io/kribbr2022/lecture-note.html#lecture-07-0804


- Reference genome 
- examples 디렉토리 생성


```r
library(genbankr)
library(Biostrings)

download.file(url="http://github.com/greendaygh/kribbr2022/raw/main/ecoli-mg1655.gb", destfile="examples/ecoli-mg1655.gb")

ecoli <- readGenBank(file.path("examples", "ecoli-mg1655.gb"))
ecoliseq <- getSeq(ecoli)

maxlen <- 1000000
ecoliseqsub <- subseq(x = ecoliseq, start = 1, end = maxlen)
writeXStringSet(x = ecoliseqsub,
                filepath = "examples/ex1/ecolisubseq.fasta")
?writeXStringSet

```

- generation of index


```r
library(Rsubread)

buildindex(basename = file.path("examples", "ex1", "ecolisingle"),
           reference = file.path("examples", "ex1", "ecolisubseq.fasta")
           )

```

- NGS fastq 파일 다운로드 
- 파일 탐색 오픈
- `D:\Rstudy\lecture12\examples\ex1` 이동 
- 주소창에서 cmd 입력
- 다음 명령 실행
- SRA DB https://greendaygh.github.io/kribbr2022/high-throughput-data.html#ngs-database


~~~
fastq-dump -X 10000 --split-files SRR11549076
~~~

- QC
- 생성된 html 파일 확인 


```r
library(Rfastp)

myreport <- rfastp(read1 = "examples/ex1/SRR11549076_1.fastq", 
       outputFastq = "examples/ex1/filterd_SRR11549076_1")

```

- mapping to reference 



```r
library(Rsubread)

myaln <- align(index = file.path("examples", "ex1", "ecolisingle"),
               readfile1 = file.path("examples", "ex1", "filterd_SRR11549076_1_R1.fastq.gz"))

```

- Sorting


```r
library(Rsamtools)

sortBam(file = file.path("examples", "ex1", "filterd_SRR11549076_1_R1.fastq.gz.subread.BAM"),
        destination = file.path("examples", "ex1", "sorted_ecolisingle"))

indexBam(file = file.path("examples", "ex1", "sorted_ecolisingle.bam"))

```


- Counting



```r
library(GenomicAlignments)
library(plyranges)
library(tidyverse)

ecolisub <- genbankr::cds(ecoli) %>% 
  plyranges::filter(end < maxlen)
class(ecolisub)
seqlengths(ecolisub) <- maxlen

mybam <- BamFile(file = file.path("examples", "ex1", "sorted_ecolisingle.bam"))
class(mybam)

mycnt <- summarizeOverlaps(features = ecolisub, 
                  reads = mybam)
class(mycnt)


assay(mycnt)
colData(mycnt)
rowRanges(mycnt)
metadata(mycnt)

```


**Class 02**

- SRA software 설치 https://greendaygh.github.io/kribbr2022/high-throughput-data.html#ngs-database
- ex2 디렉토리 만들고
- 파일탐색기 ex2 주소창에 cmd 입력
- raw data 다운로드

~~~
prefetch --option-file SRR_Acc_List.txt
~~~

- fastq 파일 생성

~~~
Fastq-dump -X 200000 --split-files SRR10009019/SRR10009019.sra
~~~



```r

no <- c(19:24)
tmps <- paste("fastq-dump -X 200000 --split-files SRR100090", no, sep="")
tmps

```

- ex2 디렉토리에서 

~~~

fastq-dump -X 200000 --split-files SRR10009019
fastq-dump -X 200000 --split-files SRR10009020
fastq-dump -X 200000 --split-files SRR10009021
fastq-dump -X 200000 --split-files SRR10009022
fastq-dump -X 200000 --split-files SRR10009023
fastq-dump -X 200000 --split-files SRR10009024

~~~


- QC


```r
library(Rfastp)


for(i in c(19:24)){
  fastq_report <- rfastp(
  read1 = file.path("examples", "ex2", paste0("SRR100090", i, "_1.fastq")),
  read2 = file.path("examples", "ex2", paste0("SRR100090", i, "_2.fastq")),
  outputFastq = file.path("examples", "ex2", paste0("filtered_SRR100090", i, ".fastq"))
)
  
}

```

- 앞에서 만든 ecolisub 파일 ex2 디렉토리에 붙여넣기 


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

buildindex(basename = file.path("examples", "ex2", "ecolimedia"), 
           reference = file.path("examples", "ex2", "ecolisubseq.fasta"))

```

- mapping


```r
library(Rsubread)


for(i in c(19:24)){
  
  alignstat <- align(
    file.path("examples", "ex2", "ecolimedia"),
    readfile1 = file.path("examples", "ex2", paste0("filtered_SRR100090", i, ".fastq_R1.fastq.gz")),
    readfile2 = file.path("examples", "ex2", paste0("filtered_SRR100090", i, ".fastq_R2.fastq.gz")),
    output_file = file.path("examples", "ex2", paste0("ecolimedia_", i, ".BAM")),
    nthreads = 6)

}

```


- sorting


```r
library(Rsamtools)



sortBam(file = file.path("examples", "ex2", "ecolimedia_19.BAM")
        , destination = file.path("examples", "ex2", "sorted_ecolimedia_19.BAM"))

sortBam(file = file.path("examples", "ex2", "ecolimedia_20.BAM")
        , destination = file.path("examples", "ex2", "sorted_ecolimedia_20.BAM"))

sortBam(file = file.path("examples", "ex2", "ecolimedia_21.BAM")
        , destination = file.path("examples", "ex2", "sorted_ecolimedia_21.BAM"))

sortBam(file = file.path("examples", "ex2", "ecolimedia_22.BAM")
        , destination = file.path("examples", "ex2", "sorted_ecolimedia_22.BAM"))

sortBam(file = file.path("examples", "ex2", "ecolimedia_23.BAM")
        , destination = file.path("examples", "ex2", "sorted_ecolimedia_23.BAM"))

sortBam(file = file.path("examples", "ex2", "ecolimedia_24.BAM")
        , destination = file.path("examples", "ex2", "sorted_ecolimedia_24.BAM"))

indexBam(files = file.path("examples", "ex2", "sorted_ecolimedia_19.BAM.bam"))
indexBam(files = file.path("examples", "ex2", "sorted_ecolimedia_20.BAM.bam"))
indexBam(files = file.path("examples", "ex2", "sorted_ecolimedia_21.BAM.bam"))
indexBam(files = file.path("examples", "ex2", "sorted_ecolimedia_22.BAM.bam"))
indexBam(files = file.path("examples", "ex2", "sorted_ecolimedia_23.BAM.bam"))
indexBam(files = file.path("examples", "ex2", "sorted_ecolimedia_24.BAM.bam"))

```


- Counting 



```r

ecolicds <- cds(ecoli)
ecolicds_sub <- ecolicds %>% 
  filter(end < maxlen)
seqlengths(ecolicds_sub) <- maxlen


filelist <- c(file.path("examples", "ex2", "sorted_ecolimedia_19.BAM.bam"),
              file.path("examples", "ex2", "sorted_ecolimedia_20.BAM.bam"),
              file.path("examples", "ex2", "sorted_ecolimedia_21.BAM.bam"),
              file.path("examples", "ex2", "sorted_ecolimedia_22.BAM.bam"),
              file.path("examples", "ex2", "sorted_ecolimedia_23.BAM.bam"),
              file.path("examples", "ex2", "sorted_ecolimedia_24.BAM.bam")
              )
mybam <- BamFileList(filelist)
myresult <- summarizeOverlaps(ecolicds_sub, mybam, ignore.strand = T)

assay(myresult)
rowData(myresult)
colData(myresult)
metadata(myresult)

```


- DESeq2


```r
library(DESeq2)

mydata <- assay(myresult)
boxplot(mydata)

```

- mean vs variance


```r

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

- DESeq2


```r
library(DESeq2)

metaData <- data.frame(Group = c("LB", "LB", "LB", "M63", "M63", "M63"), 
                       row.names = colnames(mydata))
metaData
#mygenes <- rowData(myresult)

dds <- DESeqDataSetFromMatrix(countData = mydata,
                       colData = metaData,
                       design = ~Group,
                       rowRanges = mygenes
                       )


colData(myresult)$Group <- metaData$Group
colData(myresult)
dds <- DESeqDataSetFromMatrix(myresult)
```

