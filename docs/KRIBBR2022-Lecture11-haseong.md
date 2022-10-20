## Lecture 11 (1006)

**class 01**

- install Rfastp 
- working 디렉토리에 examples 디렉토리 생성
- SRR11549087 예제 데이터는 다음 스크립트로 다운로드 

~~~
fastq-dump -X 10000 --split-files SRR11549076
~~~


```r
library(Rfastp)

download.file(url = "https://github.com/greendaygh/kribbr2022/raw/main/fastq/SRR11549087_1.fastq", destfile = "examples/SRR11549087_1.fastq")

fqfiles <- dir(path = "examples", pattern = "*.fastq")
fqfiles 

"examples/SRR11549087_1.fastq"

example_dir <- "examples"
file.path(example_dir, "SRR11549087_1.fastq")

fastq_report <- rfastp(read1 = file.path(example_dir, "SRR11549087_1.fastq"),
       outputFastq = file.path(example_dir, "filtered_SRR11549087_1.fastq"))

?rfastp

```


- 탐색기로 working direcotry 이동
- 주소창에 cmd 입력


```r
library(GEOquery)

gds <- getGEO(filename=system.file("extdata/GDS507.soft.gz",package="GEOquery"))


methods(class=class(gds))

gsm <- getGEO(filename=system.file("extdata/GSM11805.txt.gz",package="GEOquery"))

class(gsm)
gsm

methods(class=class(gsm))
Table(gsm)

```




```r


gse2553 <- getGEO('GSE2553',GSEMatrix=TRUE)
es <- gse2553[[1]]

dim(exprs(es))
fData(es)
pData(es)

```



```r
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install("airway")
```



```r
library(airway)
data(airway)


class(airway)

aw <- assay(airway)
str(aw)
methods(class=class(aw))
rowRanges(airway)
?rowRanges
colData(airway)
metadata(airway)
```





**class 02**

- airway data DEG 분석


```r
library(tidyverse)

aw_assay <- data.frame(assay(airway)[1:1000,])
aw_feature <- rowRanges(airway)[1:1000,]
aw_sample <- colData(airway)

aw_assay_t <- aw_assay %>% 
  rownames_to_column() %>% 
  pivot_longer(cols = -rowname) %>% 
  pivot_wider(names_from = rowname, values_from = value)

aw_assay2 <- aw_assay_t %>% 
  mutate(grp = aw_sample$dex) %>% 
  dplyr::select(name, grp, everything()) 


aw_assay3 <- aw_sample %>% 
  as.data.frame() %>% 
  dplyr::select(dex) %>% 
  rownames_to_column() %>% 
  left_join(aw_assay_t, by=c("rowname" =  "name")) 

```

- plot


```r

aw_assay_mean <- aw_assay2 %>% 
  group_by(grp) %>% 
  summarise(across(where(is.numeric), mean))

mydata <- aw_assay_mean %>% 
  pivot_longer(cols=-grp) 

ggplot(mydata, aes(x=name, y=value, color=grp)) +
  geom_point()

mydata2 <- aw_assay_mean %>% 
  pivot_longer(cols=-grp) %>% 
  pivot_wider(names_from = grp, values_from = value)

ggplot(mydata2, aes(x=trt, y=untrt)) +
  geom_point() + 
  scale_x_log10() +
  scale_y_log10()

```


- t-test, t.test 함수 사용 
- 아래 코드는 수정 필요 (기존 ExpressionSet 예제 포함)


```r
grp <- aw_assay3$dex

myttest <- function(x){
  z <- t.test(x[grp=="untrt"], x[grp=="trt"])
  z$p.value
}

tmp <- t.test(1:10, 21:30)
names(tmp)

testpval <- aw_assay3 %>% 
  summarise(across(where(is.numeric), myttest))

```

**class 03**

- reference genome 생성 


```r
library(genbankr)
library(Biostrings)

download.file(url="https://github.com/greendaygh/kribbr2022/raw/main/ecoli-mg1655.gb", destfile="examples/ecoli-mg1655.gb")

# ecoli <- readGenBank("examples/ecoli-mg1655.gb")

mydir <- "examples"
ecoli <- readGenBank(file.path(mydir, "ecoli-mg1655.gb"))

ecoliseq <- getSeq(ecoli)
ecolisub <- subseq(ecoliseq, 1, 10000)

writeXStringSet(ecolisub, file.path(mydir, "ecolisub.fasta"))

```

- Rsubread 패키지 사용해서 mapping


```r
if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install("Rsubread")

```



```r
library(Rsubread)

buildindex(basename = file.path(mydir, "example"),
           reference = file.path(mydir, "ecolisub.fasta"))
```


- mapping


```r
myaln <- align(index = file.path(mydir, "example"),
               readfile1 = file.path(mydir, "filtered_SRR11549087_1.fastq_R1.fastq.gz"),
               output_file = file.path(mydir, "ecoliexample.BAM"),
               nthreads = 5)


myaln

```

- sorting


```r

if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install("Rsamtools")

```



```r
library(Rsamtools)

sortBam(file = file.path(mydir, "ecoliexample.BAM"),
        destination = file.path(mydir, "sorted_ecoliexample.BAM"))

indexBam(files = file.path(mydir, "sorted_ecoliexample.BAM.bam"))

```


- genome browser (IGV, UCSC Genome Browser)
- IGV 설치
- Reference genome, sorted bam/bai file 

- Counting, GenomicAlignments::summarizeOverlaps
- counting을 위해서는 bam 파일과 annotation 파일이 필요



```r
library(GenomicAlignments)
library(plyranges)

## ecoli genbank information
class(ecoli)
methods(class="GenBankRecord")

ecolicds <- cds(ecoli)
class(ecolicds)
ecolicdssub <- ecolicds %>% 
  plyranges::filter(end < 10000)

seqlengths(ecolicdssub) <- 10000
seqlengths(ecolicdssub)


?summarizeOverlaps

mybam <- BamFile(file.path(mydir, "sorted_ecoliexample.BAM.bam"))

mycnt <- summarizeOverlaps(features = ecolicdssub,
                  reads = mybam)
class(mycnt)

assay(mycnt)
rowRanges(mycnt)
colData(mycnt)
```


















