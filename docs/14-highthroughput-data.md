# DEG analysis with Bioconductor

## High-throughput genomic data analysis


![](images/12/htanalysis.png){width=600}


## DEG analysis 

Differentially Expressed Gene 분석은 전통적 two-channel microarray나 RNA-Seq 데이터를 활용한 분석입니다. Genome reference에 fastq 파일의 read들을 맵핑하고 맵핑된 read들을 카운팅하여 해당 유전자의 발현을 정량화하고 이를 기준이 되는 발현값과 비교하여 질병이나 조건에 따라 다른 발현값을 갖는 유전자를 찾는 방법입니다. Reference 서열 정보와 발현된 mRNA 서열을 분석한 fastq 파일이 필요합니다. 다양한 분석 툴이 소개되고 있으며 진핵세포의 경우 splicing 등을 고려한 mapping 기술이 필요합니다. 


## Creating a reference genome



```r
library(genbankr)
library(Biostrings)

download.file(url="https://github.com/greendaygh/kribbr2022/raw/main/ecoli-mg1655.gb", destfile="examples/ecoli-mg1655.gb")

ecoli <- readGenBank("examples/ecoli-mg1655.gb")
ecoliseq <- getSeq(ecoli)

ecoliseqsub <- subseq(ecoliseq, 1, 100000)
names(ecoliseqsub) <- "ecolisub"
writeXStringSet(ecoliseqsub, "examples/ecolisub.fasta")

```



## Read Alignment and Mapping 

Rsubread 패키지를 활용해서 align을 수행합니다. align 전에 reference genome의 index를 생성해야 하며 지놈이 클 경우 오랜 시간이 소요될 수 있으니 주의가 필요합니다.


```r
if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install("Rsubread")
```


index 생성은 Rsubread 패키지의 buildindex 함수를 사용하며 memory 옵션을 이용해서 데이터 사이즈에 따라 할당되는 메모리를 조절할 수 있습니다. 기본은 8GB로 memory = 8000 입니다. indexSplit 옵션으로 크로모좀별로 인덱스를 분리해서 메모리 효율을 증가시킬 수 있습니다.  


```r
library(Rsubread)

buildindex(basename = file.path("examples", "ecoliexample"), 
           reference = file.path("examples", "ecolisub.fasta"))
```

# RNA-Seq alignment (Mapping)

Rsubread 패키지를 활용해서 mapping을 수행합니다. align 함수를 사용하며 splicing 여부에 따라 옵션이 조금씩 다를 수 있습니다. 



```r
library(Rsubread)

alignstat <- align(file.path("examples", "ecoliexample")
                   , readfile1 = file.path("fastq", "ftn-rep1_S10_L001_R1_001.fastq.gz")
                   , output_file = file.path(targetdir, "test.BAM")
                   , nthreads = 13)


#alignstat

```



# sorting



```r

library(Rsamtools)

sortBam(file = file.path(targetdir, "test.BAM")
        , destination = file.path(targetdir, "sort_test.BAM"))

indexBam(files = file.path(targetdir, "sort_test.BAM.bam"))

```



SAM 파일은 Sequence alignment data를 담고 있는 텍스트 파일(.txt)로 각 내용은 탭(tab)으로 분리되어 alignment, mapping 정보를 담고 있습니다. 


SAM 파일은 텍스트 파일의 문자열 형식으로 저장하여 바로 열람할 수 있으며, 이를 압축하고 색인화하여 바이너리 형식으로 변환한 것이 BAM 파일이다.

NGS의 발달로 불특정 다수의 organism에서 유전체 혹은 전사체 서열이 대량으로 시퀀싱 되고 있다. Human의 경우는 개인차에 의한 다수의 변이 정보를 밝히고 이것이 질병과 연관된 변이인지를 밝히기 위해 시퀀싱 된 reads는 유전체 서열에 다시 remapping 되기도 하고, 새로운 organism의 유전체 정보를 밝히기 위해서도 remapping이 이뤄지고 있다. SAM 파일은 '1000 genome project'를 진행하면서 공동 연구의 효율성을 위해 데이터의 공유를 표준화 하려는 방안으로 채택된 remapping의 표준 포맷이다.

Remapping을 위해 많이 이용하는 software인 BWA, Bowtie, CLCAssemblyCell 등은 모두 mapping output으로 SAM 파일을 형성한다.

SAM 파일들을 다루는데 필요한 software package로는 SAMtools 가 있다.



## Binary Alignment Map (BAM) 

시퀀싱 파일은 GEO 데이터베이스 [GSE52778](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE52778)에 있으며 패키지에도 기본으로 포함되어 있습니다. 


```r

dir <- system.file("extdata", package="airway", mustWork=TRUE)
list.files(dir)
readLines(file.path(dir, "Homo_sapiens.GRCh37.75_subset.gtf"))

```

.bam 파일은 RNA-Seq read를 포함하고 있으며 .csv 파일은 실험 디자인, .gtf 파일은 

- gff 파일 포맷과 gtf 파일 포멧 차이



```r

csvfile <- file.path(dir, "sample_table.csv")
sampleTable <- read.csv(csvfile, row.names = 1)

bamfilenames <- file.path(dir, paste0(sampleTable$Run, "_subset.bam"))

```

Rsamtools는 bam이나 sam 파일을 읽을 수 있는 패키지로 yieldSize 파라메터로 메모리 과사용을 막을 수 있습니다. 


```r
library(Rsamtools)

bamfiles <- BamFileList(bamfilenames, yieldSize = 2000)
bamfiles
class(bamfiles)
seqinfo(bamfiles[1])

```



---


<a rel="license" href="http://creativecommons.org/licenses/by-nc-nd/4.0/"><img alt="크리에이티브 커먼즈 라이선스" style="border-width:0" src="https://i.creativecommons.org/l/by-nc-nd/4.0/88x31.png" /></a><br />이 저작물은 <a rel="license" href="http://creativecommons.org/licenses/by-nc-nd/4.0/">크리에이티브 커먼즈 저작자표시-비영리-변경금지 4.0 국제 라이선스</a>에 따라 이용할 수 있습니다.

