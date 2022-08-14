# Tools for genome analysis


## genbank file


::: rmdnote
**Exercises **
 
1. `NC_045512.2`는 우한에서 발생한 코로나바이러스의 accession number임. `entrez_fetch` 함수를 사용하여 `nuccore` 데이터베이스에서 genbank 정보를 다운로드 받으시오

2. 받은 택스트를 `covid19wohan.gb`라는 파일로 저장하시오 



:::




```r
require(genbankr)

covid19 <- readGenBank("covid19wuhan.gb")
covid19
methods(class="GenBankRecord")
cds(covid19)
exons(covid19)

covid19seq <- getSeq(covid19)

```

 
## IRanges

유전체 데이터의 대부분을 차지하는 정보는 전체 지놈 서열 중 어디서 어디까지가 유전자 또는 coding sequence 이고 그 번역된 정보가 무엇인지 설명하는 정보 입니다. 즉, 일련의 feature에 대한 위치와 특성 정보를 분석하는 것이 효율적인 지놈을 분석하기 위해 필수입니다. `bioconductor` 에서는 이러한 유전체 정보를 효율적으로 분석하고 가시화 하기위한 방법들이 다양하게 개발되어 왔으며 `IRanges` 와 `GenomicRanges` 라는 패키지가 대표적으로 사용될 수 있습니다. 

IRanges는 간격을 나타내는 임의의 숫자 세트이며 지놈상에 위치한 특정 feature들의 간격이나 통계적 수치들을 효율적으로 나타내기 위해서 만들어진 패키지 입니다 [@Lawrence2013]. 임의의 feature에 대한 시작, 끝, 넓이를 나타내는 숫자들이 리스트로 이루어져 있습니다. 


```r
library(IRanges)

ir <- IRanges(start = c(1,3,5), end = c(3,5,7))
ir

ir <- IRanges(start = 1:10, width = 10:1)
ir
class(ir)
methods(class="IRanges")
?IRanges

```

IRange 객체로부터 몇몇 함수를 사용하여 정보를 추출할 수 있습니다. Rle (run-length encoding format) class 를 사용하며 Rle는 후에 좀 더 자세히 알아보겠습니다. 


```r
ir <- IRanges(start = c(1,3), end = c(4,5))
ir

start(ir)
end(ir)
width(ir)
disjointBins(ir)

ir <- IRanges(start = c(1,3,6), end = c(4,5,7))
ir
bins <- disjointBins(ir)
bins
ir2 <- disjoin(ir)
Rle(1:10, 1:10)

reduce(ir)
```


이러한 정보를 가시화하는 가장 간단한 방법은 `ggbio`라는 패키지를 사용하는 것 입니다. 



```r
library(ggbio)

autoplot(ir) 
autoplot(ir2) 

autoplot(ir) + 
  theme_bw()

autoplot(ir, aes(fill=width)) +
  theme_bw()
```


## Genomic ranges

[GenomicRanges](https://bioconductor.org/packages/release/bioc/html/GenomicRanges.html)는 지놈상의 위치정보와 Bioconductor에서 제공하는 다양한 high-throughput 정보들을 같이 표현하기 위해서 만들어진 패키지입니다. 

먼저 Rle (Run-length encoding) 개념을 알아봅니다. Rle는 런 렝스 부호화라고 하며 일종의 압축 방법입니다. 예를 들어 GATTGCCCCCCTAG 라는 서열이 있다고 하면 이를 그대로 text 파일에 저장하지 않고 GAT2GC6TAG 라고 표현함으로써 용량을 줄이는 압축의 기능을 합니다. GenomicRange는 이러한 Rle 개념을 사용하기 위해서 `Rle`라는 기본 함수를 사용합니다. 


```r
library(IRanges)

x <- "GATTGCCCCCCTAG"
y <- unlist(strsplit(x, split=""))
yrle <- Rle(y)
yrle

runLength(yrle)
runValue(yrle)
nrun(yrle)

x <- Rle(values = c(1:3), lengths = c(1:3))
x
class(x)
#methods(class="Rle")

# convert Rle to IRanges
xrange <- IRanges(start(x), end(x))
xrange

```

GRanges 함수를 이용해서 생성할 수 있으며 `browseVignettes("GenomicRanges")` 나 `methods()` 함수를 이용해서 관련된 기능을 찾아서 사용할 수 있습니다.  



```r
library(GenomicRanges)

gr <- GRanges(
    seqnames = Rle(c("chr1", "chr2", "chr1", "chr3"), c(1, 3, 2, 4)),
    ranges = IRanges(101:110, end = 111:120, names = head(letters, 10)),
    strand = Rle(strand(c("-", "+", "*", "+", "-")), c(1, 2, 2, 3, 2)),
    score = 1:10,
    GC = seq(1, 0, length=10))
gr
class(gr)

seqnames(gr)
ranges(gr)
strand(gr)

granges(gr) 
mcols(gr) #meta data

seqlengths(gr) <- c(249250621, 243199373, 198022430)
seqlengths(gr)
names(gr)

sp <- split(gr, rep(1:2, each=5))

autoplot(gr)
```



::: rmdnote
**Exercises **
 
코로나 바이러스서열들의 genbank 파일을 모두 다운로드 하고 각 바이러스별로 CDS 추출, 서열 비교, 가시화 등을 수행하시오 




:::



위 GenomicRanges 데이터를 dplyr 형태로 좀 더 쉽게 다루기 위한 패키지가 `plyragnes` 입니다. 



```r
library(plyranges)

covid19 <- readGenBank("covid19wuhan.gb")
covid19seq <- getSeq(covid19)
covid19cds <- cds(covid19)
covidcdsseq <- getSeq(covid19seq, covid19cds)

gcr <- rowSums(letterFrequency(covidcdsseq, c("G", "C"), as.prob=T))

covid19cds %>% 
  dplyr::select(gene, product) %>% 
  mutate(gc = gcr) #%>%
  #filter(grepl(pattern = "ORF", gene)) 

```


::: rmdnote
**Exercises **
 
위에서 계산된 GC 비율로 bar 그래프를 그래되 product를 라벨로 지정하여 그리시오 

:::


---


<a rel="license" href="http://creativecommons.org/licenses/by-nc-nd/4.0/"><img alt="크리에이티브 커먼즈 라이선스" style="border-width:0" src="https://i.creativecommons.org/l/by-nc-nd/4.0/88x31.png" /></a><br />이 저작물은 <a rel="license" href="http://creativecommons.org/licenses/by-nc-nd/4.0/">크리에이티브 커먼즈 저작자표시-비영리-변경금지 4.0 국제 라이선스</a>에 따라 이용할 수 있습니다.
