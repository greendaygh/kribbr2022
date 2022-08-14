# Lecture 08 Note

- 20220811

## class 01

- `rentrez` 패키지를 사용한 NCBID 데이터 다운로드 


```r
library(rentrez)

entrez_dbs()
entrez_db_summary("nuccore")

covid_paper <- entrez_search(db="pubmed", term="covid-19", retmax=40)
class(covid_paper)
methods(class="esearch")
names(covid_paper)
covid_paper$ids

katipo_search <- entrez_search(db="popset", term="Latrodectus katipo[Organism]")
class(katipo_search)
names(katipo_search)
katipo_search$ids

entsum <- entrez_summary(db="popset", id = katipo_search$ids)
entsum[[2]]$title
entsum[[2]]$article

extract_from_esummary(entsum, "title")

```

- 텍스트 형태로 다운로드 받으므로 파일에 쓰고 `readDNAString` 형태로 다시 읽어옴 


```r
library(Biostrings)

coi <- entrez_fetch(db="popset", id=katipo_search$ids, rettype="fasta")
write(coi, file = "coi.fasta")
coiseq <- readDNAStringSet(filepath = "coi.fasta")

```

- 뎅기바이러스 예제


```r

dv <- c("NC_001477", "NC_001474", "NC_001475", "NC_002640")
dvseq <- entrez_fetch(db="nuccore", id=dv, rettype = "fasta")
dvseq

write(dvseq, file="dvseq.fasta") 
dvseq <- readDNAStringSet("dvseq.fasta")  
dvseq

strwrap(as.character(dvseq[[1]]), width=10)
```

- Covid-19 예제 


```r

sresult <- entrez_search(db="popset", term="Covid-19", retmax=40)
class(sresult)
names(sresult)
sresult$ids

mysummary <- entrez_summary("popset", sresult$ids)
mysummary

#lapply
#sapply

mysummary
extract_title <- function(x){
  return(x$title)
}
mysummary[[1]]$title
extract_title(x=mysummary[[1]])


sapply(mysummary, extract_title)
tmp <- extract_from_esummary(mysummary, "title")
tmp
class(tmp)
mydata <- data.frame(title=tmp)
mydata


## sequence download

tmpstr <- entrez_fetch("popset", id=sresult$ids, rettype = "fasta")
write(tmpstr, "covid19.fasta")

covid19seq <- readDNAStringSet("covid19.fasta")
covid19seq
names(covid19seq)

```


## class 02

- `popset` 데이터베이스에서 선택한 각 아이템 (스터디별) 사용한 샘플 수 구하기
- `sapply` 및 사용자 정의 함수 `extract_statistics_count` 사용한 barplot 그리기 



```r
library(tidyverse)

mydatadf <- mydata %>% 
  rownames_to_column(var = "uid")

mysummary$`2279969783`$statistics$count[4]
mysummary$`2271163415`$statistics$count[4]
mysummary$`2252835834`$statistics$count[4]
## ..

sapply(mysummary, extract_statistics_count)


## build function
extract_statistics_count <- function(x){
  ## x == mysummary$`2279969783`
  z <- x$statistics$count[4]
  return(z)
}

cnt <- sapply(mysummary, extract_statistics_count)
cntdf <- data.frame(cnt) %>% 
  rownames_to_column("uid")

# mydata %>% 
#   rownames_to_column(var = "uid") %>% 
#   left_join(y=cnt, by="uid")

mydatadf
cntdf

mydata2 <- left_join(y=mydatadf, x=cntdf, by="uid") 


## barplot

ggplot(mydata2, aes(x=uid, y=cnt)) +
  geom_bar(stat = "identity") +
  theme(
    axis.text.x = element_text(angle=90)
  ) +
  ylab("Count")

```

- 코로나 바이러스 ID들에 대해서 정렬을 수행하고 `writePairwiseAlignments`를 이용한 가시화 


```r
covid <- data.frame(
species = c(rep("Human", 7), c("Civet", "Civet"), rep("Bat", 3), "Pangolin"),
coronavirus = c("SARS-CoV-2", "SARS-CoV-2", "SARS-CoV-1", "SARS-CoV-1", "SARS-CoV-1", "H-CoV-OC43", "MERS-CoV", "SARS-CoV", "SARS-CoV", "SL-CoV", "SL-CoV", "SL-CoV", "SL-CoV"),
isolate = c("Wuhan Hu-1", "USA-WA-1", "Urbani", "Tor2", "GD03T10013", "UK/London",  "EMC-2012", "SZ3", "Civet007", "ZXC21", "WIV16", "RaTG13", "MP789"),
year = c("2020", "2020", "2002", "2002", "2003", "2011", "2011", "2003", "2004", "2015", "2013", "2013", "2020"),
gbacc = c("NC_045512.2", "MN985325.1", "AY278741.1", "AY274119.3", "AY525636.1", "KU131570.1", "NC_019843.3", "AY304486.1", "AY572034.1", "MG772934.1", "KT444582.1", "MN996532.1", "MT084071.1"))
write.csv(covid, file = "covid_table.csv", quote = F, row.names=F)


covid19info <- read.csv("covid_table.csv")

str1 <- entrez_fetch("nuccore", covid19info$gbacc, rettype = "fasta")
str2 <- entrez_fetch("nuccore", covid19info$gbacc, rettype = "gb")
write(str1, "covid19seqs.fasta")
write(str2, "covid19seqs.gb")

library(Biostrings)

covid19seqs <- readDNAStringSet("covid19seqs.fasta")
covid19seqs

aln <- pairwiseAlignment(covid19seqs[1], covid19seqs[2])
aln

writePairwiseAlignments(aln, block.width=50)
writePairwiseAlignments(aln, file = "covidaln.txt", block.width=50)

summary(aln)

consensusMatrix(aln)[,1:10]
consensusString(aln)

```


- (참고) 문자열 나누기 코드, `sapply`와 사용자 정의 함수 `extract_first_element` 사용법 익히기 



```r
covid19seqs
names(covid19seqs)

strsplit("AT GC", split=" ")

strsplit(names(covid19seqs)[1], split=" ")
strsplit(names(covid19seqs)[1:2], split=" ")

tmpstr <- strsplit(names(covid19seqs), split=" ")
tmpstr[1:5]
sapply(tmpstr, extract_first_element)

extract_first_element <- function(x){
  return(x[1])
}
myacc <- sapply(tmpstr, extract_first_element)
mytitle <- sapply(tmpstr, function(x){
  z <- paste(x[-1], collapse=" ")
  return(z)
  })

mydf <- data.frame(myacc, mytitle)
names(covid19seqs) <- myacc
covid19seqs

aln <- pairwiseAlignment(covid19seqs[1], covid19seqs[2])
writePairwiseAlignments(aln, block.width=50)
```

## class 03

- `DECIPHER` 패키지 활용한 다중 서열 정렬 


```r

if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install("DECIPHER")

library(DECIPHER)

```

- `DECIPHER` 패키지 서열 관리 


```r

dbConn <- dbConnect(SQLite(), ":memory:")
Seqs2DB(covid19seqs, "XStringSet", dbConn, "covid19")
BrowseDB(dbConn)

coi <- readDNAStringSet("coi.fasta")
Seqs2DB(coi, "XStringSet", dbConn, "coi")
BrowseDB(dbConn)

l <- IdLengths(dbConn)
Add2DB(l, dbConn)
BrowseDB(dbConn)

covid19seq2 <- SearchDB(dbConn, identifier = "covid19")

dbDisconnect(dbConn)
BrowseDB(dbConn)
```

- `DECIPHER` 패키지 `AlignSeqs` 활용한 서열정렬 및 `BrowseSeqs` 이용한 가시화


```r

covid19sub <- subseq(covid19seq2, 1, 2000)
BrowseSeqs(covid19seq2, colWidth = 80, patterns=DNAStringSet(c("ACTG", "CSC")))


aln2 <- AlignSeqs(covid19sub)

BrowseSeqs(aln2, colWidth = 80)
BrowseSeqs(aln2, colWidth = 80, patterns=DNAStringSet(c("ACTG", "CSC")))
BrowseSeqs(aln2, colWidth = 80, patterns="-", colors="black")

```

- `RESTRICTION_ENZYMES` 활용 예제 


```r
data(RESTRICTION_ENZYMES)
RESTRICTION_ENZYMES

rsite <- RESTRICTION_ENZYMES["BsaI"]
d <- DigestDNA(rsite, covid19seq2[1])
unlist(d)

pos <- DigestDNA(rsite, covid19seq2[1], type="positions")
unlist(pos)

BrowseSeqs(covid19seq2[1], colWidth = 100, patterns=rsite)
library(stringr)
str_extract_all(rsite, "[[A-Z]]", simplify = T)


rsite2 <- paste(str_extract_all(rsite, "[[A-Z]]", simplify = T), collapse="")
rsite3 <- as.character(reverseComplement(DNAString(rsite2)))
BrowseSeqs(covid19seq2[1], colWidth = 100, patterns=c(rsite2, rsite3))

```

~~~
    아래 관련 알아보기 
    - DECIPHER AA 적용 여부
    - DECIPHER --> ggtree
~~~

- 웹에서 BLAST 실행 후 결과 데이터 분석


```r

myhit <- read.csv("F9BA23TY01R-Alignment-HitTable.csv", header = F)
myhit
myseq <- readDNAStringSet("seqdump.txt")
myseq

myaln <- AlignSeqs(myseq)
#AAStringSet(myaln)
BrowseSeqs(myaln)

myconm <- consensusMatrix(myaln)
mycons <- consensusString(aln)
dim(myconm)
myconm[1:10, 1:10]
```

## class 04

- ggtree 활용


```r
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install("ggtree")
BiocManager::install("tidytree")
BiocManager::install("genbankr")

```



```r

library(ggtree)
library(tidytree)


tree <- rtree(n = 20)
class(tree)           
methods(class="phylo")
ggtree(tree)
names(tree)

ggplot(tree) +
  geom_tree() +
  theme_tree()
  

ggtree(tree, branch.length="none")

ggtree(tree, layout="circular") +
  geom_tiplab(size=3, aes(angle=angle))

ggtree(tree, layout="circular", branch.length="none") +
  geom_tiplab(size=3, aes(angle=angle))

ggtree(tree) +
  theme_tree2()


ggtree(tree) +
  theme_tree2() +
  geom_tiplab() +
  geom_hilight(node=23, fill="steelblue", alpha=.4) 

as_tibble(tree) ## tidytree


ggtree(tree) +
  theme_tree2() +
  geom_tiplab() +
  geom_hilight(node=23, fill="steelblue", alpha=.4) 


dat <- data.frame(node=c(23, 35), type=c("AA", "BB"))

ggtree(tree) +
  theme_tree2() +
  geom_tiplab() +
  geom_hilight(dat, aes(node=node, fill=type), alpha=.4) +
  scale_fill_manual(values=c("steelblue", "darkgreen"))

```

- phylo 클래스 데이터를 tibble 형태로 변환 후 데이터 참고 


```r

as_tibble(tree) %>% 
  filter(label==c("t2", "t4"))

```

- genbank 타입 파일에서 정보 읽어오기 (다음시간 계속)


```r
library(rentrez)
require(genbankr)

mygb <- entrez_fetch("nuccore", id="NC_045512.2", rettype = "gb")
write(mygb, "co.gb")

mygbinfo <- readGenBank("co.gb")
class(mygbinfo)
methods(class="GenBankRecord")
chrseq <- getSeq(mygbinfo)
cdsrng <- cds(mygbinfo)
cdsrng
cdsseq <- getSeq(chrseq, cdsrng)
cdsseq
```



