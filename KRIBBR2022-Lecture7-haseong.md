## Lecture 07 (0804)

**class 01**


- 프로젝트 만들기 
    - Rstudio 실행 
    - File > New project > New directory > New project 
    - 디렉토리 `d:\rstudy` 지정 
    - 디렉토리 이름을 lecture7 입력 후 create project

- Rmarkdown 만들기 
    - File > New file > R markdown > 기본값 OK
    - title, output 남기고 다 지우기 
    - Ctrl + S 파일 저장


- Bioconductor 패키지 설치 


```r
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install("IRanges")
```



```r
if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install("Homo.sapiens")

```

- (참고) 코드청크에서는 커서를 코드에 두고 Ctrl + Enter 누르면 실행 


```r
library(IRanges)

vignette("IRanges")
browseVignettes("IRanges")

ir1 <- IRanges(start=1:10, width=10)
ir1

i
```


- OOP 설명


```r

df <- data.frame(x=c(1:5), y=LETTERS[1:5])
class(df)

class(df) <- c("data.frame", "myclass")
class(df)

```

- class에 따른 method 작성 예제, [function 참고](https://greendaygh.github.io/kribbr2022/r-programming.html#functions)


```r

x <- 1:10
y <- c("A", "B", "C", "A", "B")
class(x)
class(y)

mean(x)
table(y)/length(y)

mysummry <- function(z){
  if(class(z)=="integer"){
    retval <- mean(z)
  }else if(class(z)=="character"){
    retval <- table(z)/length(z)
  }
  return(retval)
}


mysummry(x)
mysummry(y)

?summary

```

- 패키지 사용법 모를 경우: class 확인, 해당 class에 사용되는 methods 확인


```r

library(Homo.sapiens)
class(Homo.sapiens)
?OrganismDb
methods(class="OrganismDb")

mygenes <- genes(Homo.sapiens)
class(mygenes)
columns(Homo.sapiens)
mygenes[1:10]

myexons <- exons(Homo.sapiens)
?GRanges
methods(class="GRanges")
gaps(myexons)

```




**class 02**

- Biostrings 설치 


```r
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install("Biostrings")

library(Biostrings)

x <- DNAString("ATGC")
x <- "ATGC"
x2 <- DNAString(x)
class(x2)


x2[2]

x <- DNAStringSet(c("ATTT", "CCCTA"))
x
?DNAStringSet

x[1]
class(x[1])

x[[1]]
class(x[[1]])

```

- 내장 변수


```r
DNA_BASES
AA_ALPHABET
```

- DNA 서열 만들기 


```r

x <- sample(DNA_BASES, 10, replace = T)
y <- paste(x, collapse = "")
DNAString(y)
y

reverseComplement(DNAString(y))
?reverseComplement

y2 <- DNAString(y)
letterFrequency(y2, c("G", "C"), as.prob=TRUE)
```

- 랜덤 DNA 생성 연습
- [pipe operator 참고](https://greendaygh.github.io/kribbr2022/tidyverse.html#pipe-operator)



```r
library(tidyverse)

x0 <- paste(sample(DNA_BASES, 10, replace = T), collapse = "")

x0 <- sample(DNA_BASES, 10, replace = T) %>%
        paste(collapse = "")

x <- paste("ATG", x0, "TAG", sep="") %>% 
      DNAString
x
```

- DNAstringset 설명


```r
x0 <- c("CTC-NACCAGTAT", "TTGA", "TACCTAGAG")
x <- DNAStringSet(x0)
x

length(x)
width(x)
nchar(x)
x
```

- 10개 DNAstring 만들기, 함수의 생성과 사용 
- [function 참고](https://greendaygh.github.io/kribbr2022/r-programming.html#functions)



```r

x0 <- sample(DNA_BASES, 30, replace = T) %>%
        paste(collapse = "") %>% 
        paste("ATG", ., "TAG", sep="") %>% 
        DNAStringSet

x1 <- sample(DNA_BASES, 30, replace = T) %>%
        paste(collapse = "") %>% 
        paste("ATG", ., "TAG", sep="") %>% 
        DNAStringSet

x2 <- sample(DNA_BASES, 30, replace = T) %>%
        paste(collapse = "") %>% 
        paste("ATG", ., "TAG", sep="") %>% 
        DNAStringSet

c(x0, x1, x2)


random_dna <- function(x){
  z <- sample(DNA_BASES, x, replace = T) %>%
        paste(collapse = "") %>% 
        paste("ATG", ., "TAG", sep="") %>% 
        DNAStringSet
  return(z)
}

random_dna(30)

```

- `apply`류 (`replicate`) 함수를 이용한 함수의 반복적인 사용
- [apply 참고](https://greendaygh.github.io/kribbr2022/data-transformation.html#apply)


```r

for(i in 1:10){
  cat(i)
}

z <- rep("", 10)
for(i in 1:10){
  z[i] <- random_dna(30)
}
DNAStringSet(z)


z <- replicate(10, random_dna(30))
z2 <- do.call(c, z)
z2
names(z2) <- paste("dna", 1:length(z2), sep="")
z2
```

- ggplot 활용 barplot 그리기 


```r
gcratio <- letterFrequency(z2, c("G", "C"), as.prob=TRUE) %>% 
  rowSums()
df <- data.frame(names=names(z2), gcratio)
df
barplot(df$gcratio)

ggplot(df, aes(x=names, y=gcratio, fill=names)) +
  geom_bar(stat="identity") +
  scale_y_continuous(limits = c(0, 1)) +
  scale_fill_brewer(palette = "green") +
  theme_bw() +
  labs(title="GC ratio") +
  xlab("DNA") + 
  ylab("GC ratio")

```

- successiveViews 활용한 서열 보기 


```r
class(z2)
z2[1]
Views(z2[[2]], start=1, end=10)
Views(z2[[2]], start=1, width=10)

successiveViews(z2[[2]], width=c(10, 10, 10, 10))
successiveViews(z2[[2]], width=rep(10, 4))
#successiveViews(z2[[2]], width=rep(length(z2[[2]]), 4))

```



```r
x <- random_dna(994)
successiveViews(x[[1]], width=rep(40, nchar(x[[1]])/40))

```


**class 03**

- DNA sequence read and write


```r
names(x) <- "mysequence"
writeXStringSet(x, filepath = "myseq.fasta", format = "fasta")

mynewseq <- readDNAStringSet(filepath = "myseq.fasta")
mynewseq

z2
writeXStringSet(z2, filepath = "z2.fasta")
readDNAStringSet(filepath = "z2.fasta")
```

- yeast chr1 ORF 추정 
- [orfinder, ncbi](https://www.ncbi.nlm.nih.gov/orffinder/)


```r
data(yeastSEQCHR1)
yeastSEQCHR1
nchar(yeastSEQCHR1)

yeast1 <- DNAStringSet(yeastSEQCHR1)
writeXStringSet(yeast1, filepath = "yeast1.fasta")

```


- yeast chr1 코돈 비율 계산


```r
yeast1orf <- readDNAStringSet(filepath = "yeast1.cds")
yeast1orf
tricodon <- trinucleotideFrequency(yeast1orf)
whole_codon <- colSums(tricodon)

class(whole_codon)

df <- data.frame(whole_codon) %>% 
  rownames_to_column() 

ggplot(df, aes(x=rowname, y=whole_codon)) +
  geom_bar(stat="identity") +
  scale_y_continuous() +
  #scale_fill_brewer(palette = "green") +
  theme_bw() +
  labs(title="Yeast chr1 codon frequency") +
  xlab("Codon") + 
  ylab("Codon frequency") +
  theme(
    axis.text.x = element_text(angle=90)
  )

```


**class 04**

- 코돈 아미노산 변환, barplot 그리기
- [data analysis with tidyverse](https://greendaygh.github.io/kribbr2022/tidyverse.html#dplyr)



```r
df
GENETIC_CODE
AMINO_ACID_CODE

GENETIC_CODE[1]
GENETIC_CODE["TTT"]
GENETIC_CODE[df$rowname]

df2 <- df %>% 
  mutate(AA1 = GENETIC_CODE[rowname]) %>% 
  mutate(AA3 = AMINO_ACID_CODE[AA1]) %>% 
  group_by(AA3) %>% 
  summarise(val=sum(whole_codon))
df2

ggplot(df2, aes(x=AA3, y=val)) +
  geom_bar(stat="identity") +
  #scale_y_continuous() +
  #scale_fill_brewer(palette = "green") +
  theme_bw() +
  labs(title="Yeast chr1 AA frequency") +
  xlab("AA") + 
  ylab("AA frequency") +
  theme(
    axis.text.x = element_text(angle=90)
  )

```

- matchpattern 사용하기 


```r
class(yeast1[[1]])
hits <- matchPattern("ATG", yeast1[[1]], min.mismatch=0, max.mismatch=0)
hits

class(z2)
hits <- vmatchPattern("ATG", z2, min.mismatch=0, max.mismatch=0)
hits
```


