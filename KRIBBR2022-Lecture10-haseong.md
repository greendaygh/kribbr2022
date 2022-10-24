## Lecture 10 (0915)

**class 01**


명령창에서 fasterq-dump 실행시 에러가 발생할 경우 아래 실행  

~~~
vdb-config --interactive
~~~

Enable Remote Access 가 (X)로 선택되어 있으면 save 후 exit. tab 키로만 움직일 수 있음. Space 키가 선택 (ok). 


- Cmd 창에서 dir 입력하면 파일 리스트가 보임

~~~
prefetch --
~~~



```r

if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install("GEOquery")

```

- GDS 


```r
library(GEOquery)
.libPaths()

gds <- getGEO(filename = system.file("extdata/GDS507.soft.gz",package="GEOquery"))
gds  
class(gds)
methods(class="GDS")
Table(gds)
Columns(gds)

```

- GSM


```r
gsm <- getGEO(filename=system.file("extdata/GSM11805.txt.gz",package="GEOquery"))
methods(class=class(gsm))
head(Meta(gsm))
Table(gsm)
Columns(gsm)
```


- GPL


```r

gpl <- getGEO(filename=system.file("extdata/GPL97.annot.gz",package="GEOquery"))
class(gpl)
methods(class="GPL")
Table(gpl)

```


- GES class


```r
gse <- getGEO(filename=system.file("extdata/GSE781_family.soft.gz",package="GEOquery"))
methods(class=class(gse))
Meta(gse)
head(GSMList(gse))
gsm <- GSMList(gse)[[1]]
Meta(gsm)
Table(gsm)
Columns(gsm)

```


**class 02**


- ExpressionSet 


```r
myexp <- Table(gds)
mysample <- Columns(gds)
mygene <- Table(gpl)

myeset <- list()
myeset$exp <- myexp
myeset$sample <- mysample
myeset$gene <- mygene

myeset

```


- Download GSE


```r
gse2553 <- getGEO('GSE2553', GSEMatrix=TRUE)
gse2553
class(gse2553)
class(gse2553[[1]])
myeset <- gse2553[[1]]

exprs(myeset)
fData(myeset)
pData(myeset)

```



```r
gds <- getGEO(filename=system.file("extdata/GDS507.soft.gz",package="GEOquery"))
class(gds)
methods(class="GDS")
GPL(gds)
eset <- GDS2eSet(gds, do.log2=TRUE)
eset
fData(eset)
```




```r
if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install("SummarizedExperiment")
```



```r
if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install("airway")
```



```r
library(airway)

data(airway)
airway

myexp <- assay(airway)
dim(myexp)

myfeature <- rowRanges(airway)
length(myfeature)
myfeature[[1]]


mysample <- colData(airway)
mysample
```

- SummarizedExpriment 생성


```r
nrows <- 200
ncols <- 6 
counts <- matrix(runif(nrows * ncols, 1, 1e4), nrows)
dim(counts)

rowRanges <- GRanges(rep(c("chr1", "chr2"), c(50, 150)),
                     IRanges(floor(runif(200, 1e5, 1e6)), width=100),
                     strand=sample(c("+", "-"), 200, TRUE),
                     feature_id=sprintf("ID%03d", 1:200))

colData <- DataFrame(Treatment=rep(c("ChIP", "Input"), 3),
                     row.names=LETTERS[1:6])

se <- SummarizedExperiment(assays = list(counts = counts),
                           rowRanges = rowRanges,
                           colData = colData
                           )
se
```


**class 03**

- 표준정규분포 모수: 


```r
library(tidyverse)

x <- seq(-4, 4, 0.01)
y <- dnorm(x, 0, 1)
?dnorm
dat <- data.frame(x, y)
ggplot(dat, aes(x=x, y=y)) +
  geom_point()
```



```r
y2 <- dnorm(x, 0, 0.5)
dat <- data.frame(x, y, y2)

dat %>% 
  pivot_longer(-x) %>% 
  ggplot(aes(x=x, y=value, color=name)) +
  geom_point()

```



```r
nsample <- 16
x <- rnorm(nsample, 0, 1)
y <- rep(1, nsample)
hist(x)
xbar <- mean(x)
dat <- data.frame(x, y)
ggplot(dat, aes(x, y)) +
  geom_point()+
  geom_point(aes(x=mean(x), y=1), color="blue", size=5, shape=15)

```



```r
nsample <- 16
x <- rnorm(nsample*2, 0, 1)
y <- c(rep(1, nsample), rep(0.9, nsample))
g <- factor(c(rep(1, nsample), rep(2, nsample)))
dat <- data.frame(x, y, g)

ggplot(dat, aes(x, y)) +
  geom_point() +
  geom_point(aes(x=mean(x[1:nsample]), y=1), colour="blue", size=5, shape=15) +
  geom_point(aes(x=mean(x[(nsample+1):length(x)]), y=0.9), colour="blue", size=5, shape=15) +
  scale_y_continuous(limits=c(0, 1.2))

```


```r
nsample <- 16
nrep <- 10

x <- rnorm(nsample*nrep, 0, 1)
tmpy <- seq(0.1, 1, length.out=nrep)
y <- rep(tmpy, each=nsample)
## ?rep
g <- factor(y)

dat <- data.frame(x, y, g)
head(dat)
## sample means
dat_mean <- dat %>% 
  group_by(g) %>% 
  summarise(mean=mean(x))
head(dat_mean)

ggplot() + 
  geom_point(data=dat, aes(x, y)) +
  geom_point(data=dat_mean, 
             aes(x=mean, y=as.numeric(as.character(g))), 
             colour="blue", 
             size=5, 
             shape=15) +
  theme_bw()
  
```

- 1000개 샘플링 해서 평균 값을 모집단의 평균값이다 라고 하는 것 보다 100개 샘플링 후 평균을 10번 반복 후 (샘플의) 평균의 분포를 가지고 모집단의 평균을 말하는 것이 더 정확하다 


```r

head(dat)
head(dat_mean)

x <- seq(-4, 4, by=0.01)
# distribution of the samples 
y <- dnorm(x, mean(dat$x), sd(dat$x))
# distribution of the sample means
y2 <- dnorm(x, mean(dat_mean$mean), sd(dat_mean$mean))
dat2 <- data.frame(x, y, y2)
dat2_long <- dat2 %>% 
  pivot_longer(cols=c(y, y2))
head(dat2_long)

ggplot(dat2_long, aes(x=x, y=value, color=name)) +
  geom_line(size=1.2) +
  labs(y="Density") +
  theme_bw() +
  scale_color_manual(name="Type", 
                     labels=c("Samples", "Sample Means"),
                     values=c("red", "blue"))
```




```r
nsample <- 100
x <- rnorm(nsample, 100, 16)
hist(x)

z <- (90 - 100)/(10/sqrt(10))
z

x <- seq(-4, 4, 0.01)
y <- dnorm(x, 0, 1)
dat <- data.frame(x, y)
ggplot(dat, aes(x=x, y=y)) +
  geom_point() +
  geom_point(aes(x=z, y=0), size=5, color="blue")

```

**class 04**


- 모집단의 분포를 가정을 하고 샘플링을 통해서 통계량을 구한 후 해당 통계량이 가정한 모집단 분포의 어디에 위치하는지를 보고 유의성 판단 


```r

x1 <- seq(-6, 6, by=0.01)
y1 <- dnorm(x1, 0, 1)
z <- data.frame(x1, y1)
pval <- round(1-pnorm(2,0,1), 4)

ggplot(z) +
  geom_line(aes(x1, y1), color="purple") +
  geom_vline(xintercept = 2) +
  geom_hline(yintercept = 0) +
  geom_area(data=filter(z, x1 > 2), 
            aes(x1, y1), 
            fill="#80008055") +
  annotate("text", x=3, y=0.1, label=pval) 
```



```r

x1 <- seq(-6, 6, by=0.01)
y1 <- dnorm(x1, 0, 1)
x2 <- seq(-1, 11, by=0.01)
y2 <- dnorm(x2, 3, 1)
z <- data.frame(x1, y1, x2, y2)
#pval <- round(1-pnorm(2,0,1), 4)
ggplot(z) +
  geom_line(aes(x1, y1), color="blue") +
  geom_line(aes(x2, y2), color="red") +
  geom_vline(xintercept = 2) +
  geom_hline(yintercept = 0) +
  geom_area(data=filter(z, x1 > 2), 
            aes(x1, y1), 
            fill="#0000ff55") +
  geom_area(data=filter(z, x2 < 2), 
            aes(x2, y2), 
            fill="#ff000055") +
  annotate("text", x=3, y=0.1, label=pval) 

```



- example 


```r
mpg <- c(11.4, 13.1, 14.7, 14.7, 15, 15.5, 15.6, 15.9, 16, 16.8)

tstat <- (mean(mpg)-17)/(sd(mpg)/sqrt(length(mpg)))
tstat

x <- seq(-5, 5, length=100)
y <- dt(x, df=length(mpg)-1)
dat <- data.frame(x, y)
ggplot(dat, aes(x, y)) +
  geom_line() +
  geom_vline(xintercept = tstat)

pt(tstat, df=9)

```


- DEG (lecture note 06, class 03) 분석 예제 
- 아래 예제는 코드를 보충해서 본문으로 옮겼습니다 


```r
gds <- getGEO(filename = system.file("extdata/GDS507.soft.gz",package="GEOquery"))
gds
eset <- GDS2eSet(gds, do.log2=TRUE)
eset

myexp <- data.frame(exprs(eset))
myfeature <- fData(eset)
mypheno <- pData(eset)

table(mypheno$disease.state)

class(myexp)

## transformation
mydat1 <- myexp %>% 
  rownames_to_column() %>% 
  pivot_longer(cols = -rowname) %>% 
  pivot_wider(names_from = rowname, values_from = value)  

mydat2 <- mypheno %>% 
  dplyr::select(sample, disease.state) %>% 
  left_join(mydat1, by = c("sample" = "name"))


mymean <- mydat2 %>% 
  group_by(disease.state) %>% 
  summarise(across(where(is.numeric), mean))



mymean2 <- mymean %>% 
  pivot_longer(-disease.state) %>% 
  pivot_wider(names_from=disease.state, values_from = value)

ggplot(mymean2, aes(x=normal, y=RCC)) +
  geom_point()

```






