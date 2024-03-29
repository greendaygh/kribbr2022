# Lecture Note

## Lecture 02 (0526)

설명이나 실습을 위해 `R 사용법 및 데이터 분석 기초 5.19(목), 5.26(목)` 강의를 진행하며 작성한 코드입니다.

**numeric vector**

-   Ctrl + Alt + i 누르면 코드청크 생성
-   커서를 해당 라인에 두구 Ctrl + Enter 누르면 해당 라인 실행


```r
2 + 2
((2-1)^2 + (1-3)^2 )^(1/2)
2 + 2; 2 - 2
```

- 기억할 단축키
    -   Ctrl + 1 : 편집창
    -   Ctrl + 2 : 콘솔창



```r
sqrt((4+3)*(2+1))
2^3 + 3^2
```

- 변수에 값 저장하기


```r
x <- 1
y <- 2
z <- x + y
```

- 변수 값 보기


```r
x
y
z
print(x)
cat(x)
```

- numeric vector 만들기


```r
x <- 1:5
x
y <- seq(1, 5, 1)
y <- seq(from=1, to=5, by=1)
y <- seq(to=5, from=1, by=1)
y <- seq(5, 1, -1)
?seq
y
```

- 연습문제 풀이


```r
odds <- seq(1, 100, by=2)
odds  
evens <- seq(2, 100, 2)
evens
```

- vector의 인덱싱 (첫번째 값의 인덱스는 1부터 시작)


```r

odds[1]
odds[1:10]
i <- 1:10
i
odds[i]
dim(odds)
length(odds)
```


```r
precip
?precip

head(precip)
str(precip)
dim(precip)

precip[1]
precip["Mobile"]
```

-   logical vector 설명
-   기본 그래픽 함수 이용하는 방법은 필요할때만 설명


```r
precip
length(precip)
plot(precip)
```

- `which` 함수 활용한 40 이상만 선택


```r
precip > 40
precip[precip > 40]
idx <- precip > 40
which(idx)
myprecip <- precip[which(idx)]
myprecip
plot(myprecip)
```

-  `which`함수 활용한 짝수 만들기 


```r
mynumers <- 1:1000
mynumers_res <- mynumers %% 2
i <- which(mynumers_res == 0)
evens <- mynumers[i]
evens
```

- 홀수 값을 저장하는 벡터 만들고 하나씩 샘플링 (`sample` 함수 사용)


```r
odds <- seq(1, 1000, 2)
length(evens)
length(odds)
?sample
mysample <- c(sample(evens, 1), sample(odds, 1))
print(mysample[1])

```

- 문자열 붙이기, `,`로 나누어진 벡터들 각각의 원소를 붙여줌 


```r
paste("X", "Y", "Z", sep="_")
paste("X", "Y", "Z", sep="")
paste("X", "Y", "Z", "X", "Y", "Z", "X", sep="")
```

- 여러 벡터에서 각각의 원소를 붙여주는 기능


```r
paste(c("X","Y"), 1:10, sep="")

gene_names <- paste("gene", 1:100, sep="")
```

-   collapse (하나의 벡터에서 해당 벡터의 원소들을 붙여주는 기능)


```r
paste(c("X", "Y"), collapse = "")
```

- 예제 (`sample()`, `paste()`, `rep()`)


```r

x <- sample(c("A", "C", "G", "T"), size=20, replace = T)
x2 <- paste(x, collapse = "")

myseq <- rep("", 5)
myseq[1] <- paste(sample(c("A", "C", "G", "T"), size=20, replace = T), collapse="")
myseq[2] <- paste(sample(c("A", "C", "G", "T"), size=20, replace = T), collapse="")
myseq[3] <- paste(sample(c("A", "C", "G", "T"), size=20, replace = T), collapse="")
myseq[4] <- paste(sample(c("A", "C", "G", "T"), size=20, replace = T), collapse="")
myseq[5] <- paste(sample(c("A", "C", "G", "T"), size=20, replace = T), collapse="")

myseq
```

-  예제 'TAAGTCT' 바코드를 각 서열의 3'에 붙여보기


```r
paste(myseq, c("TAAGTCT"), sep="")
strsplit("XYZ", split="")
```

-  `split` 함수 사용


```r
x <- strsplit("XYZ", split="")
class(x)
x
y <- unlist(x)
class(y)
y
```

- `factor` 간단한 설명


```r
x <- c("Red", "Blue", "Yellow", "Red", "Blue")
x
y <- factor(x)
y
y[1] <- "gold"

levels(y)[4] <- "gold"
y
y[1] <- "gold"
y

```

- 데이터에서 `factor`들 보기


```r
library(MASS)
Cars93
str(Cars93)
```

-  아미노산 및 각 아미노산에 해당하는 코돈 표현 예제


```r

aa <- c("Phe", "Leu", "Ser")
class(aa)
aa <- factor(aa)
class(aa)

aa <- list()
aa[[1]] <- c("UUU", "UUC")
aa
aa[[2]] <- c("UUA", "UUG", "CUA", "CUU", "CUG", "CUC")
aa
class(aa)
names(aa) <- c("Phe", "Leu")
aa
aa[[1]][1]
aa[[2]][3]
```

-  useful functions


```r
z <- sample(1:10, 100, T)
?sample
z
head(z)
sort(z)
order(z)
table(z)

```

**matrix**

-   연습문제 풀이 성적별 테이블 정렬


```r
mynum <- sample(1:100, 20, T)
mynum
score <- matrix(mynum, nrow = 10, ncol = 2)
score
myrowname <- paste("Name", 1:10, sep="")
myrowname
rownames(score) <- myrowname
score
colnames(score) <- c("Math", "Eng")

total_score <- score[,1] + score[,2]
total_score <- score[,"Math"] + score[,"Eng"]
sort(total_score, decreasing = T)
score
o <- order(total_score, decreasing = T)
o
score[o,]
```

**data.frame**


```r

math_score <- sample(1:100, 10, T)
eng_score <- sample(1:100, 10, T)
score <- data.frame(math_score, eng_score)
score
score$math_score
score$eng_score

```

-   추가 연습문제: 수학, 영어 성적을 더해서 `total_score`를 만들고 이 값을 기준으로 내림차순으로 score 데이터프레임을 정렬 하시오.


```r
total_score <- score$math_score + score$eng_score
o <- order(total_score, decreasing = T)
score[o,]

```

**list**


```r
score
class(score)
mynum
class(mynum)

z <- list()
z[[1]] <- score
z[[2]] <- mynum
z
names(z) <- c("dataframe", "numericvector")
z
z$dataframe
z$numericvector
```

-  리스트 만들때 미리 원소의 개수를 알고 있으면 그 원소의 개수에 맞게 생성해 주는 것이 좋음 `aa <- vector("list", 5)`
- 각 아미노산에 해당하는 코돈 길이가 달라도 list 형태로 저장 가능


```r
aa <- list()
aa[[1]] <- c("UUU", "UUC")
aa
aa[[2]] <- c("UUA", "UUG", "CUA", "CUU", "CUG", "CUC")
aa
class(aa)
names(aa) <- c("Phe", "Leu")
aa

as.data.frame(aa)
```

- 바람직한 데이터는 column은 변수, row는 샘플 구조의 데이터. 예를 들어서, 변수:먹이, 수명 --> 컬럼, 샘플: 마우스1, 마우스2 --> Row, 등


```r
aa <- c("Phe", "Leu")
codon <- c("UUU", "UUC", "UUA", "UUG", "CUA", "CUU", "CUG", "CUC")
data.frame(aa, codon)
aa <- c(rep("Phe", 2), rep("Leu",6))
codon <- c("UUU", "UUC", "UUA", "UUG", "CUA", "CUU", "CUG", "CUC")
data.frame(aa, codon)
```

**functions**



```r
source("myscript.R")
```

-   함수만들기

```{=html}
<!-- -->
```
    my_function_name <- function(parameter1, parameter2, ... ){
      ##any statements
      return(object)
    }


```r

## mynumers: numeric vector
mymean <- function(mynumers){
  #cat("Input numbers are", mynumers, "\n")
  numbers_mean <- sum(mynumers)/length(mynumers)
  #out <- paste("The average is ", numbers_mean, ".\n", sep="")
  #cat(out)
  return(numbers_mean)
}

```

-   함수 만든 후 loading 할 때 `{`, `}` 밖에서 Ctrl + Enter로 로딩


```r
mymean(c(1,2,3))

x <- c(1, 2, 3, 0.452, 1.474, 0.22, 0.545, 1.205, 3.55)
mymean(x)

mymean()
```

-   데이터 표준화 예제


```r
mysd <- function(x){
  numbers_sd <- sqrt(sum((x - mymean(x))^2)/(length(x)-1))  
  return(numbers_sd)
}

#x <- sample(1:100, 1000, T)
x <- rnorm(1000, 10, 5)
z <- (x - mymean(x))/mysd(x)
mymean(z)
mysd(z)
mymean(x)
mysd(x)

plot(density(x))
density(x)
```

- 코드비교


```r
x <- c(10, 20, 30)
x + 10

y <- rep(0, 3)
for(i in 1:3){
  y[i] <- x[i] + 10
}
y
```

**for**


```r
x <- 1:10
for(i in x){
  cat(i, "\n")
  flush.console()
}
```


```r
x <- 1
while(x <= 10){
  cat(x, "\n")
  flush.console()
  x <- x + 1
}

```

-  랜덤 서열만들기 예제


```r

x <- sample(c("A", "C", "G", "T"), size=20, replace = T)
x2 <- paste(x, collapse = "")

myseq <- rep("", 5)
myseq[1] <- paste(sample(c("A", "C", "G", "T"), size=20, replace = T), collapse="")
myseq[2] <- paste(sample(c("A", "C", "G", "T"), size=20, replace = T), collapse="")
myseq[3] <- paste(sample(c("A", "C", "G", "T"), size=20, replace = T), collapse="")
myseq[4] <- paste(sample(c("A", "C", "G", "T"), size=20, replace = T), collapse="")
myseq[5] <- paste(sample(c("A", "C", "G", "T"), size=20, replace = T), collapse="")


numseq <- 7
myseq <- rep("", numseq)
for(i in 1:length(myseq)){
  x <- sample(c("A", "C", "G", "T"), size=20, replace = T)
  myseq[i] <- paste(x, collapse="")
}
myseq

```

**Data transformation**

- `UsingR` 패키지의 babies 데이터셋을 적절히 변환하는 예제. `with`와 `within` 활용법 알아두기 (설명 안 함)


```r
library(UsingR)
head(babies)
str(babies)


new_babies <- within(babies, {
  gestation[gestation==999] <- NA
  dwt[dwt==999] <- NA
  smoke = factor(smoke)
  levels(smoke) = list(
    "never" = 0, 
    "smoke now" = 1, 
    "until current pregnancy" = 2,
    "once did, not now" = 3)
  })
str(new_babies)

fit <- lm(gestation~smoke, new_babies)
summary(fit) ## t-test 결과 
anova(fit)
```
