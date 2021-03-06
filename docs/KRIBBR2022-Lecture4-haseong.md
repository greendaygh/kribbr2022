# Lecture 04 Note

## Class 1

- 강의노트 주소 https://greendaygh.github.io/kribbr2022/index.html
- Rmarkdown 사용시 코드 청크 입력 Ctrl + Alt + i
- tidyverse 설치시 dbplyr 패키지 문제로 설치가 되지 않음. 우선 아래 4개 패키지 개별 설치 후 진행


```r

#install.packages("tidyverse")
install.packages("tibble")
install.packages("tidyr")
install.packages("dplyr")
install.packages("ggplot2")

```


- Shortcut Ctrl + enter: 커서 있는 라인 콘솔에서 실행 
- 인덱싱에 의한 subset과 `subset` 함수를 이용한 subset 


```r

library(UsingR)
str(babies)

head(babies)

## variables 

newbabies <- babies[,c("gestation", "dwt")]
head(newbabies)

newbabies[newbabies==999] <- NA
newbabies

df <- subset(babies, select=c(gestation, wt, dwt))
df
mean(df$gestation)
mean(df$wt)
mean(df$dwt)

apply(df, 2, mean)
?apply

```

- 인덱스 (숫자)로만 인덱싱하는 방식은 지양함 (가독률 낮음)
- `subset` 함수 사용할 수 있으며 `tidyverse` 패키지의 `dplyr::select` 사용 추천


```r
# 1st column
babies[,1]
babies[,c(1, 5)]

df <- subset(babies, select=c(id, gestation))

```

- `apply`, `lapply` 많이 씀



```r

airquality
str(airquality)
class(airquality)
grp <- factor(airquality$Month)
grp
class(grp)

airlist <- split(airquality, grp)
airlist
class(airlist)

# remove Month, Day
df <- subset(airquality, select=-c(Month, Day))
df
df2 <- split(df, grp)
length(df2)

lapply(df2, colMeans)

colMeans(df2$`5`, na.rm = T)

lapply(df2, colMeans, na.rm = T)
```

- `gse93819` 예제 데이터 두 그룹으로 나누고 각 유전자별 평균



```r

myexp <- read.csv("https://raw.githubusercontent.com/greendaygh/kribbr2022/main/examples/gse93819_expression_values.csv", header=T)

```

~~~
https://raw.githubusercontent.com/greendaygh/kribbr2022/main/examples/gse93819_expression_values.csv

위 주소 그대로 복사 후 브라우저 주소창에 입력
CTRL + A 로 모두 선택 후 CTRL+C 복사 
Rstudio > File > New file > Text file 생성 후 붙여넣기
CTRL + S "myexample1.csv"로 저장 

~~~

- 발현데이터를 둘로 나누고 각 그룹의 유전자 발현 평균 구하기 


```r

myexp <- read.csv("myexample1.csv", header=T)
str(myexp)

head(myexp)
myexp1 <- myexp[,1:10]
myexp2 <- myexp[,11:20]

myexp1mean <- apply(myexp1, 1, mean)
myexp2mean <- apply(myexp2, 1, mean)

myexpmean <- cbind(myexp1mean, myexp2mean)
myexpmean

```


- cbind 주의점: 이름이 달라도 병합이 되고 이런 문제 때문에 dplyr `join` 사용


```r

head(names(myexp1mean))
head(names(myexp2mean))
head(myexpmean)

names(myexp1mean)[1] <- "myexp"

head(names(myexp1mean))
head(names(myexp2mean))
myexpmean2 <- cbind(myexp1mean, myexp2mean)
head(myexpmean2)

plot(myexpmean)

class(myexpmean)
str(myexpmean)
df <- as.data.frame(myexpmean)
str(df)

plot(df)
mydiff <- df$myexp1mean - df$myexp2mean
mydiff
hist(mydiff, br=100)
```


- airquality long 형 변환
- pipe operator short cut: SHift + ctrl + m


```r
library(dplyr)
library(tidyr)

data(airquality)

airquality %>% head
airquality %>% str

airquality %>% 
   pivot_longer(c(Ozone, Solar.R, Wind, Temp))

myexpmeandf <- as.data.frame(myexpmean)
myexpmeandf %>% str
myexpmeandf %>% head

myexpmeandf %>% 
  pivot_longer(c(myexp1mean, myexp2mean))

```




```r
babies %>% str

newbabies <- babies %>% 
  dplyr::select(id, age, gestation, wt, dwt, smoke)
newbabies %>% str

newbabies %>% 
  filter(gestation != 999 & dwt != 999) 

```


## Class 2

- random sequence 생성 grep, grepl 실습


```r

mydna <- sample(c("A", "C", "G", "T"), 20, replace = T)
mydna <- paste(mydna, collapse = "")

n <- 100
mydna <- rep("", n) 
mydna

## case1
for(i in 1:n){
  mydna[i] <- paste(sample(c("A", "C", "G", "T"), 20, replace = T), 
                    collapse = "")
}
mydna

## case2
mydna <- rep("", n)
for(i in 1:n){
  mydna[i] <- sample(c("A", "C", "G", "T"), 20, replace = T) %>% 
    paste(collapse = "")
}
mydna

## case3
mydna <- sapply(1:n, function(x){
  sample(c("A", "C", "G", "T"), 20, replace = T) %>% 
    paste(collapse = "")
})
mydna

## return index
grep("ATG", mydna) 

grepl("ATG", mydna)

grep("Se", colnames(iris))
colnames(iris)[grep("Se", colnames(iris))]
```


- dplyr의 주요 함수들과 같이 사용되는 helper function 


```r

#newbabies %>% 
#  filter(man1 != 999 & man2 != 999.. )  x

#newbabies$gestation == 999
newbabies %>% 
  filter(!if_any(.fns = function(x){x==999}))

```

- NA 제외 후 새로 데이터 생성


```r

is.na(airquality$Ozone)
airquality$Ozone

mynewdata <- airquality %>% 
  filter(!if_any(.fns=is.na))
mynewdata

```

- arrange 기존 방법과 비교 


```r
mynewdata$Ozone
sort(mynewdata$Ozone)
i <- order(mynewdata$Ozone) # return index
i
mynewdata[i,]


mynewdata %>% 
  arrange(Ozone) 
```

- mutate 새로운 변수 추가 가능 


```r
class(myexpmean)
df <- as.data.frame(myexpmean)
class(df)
newdf <- df %>% 
  mutate(diff = sqrt((myexp1mean-myexp2mean)^2))
  
newdf$diff > mean(newdf$diff)

newdf <- df %>% 
  mutate(diff = sqrt((myexp1mean-myexp2mean)^2)) %>% 
  mutate(diffl = diff>mean(diff))

newdf %>% 
  filter(diffl==T)
```

- summarise 를 이용해서 타입별 평균과 표준편차 구하기 


```r
iris %>% str

iris %>% 
  group_by(Species) %>% 
  summarise(Sepal.Length.mean = mean(Sepal.Length),
            Sepal.Width.mean = mean(Sepal.Width))


irismean <- iris %>% 
  group_by(Species) %>% 
  summarise(across(.fns=mean))

irissd <- iris %>% 
  group_by(Species) %>% 
  summarise(across(.fns=sd))

```

- `cbind` 대신 `join` 사용 


```r
df1 <- data.frame(id=c(1,2,3,4,5,6), age=c(30, 41, 33, 56, 20, 17))
df2 <- data.frame(id=c(4,5,6,7,8,9), gender=c("f", "f", "m", "m", "f", "m"))

inner_join(df1, df2, by="id")
left_join(df1, df2, "id")
right_join(df1, df2, "id")
full_join(df1, df2, "id")

```

- ggplot mean


```r
irismean
irissd

irisplot <- irismean %>% 
  pivot_longer(-Species)

ggplot(irisplot, aes(x=Species, y=value, fill=name)) +
  geom_bar(stat="identity", position="dodge") 
```

- ggplot error bar


```r

d1 <- irismean %>% 
  pivot_longer(-Species)

d2 <- irissd %>% 
  pivot_longer(-Species)

df <- left_join(d1, d2, by=c("Species", "name"))
df
ggplot(df, aes(x=Species, y=value.x, fill=name)) +
  geom_bar(stat="identity", position="dodge") +
  geom_errorbar(aes(ymin=value.x-value.y, ymax=value.x+value.y),
                position=position_dodge(width=0.9), 
                width=0.4)

```

- `gse93819`를 4개 그룹으로 나누고 각 그룹별 평균, 표준편차 구하고 barplot 그리기
- 산포도 그리고 평균보다 2배 이상 차이나는 유전자를 골라서 색을 다르게 그리기기


```r
tmpid <- rep(c("A", "B", "C", "D"), each=5)
groupid <- paste(tmpid, 1:5, sep="")
groupid

colnames(myexp) <- groupid

data_meana <- myexp %>% 
  dplyr::select(starts_with("A")) %>% 
  apply(1, mean)
  
data_meanb <- myexp %>% 
  dplyr::select(starts_with("B")) %>% 
  apply(1, mean)
  
data_sda <- myexp %>% 
  dplyr::select(starts_with("A")) %>% 
  apply(1, sd)

data_sdb <- myexp %>% 
  dplyr::select(starts_with("B")) %>% 
  apply(1, sd)

mydata <- data.frame(data_meana,
           data_meanb,
           data_sda,
           data_sdb)

mydata

## for the rownames_to_column()
library(tibble)

mydataplot <- mydata %>% 
  rownames_to_column() %>% 
  pivot_longer(-rowname)
  
d1 <- mydataplot %>% 
  filter(grepl("data_mean",name)) %>% 
  slice(1:100) 

d2 <- mydataplot %>% 
  filter(grepl("data_sd",name)) %>% 
  slice(1:100) 

ggplot(d1, aes(x=rowname, y=value, fill=name)) +
    geom_bar(stat="identity", position="dodge") 
  

```








