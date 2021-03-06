# Lecture 03 Note

## Class 1

### file write

- package 설치 
- package 사용하기 위해서는 `library` 로 불러와야 함
- ctrl + enter 누르면 코드청크 코드 실행 (커서가 위치한)


```r
library(UsingR)
batting
?batting

str(batting)

batting$playerID
batting$HR
batting$SO

mydata <- data.frame(batting$playerID, batting$HR, batting$SO)
str(mydata)
mydata

mydata <- data.frame(playerID = batting$playerID,
                     HR = batting$HR,
                     SO = batting$SO)
mydata
```

- file writing
- Arguments 이름을 지정할 경우 순서를 바꿔도 됨


```r
write.table(x=mydata, file="mydata.txt")
?write.table

write.table(mydata, 
            file = "mydata.csv", 
            quote = F, 
            row.names = FALSE,
            sep = ",")

```

### file read

- file read


```r
myread <- read.table("mydata.csv", sep=",", header = T)
myread
str(myread)

?read.table

myread$HR
```

- 상관계수
- 회귀모형


```r
plot(myread$HR, myread$SO)
mycor <- cor(myread$HR, myread$SO)
mycor

fit <- lm(myread$HR ~ myread$SO)

plot(myread$HR, myread$SO)
abline(fit)
text(50, 170, round(mycor,2))

```


- 엑셀파일 읽기


```r
library(readxl)
#read.table("mydata.xlsx")
mydf <- read_xlsx("mydata.xlsx")
str(mydf)
class(mydf)
mydf$playerID
```

## Class 2

### apply

- 반복작업


```r
#library(UsingR)

mydata <- data.frame(playerID = batting$playerID,
                     HR = batting$HR,
                     SO = batting$SO)


mean(mydata$SO)
mean(mydata[,3])

mean(mydata$HR)
mean(mydata[,2])

mymean <- rep(0, 2)
mymean <- c(0, 0)

for(i in 1:2){
  mymean[i] <- mean(mydata[,i+1])
}
mymean

# ctrl + shift + c 를 누르면 주석
# x <- 1:10
# for(i in x){
#   cat(i, "\n")
#   flush.console()
# }

```
- apply 사용 


```r
?apply

apply(mydata[,c(2,3)], 2, mean)
```

- airquality data example


```r
?airquality
data(airquality)
str(airquality)

airquality
grp <- airquality$Month
class(grp)
grpf <- factor(grp)
airlist <- split(airquality, grpf)
?split
class(airlist)
airlist
```

- list 설명


```r

a <- 1:100
b <- 11:111
class(a)
class(b)
length(a)
length(b)
mydf <- data.frame(a, b)

mylist <- list(a=a, b=b)
mylist
mylist$a
mylist$b

mylist[[1]]
mylist$a
```


- ozone 의 평균 구하는 함수 만들기


```r
length(airlist)
airlist$`9`
class(airlist[[5]])
mean(airlist[[5]]$Ozone)
airlist[[5]]$Ozone
?mean
mean(airlist[[5]]$Ozone, na.rm=T)
```


- list 의 ozone 별 평균


```r
airlist
mymean <- c(0,0,0,0,0)
mymean[1] <- mean(airlist[[1]]$Ozone, na.rm=T)
mymean[2] <- mean(airlist[[2]]$Ozone, na.rm=T)
mymean[3] <- mean(airlist[[3]]$Ozone, na.rm=T)
mymean[4] <- mean(airlist[[4]]$Ozone, na.rm=T)
mymean[5] <- mean(airlist[[5]]$Ozone, na.rm=T)
mymean


lapply(airlist, function(x){mean(x$Ozone, na.rm=T)})


myozone <- function(x){
  z <- mean(x$Ozone, na.rm=T)
  return(z)
}
lapply(airlist, myozone)

```

## Class 3

### graphics 

- 산포도 


```r
x <- c(1:100)
y <- x*2 + rnorm(100)
myxy <- data.frame(x,y)


plot(myxy)
plot(myxy$x, myxy$y)
plot(x=myxy$x, y=myxy$y)
plot(y~x, data=myxy)
```

- histogram


```r
x <- rnorm(100)
hist(x, 
     br=20, 
     xlim=c(-3,3), 
     main="Main text", 
     xlab="X label", 
     ylab="y label")

airquality$Wind
hist(airquality$Wind, br=50)
hist(airquality$Wind, br=10)
```

- boxplot


```r
x <- rnorm(100)
class(x)
boxplot(x)


mydf <- airquality[,c(1, 2, 3, 4)]
mydf <- airquality[,c("Ozone", "Solar.R", "Wind", "Temp")]
class(mydf)
boxplot(mydf)

```

- barplot


```r
x <- sample(1:12, 200, replace = T)
x
tab_x <- table(x)


y <- sample(1:12, 200, replace = T)
tab_y <- table(y)
tab_xy <- rbind(tab_x, tab_y)
barplot(tab_xy)
barplot(tab_xy, beside = T)
barplot(tab_xy, beside = T, col=c("darkblue","red"))
barplot(tab_xy, beside = T, col=c("darkblue","red"), xlab="Month")
barplot(tab_xy, beside = T, col=c("darkblue","red"), xlab="Month", horiz=TRUE, legend.text = c("x", "y"))

```



```r
x <- rnorm(500)
hist(x, 100)
y <- 2*x + rnorm(500, mean=5, sd=1)
z <- c(x,y)
hist(z, br=100)


hist(z, br=100, probability = T)
zd <- density(z)
lines(zd)

```



```r
x <- rnorm(500)
y <- 2*x + rnorm(500, mean=5, sd=1)
myxy <- data.frame(x, y)
myxy

plot(x, y, data=myxy, xlim=c(-5, 5), ylim=c(-5, 15), pch=3)
idx <- which(x<0)
points(x[idx], y[idx], col="red")
fit <- lm(y~x)
abline(fit)


plot(y~x, data=myxy, xlim=c(-5, 5), ylim=c(-5, 15), pch=3)
idx <- which(x<0)
points(myxy[idx,], col="red")
fit <- lm(y~x, data=myxy)
abline(fit)

```



### tidyverse

- 우선 필요한 패키지 설치 및 로딩


```r
#library(tidyverse)
#install.packages("tibble")
#install.packages("dplyr")
#install.packages("tidyr")
library(tibble)
library(dplyr)
library(tidyr)
```


- code chunk shortcut CTRL + ALT + I


```r
df1 <- data.frame(x = 1:3, y = 3:1)
class(df1)
df1

df2 <- tibble(df1)
class(df2)
```



```r
airquality

myair <- airquality[1:5,]
myair_long <- pivot_longer(myair, cols = c("Ozone", "Solar.R", "Wind", "Temp"))
myair_long 
myair_long2 <- pivot_longer(myair, c(Ozone, Solar.R, Wind, Temp))
myair_long2 
myair_long3 <- pivot_longer(myair, !c(Month, Day))
myair_long3

?pivot_longer

myair_long <- pivot_longer(myair, 
                          c(Ozone, Solar.R, Wind, Temp), 
                          names_to = "Type", 
                          values_to = "Observation")

myair_long

stocks <- tibble(
  year   = c(2015, 2015, 2016, 2016),
  month  = c(   1,    2,     1,    2),
  profit = c(1.88, 0.59, 0.92, 0.17)
)

stocks
pivot_wider(stocks, names_from = year, values_from = profit)
?pivot_wider
```

A 41  O M D
A 190 S M D



- dplyr
- pipe operator 단축키 shift + ctrl + m


```r
library(dplyr)

x <- 1:100
y <- mean(x)
z <- sin(y)
sqrt(z)

sqrt(sin(mean(1:100)))

1:100 %>% mean %>% sin %>% sqrt

1:100 %>% 
  mean %>% 
  sin %>% 
  sqrt


x <- 1:5
paste(x, "a", sep="-")
x %>% 
  paste("a", sep="-") %>% 
  paste(collapse = ":")

```

- filter


```r

head(iris)

iris %>% head

iris %>% filter(Species=="setosa")
iris %>% filter(Species=="setosa" | Species=="versicolor")
iris %>% filter(Species=="setosa" & Species=="versicolor")
iris %>% 
  filter(Species=="setosa" | Species=="versicolor") %>% 
  str



iris %>% select(Species, everything()) %>% head(5)
iris %>% select(Species, everything())
iris %>% select(-Species)
iris %>% select(Petal.Length, starts_with('S'))
iris %>% select(starts_with('S'))
iris %>% select(obs = starts_with('S'))

iris
irisratio <- iris$Sepal.Length/iris$Sepal.Width
iris2 <- cbind(iris, irisratio)

iris2 <- iris %>% mutate(sepal_ratio = Sepal.Length/Sepal.Width)
head(iris2)



iris %>% summarise(m1 = mean(Sepal.Length), m2 = mean(Sepal.Width))
iris %>% 
  group_by(Species) %>% 
  summarise(mean(Sepal.Width))

```

- airquality 평균


```r

airquality %>% 
  group_by(Month) %>% 
  summarise(mean(Ozone, na.rm=T))

```

