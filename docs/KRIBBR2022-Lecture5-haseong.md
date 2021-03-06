# Lecture 05 Note

- 20220707


## Class 01

- 수업 전 확인사항 `tidyverse` 로딩 --> 권한 문제로 인한 설치 에러는 조사 후 업데이트 예정 
- 아래 `https://greendaygh.github.io/kribbr2022/` 코드 다운로드 

### basic


```r
library(tidyverse)

myexp <- read.csv("https://raw.githubusercontent.com/greendaygh/kribbr2022/main/examples/gse93819_expression_values.csv", header=T)

```


- R재설치
  - r-project.org > CRAN > korea (제일 마지막) > download R for windows > base > download R 최신버전 
  - 디렉토리 실행 설치 (디폴트 옵션)
  - Rstudio > Tools > Global options > `[Default] [64-bit] C:\Program Files\R\R-4.2.1 ` 선택 > Rstudio 재실행 
  - install tidyverse


- ggplot 사용법 


```r
library(ggplot2)
head(iris)
str(iris)
ggplot(data=iris) + 
  geom_point(mapping=aes(x=Petal.Length, y=Petal.Width))

```



```r
ggplot(data=iris, mapping=aes(x=Petal.Length, y=Petal.Width)) +
  geom_point()
```

- aes 옵션은 다음 여섯 가지 주로 사용 x, y, group, color, fill, shape 
- aes 에 쓰인 옵션은 다른 그룹의 데이터일 경우 다른 모양으로 (컬러로) 표현 한다는 의미 


```r
ggplot(iris, aes(x=Petal.Length, 
                 y=Petal.Width, 
                 color=Species, 
                 shape=Species)) + 
  geom_point(color="black")
```

### bargraph

- geom_bar 기본 설정 stat = "count"


```r
dat <- data.frame(x1=rnorm(100))
str(dat)

ggplot(dat, aes(x=x1)) +
  geom_bar()

ggplot(dat, aes(x=x1)) +
  geom_bar(stat="bin", bins=30)
```



```r
x1 <- rnorm(10)
x2 <- rnorm(10)
dat <- data.frame(x1, x2)

ggplot(dat, aes(x=x1, y=x2)) +
  geom_bar(stat="identity")
```

- geom_xx 사용하는 어떤 레이어에서건 data, aes


```r
ggplot(dat, aes(x=x1, y=x2)) +
  geom_bar(stat="identity") +
  geom_point(aes(col="red", size=5))
```

- aes 안에서 정의된 옵션은 guide 붙음 


```r
ggplot(dat, aes(x=x1, y=x2)) +
  geom_bar(stat="identity") +
  geom_point(col="red", size=5)
```

- 이산형 데이터 --> barplot 
- 연속형 --> histogram (범위 지정 필요)


```r
x1 <- as.factor(1:3)
y1 <- tabulate(sample(x1, 100, replace=T))
dat <- data.frame(x1, y1)

ggplot(dat, aes(x=x1, y=y1)) +
  geom_bar()

ggplot(dat, aes(x=x1, y=y1)) +
  geom_bar(stat="identity")

ggplot(dat, aes(x=x1, y=y1, fill=x1)) +
  geom_bar(stat="identity") +
  xlab("Category") +
  ylab("Count") +
  ggtitle("Barplot") +
  guides(fill="none")

```


## Class 02


- 그룹별로 각 유전자의 발현의 평균을 bar 그래프로 비교 
- 데이터 기본: row: 하나의 샘플, col: 하나의 변수 (특히 data.frame 형태)


```r
myexp <- read.csv("https://github.com/greendaygh/kribbr2022/raw/main/examples/gse93819_expression_values.csv", header=T)
myexp
```


### Case1 (고전적 프로그래밍 방법)

- for문 사용 (느림), 저장공간 미리 준비
- 각 데이터들마다 동일한 코드 반복, 데이터 바뀔시 재사용성 낮음


```r
group1data <- as.matrix(myexp[,c(1:5)])
##myexp[,c(6:10)]
mymean1 <- rep(0, 5)
for(i in 1:5){
  mymean1[i] <- mean(group1data[i,], na.rm=T)
}
```


### Case 2 (apply)

- 비교적 빠른 처리 가능
- 통계량 계산 후 그래프용 데이터 재구성 필요
- apply 함수는 익숙해질 필요 있음


```r
group1data <- as.matrix(myexp[,c(1:5)])
mymean1 <- apply(group1data, 2, mean)
mymean1

```

### Case 3 (tidyverse)

- 가장 지향하는 방식 
- 분석 목적에 따른 샘플, 변수, 값 구분하기
- 목적: 그룹별로 각 유전자의 발현의 평균을 bar 그래프로 비교
- 발현의 평균을 계산해야 하므로 유전자가 변수가 (컬럼) 되야함
- 발현 데이터셋 3종: 1) main expression data 2) sample metadata 3) probe metadata (gene info)


```r
mygroup <- rep(c("A", "B", "C", "D"), each=5)
mysample <- data.frame(sample_name=names(myexp), group_name=mygroup)
mysample


myexpt <- t(myexp)
myexpt <- myexp %>% 
  rownames_to_column() %>% 
  pivot_longer(cols = -rowname) %>% 
  pivot_wider(names_from = rowname, values_from = value)



myexpt[1:5, 1:5]
dim(myexpt)
class(myexpt)
myexpt <- data.frame(myexpt)

```

- 두 데이터셋 통합하기 위해 같은 변수이름 사용 (다른 변수끼리도 통합하는데 기준 변수로 사용 가능)


```r
myexp2 <- myexpt %>% 
  rownames_to_column(var = "sample_name") %>% 
  left_join(mysample, by=c("sample_name")) %>% 
  select(group_name, everything()) %>% 
  group_by(group_name)

myexp2 %>% 
  summarise(mean(X1415670_at))


is.numeric(myexp2$X1415670_at)

myexp2 %>% 
  summarise(across(everything(), mean))

mymean <- myexp2 %>% 
  summarise(across(where(is.numeric), mean))
mymean

```

- ggplot 의 입력 데이터는 가능하면 long형 tidy 데이터
- long 형 데이터는 3가지 타입의 컬럼 (id, name, value) 를 가짐 (id 는 두개 이상 가능)


```r
myplotdata <- mymean %>% 
  pivot_longer(cols = -c("group_name"))
  
ggplot(myplotdata, aes(x=name, y=value, fill=group_name)) +
  geom_bar(stat="identity", position = "dodge") +
  scale_x_discrete(limit=c("X1415670_at", "X1415671_at")) +
  ylim(-10, 2000)

```


## Class 03

- error bar 그리기


```r
myexp2 <- myexpt %>% 
  rownames_to_column(var = "sample_name") %>% 
  left_join(mysample, by=c("sample_name")) %>% 
  select(group_name, everything()) %>% 
  group_by(group_name)

mymean <- myexp2 %>% 
  summarise(across(where(is.numeric), mean))

mysd <- myexp2 %>% 
  summarise(across(where(is.numeric), sd))

mymean2 <- mymean %>% 
  pivot_longer(-group_name, values_to = "mean")
```



```r
mysd2 <- mysd %>% 
  pivot_longer(-group_name, values_to = "sd")

myplotdata <- mymean2 %>% 
  left_join(mysd2, by=c("group_name", "name")) 

myplotdata$name[1:10]

ggplot(myplotdata, aes(x=name, y=mean, fill=group_name)) +
  geom_bar(stat="identity", position = "dodge") +
  geom_errorbar(aes(ymin=mean-sd, ymax=mean+sd), 
                position = position_dodge(width=0.9), 
                width=0.4, 
                color="#666666") +
  theme_minimal() +
  scale_x_discrete(limit=myplotdata$name[1:10]) +
  ylim(-10, 3000) +
  theme(axis.text.x = element_text(angle=90))

```

- 연속형 변수는 line graph 가능


```r
x1 <- c(12, 21, 40)
x2 <- c(33, 10, 82)
dat <- data.frame(x1, x2)
ggplot(dat, aes(x=x1, y=x2)) +
  geom_line()
```

- 범주형 변수는 기본적으로 barplot


```r
x1 <- as.factor(c(1:3))
y1 <- c(33, 10, 82)
dat <- data.frame(cate=x1, count=y1)
str(dat)

ggplot(dat, aes(x=cate, y=count, group="a")) +
  geom_bar(stat="identity") +
  geom_line() +
  geom_point(size=5)

```

- smooth line


```r
ggplot(mtcars, aes(x=mpg, y=hp)) +
  geom_point() +
  geom_smooth(method = "lm", se=FALSE)

ggplot(mtcars, aes(x=mpg, y=hp)) +
  geom_point() +
  geom_smooth(span=0.5)
```

- simulation 


```r

weights <- rnorm(200, 75, 5)
heights <- weights + rnorm(200, 100, 5)
classes <- sample(c("A", "B", "C", "D"), size=length(heights), replace = T)
mydata <- data.frame(heights, weights, classes)
str(mydata)

```

- 몸무게와 키의 산포도를 그리고 반별로 어떻게 상관성이 다른지 그리고 전체적으로 어떻게 다른지 알 수 있는 그래프를 그리시오 



```r

ggplot(mydata, aes(x=heights, y=weights)) +
  geom_point(aes(color=classes)) +
  geom_smooth()
```

## Class 04

- 두 그룹의 프루브들에 대해서 산포도를 그리기 


```r
myexp2 <- myexpt %>% 
  rownames_to_column(var = "sample_name") %>% 
  left_join(mysample, by=c("sample_name")) %>% 
  select(group_name, everything()) %>% 
  group_by(group_name)

mymean <- myexp2 %>% 
  summarise(across(where(is.numeric), mean))

```

- transpose (더 효율적 방법 찾아보기) 
- 그룹별 산포도 


```r
mymeant <- t(mymean)
class(mymeant)

mymeant <- data.frame(t(mymean)[-1,])
names(mymeant) <- mymean$group_name

myplotdata <- mymeant %>% 
  dplyr::select(A, B) %>% 
  mutate(Anum = as.numeric(A), Bnum = as.numeric(B)) %>% 
  dplyr::select(Anum, Bnum)

ggplot(myplotdata, aes(x=Anum, y=Bnum)) +
  geom_point()


ggplot(myplotdata, aes(x=Anum, y=Bnum)) +
  geom_point() +
  scale_x_continuous(trans = "log") +
  scale_y_continuous(trans = "log") +
  geom_smooth(method = "lm") +
  geom_abline(intercept = 0, slope = 1, color="red", size=2) 


mymeant %>% 
  dplyr::select(A, C) %>% 
  mutate(X = as.numeric(A), Y = as.numeric(C)) %>% 
  dplyr::select(X, Y) %>% 
  ggplot(aes(x=X, y=Y)) +
  geom_point() +
  scale_x_continuous(trans = "log") +
  scale_y_continuous(trans = "log") +
  geom_smooth(method = "lm") +
  geom_abline(intercept = 0, slope = 1, color="red", size=2) 


mymeant %>% 
  dplyr::select(A, D) %>% 
  mutate(X = as.numeric(A), Y = as.numeric(D)) %>% 
  dplyr::select(X, Y) %>% 
  ggplot(aes(x=X, y=Y)) +
  geom_point() +
  scale_x_continuous(trans = "log") +
  scale_y_continuous(trans = "log") +
  geom_smooth(method = "lm") +
  geom_abline(intercept = 0, slope = 1, color="red", size=2) 

```


- facet
- R에서 공식 표현하는 방법 y = ax + b, R expression: y ~ x


```r
ggplot(iris, aes(x=Petal.Length, y=Petal.Width)) + 
  geom_point(aes(color=Species, shape=Species)) +
  facet_wrap(~Species, nrow=2)
```


- iris 예제


```r
str(iris)
mycate <- factor(sample(c(0,1), nrow(iris), replace=T))
myiris <- data.frame(iris, mycate)
str(myiris)


ggplot(myiris, aes(x=Petal.Length, y=Petal.Width)) + 
  geom_point(aes(color=Species, shape=Species)) +
  facet_grid(mycate~Species)
```

- plot 을 변수에 저장할 수 있음


```r
myplot <- ggplot(myiris, aes(x=Petal.Length, y=Petal.Width)) + 
  geom_point(aes(color=Species, shape=Species)) +
  facet_grid(mycate~Species) +
  labs(x='Four classes',
       y='Number of students',
       title='Blood type distribution',
       subtitle = 'Blood type distribution from the 200 students',
       fill='Blood Types') 

myplot + 
  theme_bw() +
  #scale_color_brewer(palette="YlGnBu")
  scale_fill_manual(values = c("orange", "skyblue", "royalblue", "blue"))
```


- A, B, C, D 서로간의 비교에 대한 산포도 여러개 켄버스에 그리기



```r


mydat <- mymeant %>% 
  slice(1:1000) %>% 
  mutate(X1 = as.numeric(A), 
         X2 = as.numeric(B),
         X3 = as.numeric(C),
         X4 = as.numeric(D)) %>% 
  dplyr::select(starts_with("X")) 

tmp1 <- mydat %>% 
  select(x=X1, y=X2) %>% 
  mutate(group_name = "AB") 

tmp2 <- mydat %>% 
  select(x=X1, y=X3) %>% 
  mutate(group_name = "AC")

tmp3 <- mydat %>% 
  select(x=X1, y=X4) %>% 
  mutate(group_name = "AD")

tmp4 <- mydat %>% 
  select(x=X2, y=X3) %>% 
  mutate(group_name = "BC")

tmp5 <- mydat %>% 
  select(x=X3, y=X4) %>% 
  mutate(group_name = "CD")


tmp1 %>% 
  bind_rows(tmp2) %>% 
  bind_rows(tmp3) %>% 
  bind_rows(tmp4) %>% 
  bind_rows(tmp5) %>% 
  ggplot(aes(x=x, y=y)) +
    geom_point() +
    scale_x_continuous(trans = "log") +
    scale_y_continuous(trans = "log") +
    geom_smooth(method = "lm") +
    geom_abline(intercept = 0, slope = 1, color="red", size=2) +
    facet_wrap(~group_name)

```




