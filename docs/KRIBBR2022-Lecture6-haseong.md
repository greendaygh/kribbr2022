# Lecture 06 Note

- 20220714

## Class 01

###  확인사항
- 1) package install (tidyverse) 가능여부 확인, 설치가 되지 않을경우 (https://cran.rstudio.com/ blablabla 접근 에러) 아래와 같이 `options`에 repos 변수 CRAN을 "http://" 로 바꿔줌



```r

options(repos = c(CRAN="http://cran.rstudio.com"))
install.packages("tidyverse")

```


- 2) dbplyr 설치 에러가 발생할 경우 R 재설치 (https://greendaygh.github.io/kribbr2022/lecture-5-note.html)
- 3) 아래 데이터 download 되는지 확인



```r
library(tidyverse)

myexp <- read.csv("https://raw.githubusercontent.com/greendaygh/kribbr2022/main/examples/gse93819_expression_values.csv", header=T)
```

- 4) 한글로된 디렉토리 이름 사용하지 않기
    - D 또는 C 드라이브에 영문 디렉토리 Rstudy 하나 만들기
    - Rstudio > File > new project > .. > "Rstudy" 디렉토리 설정 > "lecture6" 프로젝트이름 입력 
    - working directory 가 D:/Rstudy/lecture6 확인, Rmarkdown 파일 하나 생성 



- boxplot + violin 그리기


```r
library(tidyverse)
#install.packages("viridis")
library(viridis)

# create a dataset
mydata <- data.frame(name = c(rep("A", 100), rep("B", 100)),
                     value= c(rnorm(100, 10, 1), rnorm(100, 12, 1)))
mydata

ggplot(mydata, aes(x=name, y=value, fill=name)) +
  geom_violin() +
  geom_boxplot(width=0.2) 
  
mylab <- mydata %>% 
  group_by(name) %>% 
  summarise(num = n())


mylab <- paste(mylab$name, "\n n=", mylab$num, sep="")


ggplot(mydata, aes(x=name, y=value, fill=name)) +
  geom_violin() +
  geom_boxplot(width=0.2) +
  scale_x_discrete(labels = mylab) +
  scale_fill_viridis(discrete = TRUE) +
  theme(legend.position="none")
  

```

- boxplot for expression data


```r

myexp <- read.csv("https://raw.githubusercontent.com/greendaygh/kribbr2022/main/examples/gse93819_expression_values.csv", header=T)

myplotdata <- myexp %>% 
  pivot_longer(cols = everything()) 

ggplot(myplotdata, aes(x=name, y=value)) + 
  geom_boxplot()
```

- scale 패키지의 `comma` 를 사용해서 천단위 수 표현 


```r
library(scales)


ggplot(myplotdata, aes(x=name, y=value, fill=name)) + 
  geom_boxplot() +
  scale_fill_viridis(discrete = TRUE) +
  scale_y_continuous(trans="log",
                     breaks = c(0.1, 20, 1000),
                     labels = comma) +
  theme(axis.text.x = element_text(angle=90),
        legend.position = "none") +
  labs(title = "baplot", 
       subtitle = "subtitle",
       x = "xlab", 
       y = "ylab")
  

```


## Class 02

- 산포도


```r
library(UsingR)
#install.packages("UsingR")

mydata <- data.frame(mpg$displ, mpg$cty)

mpg %>% 
  dplyr::select(displ, cty) %>% 
  ggplot(aes(x=cty, y=displ)) +
  geom_point(position = "jitter")

mpg
```

- 추세선 


```r
mydata <- mpg %>% 
  dplyr::select(displ, cty, hwy) %>% 
  pivot_longer(cols = c(cty, hwy))

ggplot(mydata, aes(x=displ, y=value, color=name)) +
  geom_point(position = "jitter") +
  scale_color_manual(values = c("#ff0000", "blue"), 
                     labels = c("a", "b")) +
  geom_smooth(span=1) +
  labs(title="Cars") +
  theme(
    title = element_text(size=20),
    axis.title = element_text(size=12),
    legend.title = element_text(size=12)
  )
  

```

- `scale_size`로 point 특성 변경 



```r
mydata
ggplot(mydata, aes(x=displ, y=value, size=value, color=name)) +
  geom_point(position = "jitter") +
  scale_size(range = c(0.5, 8), guide=F) +
  scale_color_viridis(alpha=0.4, discrete = T) +
  theme_minimal()
```



```r
ggplot(mydata, aes(x=displ, y=value, size=value, fill=name)) +
  geom_point(position = "jitter", shape=21) +
  scale_size(range = c(0.5, 8), guide=F) +
  scale_fill_viridis(alpha=0.4, discrete = T) +
  theme_bw()
```

## Class 03

- expression data 두 그룹의 평균 발현양 비교 
- x: normal group mean vs. y: treatment group  mean
- plot: scatter plot
- 차이나는 유전자를 더 크게 (더 밝게) 그리기 



```r

myexp2 <- myexp[,1:10]

mygrp <- c(rep("N", 5), rep("T", 5))
samplemeta <- data.frame(sample_name = names(myexp2), group_name=mygrp)
samplemeta

```

- dplyr 패키지 (pivot_wider) 사용한 데이터 transpose


```r

tmp <- t(myexp2)
dim(tmp)
tmp[1:10, 1:10]
class(tmp)
tmp[,1]


myexp2t <- myexp2 %>% 
  rownames_to_column() %>% 
  pivot_longer(cols = -rowname) %>% 
  #pivot_wider(names_from = name, values_from = value)
  pivot_wider(names_from = rowname, values_from = value)


merdata <- samplemeta %>% 
  left_join(myexp2t, by=c("sample_name" = "name"))

```

- stat


```r
# mymean <- merdata %>% 
#   group_by(group_name) %>% 
#   summarise(across(where(is.numeric), mean))

# mysd <- merdata %>% 
#   group_by(group_name) %>% 
#   summarise(across(where(is.numeric), sd))
# 
# mymean2 <- mymean %>% 
#   pivot_longer(-group_name) 
# 
# mysd2 <- mysd %>% 
#   pivot_longer(-group_name) 
# 
# mydata <- mymean2 %>% 
#   left_join(mysd2, by=c("group_name", "name"))

```


- 두 그룹 각 유전자 발현의 평균 구하고 산포도 그리기 


```r
mymean <- merdata %>% 
  group_by(group_name) %>% 
  summarise(across(where(is.numeric), mean))

mymean2 <- mymean %>% 
  pivot_longer(-group_name) %>% 
  pivot_wider(names_from=group_name, values_from = value)


ggplot(mymean2, aes(x=N, y=T)) +
  geom_point()

```

- t-test 통한 pvalue 구하기 


```r
merdata %>% 
  group_by(group_name) %>% 
  summarise(across(where(is.numeric), mean))

mean(c(1:10))

tmp <- t.test(1:10, 2:11)
class(tmp)
#?help
names(tmp)
#attributes(tmp)
#attributes(tmp)$names
tmp$statistic
tmp$p.value
```




```r
mytest <- function(x){
  t.test(x[1:5], x[6:10])$p.value
}

merdata %>% 
  dplyr::select(-c("sample_name", "group_name")) %>% 
  apply(2, function(x){t.test(x[1:5], x[6:10])$p.value})
  

mypvalue <- merdata %>% 
  dplyr::select(-c("sample_name", "group_name")) %>% 
  apply(2, mytest)

myplotdata <- data.frame(mypvalue) %>% 
  rownames_to_column() %>% 
  left_join(mymean2, by = c("rowname" = "name"))
  
```



```r
myplotdata %>% 
  mutate(mlogp = -log(mypvalue)) %>% 
  ggplot(aes(x=N, y=T, size=mlogp)) +
  geom_point(shape=21, fill="black", alpha=0.2) +
  scale_x_continuous(trans = "log", breaks = c(1, 20, 400, 8000)) +
  scale_y_continuous(trans = "log", breaks = c(1, 20, 400, 8000))



myplotdata %>% 
  mutate(mlogp = -log(mypvalue)) %>% 
  ggplot(aes(x=N, y=T, size=mlogp, fill=mlogp)) +
  geom_point(shape=21) +
  scale_x_continuous(trans = "log", breaks = c(1, 20, 400, 8000)) +
  scale_y_continuous(trans = "log", breaks = c(1, 20, 400, 8000)) +
  scale_fill_viridis(alpha=0.2)
```


- 에러바 있는 barplot 


```r
mymean <- merdata %>%
 group_by(group_name) %>%
 summarise(across(where(is.numeric), mean)) %>%
 pivot_longer(-group_name)

mysd <- merdata %>%
 group_by(group_name) %>%
 summarise(across(where(is.numeric), sd)) %>%
 pivot_longer(-group_name)

mydata <- mymean %>%
 left_join(mysd, by=c("group_name", "name"))


ggplot(mydata, aes(x=name, y=value.x, fill=group_name)) +
  geom_bar(stat="identity", position="dodge", color="black") +
  scale_x_discrete(limits = c("1415670_at", "1415671_at")) +
  scale_y_continuous(limits = c(0, 3000)) +
  geom_errorbar(aes(ymin=value.x, ymax=value.x+value.y),
                position = position_dodge(width=0.9), 
                width = 0.2) +
  scale_fill_viridis(discrete = T, option = "D")
  

```

- 옆으로 눕히기 



```r
ggplot(mydata, aes(x=name, y=value.x, fill=group_name)) +
  geom_bar(stat="identity", position="dodge", color="black") +
  scale_x_discrete(limits = c("1415670_at", "1415671_at")) +
  scale_y_continuous(limits = c(0, 3000)) +
  geom_errorbar(aes(ymin=value.x, ymax=value.x+value.y),
                position = position_dodge(width=0.9), 
                width = 0.2) +
  scale_fill_viridis(discrete = T, option = "D") +
  coord_flip() 

```


```r

ggplot(mydata, aes(x=name, y=value.x, fill=group_name)) +
  geom_bar(stat="identity", position="dodge", color="black") +
  scale_x_discrete(limits = c("1415670_at", "1415671_at")) +
  scale_y_continuous(limits = c(0, 3000)) +
  geom_errorbar(aes(ymin=value.x, ymax=value.x+value.y),
                position = position_dodge(width=0.9), 
                width = 0.2) +
  scale_fill_viridis(discrete = T, option = "D") +
  coord_polar(start = 0)  

```

## Class 04

- heat map 


```r

val <- mpg$class
num <- 10
df <- expand.grid(y = 1:num, x = 1:num)
df
categ_table <- round(table(val) * ((num*num)/(length(val))))
categ_table
df$category <- factor(rep(names(categ_table), categ_table))  
df

ggplot(df, aes(x=x, y=y, fill=category)) + 
  geom_tile(color="black", size=0.5)
```

- myexp 데이터로 heatmap 


```r
myexp2 <- myexp %>% 
  rownames_to_column() %>% 
  slice(1:30) %>% 
  pivot_longer(-rowname, names_to = "sample_name") 
  



ggplot(myexp2, aes(x=rowname, y=sample_name, fill=value)) +
  geom_tile(color="gray", size=0.2) + 
  theme(
    axis.text.x = element_text(angle=90)
  ) +
  scale_fill_viridis_c()
  
  

```


- density plot


```r

myexp %>% 
  dplyr::select("GSM2462948", "GSM2462949") %>%
  rownames_to_column() %>% 
  pivot_longer(-rowname) %>% 
  ggplot(aes(x=value, fill=name)) + 
    geom_histogram() +
    scale_x_continuous(trans="log") 

```



```r
myplotdata <- myexp %>% 
  dplyr::select("GSM2462948", "GSM2462949") %>%
  rownames_to_column() %>% 
  pivot_longer(-rowname) 

ggplot(myplotdata, aes(x=value)) +
  geom_density(aes(fill=factor(name)), alpha=0.4) +
  scale_x_continuous(trans = "log") +
  scale_fill_viridis(discrete = T) +
  theme_bw()

#install.packages("ggridges")

```


- ggridges 


```r
library(ggridges)

ggplot(myplotdata, aes(x=value)) +
  geom_density_ridges(aes(fill=factor(name)), alpha=0.4) +
  scale_x_continuous(trans = "log") +
  scale_fill_viridis(discrete = T) +
  theme_bw()

```






