# R/Rstudio basics

## What is R / Rstudio

![](images/01/r.jpg)

R은 통계나 생물통계, 유전학을 연구하는 사람들 사이에서 널리 사용되는 오픈소스 프로그래밍 언어 입니다. Bell Lab에서 개발한 S 언어에서 유래했으며 많은 라이브러리 (다른 사람들이 만들어 놓은 코드)가 있어서 쉽게 가져다 사용할 수 있습니다. R은 복잡한 수식이나 통계 알고리즘을 간단히 구현하고 사용할 수 있으며 C, C++, Python 등 다른 언어들과의 병행 사용도 가능합니다. R은 IEEE에서 조사하는 Top programming languages에서 2018년 7위, 2019년 5위, 2020년 6위, 2021년 7위로 꾸준히 높은 사용자를 확보하며 빅데이터, AI 시대의 주요한 프로그래밍 언어로 사용되고 있습니다.

![<https://spectrum.ieee.org/top-programming-languages/>](images/01/toplanguage2021.png){width="450px"}

R은 데이터를 통계분석에 널리 사용되는데 이는 데이터를 눈으로 확인하기 위한 visualization 이나 벡터 연산 등의 강력한 기능 때문에 점점 더 많은 사람들이 사용하고 있습니다. 기존에는 속도나 확장성이 다른 언어들에 비해 단점으로 지적되었으나 R 언어의 계속적인 개발과 업데이트로 이러한 단점들이 빠르게 보완되고 있습니다. R 사용을 위해서는 R 언어의 코어 프로그램을 먼저 설치하고 그 다음 R 언어용 IDE(Integrated Development Environment)인 RStudio 설치가 필요합니다.

![](images/01/rstudio.png)

Rstudio는 R 언어를 위한 오픈소스 기반 통합개발환경(IDE)으로 R 프로그래밍을 위한 편리한 기능들을 제공해 줍니다. R언어가 주목을 받고 두터운 사용자 층을 확보할 수 있게된 핵심 동력이 Rstudio 입니다. 자체적으로 최고수준의 오픈소스 개발팀이 있으며 `tidyverse`, \``,`shiny\` 등의 데이터 분석 관련 주요 패키지를 개발하였고 정기적으로 conference 개최를 하면서 기술 보급의 핵심 역할을 하고 있습니다.

![<https://www.rstudio.com/>](images/01/rstudiobottom.png){width="550px"}

## R / Rstudio Installation

### R 설치

-   R 사이트에 접속 후 (<https://www.r-project.org/>) 좌측 메뉴 상단에 위치한 CRAN 클릭.
-   미러 사이트 목록에서 Korea의 아무 사이트나 들어감
-   Download R for Windows를 클릭 후 base 링크 들어가서
-   Download R x.x.x for Windows 링크 클릭으로 실행 프로그램 다운로드
-   로컬 컴퓨터에 Download 된 R-x.x.x-win.exe 를 실행 (2022.5 현재 R 버전은 4.2.0).
-   설치 프로그램의 지시에 따라 R 언어 소프트웨어 설치를 완료

### Rstudio 설치

-   사이트에 접속 (<https://www.rstudio.com/>), 상단의 Products \> RStudio 클릭
-   RStudio Desktop 선택
-   Download RStudio Desktop 클릭
-   RStudio Desktop Free 버전의 Download를 선택하고
-   Download RStudio for Windows 클릭, 다운로드
-   로컬 컴퓨터에 다운로드된 RStudio-x.x.x.exe 실행 (2022.5 현재 RStudio Desktop 2022.02.2+485)
-   설치 가이드에 따라 설치 완료

## Rstudio interface

![](images/01/01-11.PNG){width="650px"}

-   기본 화면에서 좌측 상단의 공간은 코드편집창, 좌측 하단은 콘솔창
-   각 위치를 기호에 따라서 바꿀 수 있음 (View --\> Pane)

### Keyboard shortcuts

-   참고사이트
    -   <https://support.rstudio.com/hc/en-us/articles/200711853-Keyboard-Shortcuts>
    -   Tools --\> Keyboard shortcut Quick Reference (`Alt + Shift + K`)
-   코드편집창 이동 (`Ctrl + 1`) 콘솔창 이동(`Ctrl + 2`)
-   한 줄 실행 (`Ctrl + Enter`)
-   저장 (`Ctrl + S`)
-   주석처리 (`Ctrl + Shift + C`)
    -   또는 `#`으로 시작하는 라인
-   텝 이동 (`Ctrl + F11`, `Ctrl + F12`)
-   코드편집창 확대 (`Shift + Ctrl + 1`) 콘솔창 확대 (`Shift + Ctrl + 2`)
-   컬럼 편집 (`Alt + 마우스 드레그`)

::: rmdnote
### Exercise

1.  코드편집창에서 다음을 입력/실행하고 단축키를 사용하여 주석을 넣으시오

![](images/01/01-14.PNG)

-   단축키 `Ctrl + enter`로 코드 실행
-   단축키 `Ctrl + 2`로 커서 콘솔창으로 이동
-   `x`값 `x+y`값 확인
-   단축키 `Ctrl + 1`로 코드편집창 이동
-   단축키 `Ctrl + Shift + C` 사용


```r
# x <- 10
# y <- 20
```
:::

### Set a project

프로젝트를 만들어서 사용할 경우 파일이나 디렉토리, 내용 등을 쉽게 구분하여 사용 가능합니다. 아래와 같이 임의의 디렉토리에 `kribbR` 이라는 디렉토리를 생성하고 `lecture1` 프로젝트를 만듭니다. 

> File \> New Project \> New Directory \> New Project \> "kribbR" \> Create Project


시작할 때는 해당 디렉토리의 `xxx.Rproj` 파일을 클릭합니다. Rstudio 오른쪽 상단 프로젝트 선택을 통해서 빠르게 다른 프로젝트의 작업공간으로 이동할 수 있습니다. 


### Terminology

-   Session: R 언어 실행 환경
-   Console: 명령어 입력하는 창
-   Code: R 프로그래밍 변수/제어문 모음
-   Object types:
    -   vector: 값들의 모임 combine function `c()` EX: c(6, 11, 13, 31, 90, 92)
    -   matrix: 2D 형태 값들의 모임
    -   array: 1D, 2D, 3D, ... 형태 값들의 모임
    -   factor: 범주형 데이터
    -   data frame: 2D 형태 값들의 모임 (다른 타입 값 가능)
    -   list:
    -   function: 특정 기능 수행, [함수이름, 입력값 (arguments), 출력값 (return)] 으로 구성
-   Data (value) types:
    -   Integers
    -   doubles/numerics
    -   logicals
    -   characters.
-   Conditionals (조건, 제어):
    -   `if`, `==`, `&` (AND), `|` (OR) Ex: `(2 + 1 == 3) & (2 + 1 == 4)`
    -   `for`, `while`: 반복 수

## Getting help

R은 방대한 양의 도움말 데이터를 제공하며 다음과 같은 명령어로 특정 함수의 도움말과 예제를 찾아볼 수 있습니다.


```r
help("mean")
?mean
example("mean")
help.search("mean")
??mean
help(package="MASS")
```

또한 <https://www.rstudio.com/resources/cheatsheets/> 에서는 다양한 R언어의 기능을 한 눈에 알아볼 수 있게 만든 cheatsheet 형태의 문서를 참고할 수 있습니다. 

![](images/01/base-r_1.png){width="650px"} 
![](images/01/base-r_2.png){width="650px"}

## R packages and Dataset

R 패키지는 함수와 데이터셋의 묶음으로 다른 사람들이 만들어 놓은 코드나 기능을 가져와서 사용하므로써 코드 작성의 수고로움을 줄이고 편리하고 검증된 함수(기능)를 빠르게 도입하여 사용할 수 있다는 장점이 있습니다. 예를 들어 `sd()` 함수는 `stats` package에서 제공하는 함수로써 표준편차 계산을 위한 별도의 함수를 만들어서 사용할 필요가 없이 바로 (stats 패키지는 R 기본 패키지로) 별도 설치 없이 바로 사용 가능합니다.  


패키지를 구할 수 있는 가장 대표적인 사이트는 (저장소는) The Comprehensive R Archive Network (CRAN) <http://cran.r-project.org/web/views/> 이며 우리가 사용할 Bioconductor - <http://www.bioconductor.org/packages/release/bioc/> 저장소도 체계적으로 관리되고 많은 사람들이 활용하고 있는 사이트 중 하나입니다. 

![](images/01/01-18.png){width="350px"}

-   UsingR package installation

![](images/01/01-19.png){width="350px"}


```r
install.packages("UsingR")
```

-   UsingR package loading


```r
library(UsingR)
help(package="UsingR")
```

일반적으로 패키지 안에 관련된 데이터도 같이 저장되어 있으며 `data()` 함수를 이용해서 패키지 데이터를 사용자 작업공간에 복사해서 사용 가능합니다. 


```r
head(rivers)
length(rivers)
class(rivers)
data(rivers)
data(package="UsingR")
library(HistData)
head(Cavendish)
str(Cavendish)
head(Cavendish$density2)
```

## R programming

### Console calculator


```r
2 + 2
((2 – 1)^2 + (1 – 3)^2 )^(1/2)
2 + 2; 2 - 2
```

### Exercise

다음 공식들을 계산하는 R 코드를 작성하시오

$$ \sqrt{(4+3)(2+1)} $$

$$ 2^3 + 3^2 $$

$$ \frac{0.25 - 0.2}{\sqrt{0.2 (1-0.2)/100}}$$

### Variables and values

-   프로그래밍 언어의 공통적 개념 `변수`, `함수`, `자료형`, `조건문`, `반복문`\
-   Assignment operator ( `<-` OR `=` )
    -   Valid object name `<-` value
    -   단축키: `Alt + -` (the minus sign)


```r
x <- 2
y <- x^2 – 2*x + 1
y
x <- "two"  
some_data <- 9.8
```

-   내장 변수 Built-in variables


```r
pi
```

-   변수이름 작명법
    -   문자, 숫자, "\_", "." 등으로 구성
    -   대소문자 구분
    -   가독성, 의미있는 변수 이름
    -   길이 제한 없음


```r
i_use_snake_case <- 1
otherPeopleUseCamelCase <- 2
some.people.use.periods <- 3
And_aFew.People_RENOUNCEconvention <- 4
```

-   자동 완성 기능 (Tab completion) in RStudio
-   이전 명령은 콘솔에서 위 아래 화살표

### Exercise

1.  변수 `x`에 1, 3, 5, 7, 9를, 변수 `y`에 2, 4, 6, 8, 10을 저장하는 코드를 작성하시오

2.  앞서 변수 `x`와 `y`를 더한 값을 `z`에 저장하는 코드를 작성하시오

3.  변수 x에 "hello world!" 를 저장하고 x의 값을 출력하는 코드를 작성하시오

### Functions

함수(Function)란 사용자가 원하는 기능을 수행하는 코드의 모음으로서 반복적으로 쉽게 사용할 수 있도록 만들어 놓은 코드 입니다. 특정 데이터를 입력으로 받아 원하는 기능을 수행한 후 결과 데이터를 반환하는 구조를 가집니다. 함수는 일반적으로 다음과 같은 포멧으로 구현할 수 있습니다.


```r
my_function_name <- function(parameter1, parameter2, ... ){
  ##any statements
  return(object)
}
```

예를 들어 다음과 같은 `my_sine` 함수를 만들 수 있으며 parameter (매개변수)는 `x`이고 `y`는 반환값을 저장하는 지역변수 입니다.


```r
my_sine <- function(x){
	y <- sin(x)
	return(y)
}
```

-   내장 함수 (Built-in functions)


```r
x <- pi
sin(x)
sqrt(x)
log(x)
log(x, 10)
x <- c(10, 20, 30)
x + x
mean(x)
sum(x)/length(x)
```

### Exercise

1.  변수 `x`에 1, 3, 5, 7, 9를, 변수 `y`에 2, 4, 6, 8, 10을 저장하는 코드를 작성하시오

2.  `x`와 `y`를 더한 값을 `z`에 저장하는 코드를 작성하시오

3.  `mysum` 이라는 이름의 함수를 작성하되 두 변수를 입력으로 받아 더한 후 결과를 반환하는 코드를 작성하시오

4.  `mymean` 이라는 이름의 함수를 작성하되 두 변수를 입력으로 받아 평균을 구한 후 결과를 반환하는 코드를 작성하시오

<a rel="license" href="http://creativecommons.org/licenses/by-nc-nd/4.0/"><img src="https://i.creativecommons.org/l/by-nc-nd/4.0/88x31.png" alt="크리에이티브 커먼즈 라이선스" style="border-width:0"/></a><br />이 저작물은 <a rel="license" href="http://creativecommons.org/licenses/by-nc-nd/4.0/">크리에이티브 커먼즈 저작자표시-비영리-변경금지 4.0 국제 라이선스</a>에 따라 이용할 수 있습니다.
