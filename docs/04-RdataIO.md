# Data I/O {#rdataio}


일반적인 데이터 분석은 데이터 전처리(변환), 가시화, 모델링(통계분석)의 반복적인 수행으로 진행될 수 있습니다. R에서는 `data.frame` 형식의 데이터 타입이 주로 사용되며 (최근 `tibble`형식) 따라서 `data.frame` 기반의 데이터를 다루기 위한 다양한 함수를 익힐 필요가 있습니다. 이번 강의에서는 `data.frame` 데이터를 읽거나 쓰는 함수들과 함께 데이터 전처리를 (변환) 위한 함수들을 배워보겠습니다. 

앞에서 배웠던 데이터를 저장하는 object의 종류를 먼저 간략히 정리해 봅니다.

-   Vectors - 같은 타입의 데이터를 (Numeric, character, factor, ...) 저장한 오브젝트 타입으로 인덱스는 `[`, `]` 사용.
-   Lists - 여러개의 `vector`를 원소로 가질 수 있으며 이 원소들은 문자나 숫자 어떤 데이터 타입도 가능하고 길이가 달라도 됨. list의 인덱싱에서 `[` `]`는 리스트를 반환하고 `[[` `]]`는 vector를 반환함.
-   Matrices - 같은 타입의 데이터로 채워진 2차원 행렬이며 인덱스는 `[i, j]` 형태로 i는 row, j는 column 을 나타냄. 메트릭스의 생성은 `matrix` 명령어를 사용하며 왼쪽부터 column 값을 모두 채우고 다음 컬럼 값을 채워 나가는 것이 기본 설정이며 `byrow=T` 를 통해 row를 먼저 채울수도 있음. row와 column 이름은 `rownames`와 `colnames`로 설정이 가능하며 `rbind`와 `cbind`로 두 행렬 또는 행렬과 백터를 연결할 수 있음 ( `rbind`와 `cbind`의 경우 행렬이 커지면 컴퓨터 리소스 많이 사용함)
-   data.frame - `list`와 `matrix`의 특성을 모두 갖는 오브젝트 타입으로 `list`와 같이 다른 타입의 `vector`형 변수 여러개가 컬럼에 붙어서 `matrix` 형태로 구성됨. 단, `list`와는 다르게 각 변수의 길이가 (row의 길이) 같아야 함. `$` 기호로 각 변수들을 인덱싱(접근) 할 수 있고 matrix와 같이 `[i,j]` 형태의 인덱싱도 가능.



## Loading data into R

데이터 분석을 위해서 가장 먼저 할 일은 데이터를 R로 읽어들이는 것 입니다.  [Bioinformatics Data Skills by Vince Buffalo](http://2.droppdf.com/files/5aTvl/bioinformatics-data-skills.pdf)의 Chapter 8에서 소개한 데이터 중 일부인 [Dataset_S1_sub.txt](Dataset_S1_sub.txt)를 이용하겠습니다. 대부분의 텍스트 파일은 아래와 같이 `csv` 또는 `txt` 파일로 저장하여 메모장으로 열어 확인할 수 있으며 읽어올 때 구분자 (sep 파라메터) 나 header를 (header 파라메터) 읽을지 등을 옵션으로 지정할 수 있습니다.


```r
dat <- read.csv("Dataset_S1_sub.txt")
head(dat)
```

Dataset_S1_sub.txt 파일을 열어보면 다음과 같이 header와 ","로 구분되어 있는 것을 볼 수 있습니다. `read.csv` 함수의 도움말을 보면 이 함수의 파라메터 head와 sep이 기본값으로 `T`와 `,`로 되어 있는 것을 볼 수 있습니다. `read.csv` 외에도 `read.table`, `read.delim` 등의 함수를 이용해서 택스트 파일을 읽어올 수 있습니다.


```r
str(dat)
```

참고로 위 데이터는 wide형 데이터 입니다. wide형 데이터는 일반적인 excel에서의 데이터 형태로 column은 변수, row는 샘플이 저장됩니다. 만약 새로운 변수가 추가 되면 column 오른쪽에 붙어 wide하게 확장되고 데이터(샘플)이 추가되면 아래에 붙어서 row가 추가 됩니다.

![출처: Bioinformatics Data Skills by Vince Buffalo](images/04/wide.JPG){width="548"}

반면 long 형 데이터는 아래와 같이 일반적으로 3개의 컬럼을 갖습니다. 이 경우 변수든 샘플이든 새로운 데이터가 추가될 경우 아래로 확장됩니다. wide형과 long형에 대한 추가 설명은 다음 강의에서 진행하도록 하겠습니다. 

![출처: Bioinformatics Data Skills by Vince Buffalo](images/04/long.JPG){width="552"}

## writing data into a text file

읽어들이거나 분석한 결과에 대한 데이터는 `write.table` 또는 `write.csv` 함수를 사용하여 텍스트 파일의 형태로 저장할 수 있습니다. 이 경우 알아둘 파라메터는 `quote`, `row.names`, `col.names`, `sep` 등이 있습니다.


```r

write.table(dat, file="table_write.txt")
write.table(dat, file="table_write.txt", quote=F)
write.table(dat, file="table_write.txt", quote=F, row.names=F)
write.table(dat, file="table_write.txt", quote=F, row.names=F, sep=",")
write.table(dat, file="table_write.csv", quote=F, row.names=F, sep=",")

```

## Excel file read

텍스트 파일 외에 엑셀파일은 `readxl` 이라는 R 패키지를 활용하여 읽거나 쓸 수 있습니다. 패키지는 다음과 같은 방법으로 설치할 수 있으며 `read_excel` 이라는 함수를 사용해서 데이터를 읽어들일 수 있습니다.


```r
install.packages("readxl")
library(readxl)
```

실습 파일은 형광 세포를 배양하여 형광리더기를 이용해 얻어진 실제 데이터이며 [Exp_data.xls](Exp_data.xls) 에서 다운로드 받을 수 있습니다. `read_excel` 함수를 이용하여 파일의 내용을 읽어오면 기본 자료형이 tibble 입니다. tibble은 최근 많이 쓰이는 R object로 data.frame과 유사하나 입력값의 type, name, rowname을 임으로 바꿀 수 없다는 점이 다릅니다.


```r

dat <- read_excel("Exp_data.xls", sheet=1, skip = 0, col_names=T)

```

엑셀파일에는 두 종류의 ($OD600_{nm}$, Fluorescence) 데이터가 저장되어 있습니다. 첫 번째 sheet에는 다음처럼 wide형 데이터가 저장되어 있습니다.

![](images/04/excelfile01.PNG)



------------------------------------------------------------------------

<a rel="license" href="http://creativecommons.org/licenses/by-nc-nd/4.0/"><img src="https://i.creativecommons.org/l/by-nc-nd/4.0/88x31.png" alt="크리에이티브 커먼즈 라이선스" style="border-width:0"/></a><br />이 저작물은 <a rel="license" href="http://creativecommons.org/licenses/by-nc-nd/4.0/">크리에이티브 커먼즈 저작자표시-비영리-변경금지 4.0 국제 라이선스</a>에 따라 이용할 수 있습니다.
