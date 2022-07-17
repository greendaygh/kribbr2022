--- 
title: 2022년 한국생명공학연구원 연구데이터 분석과정 R
#subtitle: R프로그래밍 1차
author: 합성생물학전문연구소 김하성
date: "2022-07-18"
output:
  pdf_document:
    toc: yes
    toc_depth: '2'
    number_sections: yes
    latex_engine: xelatex
  html_document:
    toc: yes
    toc_float: yes
    toc_depth: 2
    number_sections: yes
    df_print: paged
documentclass: book
bibliography:
- book.bib
- packages.bib
description: Lecture
biblio-style: apalike
link-citations: yes
csl: chicago-fullnote-bibliography.csl
site: bookdown::bookdown_site
mainfont : NanumGothic


---


# Introduction {#introduction}

## 강의 개요 {#Information}
- 목표: 생물 데이터 분석을 위한 R 사용법과 (Rstudio, Tidyverse, Bioconductor 포함) 프로그래밍 기술을 습득함
- 장소: 코빅 3층 전산교육장(1304호)
- 강사: 한국생명공학연구원 합성생물학전문연구단 김하성
- 연락처: 042-860-4372, haseong@kribb.re.kr 
- 강의자료: https://greendaygh.github.io/kribbr2022/


## 강의 계획 {#Schedule}

1. R 사용법 및 데이터 분석 기초 	5.19(목), 5.26(목)
2. R/Tidyverse 데이터 분석 중급 	6.9(목), 6.16(목)
3. R/Tidyverse 활용 데이터 가시화 	7.7(목), 7.14(목)
4. R/Bioconductor 활용한 바이오데이터 분석 기초	8.4(목), 8.11(목)
5. R/Bioconductor 활용한 NGS 데이터 분석 기초	9.1(목), 9.15(목)
6. R/Bioconductor 활용한 NGS 데이터 분석 및 Workflow 	10.6(목), 10.13(목)


## 참고 자료 {#References}

- [R 홈페이지](https://www.r-project.org/)
- [Rstudio 홈페이지](https://www.rstudio.com/)
- [Bioconductor](https://www.bioconductor.org/)
- [R 기본 문서들](https://cran.r-project.org/manuals.html) 
- [R ebooks](https://bookdown.org/)
- [Cheat Sheets](https://www.rstudio.com/resources/cheatsheets/)
- [RStudio Webinars](https://resources.rstudio.com/)
- [Shiny](http://shiny.rstudio.com/tutorial/)
- [Hadley github](https://github.com/hadley)
- [R for Data Science](https://r4ds.had.co.nz) 
- Using R for Introductory Statistics by John Verzani
  - Free version of [1st Edition](https://cran.r-project.org/doc/contrib/Verzani-SimpleR.pdf)
  - [Second edition](https://www.crcpress.com/Using-R-for-Introductory-Statistics-Second-Edition/Verzani/p/book/9781466590731)
- [Bioinformatics Data Skills](http://2.droppdf.com/files/5aTvl/bioinformatics-data-skills.pdf) by Vince Buffalo
- [Introductory Statistics with R](http://www.academia.dk/BiologiskAntropologi/Epidemiologi/PDF/Introductory_Statistics_with_R__2nd_ed.pdf) by Dalgaard
- 일반통계학 (영지문화사, 김우철 외)
