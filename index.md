--- 
title: 2022년 한국생명공학연구원 연구데이터 분석과정 R
#subtitle: R프로그래밍 1차
author: 합성생물학전문연구소 김하성
date: "2022-10-24"
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
- 강의관련 게시판: https://github.com/greendaygh/kribbr2022/issues


## 강의 계획 {#Schedule}

1. (05.19, 13:30~17:30) [R/Rstudio] 사용법 및 데이터 구조
2. (05.26, 13:30~17:30) [R/Rstudio] 프로그래밍 기초 
3. (06.09, 13:30~17:30) [R/Tidyverse] 데이터 분석 기초 
4. (06.16, 13:30~17:30) [R/Tidyverse] 데이터 분석 중급  
5. (07.07, 13:30~17:30) [R/Tidyverse] 데이터 가시화 
6. (07.14, 13:30~17:30) [R/Tidyverse] 데이터 가시화 활용 
7. (08.04, 13:30~17:30) [R/Bioconductor] 바이오 데이터 분석 기초
8. (08.11, 13:30~17:30) [R/Bioconductor] 서열 비교 및 계통 분석 
9. (09.01, 13:30~17:30) [R/Bioconductor] 지놈 스케일 바이오 데이터 분석 
10. (09.15, 13:30~17:30) [R/Bioconductor] NGS 데이터 소개 및 통계 기초
11. (10.06, 13:30~17:30) [R/Bioconductor] NGS 기반 Differentially Expressed Genes 분석 
12. (10.20, 13:30~17:30) [R/Bioconductor] NGS 기반 Gene Set Enrichment Analysis 


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
- [Modern Statistics for Modern Biology](http://web.stanford.edu/class/bios221/book/index.html)
- [bios221](https://web.stanford.edu/class/bios221/labs/)
- [Annotation Workshop 2021](https://jmacdon.github.io/Bioc2021Anno/articles/AnnotationWorkshop.html#summarizedexperiment-objects-1)
- [CSAMA 2022](https://www.bioconductor.org/help/course-materials/2022/CSAMA/)
- [RNA-Seq CSAMA 2022](https://www.bioconductor.org/help/course-materials/2022/CSAMA/lab/2-tuesday/lab-03-rnaseq/rnaseqGene_CSAMA2022.html)
- [Annotation_Resources 2015](https://bioconductor.org/help/course-materials/2015/BioC2015/Annotation_Resources.html)
- [Sequence analysis](http://bioconductor.org/help/course-materials/2015/LearnBioconductorFeb2015/A01.3_BioconductorForSequenceAnalysis.html)
- 일반통계학 (영지문화사, 김우철 외)
