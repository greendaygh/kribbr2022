
# Tools for sequences

## Sequences from NCBI

전세계 연구자들이 서열 데이터를 분석하는데 가장 많이 이용하는 사이트 중 하나가 NCBI 이며 따라서 NCBI에서는 연구자들이 데이터베이스에 접근하기위한 편리한 방법을 제공하고 있고 그 중 하나가 Entrez 입니다. 

R에서도 Entrez 기능을 도입한 package들이 제공되고 있으며 그 중 하나가 `rentrez` 입니다. https://www.ncbi.nlm.nih.gov/books/NBK25500/ 이 곳의 Downloading Full Records 를 참고하시면 좋습니다. Entrez는 대략적으로 다음 9개의 유틸리티를 제공합니다. 


> EInfo (database statistics)  
> ESearch (text searches)  
> EPost (UID uploads)  
> ESummary (document summary downloads)  
> EFetch (data record downloads)  
> ELink (Entrez links)  
> EGQuery (global query)  
> ESpell (spelling suggestions)  
> ECitMatch (batch citation searching in PubMed)  

이 중 `ESerach`, `EPost`, `ESummary`, `EFetch` 등이 많이 사용하는 유틸이며 정보를 다운로드 받을 경우는 `EFetch` 를 주로 사용하게 됩니다. rentrez 는 위와 같은 NCBI Eutils API를 활용하여 R 환경에서 탐색이나 다운로드 등 NCBI 데이터베이스와 상호작용이 용이하도록 만들어 놓은 tool 입니다. [rentrez landing page](https://cran.r-project.org/web/packages/rentrez/vignettes/rentrez_tutorial.html) `entrez_dbs`명령은 NCBI에서 제공하는 데이터베이스의 리스트를 볼 수 있으며 특정 DB에 대한 설명은 `entrez_db_summary`를 사용하면 되겠습니다. `entrez_search`는 각종 키워드를 사용한 검색 기능을 제공합니다. 


```r
library(rentrez)
require(Biostrings)

entrez_dbs()
entrez_db_summary("nuccore")

covid_paper <- entrez_search(db="pubmed", term="covid19")
covid_paper$ids

names(covid_paper)
covid_paper$ids


covid_link <- entrez_link(db="all", id=covid_paper$ids, dbfrom="pubmed")
names(covid_link)
names(covid_link$links)
head(covid_link$links$pubmed_pubmed)

```

`entrez_search`에서 검색어를 입력하는 방식은 [이곳](https://cran.r-project.org/web/packages/rentrez/vignettes/rentrez_tutorial.html#building-search-terms)을 참고하세요. 검색으로 찾아진 특정 오브젝트(객체)에 대한 내용은 `entrez_summary` 함수를 사용하여 조회할 수 있으며 `extract_from_esummary`로 조회된 아이템들에 대한 정보를 추출할 수 있습니다. 특정 id에 대한 서열 등 다양한 타입의 데이터를 실제로 다운로드 받는 기능은 `entrez_fetch` 함수가 제공하고 있습니다.  `entrez_fetch` 함수의 `rettype` 옵션에서 지원하는 데이터 타입을 다운로드 받을 수 있으며  rettype (return type)의 자세한 정보는 [Eutils table](https://www.ncbi.nlm.nih.gov/books/NBK25499/table/chapter4.T._valid_values_of__retmode_and/) 또는 [NCBI Eutils](https://www.ncbi.nlm.nih.gov/books/NBK25499/) 페이지를 참고하시기 바랍니다. 



```r
# popset database is a collection of related DNA sequences derived from population
katipo_search <- entrez_search(db="popset", term="Latrodectus katipo[Organism]")
katipo_search$ids

katipo_summs <- entrez_summary(db="popset", id=katipo_search$ids)
names(katipo_summs)
katipo_summs$`41350664`
class(katipo_summs)
methods(class="esummary_list")

titles <- extract_from_esummary(katipo_summs, "title")
unname(titles)

print(katipo_summs)
katipo_summs$`1790798044`$gi


COI_ids <- katipo_search$ids[c(2,6)]
trnL_ids <- katipo_search$ids[4]
COI <- entrez_fetch(db="popset", id=COI_ids, rettype="fasta")
trnL <- entrez_fetch(db="popset", id=trnL_ids, rettype="fasta")

write(COI, "COI.fasta")
write(trnL, "trnl.fasta")

#library(Biostrings)
coi <- readDNAStringSet("COI.fasta")
trnl <- readDNAStringSet("trnl.fasta")
```


::: rmdnote
**Exercises **


뎅기바이러스 서열 4종에 대한 NCBI의 accession 번호가 다음과 같음 NC_001477, NC_001474, NC_001475, NC_002640 해당 DNA 서열을 fasta 형식으로 `nuccore` 데이터베이스에서 다운로드 하시오. (참고로 `strwrap` 함수 사용법을 익혀두면 좋습니다)




:::


::: rmdnote
**Exercises **

1. popset 데이터베이스에서 "Covid-19" 단어가 들어간 유전자 40개를 찾고 (`entrez_search`에서 `retmax=40` 옵션 사용) 이들의 요약 정보 중 title 속성을 출력하시오 (`entrez_summary`와 `extract_from_esummary` 함수 사용). 






2. 위 결과에서 찾아진 유전자들 각각이 몇 개의 서열 샘플에 (population) 대해서 연구된 것인지 각각의 서열을 fasta 형태로 다운로드 받고 샘플의 개수에 대한 `barplot`을 그리시오

- `summary_record` 결과를 받아서 `extract_from_esummary`로 title을 추출 후 `data.frame`으로 변환 
- `tidyverse`의 `rownames_to_column()` 함수로 uid 정보 변수로 변환, mydata 이름으로  저장
- `entrez_fetch` 함수로 모든 uid에 대한 샘플 서열 `fasta` 파일 다운로드 후 파일 저장 (`write`함수 사용)
- `readDNAStringSet` 함수로 읽은 후 앞서 title 정보 비교를 통해서 앞서 mydata 와 병합
- 각 uid 별로 몇 개의 서열 샘플이 있는지 정보를 추출 후 barplot 그리기 





:::


::: rmdnote
**Exercises **

[Comparative sequence analysis of SARS-CoV-2 suggests its high transmissibility and pathogenicity](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC7938774/) 논문을 참고하여 COVID-19 서열의 NCBI accession 번호를 찾고 `nuccore` 데이터베이스에서 `fasta` 포멧과 `genbank` 포멧의 정보를 다운로드 하겠습니다. 데이터는 "covid_table.csv" 파일에 저장되어 있습니다.  


```r
covid <- data.frame(
species = c(rep("Human", 7), c("Civet", "Civet"), rep("Bat", 3), "Pangolin"),
coronavirus = c("SARS-CoV-2", "SARS-CoV-2", "SARS-CoV-1", "SARS-CoV-1", "SARS-CoV-1", "H-CoV-OC43", "MERS-CoV", "SARS-CoV", "SARS-CoV", "SL-CoV", "SL-CoV", "SL-CoV", "SL-CoV"),
isolate = c("Wuhan Hu-1", "USA-WA-1", "Urbani", "Tor2", "GD03T10013", "UK/London",	"EMC-2012", "SZ3", "Civet007", "ZXC21",	"WIV16", "RaTG13", "MP789"),
year = c("2020", "2020", "2002", "2002", "2003", "2011", "2011", "2003", "2004", "2015", "2013", "2013", "2020"),
gbacc = c("NC_045512.2", "MN985325.1", "AY278741.1", "AY274119.3", "AY525636.1", "KU131570.1", "NC_019843.3", "AY304486.1", "AY572034.1", "MG772934.1", "KT444582.1", "MN996532.1", "MT084071.1"))
write.csv(covid, file = "covid_table.csv", quote = F, row.names=F)
```




```r
require(kableExtra)
#> Loading required package: kableExtra
covid19 <- read.csv("covid_table.csv")
kable_classic(kable(covid19))
```

<table class=" lightable-classic" style='font-family: "Arial Narrow", "Source Sans Pro", sans-serif; margin-left: auto; margin-right: auto;'>
 <thead>
  <tr>
   <th style="text-align:left;"> species </th>
   <th style="text-align:left;"> coronavirus </th>
   <th style="text-align:left;"> isolate </th>
   <th style="text-align:right;"> year </th>
   <th style="text-align:left;"> gbacc </th>
  </tr>
 </thead>
<tbody>
  <tr>
   <td style="text-align:left;"> Human </td>
   <td style="text-align:left;"> SARS-CoV-2 </td>
   <td style="text-align:left;"> Wuhan Hu-1 </td>
   <td style="text-align:right;"> 2020 </td>
   <td style="text-align:left;"> NC_045512.2 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> Human </td>
   <td style="text-align:left;"> SARS-CoV-2 </td>
   <td style="text-align:left;"> USA-WA-1 </td>
   <td style="text-align:right;"> 2020 </td>
   <td style="text-align:left;"> MN985325.1 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> Human </td>
   <td style="text-align:left;"> SARS-CoV-1 </td>
   <td style="text-align:left;"> Urbani </td>
   <td style="text-align:right;"> 2002 </td>
   <td style="text-align:left;"> AY278741.1 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> Human </td>
   <td style="text-align:left;"> SARS-CoV-1 </td>
   <td style="text-align:left;"> Tor2 </td>
   <td style="text-align:right;"> 2002 </td>
   <td style="text-align:left;"> AY274119.3 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> Human </td>
   <td style="text-align:left;"> SARS-CoV-1 </td>
   <td style="text-align:left;"> GD03T10013 </td>
   <td style="text-align:right;"> 2003 </td>
   <td style="text-align:left;"> AY525636.1 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> Human </td>
   <td style="text-align:left;"> H-CoV-OC43 </td>
   <td style="text-align:left;"> UK/London </td>
   <td style="text-align:right;"> 2011 </td>
   <td style="text-align:left;"> KU131570.1 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> Human </td>
   <td style="text-align:left;"> MERS-CoV </td>
   <td style="text-align:left;"> EMC-2012 </td>
   <td style="text-align:right;"> 2011 </td>
   <td style="text-align:left;"> NC_019843.3 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> Civet </td>
   <td style="text-align:left;"> SARS-CoV </td>
   <td style="text-align:left;"> SZ3 </td>
   <td style="text-align:right;"> 2003 </td>
   <td style="text-align:left;"> AY304486.1 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> Civet </td>
   <td style="text-align:left;"> SARS-CoV </td>
   <td style="text-align:left;"> Civet007 </td>
   <td style="text-align:right;"> 2004 </td>
   <td style="text-align:left;"> AY572034.1 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> Bat </td>
   <td style="text-align:left;"> SL-CoV </td>
   <td style="text-align:left;"> ZXC21 </td>
   <td style="text-align:right;"> 2015 </td>
   <td style="text-align:left;"> MG772934.1 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> Bat </td>
   <td style="text-align:left;"> SL-CoV </td>
   <td style="text-align:left;"> WIV16 </td>
   <td style="text-align:right;"> 2013 </td>
   <td style="text-align:left;"> KT444582.1 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> Bat </td>
   <td style="text-align:left;"> SL-CoV </td>
   <td style="text-align:left;"> RaTG13 </td>
   <td style="text-align:right;"> 2013 </td>
   <td style="text-align:left;"> MN996532.1 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> Pangolin </td>
   <td style="text-align:left;"> SL-CoV </td>
   <td style="text-align:left;"> MP789 </td>
   <td style="text-align:right;"> 2020 </td>
   <td style="text-align:left;"> MT084071.1 </td>
  </tr>
</tbody>
</table>

:::






## Align two sequences

Biostrings 패키지에는 다음과 같이 local, global alignment를 수행할 수 있는 함수를 제공하고 있습니다. 첫 번째 파라메터는 pattern이며 두 번째는 subject 로서 pattern은 query로서 해당 서열이 subject (target)에 있는지를 보는 것과 같습니다. 


```r
covid19

aln <- pairwiseAlignment(covid19[[1]], covid19[[2]])
alnseqs <- c(alignedPattern(aln), alignedSubject(aln))
class(aln)
class(alnseqs)
methods(class="PairwiseAlignmentsSingleSubject")

writePairwiseAlignments(aln, Matrix="BLOSUM62", block.width=10)

Views(aln)

```

## Multiple sequence alignment

Multiple sequence alignment(MSA) tool은 서열 데이터의 양과 계산량의 문제로 linux 기반 commandline 프로그램들이 많습니다. 대표적으로 [CLUSTAL-Omega](https://www.ebi.ac.uk/Tools/msa/clustalo/), [MUSCLE](https://www.ebi.ac.uk/Tools/msa/muscle/). window 기반 환경에서는 docker 등을 활용해서 관련 분석을 수행할 수 있습니다. 본 강의에서는 `DECIPHER` 패키지를 활용합니다. 


[DECIPHER](https://www.bioconductor.org/packages/release/bioc/html/DECIPHER.html) 패키지는 서열 alignment나 primer design 등을 수행할 수 있는 패키지로 다음과 같이 별도 메모리에 서열을 저장하고 빠르게 alignment를 수행할 수 있어서 중소 규모의 서열에 대한 분석으로 유용하게 사용될 수 있습니다. 



```r
library(DECIPHER)
dbConn <- dbConnect(SQLite(), ":memory:")
Seqs2DB(coi, "XStringSet", dbConn, "coi")
BrowseDB(dbConn)

l <- IdLengths(dbConn)
Add2DB(l, dbConn)
BrowseDB(dbConn)

Seqs2DB(trnl, "XStringSet", dbConn, "trnl")
BrowseDB(dbConn)

## extract sequences
dna <- SearchDB(dbConn, identifier="coi")
BrowseSeqs(dna)

dbDisconnect(dbConn)
```

앞서 다운로드 받은 Latrodectus katipo 서열 데이터를 비교하는 코드입니다. 


```r
coi <- readDNAStringSet("COI.fasta")
BrowseSeqs(coi)
alignedcoi <- AlignSeqs(coi)
BrowseSeqs(alignedcoi)
class(alignedcoi)

conseq <- ConsensusSequence(alignedcoi)
IUPAC_CODE_MAP

```



```r
BrowseSeqs(alnseqs)
BrowseSeqs(alnseqs, colWidth=200)
BrowseSeqs(alnseqs, colWidth=200, patterns = "TCCTGCCCGGGGCCT")
```





## Phylogenetic trees with clustering

DECIPHER 패키지에는 XStringSet 서열의 거리를 계산해주는 `DistanceMatrix` 함수가 있습니다. 이 함수를 이용하면 역시 같은 패키지에서 제공하는 `IdClusters` 함수를 이용해서 유사한 서열끼리 묶어주는 tree 를 만들 수 있습니다. 



```r
dm <- DistanceMatrix(alignedcoi)
class(dm)
dim(dm)
dm[1:2,1:2]

tree <- IdClusters(dm, cutoff=10, method="NJ", showPlot=TRUE, type="dendrogram")
class(tree)
methods(class="dendrogram")
plot(tree)

```

최근에는 ggplot 형태의 [ggtree](https://yulab-smu.top/treedata-book/chapter12.html), [reference](https://www.molecularecologist.com/2017/02/08/phylogenetic-trees-in-r-using-ggtree/)를 사용하면 쉽게 tree를 그릴 수 있습니다. 







