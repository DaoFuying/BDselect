# BinomialDistribution
A feature ranking tool based on binomial distribution

Binomial distribution is one of the wonderful feature selection techniques that have been successfully applied in many works. For the convenience of researchers, a feature ranking tool based on binomial distribution was developed.

#### 1.Installation：
*Download*： [github](https://github.com/)
environment(recommend python3.8):  
  ```
  pip3 install -r requirements.txt 
  ```  
#### 2. parameters:
 |parameters|description|

|-k |the number of kmer|
|-t |type of sequence(require:DNA ,RNA or protein)|   
|-f |the path of sequencefile|  
|——————————————————|————————————————| 

#### 3.Example
 ```
python3 BinomialDistribution.py -k 2 -t DNA -f C:\Users\DNASeqfile
or
BinomialDistribution.exe -k 2 -t DNA -f C:\Users\DNASeqfile
 ```
#### 4.Result
|file name|description|
|feature.csv | the csv feature file before sorting by BD|
|BDSortedfeature.csv | the csv feature file after sorting by BD|
|BinomialDistributionCL.txt | rank file|
|clresult.txt | the detail of CL values|


contact fydao@std.uestc.edu.cn
