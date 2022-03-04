# Some useful scripts for parsing files

##  Replace NAs: Analyze NAs in VAST-TOOLS INCLUSION final Table: 

* We count how many NAs contain a quality score for every PSI value computed. If they have more than 3, we replace this PSI value with NAnew3. 

* Input is the default vast-tools table (INCLUSION-...)
-Output: desired output, script 1: only replacing NAnew3

* If INCLUSION table was produced using last version of vast-tools, we can replace Intron Retention PSI value with NAI  if the adjusted pvalue for the "read umbalance" is below 0.05. 

* Input is the default vast-tools table (INCLUSION-...)
* Output: desired output, script 1: replacing NAnew3+NAI

## How to run:
- 1.- Download scripts in your desired folder.
- 2.- Load libraries:

library(data.table)
library(stringr)

- 3.- Source the function you want to use:

- source("replaceNA.R") # only NewNAs
- source("replaceNAI.R") # NewNAs + NAIs

* Define the input / output files:

- INCFile<-"test.tab"
- OUTFile<-"test_NA_NewN3.tab"
- OUTFile2<-"test_NA_NewN3_NAIs.tab"

- replaceNA(INCFile, OUTFile)
- replaceNAI(INCFile, OUTFile2)

