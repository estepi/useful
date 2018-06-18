# useful
Some useful scripts for parsing files

Analyze NAs in VAST-TOOLS INCLUSION final Table: 
-We count how many NAs contain a quality score for every PSI value computed. If they have more than 3, we replace this PSI value with NAnew3
(step1)
- We also replace Intron Retention PSI value with NAI  if the adjusted pvalue for the "read umbalance" is below 0.05
(step2)
-Input is the default vast-tools table (INCLUSION-...)
-Output: desired output, step1: only replacing NAnew3, step2: replacing NAnew3+NAI


