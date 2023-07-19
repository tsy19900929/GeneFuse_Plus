# GeneFuse_Plus
  * identical to output from GeneFuse 0.8.0
  * more flexible for high positive rate
  * [fix stuck from multi thread lock](https://github.com/OpenGene/GeneFuse/issues/30)


# New features
  -f, &emsp; list of csv files (string) [e.g.](https://raw.githubusercontent.com/tsy19900929/GeneFuse_Plus/master/csv.list)      
  -h, &emsp;file name to store HTML report; 1st csv file -> _1.html, 2nd csv file -> _2.html ...  (string [=genefuse.html])  
  -t, &emsp; worker thread number, default is 8 (int [=8])

# Issue
  `"version GLIBC_2.33 not found"`    
  [download source code then `make`] or [google] 
