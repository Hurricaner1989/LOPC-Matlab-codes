LOPC Matlab package
========

07/18/16
--------
I reorganize the Data folder and make a few changes on the demo file. 

* Note falseDiscoveryRate.m depends on mafdr function from bioinformatics toolbox and firstParCorr.m, secParCorr.m and lopc.m depend on normcdf function from statistics and machine learning toolbox. Thanks for Dr. Chris Gaiteri to point this out. 
* For Matlab users without access to these two toolboxes, I'm working on rebuilding the codes to eliminate the dependences. More progress will come up soon.

07/12/16
--------
I make a few changes on the demo file. 

07/11/16 
--------
I submit version 1.0 which contains four Matlab functions (i.e., falseDiscoveryRate.m, firstParCorr.m, secParCorr.m, lopc.m) and one demo (i.e., demo.m and files in Data folder). 

* lopc.m implements our proposed LOPC algorithm in Methods paper.
* firstParCorr.m and secParCorr.m are traditional methods to calculate low order partial correlation up to first and second order, respectively. They are used for comparisons in demo. 
* falseDiscoveryRate.m implements a function to calculate adjusted p-value from a p-value matrix.




