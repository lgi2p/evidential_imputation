# Imputation_credibiliste
This repo contains 3 Evidential K Nearest Neighbor versions
All versions use Dempster Shafer Framework and Dempster Rule to make prediction


Evidential K Nearest Neighbor Denoeux 1998 that use distance

Evidential K Nearest Neighbor Denoeux 1998 that use distance and labels uncertainties

Evidential K Nearest Neighbor that use only labels uncertainties 



Also we compare different imputation methods for time series 3 of them use uncertainty models 


EVID_IMP : Imputation method that use for each missing value, the first known value after (futur) and the first known value before (past)  

EVID_IMP_POND : Same method as EVID_IMP but that give a weights to the past and futur values depending on their distance from the missing value

Tree_IMP : Imputation that use Decision Tree for imputation 




