# matlab scripts for SCZ_treatment project

## 1. run fmriprep on server
```
S_run_fmriprep.sh
```

## 2. extract freesurfer(FS) output 
```
S_FS_table.m
```

## 3. QC freesurfer output
check freesurfer output using [ENIGMA-QC methods](https://enigma.ini.usc.edu/protocols/imaging-protocols/)


## 4. Run fitlme model to calculate group and group-by-wave effects
```
S_FS_region_fitlme.m  
S_FS_region_fitlme_subgroup.m 
```
*separate for medicated and unmedicated patients at recruitment*

## 5. Gene expression of gene set of interests (GOI) using GAMBA
```
S_GE_GOI.m
```

## 6. Brain changes correlated with GOI
```
S_FS_corr_gene_expression.m
S_GE_spin.m
```

## 7. Peripheral inflammatory markers
```
S_blood_measure.m
```
