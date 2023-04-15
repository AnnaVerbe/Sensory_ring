# Sensor_ring
You will find here the data files, matlab script and simulink file link to the modelisation of the sensory fusion in the hoverfly righting reflex. Simply run RA_model.m under matlab to simulate for the 5 experimental conditions. Note that the file requires the simulink toolbox.   

>  Linked article: Verbe, A., Martinez, D. & Viollet, S. Sensory fusion in the hoverfly righting reflex. Sci Rep 13, 6138 (2023). doi: https://doi.org/10.1038/s41598-023-33302-z


## Experimental data: 
```
Roll angle through time for all conditions and trials: see Tab_alltrials.xls
```

``` 
Mean roll angle (column 3) through time (column 2) for all the conditions (column 1): see Tab_Mean.xls
```
Conditions and corresponding numbers: 

* P A+ Vt: 1
* P A- Vb: 2
* P A+ Vb: 3
* P A+ Vdark: 4
* P A- Vdark: 5

``` 
Data from the model: see Tab_dataModel.xlsx
```

## RA-based model

```
File: RA_model.m 
```

> Ring attractor network with sigma-pi units, as described in the SI of (Verbe, Martinez & Viollet)

## Matlab/simulink-based model:
``` 
File: model_mouche_integrateur_V2.slxc
```
>  Dynamic model of the righting reflex. Adapted from: Anna Verbe, Léandre P. Varennes, Jean-Louis Vercher, Stéphane Viollet; How do hoverflies use their righting reflex?. *J Exp Biol* 1 July 2020; 223 (13): jeb215327. doi: https://doi.org/10.1242/jeb.215327






