<<<<<<< HEAD
# HTRgene
=======
# memoryPyaffy
>>>>>>> d992c9b14ffef796a8a09cd4d98362574e54f37c
This program was memory-efficient version of pyaffy

## Installation
To download all the examples, simply clone this repository:
```
<<<<<<< HEAD
git clone https://github.com/hongryulahn/HTRgene
python pyaffy-0.3.2/setup.py install 
=======
git clone https://github.com/hongryulahn/memoryPyaffy
python pyaffy-0.3.2/setup.py install
>>>>>>> d992c9b14ffef796a8a09cd4d98362574e54f37c
```

## Requirement
Prepare proper a custom cdf file from
**http://brainarray.mbni.med.umich.edu/Brainarray/Database/CustomCDF/CDF_download.asp**

## Making expression table
```
python memoryRMA.py ATH1121501_At_TAIRG.cdf raw_data_cold exp.cold.raw
<<<<<<< HEAD
=======
python memoryRMA.py ATH1121501_At_TAIRG.cdf raw_data_heat exp.heat.raw
>>>>>>> d992c9b14ffef796a8a09cd4d98362574e54f37c
```

## Modifying labels
```
<<<<<<< HEAD
python modifyLabel_removeAt.py exp.cold.raw sample2label.txt exp.cold 
```
=======
python modifyLabel_removeAt.py exp.cold.raw sample2label.cold.txt exp.cold
python modifyLabel_removeAt.py exp.cold.raw sample2label.heat.txt exp.heat
>>>>>>> d992c9b14ffef796a8a09cd4d98362574e54f37c
