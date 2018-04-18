# memoryPyaffy
This program was memory-efficient version of pyaffy

## Installation
To download all the examples, simply clone this repository:
```
git clone https://github.com/hongryulahn/memoryPyaffy
python pyaffy-0.3.2/setup.py install
```

## Requirement
Prepare proper a custom cdf file from
**http://brainarray.mbni.med.umich.edu/Brainarray/Database/CustomCDF/CDF_download.asp**

## Making expression table
```
python memoryRMA.py ATH1121501_At_TAIRG.cdf raw_data_cold exp.cold.raw
python memoryRMA.py ATH1121501_At_TAIRG.cdf raw_data_heat exp.heat.raw
```

## Modifying labels
```
python modifyLabel_removeAt.py exp.cold.raw sample2label.cold.txt exp.cold
python modifyLabel_removeAt.py exp.cold.raw sample2label.heat.txt exp.heat
