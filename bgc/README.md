# BGC models

model name, module name, directory name,

if you use hyphens '-' in your model name,
the name will be internally changed to '_'
for instance: I-Cs, I_Cs

First comparison:

`metos3d4py`, or better `pymetos3d`, 288 time steps

```
Time:
  j: 00287, tj: 0.099653, delta: 0.419, s: 121.445                                          
```

`metos3d`, C implementation, 2880 time steps

```
jpicau@jserverpro:build> ./metos3d-simpack-I-Cs.exe test.I-Cs.option.txt 
       1.503s Metos3DInitWithFilePath                           filePath: test.I-Cs.option.txt                      
     100.791s 0000 Spinup Function norm 2.858675927513e+01 6.605300228544e+07
     102.301s Metos3DSolver
     102.428s Metos3DFinal
```

thus, the new implementation is ~12 times slower,


