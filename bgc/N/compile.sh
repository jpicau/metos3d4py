# debug
#f2py -c model.F90 insolation.F90 -m N -DF2PY_REPORT_ON_ARRAY_COPY=1
f2py -c model.F90 insolation.F90 -m N --debug-capi
#f2py -c model.F90 insolation.F90 -m N

