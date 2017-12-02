
from numpy import f2py
fsource = open("model.F90").read()
s = f2py.compile(fsource, modulename="I_Cs", source_fn="model.F90", extension=".F90", verbose=False, extra_args="--quiet")
