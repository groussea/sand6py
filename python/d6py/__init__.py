

#IMPORT d6_python

try:
    from d6py.d6_python import *
except:
    print('d6_python not build with pybind11 or has not been found')
    
from d6py.Tools import *