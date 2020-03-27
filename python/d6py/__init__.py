

#IMPORT d6_python

try:
    import d6py.d6_python_3D
except:
    print('d6_python_3D not build with pybind11 or has not been found')

try:
    import d6py.d6_python_2D
except:
    print('d6_python_2D not build with pybind11 or has not been found')


from d6py.Tools import *