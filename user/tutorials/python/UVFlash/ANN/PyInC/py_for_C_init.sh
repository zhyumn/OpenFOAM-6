export PYTHON_LIB_PATH=/scratch.nike/srini237/anaconda3/pkgs/python-3.9.12-h12debd9_0/lib
#export NUMPY_LIB_PATH=/scratch.nike/srini237/anaconda3/pkgs/python-3.9.12-h12debd9_0/lib
export LD_LIBRARY_PATH=$PYTHON_LIB_PATH:$LD_LIBRARY_PATH
export LIBRARY_PATH=$PYTHON_LIB_PATH:$LIBRARY_PATH
export NUMPY_INCLUDE_PATH=/scratch.nike/srini237/anaconda3/lib/python3.9/site-packages/numpy/core/include/
export PYTHON_INCLUDE_PATH=/scratch.nike/srini237/anaconda3/include/python3.9/
export ONNX_LIB_PATH=/scratch.nike/srini237/OpenFOAM/temp/OpenFOAM-6/user/tutorials/python/UVFlash/ANN/onnx_runtime/onnxruntime-linux-x64-1.14.1/lib
export LD_LIBRARY_PATH=$ONNX_LIB_PATH:$LD_LIBRARY_PATH
export LIBRARY_PATH=$ONNX_LIB_PATH:$LIBRARY_PATH

##compiling add /scratch.nike/srini237/anaconda3/pkgs/python-3.9.12-h12debd9_0/lib/libpython3.9.so.1.0

## to check paths 
# python
# import sysconfig
# sysconfig.get_paths()