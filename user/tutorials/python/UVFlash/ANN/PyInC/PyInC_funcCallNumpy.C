#define NPY_NO_DEPRECATED_API NPY_1_7_API_VERSION
#include <stdlib.h>
#include <stdio.h>
#include </scratch.nike/srini237/anaconda3/include/python3.9/Python.h>
#include </scratch.nike/srini237/anaconda3/lib/python3.9/site-packages/numpy/core/include/numpy/arrayobject.h>
#include <iostream>


int main(){

    setenv("PYTHONPATH",".",1);
    
    
    Py_Initialize();

    //init_numpy();

    double x = 2.2;
    std::cout<<"here 1 \n";

    
    PyObject *pName = PyUnicode_FromString((char*)"mymoduleNumpy");
    PyObject *pModule = PyImport_Import(pName);
    std::cout<<"here 2 \n";

    npy_intp dims[] = {10};
    double data[] = {0.0, 0.1, 0.2, 0.3, 0.4,
                     0.5, 0.6, 0.7, 0.8, 0.9};
    PyObject *pArg = PyArray_SimpleNewFromData(1, dims, NPY_DOUBLE, data);
    
    PyObject *pFunc = PyObject_GetAttrString(pModule,(char*)"myfunc");
    std::cout<<"here 3 \n";
    //PyObject *pArg = PyTuple_Pack(1,PyFloat_FromDouble(x));
    PyObject *pResult = PyObject_CallObject(pFunc,pArg);
    PyArrayObject *pArray = (PyArrayObject*)pResult;
    std::cout<<"here 4 \n";
    npy_intp *pArrayDims = PyArray_DIMS(pArray);
    double result = PyFloat_AsDouble(pResult);
    std::cout<<"here 5 \n";

    std::cout<<result<<"\n";

    // Clean up
    Py_DECREF(pName);
    Py_DECREF(pArg);
    Py_DECREF(pResult);
    Py_DECREF(pFunc);
    Py_DECREF(pModule);

    // Shut down the Python interpreter
    Py_Finalize();

    return 0;


}