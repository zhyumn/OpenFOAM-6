#include <stdlib.h>
#include <stdio.h>
#include </scratch.nike/srini237/anaconda3/include/python3.9/Python.h>
#include <iostream>
#include <time.h>

int main(int argc, char *argv[]){

    clock_t start,end,duration;
    setenv("PYTHONPATH",".",1);
    // PySys_SetArgv(argc, (wchar_t**)argv);
    // PyRun_SimpleString("import tensorflow as tf\n"
    //                    "print(tf.__version__)\n");

    
    Py_Initialize();
    PyErr_Occurred();
    PyErr_Print();
    

    double x = 2.2;
    std::cout<<"here 1 \n";

    
    PyObject *pName = PyUnicode_FromString((char*)"predictANN");
    PyObject *pModule = PyImport_Import(pName);
    std::cout<<"here 2 \n";

    //PyObject *pDict = PyModule_GetDict(pModule);    
    
    PyObject *pFunc = PyObject_GetAttrString(pModule,(char*)"model_pred");
    std::cout<<"here 3 \n";
    PyObject *pArg = PyTuple_Pack(1,PyFloat_FromDouble(x));
    // start = clock();
    PyObject *pResult = PyObject_CallObject(pFunc,pArg);
    std::cout<<"here 4 \n";
    double result = PyFloat_AsDouble(pResult);
    // end = clock();
    // duration = (end-start);
    // std::cout<<long(duration)<<"\t"<<CLOCKS_PER_SEC<<std::endl;
    std::cout<<"here 5 \n";

    std::cout<<result<<"\n";

    // // Clean up
    // Py_DECREF(pName);
    // Py_DECREF(pArg);
    // Py_DECREF(pResult);
    // Py_DECREF(pFunc);
    // Py_DECREF(pModule);

    // // Shut down the Python interpreter
    // Py_Finalize();

    return 0;


}