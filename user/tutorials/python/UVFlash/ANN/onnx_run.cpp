#include <iostream>
#include "/scratch.nike/srini237/OpenFOAM/temp/OpenFOAM-6/user/tutorials/python/UVFlash/ANN/onnx_runtime/onnxruntime-linux-x64-1.14.1/include/onnxruntime_cxx_api.h"
#include <vector>

    //g++ onnx_run.cpp -I/scratch.nike/srini237/OpenFOAM/temp/OpenFOAM-6/user/tutorials/python/UVFlash/ANN/onnx_runtime/onnxruntime-linux-x64-1.4.1/include/ /scratch.nike/srini237/OpenFOAM/temp/OpenFOAM-6/user/tutorials/python/UVFlash/ANN/onnx_runtime/onnxruntime-linux-x64-1.14.1/lib/libonnxruntime.so.1.14.1

int main(int argc, char* argv[]){
    Ort::Env env; //(ORT_LOGGING_LEVEL_WARNING, "example");

    Ort::SessionOptions session_options;
    session_options.SetGraphOptimizationLevel(GraphOptimizationLevel::ORT_ENABLE_ALL);
    Ort::Session session(env, "/scratch.nike/srini237/OpenFOAM/temp/OpenFOAM-6/user/tutorials/python/UVFlash/ANN/OnnxModel2.onnx", session_options);
    Ort::AllocatorWithDefaultOptions allocator;

    //std::vector<std::string> input_names = session.GetInputName(0,allocator);
 

    std::cout << "Number of inputs = " << session.GetInputCount() << std::endl;
    std::cout << "Number of outputs = " << session.GetOutputCount() << std::endl;
    //Ort::AllocatorWithDefaultOptions allocator;
    //auto* inputname = session.GetInputName(0);
    auto* inputname = "input";
    std::cout << "Input name : "<<inputname << std::endl;
    auto* outputname0 = "T";//session.GetOutputName(0,env);
    auto* outputname1 = "p";//session.GetOutputName(1);
    auto* outputname2 = "phi";//session.GetOutputName(2);
    auto* outputname3 = "c";//session.GetOutputName(3);
    std::cout << "Output names : "<< outputname0<<outputname1<<outputname2<<outputname3 << std::endl;

    auto inputshape = session.GetInputTypeInfo(0).GetTensorTypeAndShapeInfo().GetShape();
   std::vector<float> inputValues = {-4712.07,0.000766946,0.7156593087182554,0.2843406912817447};
    //std::vector<std::vector<float>> inputValues = {{-5650.61,0.0012,0.72,0.28},{18526.4824219,0.0010647943709,0.965139746666,0.034860227257}};
    
    int batchsize = 1;
    inputshape[0] = batchsize;

    // for(int i=0;i<inputshape.size();i++){
    //     std::cout<<inputshape[i]<<std::endl;
    // }
    
    auto memoryinfo = Ort::MemoryInfo::CreateCpu(OrtDeviceAllocator,OrtMemTypeCPU);
    auto inputOnnxTensor = Ort::Value::CreateTensor<float>(memoryinfo,inputValues.data(),inputValues.size(),inputshape.data(),inputshape.size());
    
    // std::vector<const char*> inputNames;
    // inputNames.push_back(inputname);
    //const char* const * inputNames = {inputname};
    const char* inputNames[] = {"input"};
    const char* outputNames[] = {"T","P","phi","c"};

    std::cout<<"here"<<std::endl;
    
    // std::vector<std::string> outputNames;
    // outputNames.push_back(outputname0);
    // outputNames.push_back(outputname1);
    // outputNames.push_back(outputname2);
    // outputNames.push_back(outputname3);

    auto outputValues = session.Run(Ort::RunOptions(nullptr),inputNames,&inputOnnxTensor,1*batchsize,outputNames,4*batchsize);
    auto& T = outputValues[0];
    auto& P = outputValues[1];
    auto& phi = outputValues[2];
    auto& c = outputValues[3];

    const auto* floatT = T.GetTensorMutableData<float>();
    const auto* floatP = P.GetTensorMutableData<float>();
    const auto* floatphi = phi.GetTensorMutableData<float>();
    const auto* floatc = c.GetTensorMutableData<float>();

    std::cout<<*floatT<<std::endl;
    std::cout<<*floatP<<std::endl;
    std::cout<<*floatphi<<std::endl;
    std::cout<<*floatc<<std::endl;

    return 0;
}

