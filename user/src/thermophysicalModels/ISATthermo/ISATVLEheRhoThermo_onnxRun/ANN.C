//ANN stuff
//namespace std
template <class BasicPsiThermo, class MixtureType>
Ort::Session Foam::ISATVLEheRhoThermo<BasicPsiThermo, MixtureType>::ANN_load(){
    Ort::Env env;
    Ort::SessionOptions session_options;
    session_options.SetGraphOptimizationLevel(GraphOptimizationLevel::ORT_ENABLE_ALL);
    Ort::Session session_(env,"/scratch.nike/srini237/OpenFOAM/temp/OpenFOAM-6/user/tutorials/python/UVFlash/ANN/OnnxModel2.onnx", session_options); 
    return session_; 
}
// template <class BasicPsiThermo, class MixtureType>
// void Foam::ISATVLEheRhoThermo<BasicPsiThermo, MixtureType>::ANN_load_BC(char* saved_model_dir_BC,char* tags_BC,int ntags_BC){
//     //ANN information and initialization
//     // char* saved_model_dir = "./ANN/test_model_64_64_32_ind_pbigger_1e-3/"; // Path of the model
//     // char* tags = "serve"; // default model serving tag; can change in future
//     // int ntags = 1;
//     std::cout<<"here 1"<<std::endl;
//     Session_BC = TF_LoadSessionFromSavedModel(SessionOpts_BC, RunOpts_BC, saved_model_dir_BC, &tags_BC, ntags_BC, Graph_BC, NULL, Status_BC);
//     if(TF_GetCode(Status_BC) == TF_OK)
//     {
//         printf("TF_LoadSessionFromSavedModel 2 OK\n");
//     }
//     else
//     {
//         printf("%s",TF_Message(Status));
//     }
    
// }

template <class BasicPsiThermo, class MixtureType>
std::tuple<scalar,scalar,scalar,scalar> Foam::ISATVLEheRhoThermo<BasicPsiThermo, MixtureType>::ANN_predict(float u, float v, float x1, float x2){
    //const char* model_dir = "/scratch.nike/srini237/OpenFOAM/temp/OpenFOAM-6/user/tutorials/python/UVFlash/ANN/OnnxModel.onnx";
    //std::cout<<"here inside predict"<<std::endl;
    //Ort::Session session(env, "/scratch.nike/srini237/OpenFOAM/temp/OpenFOAM-6/user/tutorials/python/UVFlash/ANN/OnnxModel.onnx", session_options);
    auto inputshape = session_UVX.GetInputTypeInfo(0).GetTensorTypeAndShapeInfo().GetShape();
    std::vector<float> inputValues = {u,v,x1,x2};
    int batchsize = 1;
    inputshape[0] = batchsize;
    auto memoryinfo = Ort::MemoryInfo::CreateCpu(OrtDeviceAllocator,OrtMemTypeCPU); 
    auto inputOnnxTensor = Ort::Value::CreateTensor<float>(memoryinfo,inputValues.data(),inputValues.size(),inputshape.data(),inputshape.size());
    const char* inputNames[] = {"input"};
    const char* outputNames[] = {"T","P","phi","c"};
    //std::cout<<"here before run"<<std::endl;
    auto outputValues = session_UVX.Run(Ort::RunOptions(nullptr),inputNames,&inputOnnxTensor,1,outputNames,4);
    auto& T = outputValues[0];
    auto& P = outputValues[1];
    auto& phi = outputValues[2];
    auto& c = outputValues[3];

    float* floatT = T.GetTensorMutableData<float>();
    float* floatP = P.GetTensorMutableData<float>();
    float* floatphi = phi.GetTensorMutableData<float>();
    float* floatc = c.GetTensorMutableData<float>();
    // std::cout<<u<<"\t"<<v<<"\t"<<x1<<"\t"<<x2<<std::endl;
    // std::cout<<*floatT<<"\t"<<*floatP<<"\t"<<*floatphi<<"\t"<<*floatc<<std::endl;

    return std::make_tuple(*floatT,*floatP,*floatphi,*floatc);
}

template <class BasicPsiThermo, class MixtureType>
std::tuple<scalar,scalar,scalar,scalar> Foam::ISATVLEheRhoThermo<BasicPsiThermo, MixtureType>::ANN_predict_TPX(float T, float p, float x1, float x2){
    //const char* model_dir = "/scratch.nike/srini237/OpenFOAM/temp/OpenFOAM-6/user/tutorials/python/UVFlash/ANN/OnnxModel.onnx";
    //std::cout<<"here inside predict"<<std::endl;
    //Ort::Session session(env, "/scratch.nike/srini237/OpenFOAM/temp/OpenFOAM-6/user/tutorials/python/UVFlash/ANN/OnnxModel.onnx", session_options);
    auto inputshape = session_TPX.GetInputTypeInfo(0).GetTensorTypeAndShapeInfo().GetShape();
    std::vector<float> inputValues = {T,p,x1,x2};
    int batchsize = 1;
    inputshape[0] = batchsize;
    auto memoryinfo = Ort::MemoryInfo::CreateCpu(OrtDeviceAllocator,OrtMemTypeCPU); 
    auto inputOnnxTensor = Ort::Value::CreateTensor<float>(memoryinfo,inputValues.data(),inputValues.size(),inputshape.data(),inputshape.size());
    const char* inputNames[] = {"input"};
    const char* outputNames[] = {"u","rho","phi","c"};
    //std::cout<<"here before run"<<std::endl;
    auto outputValues = session_TPX.Run(Ort::RunOptions(nullptr),inputNames,&inputOnnxTensor,1,outputNames,4);
    auto& u = outputValues[0];
    auto& rho = outputValues[1];
    auto& phi = outputValues[2];
    auto& c = outputValues[3];

    float* floatu = u.GetTensorMutableData<float>();
    float* floatrho = rho.GetTensorMutableData<float>();
    float* floatphi = phi.GetTensorMutableData<float>();
    float* floatc = c.GetTensorMutableData<float>();
    // std::cout<<u<<"\t"<<v<<"\t"<<x1<<"\t"<<x2<<std::endl;
    // std::cout<<*floatT<<"\t"<<*floatP<<"\t"<<*floatphi<<"\t"<<*floatc<<std::endl;

    return std::make_tuple(*floatu,*floatrho,*floatphi,*floatc);
}

// template <class BasicPsiThermo, class MixtureType>
// std::tuple<scalar,scalar,scalar,scalar> Foam::ISATVLEheRhoThermo<BasicPsiThermo, MixtureType>::ANN_predict_BC(float T, float p, float x1, float x2){
//     int NumInputs = 1;
//     TF_Output* Input = (TF_Output*)malloc(sizeof(TF_Output) * NumInputs);

//     TF_Output t0 = {TF_GraphOperationByName(Graph_BC, "serving_default_input"), 0};
    
//     Input[0] = t0;

//     // if(t0.oper == NULL)
//     //     printf("ERROR: Failed TF_GraphOperationByName serving_default_input\n");
//     // else
//     //     printf("TF_GraphOperationByName serving_default_input is OK\n");
    

//     int NumOutputs = 4;
//     TF_Output* Output = (TF_Output*)malloc(sizeof(TF_Output) * NumOutputs);

//     TF_Output t4 = {TF_GraphOperationByName(Graph_BC, "StatefulPartitionedCall"), 3}; 
//     TF_Output t5 = {TF_GraphOperationByName(Graph_BC, "StatefulPartitionedCall"), 2};
//     TF_Output t6 = {TF_GraphOperationByName(Graph_BC, "StatefulPartitionedCall"), 1};
//     TF_Output t7 = {TF_GraphOperationByName(Graph_BC, "StatefulPartitionedCall"), 0};

//     Output[0] = t4;
//     Output[1] = t5;
//     Output[2] = t6;
//     Output[3] = t7;

//     // if(t4.oper == NULL)
//     //     printf("ERROR: Failed TF_GraphOperationByName StatefulPartitionedCall\n");
//     // else        
//     //     printf("TF_GraphOperationByName StatefulPartitionedCall is OK\n");
    

//     TF_Tensor** InputValues = (TF_Tensor**)malloc(sizeof(TF_Tensor*)*NumInputs);
//     TF_Tensor** OutputValues = (TF_Tensor**)malloc(sizeof(TF_Tensor*)*NumOutputs);

//     int ndims = 2;
//     int64_t dims[] = {1,4};
//     float data[] = {T,p,x1,x2}; 
//     int ndata = sizeof(float)*1*4;

//     TF_Tensor* tensor = TF_NewTensor(TF_FLOAT, dims, ndims, data, ndata, &NoOpDeallocator, 0);

//     // if (tensor != NULL)
//     // {
//     //     printf("TF_NewTensor is OK\n");
//     // }
//     // else
//     //     printf("ERROR: Failed TF_NewTensor\n");
    
//     InputValues[0] = tensor;	
    
//     TF_SessionRun(Session_BC, NULL, Input, InputValues, NumInputs, Output, OutputValues, NumOutputs, NULL, 0, NULL , Status_BC);
    
//     // if(TF_GetCode(Status) == TF_OK)
//     // {
//     //     printf("Session is OK\n");
//     // }
//     // else
//     // {
//     //     printf("%s",TF_Message(Status));
//     // }

//     //std::cout<<"u = "<<data[0]<<"\t v = "<<data[1]<<"\t n2 mole fraction = "<<data[2]<<std::endl;
//     void* buff1 = TF_TensorData(OutputValues[0]);
//     void* buff2 = TF_TensorData(OutputValues[1]);
//     void* buff3 = TF_TensorData(OutputValues[2]);
//     void* buff4 = TF_TensorData(OutputValues[3]);
//     float* op1 = (float*)buff1;
//     float* op2 = (float*)buff2;
//     float* op3 = (float*)buff3;
//     float* op4 = (float*)buff4;
//     //std::cout<<"T \t"<<op1[0]<<"\t p \t"<<op2[0]<<"\t phi \t"<<op3[0]<<"\t c \t"<<op4[0]<<std::endl;

//     return std::make_tuple(op1[0],op2[0],op3[0],op4[0]);
// }