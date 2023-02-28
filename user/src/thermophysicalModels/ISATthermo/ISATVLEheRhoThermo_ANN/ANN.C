//ANN stuff
//namespace std
template <class BasicPsiThermo, class MixtureType>
void Foam::ISATVLEheRhoThermo<BasicPsiThermo, MixtureType>::ANN_load(char* saved_model_dir,char* tags,int ntags){
    //ANN information and initialization
    // char* saved_model_dir = "./ANN/test_model_64_64_32_ind_pbigger_1e-3/"; // Path of the model
    // char* tags = "serve"; // default model serving tag; can change in future
    // int ntags = 1;
    Session = TF_LoadSessionFromSavedModel(SessionOpts, RunOpts, saved_model_dir, &tags, ntags, Graph, NULL, Status);
    // if(TF_GetCode(Status) == TF_OK)
    // {
    //     printf("TF_LoadSessionFromSavedModel OK\n");
    // }
    // else
    // {
    //     printf("%s",TF_Message(Status));
    // }
}

template <class BasicPsiThermo, class MixtureType>
std::tuple<scalar,scalar,scalar,scalar> Foam::ISATVLEheRhoThermo<BasicPsiThermo, MixtureType>::ANN_predict(float u, float v, float x1, float x2){
    int NumInputs = 1;
    TF_Output* Input = (TF_Output*)malloc(sizeof(TF_Output) * NumInputs);

    TF_Output t0 = {TF_GraphOperationByName(Graph, "serving_default_input"), 0};
    
    Input[0] = t0;

    // if(t0.oper == NULL)
    //     printf("ERROR: Failed TF_GraphOperationByName serving_default_input\n");
    // else
    //     printf("TF_GraphOperationByName serving_default_input is OK\n");
    

    int NumOutputs = 4;
    TF_Output* Output = (TF_Output*)malloc(sizeof(TF_Output) * NumOutputs);

    TF_Output t4 = {TF_GraphOperationByName(Graph, "StatefulPartitionedCall"), 1}; 
    TF_Output t5 = {TF_GraphOperationByName(Graph, "StatefulPartitionedCall"), 0};
    TF_Output t6 = {TF_GraphOperationByName(Graph, "StatefulPartitionedCall"), 3};
    TF_Output t7 = {TF_GraphOperationByName(Graph, "StatefulPartitionedCall"), 2};

    Output[0] = t4;
    Output[1] = t5;
    Output[2] = t6;
    Output[3] = t7;

    // if(t4.oper == NULL)
    //     printf("ERROR: Failed TF_GraphOperationByName StatefulPartitionedCall\n");
    // else        
    //     printf("TF_GraphOperationByName StatefulPartitionedCall is OK\n");
    

    TF_Tensor** InputValues = (TF_Tensor**)malloc(sizeof(TF_Tensor*)*NumInputs);
    TF_Tensor** OutputValues = (TF_Tensor**)malloc(sizeof(TF_Tensor*)*NumOutputs);

    int ndims = 2;
    int64_t dims[] = {1,4};
    float data[] = {u,v,x1,x2}; 
    int ndata = sizeof(float)*1*4;

    TF_Tensor* tensor = TF_NewTensor(TF_FLOAT, dims, ndims, data, ndata, &NoOpDeallocator, 0);

    // if (tensor != NULL)
    // {
    //     printf("TF_NewTensor is OK\n");
    // }
    // else
    //     printf("ERROR: Failed TF_NewTensor\n");
    
    InputValues[0] = tensor;	
    
    TF_SessionRun(Session, NULL, Input, InputValues, NumInputs, Output, OutputValues, NumOutputs, NULL, 0, NULL , Status);
    
    // if(TF_GetCode(Status) == TF_OK)
    // {
    //     printf("Session is OK\n");
    // }
    // else
    // {
    //     printf("%s",TF_Message(Status));
    // }

    //std::cout<<"u = "<<data[0]<<"\t v = "<<data[1]<<"\t n2 mole fraction = "<<data[2]<<std::endl;
    void* buff1 = TF_TensorData(OutputValues[0]);
    void* buff2 = TF_TensorData(OutputValues[1]);
    void* buff3 = TF_TensorData(OutputValues[2]);
    void* buff4 = TF_TensorData(OutputValues[3]);
    float* op1 = (float*)buff1;
    float* op2 = (float*)buff2;
    float* op3 = (float*)buff3;
    float* op4 = (float*)buff4;
    //std::cout<<"T \t"<<op1[0]<<"\t p \t"<<op2[0]<<"\t phi \t"<<op3[0]<<"\t c \t"<<op4[0]<<std::endl;

    return std::make_tuple(op1[0],op2[0],op3[0],op4[0]);
}
