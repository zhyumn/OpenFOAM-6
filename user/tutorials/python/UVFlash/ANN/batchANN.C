#include <stdio.h>
#include <stdlib.h>
#include <tensorflow/c/c_api.h>

int main() {
    // Set up the input tensor shapes
    const int batch_size = 5;
    const int input_size = 4;
    const int output_size = 4;
    const int64_t input_dims[2] = {batch_size, input_size};
    const int64_t output_dims[2] = {batch_size, output_size};

    // Create the input and output tensors
    TF_Tensor* input_tensor = TF_AllocateTensor(TF_FLOAT, input_dims, 2, input_size * batch_size * sizeof(float));
    TF_Tensor* output_tensor = TF_AllocateTensor(TF_FLOAT, output_dims, 2, output_size * batch_size * sizeof(float));

    // Load the saved model
    const char* saved_model_dir = "./test_model_64_64_32_ind_pbigger_1e-3"; // Replace with your own saved model directory
    TF_Status* status = TF_NewStatus();
    TF_SessionOptions* session_opts = TF_NewSessionOptions();
    TF_Session* session = TF_LoadSessionFromSavedModel(session_opts, NULL, saved_model_dir, NULL, 0, NULL, status);
    if (TF_GetCode(status) != TF_OK) {
        fprintf(stderr, "Failed to load saved model: %s\n", TF_Message(status));
        return 1;
    }

    // Create the input and output placeholders
    TF_Output input_placeholder = {TF_GraphOperationByName(TF_SessionGraph(session), "input"), 0};
    TF_Output output_placeholder = {TF_GraphOperationByName(TF_SessionGraph(session), "output"), 0};

    // Run the model
    TF_SessionRun(session, NULL,
                  &input_placeholder, (TF_Tensor**) &input_tensor, 1,
                  &output_placeholder, (TF_Tensor**) &output_tensor, 1,
                  NULL, 0, NULL, status);

    // Print the results
    float* output_data = (float*) TF_TensorData(output_tensor);
    for (int i = 0; i < batch_size * output_size; i++) {
        printf("%f ", output_data[i]);
    }
    printf("\n");

    // Clean up
    TF_CloseSession(session, status);
    TF_DeleteSession(session, status);
    TF_DeleteSessionOptions(session_opts);
    TF_DeleteTensor(input_tensor);
    TF_DeleteTensor(output_tensor);
    TF_DeleteStatus(status);

    return 0;
}
