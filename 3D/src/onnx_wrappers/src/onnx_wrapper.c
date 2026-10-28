#include <stdio.h>
#include <stdlib.h>
#include <onnxruntime_c_api.h>

static const OrtApi* g_ort = NULL;
OrtEnv* env;
OrtSessionOptions* session_options;
OrtSession* session;
OrtMemoryInfo* memory_info;

// Helper to check for ORT errors
void check_status(OrtStatus* status) {
    if (status != NULL) {
        const char* msg = g_ort->GetErrorMessage(status);
        fprintf(stderr, "ONNX Runtime Error: %s\n", msg);
        g_ort->ReleaseStatus(status);
        exit(1);
    }
}

// Function callable from Fortran
void init_onnx_cnn(
    const char* model_path
) {
    // 1. Initialize API table
    if (g_ort == NULL) {
        g_ort = OrtGetApiBase()->GetApi(ORT_API_VERSION);
    }

    // 2. Create Environment and Session Options
    check_status(g_ort->CreateEnv(ORT_LOGGING_LEVEL_WARNING, "FortranONNX", &env));

    check_status(g_ort->CreateSessionOptions(&session_options));

    // 3. Disable thread pool oversubscription across MPI ranks
    check_status(g_ort->SetIntraOpNumThreads(session_options, 16));
    check_status(g_ort->SetInterOpNumThreads(session_options, 1));

    // 4. Create Session (load model)
    check_status(g_ort->CreateSession(env, model_path, session_options, &session));

    // 5. Memory Allocator Info
    check_status(g_ort->CreateCpuMemoryInfo(OrtArenaAllocator, OrtMemTypeDefault, &memory_info));
}

void run_onnx_cnn(
    const float* input_data,  // Pointer to flat input array
    int64_t* input_dims,     
    size_t num_dims,
    float* output_data,       // Buffer to store output results
    size_t output_size
) {
    // 1. Wrap input array into ORT Tensor
    size_t input_tensor_size = 1;
    for (size_t i = 0; i < num_dims; i++) {
        input_tensor_size *= input_dims[i];
    }

    OrtValue* input_tensor = NULL;
    check_status(g_ort->CreateTensorWithDataAsOrtValue(
        memory_info,
        (void*)input_data,
        input_tensor_size * sizeof(float),
        input_dims,
        num_dims,
        ONNX_TENSOR_ELEMENT_DATA_TYPE_FLOAT,
        &input_tensor
    ));

    // Define Input and Output Tensor Names (Match your exported ONNX model input/output names)
//    const char* input_node_names[] = {"input"};   // Change to your model's input layer name
//    const char* output_node_names[] = {"output"}; // Change to your model's output layer name

    const char* input_node_names[] = {"onnx::Reshape_0"};   // Change to your model's input layer name
    const char* output_node_names[] = {"257"}; // Change to your model's output layer name

    // 2. Run Inference
    OrtValue* output_tensor = NULL;
    check_status(g_ort->Run(
        session,
        NULL,
        input_node_names,
        (const OrtValue* const*)&input_tensor,
        1,
        output_node_names,
        1,
        &output_tensor
    ));

    // 3. Copy output data back to Fortran array buffer
    float* float_out;
    check_status(g_ort->GetTensorMutableData(output_tensor, (void**)&float_out));
    for (size_t i = 0; i < output_size; i++) {
        output_data[i] = float_out[i];
    }
    // Cleanup
    g_ort->ReleaseValue(output_tensor);
    g_ort->ReleaseValue(input_tensor);
}

void close_onnx_cnn()
{
    g_ort->ReleaseMemoryInfo(memory_info);
    g_ort->ReleaseSession(session);
    g_ort->ReleaseSessionOptions(session_options);
    g_ort->ReleaseEnv(env);
}

