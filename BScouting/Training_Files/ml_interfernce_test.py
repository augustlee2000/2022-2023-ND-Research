import onnxruntime

onnx_model_path = "very_bad_btagger.onnx"
session = onnxruntime.InferenceSession(onnx_model_path)


# Define values for N and n_pf
N = 1
n_pf = 10

# Create a tensor with the specified shape
tensor_shape = (N, 2, n_pf)
pf_points = np.zeros(tensor_shape, dtype=np.float32)

tensor_shape = (N, 16, n_pf)
pf_features = np.zeros(tensor_shape, dtype=np.float32)

tensor_shape = (N, 1, n_pf)
pf_mask = np.ones(tensor_shape, dtype=np.float32)


inputs = {
    'pf_points': pf_points,
    'pf_features': pf_features,
    'pf_mask': pf_mask, 
}
# Run inference
output = session.run(["softmax"],inputs)

# Assuming you only have one output named "softmax"
softmax_output = output[0]

# Print or inspect the softmax output
print("Softmax output shape:", softmax_output.shape)
print("Softmax output values:")
print(softmax_output)

