import numpy as np
# import tensorflow as tf
import onnx
# import os
# import tf2onnx
import onnxruntime

def model_pred():
# reconstruct model
# print(tf.__version__)
# recons = tf.keras.models.load_model("test_model_64_64_32_1e-3")
    inputs = np.array([[-5650.61,0.0011777777777777778,0.7156593087182554,0.2843406912817447],[-5450.61,0.0013777777777777778,0.7,0.3]],dtype=np.float32)
    # print(inputs.shape)
    # print(recons.predict(inputs))

    # conv_model,_ = tf2onnx.convert.from_keras(recons)

    # onnx.save(conv_model,'test_model_64_64_32_1e-3.onnx')

    session = onnxruntime.InferenceSession('test_model_64_64_32_1e-3.onnx')
    print(session.get_inputs()[0].name)
    print(session.get_inputs()[0].shape)
    print(session.get_outputs()[0].name)
    print(session.get_outputs()[0].shape)
    print(session.get_outputs()[1].name)
    print(session.get_outputs()[2].name)
    print(session.get_outputs()[3].name)
    input_name = session.get_inputs()[0].name
    output_name = [session.get_outputs()[0].name, session.get_outputs()[1].name, session.get_outputs()[2].name, session.get_outputs()[3].name]
    T,p,phi,c = session.run(output_name,{input_name: inputs})
    print(T,p,phi,c)
