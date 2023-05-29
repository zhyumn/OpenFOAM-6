import numpy as np
import tensorflow as tf
import onnx
#import os
#import tf2onnx
import onnxruntime
# reconstruct model
# print(tf.__version__)
print(onnxruntime.__version__)

recons = tf.keras.models.load_model("ANN_model")
   
inputs = np.array([[-5333.63574219,0.00106602546293,0.574138760567,0.425861209631]],dtype=np.float32)
print(inputs.shape)
print(recons.predict(inputs))

#conv_model,_ = tf2onnx.convert.from_keras(recons)

#onnx.save(conv_model,'OnnxModel2.onnx')

session = onnxruntime.InferenceSession('OnnxModel2.onnx')
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
