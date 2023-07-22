import numpy as np
import tensorflow as tf
import onnx
import os
import tf2onnx
import onnxruntime
# reconstruct model
# print(tf.__version__)

print(onnxruntime.__version__)

recons = tf.keras.models.load_model("TPX")

u_in = 25600.0
v_in = 0.0042
x1_in = 0.9651397724218261
x2_in = 0.034860227578173995
   
inputs = np.array([[u_in,v_in,x1_in,x2_in]],dtype=np.float32)
print(inputs.shape)
# print(recons.predict(inputs))

conv_model,_ = tf2onnx.convert.from_keras(recons)

onnx.save(conv_model,'TPX.onnx')

session = onnxruntime.InferenceSession('OnnxModel2.onnx')
# print(session.get_inputs()[0].name)
# print(session.get_inputs()[0].shape)
# print(session.get_outputs()[0].name)
# print(session.get_outputs()[0].shape)
# print(session.get_outputs()[1].name)
# print(session.get_outputs()[2].name)
# print(session.get_outputs()[3].name)
input_name = session.get_inputs()[0].name
output_name = [session.get_outputs()[0].name, session.get_outputs()[1].name, session.get_outputs()[2].name, session.get_outputs()[3].name]
T,p,phi,c = session.run(output_name,{input_name: inputs})
print(T,p,phi,c)





# sessionTPX = onnxruntime.InferenceSession('TPX.onnx')
# inputs_TPX = np.array([[T[0,0],p[0,0],x1_in,x2_in]],dtype=np.float32)
# input_nameTPX = sessionTPX.get_inputs()[0].name
# output_nameTPX = [sessionTPX.get_outputs()[0].name, sessionTPX.get_outputs()[1].name, sessionTPX.get_outputs()[2].name, sessionTPX.get_outputs()[3].name]
# u,rho,phi2,c2 = sessionTPX.run(output_nameTPX,{input_nameTPX: inputs_TPX})
# # outputs = session.run(output_name,{input_name: inputs_TPX})
# # u = outputs[0,0]
# # rho = outputs[0,1]
# # phi2 = outputs[0,2]
# # c2 = outputs[0,3]
# v = 1/rho

# print(inputs)
# print(u,v,phi2,c2)
# inputs = np.array([[u[0,0],v[0,0],x1_in,x2_in]],dtype=np.float32)
# T,p,phi,c = session.run(output_name,{input_name: inputs})
# print(T,p,phi,c)




