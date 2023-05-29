import numpy as np
import tensorflow as tf
from tensorflow.keras import Sequential,Input,Model
from tensorflow.keras.layers import Dense
from sklearn.model_selection import train_test_split


data1 = np.loadtxt("data_gen_IC_shock_N2.txt",delimiter='\t',dtype=float)
data2 = np.loadtxt("data_gen_IC_shock_new.txt",delimiter='\t',dtype=float)
data3 = np.loadtxt("data_gen_IC_droplet_new.txt",delimiter='\t',dtype=float)

data = np.vstack((data1,data2,data3))

train,test = train_test_split(data, test_size=0.1, random_state = 1)
train,val = train_test_split(train, test_size=0.2, random_state = 1)

# print(train,"\n")
# print(val,'\n')
# print(test,'\n')
# print(train.shape)

print(train[:,4].shape)
Xtrain = np.vstack((train[:,4],train[:,5],train[:,2],train[:,3])).T
print(Xtrain)
print(Xtrain.shape)
utrain = train[:,0]
rhotrain = train[:,1]
vftrain = train[:,6]
ctrain = train[:,7]

Xval = np.vstack((val[:,4],val[:,5],val[:,2],val[:,3])).T
uval = val[:,0]
rhoval = val[:,1]
vfval = val[:,6]
cval = val[:,7]

Xtest = np.vstack((test[:,4],test[:,5],test[:,2],test[:,3])).T
utest = test[:,0]
rhotest = test[:,1]
vftest = test[:,6]
ctest = test[:,7]

# print(Xtrain,"\n")
# print(Xval,'\n')
# print(Xtest,'\n')

num1 = 64
num2 = 64
num3 = 32

input_shape = (4,) 
initializer = tf.keras.initializers.RandomNormal(stddev=0.1)
normalizer = tf.keras.layers.Normalization(axis = -1)
normalizer.adapt(Xtrain)

#model definition
inputs = Input(shape=(4,),name ='input')
l1 = normalizer(inputs)

# l2 = Dense(num1,activation='relu',name='l2',kernel_initializer=initializer,bias_initializer=tf.keras.initializers.Zeros())(l1)
# l3 = Dense(num1,activation='relu',name='l3',kernel_initializer=initializer,bias_initializer=tf.keras.initializers.Zeros())(l2)
# x = Dense(num1,activation='relu',name='x',kernel_initializer=initializer,bias_initializer=tf.keras.initializers.Zeros())(l3)


Tx1 = Dense(num1,activation='relu',name='T1',kernel_initializer=initializer,bias_initializer=tf.keras.initializers.Zeros())(l1)
Tx2 = Dense(num2,activation='relu',name='T2',kernel_initializer=initializer,bias_initializer=tf.keras.initializers.Zeros())(Tx1)
Tx3 = Dense(num2,activation='relu',name='T3',kernel_initializer=initializer,bias_initializer=tf.keras.initializers.Zeros())(Tx2)
Tx4 = Dense(num2,activation='relu',name='T4',kernel_initializer=initializer,bias_initializer=tf.keras.initializers.Zeros())(Tx3)
Tx5 = Dense(num2,activation='relu',name='T5',kernel_initializer=initializer,bias_initializer=tf.keras.initializers.Zeros())(Tx4)
Tx6 = Dense(num3,activation='relu',name='T6',kernel_initializer=initializer,bias_initializer=tf.keras.initializers.Zeros())(Tx5)

px1 = Dense(num1,activation='relu',name='p1',kernel_initializer=initializer,bias_initializer=tf.keras.initializers.Zeros())(l1)
px2 = Dense(num2,activation='relu',name='p2',kernel_initializer=initializer,bias_initializer=tf.keras.initializers.Zeros())(px1)
px3 = Dense(num2,activation='relu',name='p3',kernel_initializer=initializer,bias_initializer=tf.keras.initializers.Zeros())(px2)
px4 = Dense(num2,activation='relu',name='p4',kernel_initializer=initializer,bias_initializer=tf.keras.initializers.Zeros())(px3)
px5 = Dense(num2,activation='relu',name='p5',kernel_initializer=initializer,bias_initializer=tf.keras.initializers.Zeros())(px4)
px6 = Dense(num3,activation='relu',name='p6',kernel_initializer=initializer,bias_initializer=tf.keras.initializers.Zeros())(px5)

vfx1 = Dense(num1,activation='relu',name='vf1',kernel_initializer=initializer,bias_initializer=tf.keras.initializers.Zeros())(l1)
vfx2 = Dense(num2,activation='relu',name='vf2',kernel_initializer=initializer,bias_initializer=tf.keras.initializers.Zeros())(vfx1)
vfx3 = Dense(num3,activation='relu',name='vf3',kernel_initializer=initializer,bias_initializer=tf.keras.initializers.Zeros())(vfx2)
vfx4 = Dense(num3,activation='relu',name='vf4',kernel_initializer=initializer,bias_initializer=tf.keras.initializers.Zeros())(vfx3)
vfx5 = Dense(num3,activation='relu',name='vf5',kernel_initializer=initializer,bias_initializer=tf.keras.initializers.Zeros())(vfx3)

cx1 = Dense(num1,activation='relu',name='c1',kernel_initializer=initializer,bias_initializer=tf.keras.initializers.Zeros())(l1)
cx2 = Dense(num2,activation='relu',name='c2',kernel_initializer=initializer,bias_initializer=tf.keras.initializers.Zeros())(cx1)
cx3 = Dense(num3,activation='relu',name='c3',kernel_initializer=initializer,bias_initializer=tf.keras.initializers.Zeros())(cx2)
cx4 = Dense(num3,activation='relu',name='c4',kernel_initializer=initializer,bias_initializer=tf.keras.initializers.Zeros())(cx3)
cx5 = Dense(num3,activation='relu',name='c5',kernel_initializer=initializer,bias_initializer=tf.keras.initializers.Zeros())(cx4)


# x4 = Dense(num1,activation='relu',name='4',kernel_initializer=initializer)(x)
# x5 = Dense(num2,activation='relu',name='5',kernel_initializer=initializer)(x4)
# x6 = Dense(num3,activation='relu',name='6',kernel_initializer=initializer)(x5)

# x7 = Dense(num1,activation='relu',name='7',kernel_initializer=initializer)(x)
# x8 = Dense(num2,activation='relu',name='8',kernel_initializer=initializer)(x7)
# x9 = Dense(num3,activation='relu',name='9',kernel_initializer=initializer)(x8)

# x10 = Dense(num1,activation='relu',name='10',kernel_initializer=initializer)(x)
# x11 = Dense(num2,activation='relu',name='11',kernel_initializer=initializer)(x10)
# x12 = Dense(num3,activation='relu',name='12',kernel_initializer=initializer)(x11)

output1 = Dense(1,activation='relu',name='u',kernel_initializer=initializer)(Tx6)
output2 = Dense(1,activation='relu',name='rho',kernel_initializer=initializer)(px6)
output3 = Dense(1,activation='relu',name='phi',kernel_initializer=initializer)(vfx5)
output4 = Dense(1,activation='relu',name='c',kernel_initializer=initializer)(cx5)

model = Model(inputs=inputs, outputs=[output1,output2,output3,output4])
model.compile(optimizer=tf.keras.optimizers.Adam(learning_rate=1e-4), loss = {'u':'mae','rho':'mae','phi':'mae','c':'mae'} ) #loss=tf.keras.losses.BinaryCrossentropy(from_logits=False), metrics=[tf.keras.metrics.MeanAbsoluteError(),])

model.summary()

# Xtest = np.loadtxt("Xtest.csv",delimiter=',',dtype=float)
# Ytest = np.loadtxt("Ytest.csv",delimiter=',',dtype=float)

# #Ytest = Bcsv.copy()
# #Ytest = np.array(Ytest)
# Xliq_t= Ytest[:,0:2]
# Xgas_t = Ytest[:,2:4]
# phi_t = Ytest[:,4]

model.fit(Xtrain,[utrain,rhotrain,vftrain,ctrain],epochs = 200, validation_data=(Xval,[uval,rhoval,vfval,cval]))

# #print(model.predict(Xtest) - Ytest)
print("model testing")
print(model.evaluate(Xtest,[utest,rhotest,vftest,ctest]))
model.save("TPX")

# #test reconstruct model
# #recons = tf.keras.models.load_model("my_model")

# #print(recons.predict(Xtest))

# #converting keras model to concrete function
# #full_model = tf.function(lambda x: model(x))
# #full_model = full_model.get_concrete_function(tf.TensorSpec(model.inputs[0].shape,model.inputs[0].dtype))

# #get frozen concretefunction
# #frozen_func = convert_variables_to_constants(full_model)
# #frozen_func.graph.as_graph_def()
# #save to disk
# #tf.train.write_graph(graph_or_graph_def=frozen_func.graph,logdir=frozen_out_path,name=f"{frozen_graph_filename}.pb",as_text=False)

# #print(loss)
# #tf.keras.metrics.mean_absolute_percentage_error(Ytest, Ypred)"""
