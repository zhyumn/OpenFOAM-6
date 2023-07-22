import numpy as np
import tensorflow as tf
from tensorflow.keras import Sequential,Input,Model
from tensorflow.keras.layers import Dense, Add, ReLU
from sklearn.model_selection import train_test_split


data1 = np.loadtxt("data_gen_IC_shock_N2.txt",delimiter='\t',dtype=float)
data2 = np.loadtxt("data_gen_IC_shock_new.txt",delimiter='\t',dtype=float)
data3 = np.loadtxt("data_gen_IC_droplet_new.txt",delimiter='\t',dtype=float)

data = np.vstack((data1,data2,data3))

# mu = np.mean(data,axis=0)
# sigma = np.std(data,axis=0)

# mu_input = np.array([mu[4],mu[5],mu[2],mu[3]])
# sigma_input = np.array([sigma[4],sigma[5],sigma[2],sigma[3]])

# mu_output = np.array([mu[0],mu[1],mu[6],mu[7]])
# sigma_output = np.array([sigma[0],sigma[1],sigma[6],sigma[7]])

# norm_vals = np.array((mu, sigma))
# np.savetxt('norm_vals.dat', norm_vals, delimiter='\t')

# data = (data - mu) / sigma

train,test = train_test_split(data, test_size=0.1, random_state = 13242)
train,val = train_test_split(train, test_size=0.2, random_state = 1423)

# print(train,"\n")
# print(val,'\n')
# print(test,'\n')
# print(train.shape)

print(train[:,4].shape)
Xtrain = np.vstack((train[:,4],train[:,5],train[:,2],train[:,3])).T
print(Xtrain)
print(Xtrain.shape)
utrain = train[:,0]
vtrain = 1.0/train[:,1]
vftrain = train[:,6]
ctrain = train[:,7]

Xval = np.vstack((val[:,4],val[:,5],val[:,2],val[:,3])).T
uval = val[:,0]
vval = 1.0/val[:,1]
vfval = val[:,6]
cval = val[:,7]

Xtest = np.vstack((test[:,4],test[:,5],test[:,2],test[:,3])).T
utest = test[:,0]
vtest = 1.0/test[:,1]
vftest = test[:,6]
ctest = test[:,7]

initializer = tf.keras.initializers.RandomNormal(stddev=0.1, seed=1342)

def identityBlock(inputLayer, units):
    x = Dense(units,activation='relu', kernel_initializer=initializer, bias_initializer=initializer)(inputLayer)
    # x = BatchNormalization()(x)
    # x = tf.keras.activations.relu(x)
    x = Dense(units, activation='relu', kernel_initializer=initializer, bias_initializer=initializer)(x)
    # x = BatchNormalization()(x)
    # x = tf.keras.activations.relu(x)
    x = Dense(inputLayer.shape[1], kernel_initializer=initializer, bias_initializer=initializer)(x)
    # x = BatchNormalization()(x)
    add = Add()([x, inputLayer])
    # add = ReLU()(add)    
    return add

def denseBlock(inputLayer, units):
    x = Dense(units, activation='relu', kernel_initializer=initializer, bias_initializer=initializer)(inputLayer)
    # x = BatchNormalization()(x)
    # x = tf.keras.activations.relu(x)
    x = Dense(units, activation='relu', kernel_initializer=initializer, bias_initializer=initializer)(x)
    # x = BatchNormalization()(x)
    # x = tf.keras.activations.relu(x)
    x = Dense(units, kernel_initializer=initializer, bias_initializer=initializer)(x)
    # x = BatchNormalization()(x)
    y = Dense(units, kernel_initializer=initializer, bias_initializer=initializer)(inputLayer)
    # y = BatchNormalization()(y)
    add = Add()([x, y])
    # add = ReLU()(add)
    return add

def resBlock(inputLayer, units):
    l1 = denseBlock(inputLayer, units)
    l2 = identityBlock(l1, units)
    l3 = identityBlock(l2, units)
    return l3

def custom_norm(inputLayer,mean,sigma):
    x = tf.math.subtract(inputLayer,mean)
    x = tf.math.divide(x,sigma)
    return x

def custom_denorm(inputLayer,mean,sigma):
    x = tf.math.multiply(inputLayer,sigma)
    x = tf.math.add(x,mean)
    return x


# MODEL
input_shape = (4,)
inputs = Input(shape=(4,),name ='input')
normalizer = tf.keras.layers.Normalization(axis = -1)
normalizer.adapt(Xtrain)
l1 = normalizer(inputs)

units = 32
x = resBlock(inputs, units)
for i in range(3):
    x = resBlock(x, units)

# x = BatchNormalization()(x)
output1 = Dense(1,activation='relu',name='u',kernel_initializer=initializer)(x)
output2 = Dense(1,activation='relu',name='v',kernel_initializer=initializer)(x)
output3 = Dense(1,activation='relu',name='phi',kernel_initializer=initializer)(x)
output4 = Dense(1,activation='relu',name='c',kernel_initializer=initializer)(x)

model = Model(inputs=inputs, outputs=[output1,output2,output3,output4])
model.compile(optimizer=tf.keras.optimizers.Adam(), loss = {'u':'mae','v':'mae','phi':'mae','c':'mae'} ) #loss=tf.keras.losses.BinaryCrossentropy(from_logits=False), metrics=[tf.keras.metrics.MeanAbsoluteError(),])

model.summary()

model.fit(Xtrain,[utrain,vtrain,vftrain,ctrain],epochs = 200, batch_size = 64, validation_data=(Xval,[uval,vval,vfval,cval]))

# #print(model.predict(Xtest) - Ytest)
print("model testing")
print(model.evaluate(Xtest,[utest,vtest,vftest,ctest]))
model.save("TPX_resnet")
