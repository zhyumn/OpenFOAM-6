import numpy as np
import tensorflow as tf

# reconstruct model
recons = tf.keras.models.load_model("test_model_64_64_32_1e-3")
input = [-5650.61,0.0011777777777777778,0.7156593087182554,0.2843406912817447]
print(recons.predict(input))