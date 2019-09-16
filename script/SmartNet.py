#==================================================================================================
# Title      : SmartNet
# Description: Smart Predictiing and Designing of Pore-Network
# Author     : Zhenyu Xu; westlark@outlook.com
# Start Time : 2019.09.01
# License    : Apache 2.0
#==================================================================================================

#============================ Module Importation ==================================================
from __future__ import absolute_import, division, print_function, unicode_literals

import tensorflow as tf
from tensorflow import keras

print(tf.__version__)

import numpy as np

import os
#==================================================================================================

#============================ Basic Functions =====================================================
def create_model():
	model = tf.keras.models.Sequential([
		keras.layers.Dense(512, activation=tf.keras.activations.relu, input_shape=(784,)),
		keras.layers.Dropout(0.2),
		keras.layers.Dense(10, activation=tf.keras.activations.softmax)])

	model.compile(optimizer=tf.keras.optimizers.Adam(),
                loss=tf.keras.losses.sparse_categorical_crossentropy,
                metrics=['accuracy'])

	return model

checkpoint_path = "training_1/cp.ckpt"
checkpoint_dir = os.path.dirname(checkpoint_path)

# Create checkpoint callback
cp_callback = tf.keras.callbacks.ModelCheckpoint(checkpoint_path,
                                                 save_weights_only=True,
                                                 verbose=1)

model = create_model()

model.fit(train_images, train_labels,  epochs = 10,
          validation_data = (test_images,test_labels),
          callbacks = [cp_callback])  # pass callback to training

# This may generate warnings related to saving the state of the optimizer.
# These warnings (and similar warnings throughout this notebook)
# are in place to discourage outdated usage, and can be ignored.


model = create_model()

loss, acc = model.evaluate(test_images, test_labels)
print("Untrained model, accuracy: {:5.2f}%".format(100*acc))

model.load_weights(checkpoint_path)
loss,acc = model.evaluate(test_images, test_labels)
print("Restored model, accuracy: {:5.2f}%".format(100*acc))