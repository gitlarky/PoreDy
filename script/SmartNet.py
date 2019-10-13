#==================================================================================================
# Title      : SmartNet
# Description: Smart Predictiing and Designing of Pore-Network
# Author     : Zhenyu Xu; westlark@outlook.com
# Start Time : 2019.09.01
# License    : Apache 2.0
#==================================================================================================

#============================ Module Importation ==================================================
# from __future__ import absolute_import, division, print_function, unicode_literals

# import tensorflow as tf
# from tensorflow import keras

# print(tf.__version__)

# import numpy as np

# import os
#==================================================================================================
#============================ Create Pore-Network Samples =========================================
def WriteAT(File='', Matrix):
	return True
def CreatePoreNetworkSamples(Nx=0, Ny=0, Folder=''):
	# Determine What and How Many is each type of Throat ------------------------------------------
	TP         =[1e-6, 2e-6, 3e-6, 4e-6, 5e-6, 6e-6, 7e-6, 8e-6, 9e-6] # Parameter is Diameter here
	TypeN      =len(TP)
	print('Throat Parameter   :', TP)
	ShoTNum    =Ny*(Nx+1)+Ny*Nx
	LonTNum    =Ny*Nx
	TotTNum    =ShoTNum+LonTNum
	AvgShoTNum =ShoTNum//TypeN
	RemShoTNum =ShoTNum%TypeN
	STN        =[AvgShoTNum    for i in range(TypeN)]
	STN[4]    +=RemShoTNum
	print('Short Throat Number:', STN)
	AvgLonTNum =LonTNum//TypeN
	RemLonTNum =LonTNum%TypeN
	LTN        =[AvgLonTNum    for i in range(TypeN)]
	LTN[4]    +=RemLonTNum
	print('Long  Throat Number:', LTN)
	TN         =[STN[i]+LTN[i] for i in range(TypeN)]
	print('Total Throat Number:', TN)

	SampleIndex=0

	# for offset in range(Nx)

	# with open(CaseName+'.at', 'w') as at:
	# 	for j in range(Ny*2):
	# 		for i in range(Nx*2):
	# 			if(i%2==0 and j%2==0):
	# 				Diameter=0
	# 				PolyN   =0
	# 			elif():
	# 				Diameter=
	# 				PolyN   =5
	# 			at.write('% 9d\t% 9.6e\t' % (Diameter, PolyN))
	# 		at.write('\n')






# #============================ Basic Functions =====================================================
# def create_model():
# 	model = tf.keras.models.Sequential([
# 		keras.layers.Dense(512, activation=tf.keras.activations.relu, input_shape=(784,)),
# 		keras.layers.Dropout(0.2),
# 		keras.layers.Dense(10, activation=tf.keras.activations.softmax)])

# 	model.compile(optimizer=tf.keras.optimizers.Adam(),
#                 loss=tf.keras.losses.sparse_categorical_crossentropy,
#                 metrics=['accuracy'])

# 	return model

# checkpoint_path = "training_1/cp.ckpt"
# checkpoint_dir = os.path.dirname(checkpoint_path)

# # Create checkpoint callback
# cp_callback = tf.keras.callbacks.ModelCheckpoint(checkpoint_path,
#                                                  save_weights_only=True,
#                                                  verbose=1)

# model = create_model()

# model.fit(train_images, train_labels,  epochs = 10,
#           validation_data = (test_images,test_labels),
#           callbacks = [cp_callback])  # pass callback to training

# # This may generate warnings related to saving the state of the optimizer.
# # These warnings (and similar warnings throughout this notebook)
# # are in place to discourage outdated usage, and can be ignored.


# model = create_model()

# loss, acc = model.evaluate(test_images, test_labels)
# print("Untrained model, accuracy: {:5.2f}%".format(100*acc))

# model.load_weights(checkpoint_path)
# loss,acc = model.evaluate(test_images, test_labels)
# print("Restored model, accuracy: {:5.2f}%".format(100*acc))


#============================ Main Program ========================================================
CreatePoreNetworkSamples(Nx=20, Ny=20, Folder='/home/xu/work/PoreNetwork2020Samples')
CreatePoreNetworkSamples(Nx=10, Ny=10, Folder='/home/xu/work/PoreNetwork1010Samples')
CreatePoreNetworkSamples(Nx=40, Ny=40, Folder='/home/xu/work/PoreNetwork4040Samples')