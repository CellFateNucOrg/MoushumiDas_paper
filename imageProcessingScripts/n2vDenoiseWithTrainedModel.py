#!/usr/bin/env python3

from n2v.models import N2VConfig, N2V
import numpy as np
#from nis2pyr.convertor import convert_nd2_to_pyramidal_ome_tiff
from csbdeep.utils import plot_history
from csbdeep.io import save_tiff_imagej_compatible
from n2v.utils.n2v_utils import manipulate_val_data
from n2v.internals.N2V_DataGenerator import N2V_DataGenerator
from matplotlib import pyplot as plt
import urllib
#import nd2
from nd2reader import ND2Reader
import os
import re
import zipfile
import argparse
from tifffile import imread, imwrite
os.environ['TF_XLA_FLAGS'] = '--tf_xla_enable_xla_devices'

import ssl


def n2vDenoiseWithTrainedModel(model_dir, model_name, image_dir, image_name, channel_name, plotInputPredict=False):
	'''
	Take in name and location of a n2v trained model and a single channel image
	and denoise
	'''
	print("\n\nDenoising image...\n\n")
	# A previously trained model is loaded by creating a new N2V-object without providing a 'config'.
	model = N2V(config=None, name=model_name, basedir=model_dir)

	# We load the data we want to process.
	split_image_dir=os.path.join(image_dir, 'n2v_denoise/TMP' + channel_name)
	split_image_name=re.sub("\.nd2", "_" + channel_name + ".tif", image_name)
	img = imread(os.path.join(split_image_dir,split_image_name))

	# Here we process the data.
	# The 'n_tiles' parameter can be used if images are too big for the GPU memory.
	# If we do not provide the 'n_tiles' parameter the system will automatically try to find an appropriate tiling.
	pred = model.predict(img, axes='ZYX', n_tiles=(1,4,4))
	
	# output flattened pngs of input and predicted images
	if(plotInputPredict):
		inputpred_path=os.path.join(image_dir,'n2v_denoise/inputPred/')
		if not os.path.isdir(inputpred_path):
			os.makedirs(inputpred_path, exist_ok=True)
		
		# Let's look at the results.
		plt.figure(figsize=(30,30))
		
		# We show the noisy input...
		plt.subplot(1,2,1)
		plt.imshow(np.max(img[...],axis=0),
				   cmap='gray',
				   vmin=np.percentile(img,0.1),
				   vmax=np.percentile(img,99.9)
				  )
		plt.title('Input');
		# and the result.
		plt.subplot(1,2,2)
		plt.imshow(np.max(pred[...],axis=0), 
				   cmap='gray',
				   vmin=np.percentile(pred,0.1),
				   vmax=np.percentile(pred,99.9)
				  )
		plt.title('Prediction');
		plt.savefig(os.path.join(inputpred_path,"InputvPredict_"+re.sub("\.nd2","_"+channel_name+"_n2v.tif",image_name)))
	
	# output denoised images
	denoised_path=os.path.join(image_dir,'n2v_denoise/denoised/')
	if not os.path.isdir(denoised_path):
			os.makedirs(denoised_path, exist_ok=True)
	imwrite(os.path.join(denoised_path, re.sub("\.nd2", "_"+channel_name+"_n2v.tif", image_name)), pred)
	return None



def splitChannels(image_dir, image_name, channel_list):
	'''
	Split nd2 image into separate color channels and save each to a temporary 
	directory named after the channel to be used as input for denoising
	'''
	print("\n\nSplitting images and converting them to tiff\n\n")
	img=ND2Reader(os.path.join(image_dir,image_name))
	print("Image shape is: "+str(img.shape))
	if(len(img.shape)>3):
		print("splitting channels")
		img.bundle_axes = 'zyx'
		img.iter_axes = 'c'
		for i in range(len(channel_list)): 
			tempChannel_path=os.path.join(image_dir, 'n2v_denoise/TMP' + channel_list[i])
			if not os.path.isdir(tempChannel_path):
				os.makedirs(tempChannel_path, exist_ok=True)
			imwrite(os.path.join(tempChannel_path, re.sub("\.nd2" ,  "_" + channel_list[i] + ".tif", image_name)), img[i])
	else:
		print("saving directly as tiff")
		tempChannel_path=os.path.join(image_dir, 'n2v_denoise/TMP' + channel_list[0])
		if not os.path.isdir(tempChannel_path):
			os.makedirs(tempChannel_path, exist_ok=True)
		imwrite(os.path.join(tempChannel_path, re.sub("\.nd2" ,  "_" + channel_list[0] + ".tif", image_name)), img)
	return None



def get_args():
	'''get command line arguments'''
	parser = argparse.ArgumentParser()
	parser.add_argument("-d","--image_dir", help="Path to directory containing nd2 images", type=str, required=True)
	parser.add_argument("-i","--image_name", help="Name of nd2 image to be processed", type=str, required=True)
	parser.add_argument("-m","--model_dir", help="Path to directory containing the models", type=str, required=True)
	parser.add_argument("-n","--model_base_name", help="Base name of the models needed for denoising",default='n2v_3D_CREST_')
	parser.add_argument("-c","--channel_list", action="extend", nargs="+", type=str, help="Ordered list of channel names",required=True,default=[])
	return parser.parse_args()


def main():
	'''main'''
	args = get_args()
	image_dir=args.image_dir
	image_name=args.image_name
	model_dir=args.model_dir
	model_base_name=args.model_base_name
	channel_list=args.channel_list
	print("channels are:")
	print(*channel_list)
	# split channels and convert to tiff 
	splitChannels(image_dir, image_name, channel_list)
	# cycle through and denoise
	for channel_name in channel_list:
		model_name=model_base_name+channel_name
		n2vDenoiseWithTrainedModel(model_dir, model_name, image_dir, image_name, channel_name ,plotInputPredict=True)
		split_image_dir=os.path.join(image_dir, 'n2v_denoise/TMP' + channel_name)
		split_image_name=re.sub("\.nd2", "_" + channel_name + ".tif", image_name)
		os.remove(os.path.join(split_image_dir,split_image_name))
	return None

if __name__=='__main__':
	main()

# #move to the folder with our data
# channel_color='red'
# os.chdir(workDir)
# model_name = 'n2v_3D_CREST_green'
# model_dir = '/mnt/external.data/MeisterLab/jsemple/'+'/n2v_denoise/training/models/'
# image_name='C1-SAMS-3_L1s_0014_green.tif'
# channel_list=["green","red"]
#channel_list=["green","red"]
