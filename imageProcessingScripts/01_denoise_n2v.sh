#!/bin/bash
#SBATCH --time=0-05:00:00
#SBATCH --mem 32GB
#SBATCH --gres=gpu:1

IMAGE_DIR=/mnt/external.data/MeisterLab/jsemple/microscopy/20231213_941-9_SMC1GFP_HS
MODEL_DIR=/mnt/external.data/MeisterLab/jsemple/microscopy/20231213_941-9_SMC1GFP_HS/n2v_denoise/training/models
MODEL_BASE_NAME=n2v_3D_CREST_
CHANNEL_LIST=(green)

source $HOME/miniforge3/bin/activate n2v

# use some expression to get list of files to process 
FILE_LIST=(`ls *.nd2 | grep -v _bf.nd2`)

for IMAGE_NAME in ${FILE_LIST[@]}
do
  echo "processing " $IMAGE_NAME
  python ./n2vDenoiseWithTrainedModel.py -d $IMAGE_DIR -i $IMAGE_NAME -m $MODEL_DIR -n $MODEL_BASE_NAME -c ${CHANNEL_LIST[@]} 
done
