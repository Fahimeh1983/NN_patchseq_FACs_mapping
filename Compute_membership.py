import itertools
import argparse
import numpy as np
import pandas as pd
import tensorflow as tf
from pathlib import Path
import Ttools.utils as utils
from tensorflow import keras
from tensorflow.keras import regularizers
from sklearn.metrics import confusion_matrix
from sklearn.model_selection import train_test_split
from tensorflow.keras.layers import BatchNormalization




parser = argparse.ArgumentParser()
parser.add_argument("--data_dir",   default='Manuscript_patchseq_2019/Taxonomy',   help='parent dir')
parser.add_argument("--keras_model_dir", default='Keras_models',   help='model dir')
parser.add_argument("--output_dir", default='Manuscript_patchseq_2019/Taxonomy/',   help='model dir')
parser.add_argument("--model_id", default='TEMP',   help='model dir')

def main(data_dir = 'Manuscript_patchseq_2019', keras_model_dir = 'Keras_models', 
         output_dir = "Manuscript_patchseq_2019/Taxonomy/", model_id = "TEMP" ):

    print("Loading data")
    query_dat_norm = pd.read_csv(data_dir + "/query_dat_norm.csv", index_col= 0)
    print("Load:" + data_dir + "/query_dat_norm.csv")
    FACS_dat_norm = pd.read_csv(data_dir + "/FACs_norm.csv", index_col= 0)
    print("Load:" + data_dir + "/FACS_norm.csv")
    query_dat_norm = query_dat_norm.T

    patchseq_cells = query_dat_norm.index.tolist()
    FACS_cells = FACS_dat_norm.columns.tolist()
    V1_cl = pd.read_csv(data_dir + "select_cl.csv")['x'].tolist()

    FACS_membership = pd.DataFrame(0, index=FACS_cells, columns=V1_cl)
    Patchseq_membership = pd.DataFrame(0, index=patchseq_cells, columns=V1_cl)

    for run_id in range(100):
         model_name = model_id + str(run_id + 1) + ".h5"
         facs_prediction_file =  "facs_membership_" + str(run_id+1) + ".csv"
         print(model_name)
   
         model = keras.models.load_model(keras_model_dir + model_name)

         facs_memb = pd.read_csv(keras_model_dir + facs_prediction_file, index_col=0)
         FACS_membership = FACS_membership + facs_memb.loc[FACS_cells][V1_cl]
         patchseq_memb = model.predict(query_dat_norm)
         patchseq_memb = pd.DataFrame(patchseq_memb, index= patchseq_cells, columns=V1_cl)
         Patchseq_membership.loc[patchseq_cells] = Patchseq_membership.loc[patchseq_cells] + patchseq_memb.loc[patchseq_cells]

    temp = FACS_membership/FACS_membership.sum(axis=1)[:,None]
    temp.to_csv(output_dir + "/NN_imb_1000epoch_500batch_FACS_membership.csv")

    temp = Patchseq_membership/Patchseq_membership.sum(axis=1)[:,None]
    temp.to_csv(output_dir + "/NN_imb_1000epoch_500batch_patchseq_membership.csv")

if __name__ == "__main__":
    args = parser.parse_args()
    main(**vars(args))
