#######################################################
### Fahimeh Baftizadeh: ###############################
#######################################################

# How to run this code from the command line:
# Make sure, you have a directory for output such as "new_keras_models/"
#python -m Keras_imbalanced_trainer --batch_size 500 --epochs 10 --run_iter 0 --n_features 4020 --n_celltypes 93 --n_hidden 10 --n_test_cells 10 --dropout 0.6

#######################################################
### Loading libraries: ################################
#######################################################
import csv
import timeit
import random
import argparse
import datetime
import itertools
import numpy as np
import pandas as pd
import tensorflow as tf

from pathlib import Path
from scipy.sparse import issparse

from sklearn.base import clone
from sklearn.utils import safe_indexing
from sklearn.metrics import confusion_matrix
from sklearn.utils import check_random_state
from sklearn.utils.testing import set_random_state
from sklearn.model_selection import train_test_split

from keras import regularizers
from keras.optimizers import Adam
from keras.utils import Sequence
from keras.models import Sequential
from keras.callbacks import CSVLogger
from keras.layers import Activation, Dense, Dropout, BatchNormalization, Flatten

from imblearn.utils import Substitution
from imblearn.utils._docstring import _random_state_docstring
from imblearn.tensorflow import balanced_batch_generator as tf_bbg
from imblearn.keras import BalancedBatchGenerator

#######################################################
### Reading inputs from the command line: #############
#######################################################

print(tf.__version__)
parser = argparse.ArgumentParser()
parser.add_argument("--batch_size", default=500,         type=int,    help="Batch size")
parser.add_argument("--epochs",     default=10,          type=int,    help="Number of training epochs")
parser.add_argument("--run_iter",   default=0,           type=int,    help="Run iteration")
parser.add_argument("--data_dir",   default='Manuscript_patchseq_2019',   help='parent dir')
parser.add_argument("--output_dir", default='Keras_models',   help='model dir')
parser.add_argument("--n_features", default=4020,          type=int,    help="Number of features")
parser.add_argument("--n_celltypes",default=93,          type=int,    help="Number of celltypes or labels")
parser.add_argument("--n_hidden",   default=10,          type=int,    help="size of hidden layer")
parser.add_argument("--n_test_cells",   default=10,          type=int,    help="N test cells")
parser.add_argument("--dropout",   default=0.6,          type=float,    help="drop out")


#######################################################
### Functions: ########################################
#######################################################

def read_data(path):
    file_path = path
    data = pd.read_csv(file_path)
    colnames = list(data['Unnamed: 0'])
    data = data.T
    data.columns =  colnames
    data = data.drop(axis=0, labels='Unnamed: 0')
    return data

def read_labels(path, select_cl):
    file_path = path
    labels = pd.read_csv(file_path)
    labels.index = list(labels['Unnamed: 0'])
    labels = labels.drop(axis=1, labels='Unnamed: 0')
    new_factors = np.arange(len(select_cl)).tolist()
    cls = select_cl
    refactor_cls = []
    for items in labels.cl:
        if items in cls: 
            index = cls.index(items)
            refactor_cls = refactor_cls + [new_factors[index]]
        else:
            refactor_cls = refactor_cls + [np.nan]
    labels["old_factor_cl"] = labels["factor_cl"]
    labels["factor_cl"] = refactor_cls 
    return labels

def split_data_intwo(data, labels, test_size, cvset):
    train_data, test_data, train_labels, test_labels = train_test_split(data, labels,
                                              test_size = test_size,
                                              random_state = cvset)

    return train_data, test_data, train_labels, test_labels

def make_model(n_features, HL, FL, dropout):
    model = Sequential()
    model.add(Dropout(dropout, input_shape=(n_features,)))
    model.add(Dense(HL, 
                    kernel_regularizer= regularizers.l2(0.01),
                   bias_regularizer=regularizers.l2(0.01)))
    model.add(BatchNormalization())
    model.add(Activation('relu'))
    model.add(Dense(FL, activation=tf.nn.softmax))
    
    model.compile(optimizer='adam', loss='sparse_categorical_crossentropy', metrics=['accuracy'])
    return model

#######################################################
### Keras trainer: ####################################
#######################################################

def main(batch_size=500, epochs=10, run_iter=0,  data_dir='/Manuscript_Patchseq_2019', 
         output_dir = "/new_Keras_models", n_features = 4020, n_celltypes = 93, 
         n_hidden=10, dropout = 0.6, n_test_cells = 1464):
   

    model_id = str(epochs) + '_' + str(batch_size)    
    fileid = model_id + '_ri_' + str(run_iter)
    facs_output_id = "facs_membership_" + str(run_iter)
    cells_output_id = model_id + '_testcells_' + str(run_iter)
    V1_cl = pd.read_csv(data_dir + "/select_cl.csv")['x'].tolist()

    FACs_data = read_data(data_dir + "/FACs_norm.csv")
    FACs_labels = read_labels(data_dir + "/FACs_correct_labels.csv", V1_cl)
    FACs_labels = FACs_labels['factor_cl']
    FACs_cells = FACs_data.index.tolist()
    FACs_membership = pd.DataFrame(0, index=FACs_cells, columns=V1_cl)
    
    train_data , test_data, train_labels, test_labels = split_data_intwo(FACs_data, FACs_labels, 
                                                                         test_size = n_test_cells, cvset = random.randint(0, 10000))
    test_cells = test_data.index.tolist()

    start_time = timeit.default_timer()
    model = make_model(n_features, n_hidden, n_celltypes, dropout)
    results = model.fit(train_data, train_labels, epochs=epochs, batch_size=batch_size,  verbose=0)
    
    facs_memb = model.predict(test_data)
    facs_memb = pd.DataFrame(facs_memb, index= test_cells, columns=V1_cl)
    FACs_membership.loc[test_cells] = FACs_membership.loc[test_cells] + facs_memb.loc[test_cells]
    
    print(datetime.datetime.now())

    elapsed = timeit.default_timer() - start_time
    score, acc = model.evaluate(test_data, test_labels,
                       batch_size=batch_size, verbose=0)

    print('Test accuracy:', acc)
    print('-------------------------------')
    print('Training time:', elapsed)
    print('-------------------------------')

    FACs_membership.to_csv(output_dir + facs_output_id + '.csv')
    model.save(output_dir + fileid + '.h5')
    with open(output_dir + cells_output_id + '.csv', 'w') as myfile:
    	wr = csv.writer(myfile)
    	wr.writerow(test_cells)

if __name__ == "__main__":
    args = parser.parse_args()
    main(**vars(args))
