import numpy as np
import pandas as pd
import pickle as pk
from sklearn.model_selection import cross_val_predict
from sklearn.preprocessing import StandardScaler
from sklearn import svm
from sklearn.metrics import roc_curve, roc_auc_score, accuracy_score, precision_score, confusion_matrix, recall_score, f1_score, auc, matthews_corrcoef
from lightgbm import LGBMClassifier

# set up input (data) and output path
output_path = './IndividualCysPredictionProb/'

test_data_path = './SingleCysWindowedFile/'
list_file_path = '../Input/id_list.txt'

modelfile = 'singlecys_collectindivcysprob_ws31.model'

scalerfile = 'singlecys_collectindivcysprob_ws31.scaler'

# load the model from disk
clf = pk.load(open(modelfile, 'rb'))
scaler = pk.load(open(scalerfile, 'rb'))

print('started predicting protein wise single cysteine bonding probability')

with open(list_file_path, "r") as f:
	for id in f:
		id = id.strip()
		test_protein_path = test_data_path+id+'_singlecys_W_31.txt'
		print(test_protein_path)
		# read the test data file
		test_df = pd.read_csv(test_protein_path, header=None)

		test = test_df.to_numpy()
		y_test = test[:,0]
		X_test = test[:,1:]
		X_test = scaler.transform(X_test) # the scaler instance is used on test data to transform it the same way it did on the training set
		
		y_pred = clf.predict(X_test)
		y_pred_prob = clf.predict_proba(X_test)
		y_pred = np.column_stack([y_pred, y_pred_prob])
		
		#save the prediction output of test data		
		output_file = output_path+id+'.predict'
		out_file = open(output_file, 'wb')
		np.savetxt(fname=output_file, X=y_pred, fmt='%d %0.4f %0.4f', header='predClass, probN, probB', comments='')

