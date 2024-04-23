# coding: utf-8

import numpy as np
import pandas as pd
import pickle as pk
import math
from sklearn.ensemble import BaggingClassifier
from sklearn.preprocessing import StandardScaler
from sklearn.svm import SVC
from sklearn.metrics import roc_curve, roc_auc_score, accuracy_score, precision_score, confusion_matrix, recall_score, f1_score, auc, matthews_corrcoef
import joblib

modelfile = 'cyspair_incsingcysprob_5folddata_ws13.model'
# load the model from disk
loaded_model = joblib.load(modelfile)
output_path = './Predictions_IncSingCysProb_5folddata/'
list_file_path = '../Input/id_list.txt'
test_data_path = './CysPairWindowedFiles_IncCysProb/'

with open(list_file_path, "r") as f:
	for id in f:
		id = id.strip()
		test_protein_path = test_data_path+id+'_cyspair_incsingcysprob_ws_13.txt'
		print(test_protein_path)
		# read the test data file
		test_df = pd.read_csv(test_protein_path, header=None)

		test = test_df.to_numpy()
		y_test = test[:,0]
		X_test = test[:,1:]
		
		y_pred = loaded_model.predict(X_test)
		y_pred_prob = loaded_model.predict_proba(X_test)
		y_pred = np.column_stack([y_pred, y_pred_prob])
		
		#save the prediction output of test data		
		output_file = output_path+id+'.predict'
		out_file = open(output_file, 'wb')
		np.savetxt(fname=output_file, X=y_pred, fmt='%d %0.4f %0.4f', header='predClass, probN, probB', comments='')

with open(list_file_path, "r") as f:
	for id in f:
		id = id.strip()
		file1=open("../Output/"+id+'.predict',"w")
		file2=open('./IndividualCysPredictionProb/'+id+'.predict',"r")
		file3=open(output_path+id+'.predict',"r")
		# read content from first file
		file1.write("IndividualCysPredictionProb \n")
		for line in file2:			
			file1.write(line)

		file1.write("CysPairPredictionProb \n")
		for line in file3:			
			file1.write(line)

		file1.close()
		file2.close()
		file3.close()









