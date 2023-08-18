/*
Author: Sumaiya Iqbal
Part of balancedSSP software
---- predict secondary structure by svm model
process output into foramtted file (.SSp)
*/
#define _CRT_SECURE_NO_DEPRECATE
#define MAX_LINE_SIZE 30000 //Maximum line size that can be read from the input file
#define WINDOW_SIZE 21															// window size
#define MAX_FEATURE_COUNT 2000
#define NO_OF_FEATURES 52

#include<stdio.h>
#include<stdlib.h>
#include<string.h>
#include <malloc.h>
#include<ctype.h>


FILE *nfp;
FILE *allf;
FILE *svm;
FILE *pred;
FILE *pred_procss;
FILE *libpath;
FILE *fasta;

char fsequence[MAX_LINE_SIZE];
char fastaFileName[MAX_LINE_SIZE];
char make_directory_command[MAX_LINE_SIZE];
char libSVM_d[MAX_LINE_SIZE];
char prediction_output_file_name[MAX_LINE_SIZE];
char prediction_output_file_name_procss[MAX_LINE_SIZE];
char output_line[MAX_LINE_SIZE];
char normalOuputFile[] = "../Output/log/ss_svm_prediction.txt";

char rline[MAX_LINE_SIZE];
char svm_scale_command[MAX_LINE_SIZE];
char svm_predict_command[MAX_LINE_SIZE];
char wline[MAX_LINE_SIZE];
char id[MAX_LINE_SIZE];

int lineLength = 0;
int seqNumber = 0;
int feature_count = 0;
int seqLength = 0;
int track = 0;
int spc = 0;
int constl = 0;

/*Start of SUB-ROUTINE: substring - It returns a pointer to the substring */
char *substring(char *string, int position, int length)
{
	char *pointer;
	int c;

	pointer = (char *)malloc(length + 1);

	if (pointer == NULL)
	{
		printf("Unable to allocate memory.\n");
		exit(EXIT_FAILURE);
	}

	for (c = 0; c < position - 1; c++)
		string++;

	for (c = 0; c < length; c++)
	{
		*(pointer + c) = *string;
		string++;
	}

	*(pointer + c) = '\0';

	return pointer;
}
/*End of SUB-ROUTINE: substring - It returns a pointer to the substring */

/* SUB-ROUTINE: Trim */
char* ltrim(char *s)
{
	while (isspace(*s)) s++;
	return s;
}

char *rtrim(char *s)
{
	char* back = s + strlen(s);
	while (isspace(*--back));
	*(back + 1) = '\0';
	return s;
}

char *trim(char *s)
{
	return rtrim(ltrim(s));
}

int main(int argc, char *argv[])
{
	
	//==========================================================================================================================
	// collect id and path of libsvm as argument
	strcpy(id, argv[1]);
	strcpy(libSVM_d, argv[2]);
	//==========================================================================================================================
	
	//==========================================================================================================================
	nfp = fopen(normalOuputFile, "wb+");
	if (nfp == NULL) {
		fprintf(stderr, "Can't create normal output file!\n");
		exit(1);
	}
	fprintf(nfp, "%s\n", id);
	//==========================================================================================================================

	//==========================================================================================================================
	// Prepare fasta file name, open and collect fasta sequence for formatted output
	strcpy(fastaFileName, "../Features/");
	strcat(fastaFileName, id);
	strcat(fastaFileName, "/");
	strcat(fastaFileName, id);
	strcat(fastaFileName, ".fasta");

	fasta = fopen(fastaFileName, "r");
	if (fasta == NULL) {
		fprintf(stderr, "Can't open fasta File!\n");
		exit(1);
	}

	fgets(fsequence, sizeof fsequence, fasta);	// skip header
	fgets(fsequence, sizeof fsequence, fasta);	// collect sequence
	int tr = 0;
	seqLength = 0;
	while ((fsequence[tr] >= 'A') && (fsequence[tr] <= 'Z'))
	{
		seqLength++;
		tr++;
	}
	//==========================================================================================================================

	//==========================================================================================================================
	// scale & predict
	// prepare svm-scale command (SL, W-21)
	svm_scale_command[0] = '\0';
	strcpy(svm_scale_command, libSVM_d);
	strcat(svm_scale_command, "/svm-scale -r ");
	strcat(svm_scale_command, "../Models/SS_SVM/scale_range_libsvm ");
	strcat(svm_scale_command, "../Features/");
	strcat(svm_scale_command, id);
	strcat(svm_scale_command, "/");
	strcat(svm_scale_command, id);
	strcat(svm_scale_command, ".initialSS.input > ../Features/");
	strcat(svm_scale_command, id);
	strcat(svm_scale_command, "/");
	strcat(svm_scale_command, id);
	strcat(svm_scale_command, ".svm.ss.scale");
	system(svm_scale_command);
	//==========================================================================================================================

	//==========================================================================================================================
	// prepare svm-predict command and execute (SL, W-21)
	svm_predict_command[0] = '\0';
	strcpy(svm_predict_command, libSVM_d);
	strcat(svm_predict_command, "/svm-predict -b 1 ");
	strcat(svm_predict_command, "../Features/");
	strcat(svm_predict_command, id);
	strcat(svm_predict_command, "/");
	strcat(svm_predict_command, id);
	strcat(svm_predict_command, ".svm.ss.scale ");
	strcat(svm_predict_command, "../Models/SS_SVM/SSD_TR1001.libsvm.scale.model ");
	strcat(svm_predict_command, "../Features/");
	strcat(svm_predict_command, id);
	strcat(svm_predict_command, "/");
	strcat(svm_predict_command, id);
	strcat(svm_predict_command, ".svm.ss.predict");
	strcat(svm_predict_command, " >> ../Output/log/log_");
	strcat(svm_predict_command, id);
	strcat(svm_predict_command, "_svm_prediction");
	fprintf(nfp, "SVM PREDICT Command: %s\n", svm_predict_command);
	system(svm_predict_command);
	//==========================================================================================================================
	
	//==========================================================================================================================
	// Process Output
	// open predicted output file (annotation + both probability)
	prediction_output_file_name[0] = '\0';
	strcpy(prediction_output_file_name, "../Features/");
	strcat(prediction_output_file_name, id);
	strcat(prediction_output_file_name, "/");
	strcat(prediction_output_file_name, id);
	strcat(prediction_output_file_name, ".svm.ss.predict");

	pred = fopen(prediction_output_file_name, "r");
	if (prediction_output_file_name == NULL) {
		fprintf(stderr, "Can't open prediction output file!\n");
		exit(1);
	}
	fgets(output_line, sizeof output_line, pred);	// skip header
	//==========================================================================================================================

	//==========================================================================================================================
	// prepare processed output file
	prediction_output_file_name_procss[0] = '\0';
	strcpy(prediction_output_file_name_procss, "../Output/prediction/");
	strcat(prediction_output_file_name_procss, id);
	strcat(prediction_output_file_name_procss, "/SS/");
	strcat(prediction_output_file_name_procss, id);
	strcat(prediction_output_file_name_procss, ".SSp");

	pred_procss = fopen(prediction_output_file_name_procss, "wb+");
	if (prediction_output_file_name_procss == NULL) {
		fprintf(stderr, "Can't open dpredict output file!\n");
		exit(1);
	}
	fprintf(pred_procss, "Secondary Structure Prediction Output (Part of REGAd3p ASA Prediction)\n");
	fprintf(pred_procss, "TARGET: %s, ", id);
	fprintf(pred_procss, "Length: %d\n", seqLength);
	fprintf(pred_procss, "\n");
	fprintf(pred_procss, "SR#  AA  SSp  pCoil          pBeta          pHelix\n");
	//==========================================================================================================================

	//==========================================================================================================================
	track = 0;
	while (fgets(output_line, sizeof output_line, pred) != NULL)
	{
		//==========================================================================================================================
		// write serial number
		sprintf(wline, "%d", track + 1);
		spc = 0;
		constl = strlen(wline);
		while (spc < 5 - constl)
		{
			strcat(wline, " ");
			spc++;
		}
		//==========================================================================================================================

		//==========================================================================================================================
		// write amino acid
		char* fresidue1 = substring(fsequence, track + 1, 1);		// residues	
		strcat(wline, fresidue1);
		strcat(wline, "   ");
		track++;
		//==========================================================================================================================

		//==========================================================================================================================
		if (output_line[0] == '1')						// Write annotation
		{
			strcat(wline, "H");
		}
		else if (output_line[0] == '2')
		{
			strcat(wline, "E");
		}
		else if (output_line[0] == '4')
		{
			strcat(wline, "C");
		}
		strcat(wline, "    ");
		//==========================================================================================================================

		//==========================================================================================================================
		char *tdp;
		tdp = strtok(output_line, " ");
		int token_track = 0;
		while (tdp != NULL)
		{
			if (token_track == 1)
			{
				strcat(wline, tdp);
				spc = 0;
				constl = strlen(wline);
				while (spc < 29 - constl)
				{
					strcat(wline, " ");
					spc++;
				}
			}
			if (token_track == 2)
			{
				strcat(wline, tdp);
				spc = 0;
				constl = strlen(wline);
				while (spc < 44 - constl)
				{
					strcat(wline, " ");
					spc++;
				}
			}
			if (token_track == 3)
			{
				strcat(wline, tdp);
				strcat(wline, " ");
				int w_l = strlen(wline);
				wline[w_l] = '\0';
			}
			tdp = strtok(NULL, " ");
			token_track++;
		}
		fprintf(pred_procss, "%s\n", trim(wline));				// write annotation and both probability
	}
	fprintf(pred_procss, "END\n");
	fclose(pred);
	fclose(pred_procss);

	fclose(nfp);
}