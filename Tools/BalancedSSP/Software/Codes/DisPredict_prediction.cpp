/*
Athor: Sumaiya Iqbal
Run prediction by DisPredict by executing libSVM
Part of BalancedSSP software
*/

#define _CRT_SECURE_NO_DEPRECATE
#define MAX_LINE_SIZE 30000 //Maximum line size that can be read from the input file
#define WINDOW_SIZE 21															// window size
#define MAX_FEATURE_COUNT 2000
#define NO_OF_FEATURES 56

#include<stdio.h>
#include<stdlib.h>
#include<string.h>
#include <malloc.h>
#include<ctype.h>

FILE *mvalues;
FILE *bvalues;
FILE *list;
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
char normalOuputFile[] = "../Output/log/log_DisPredict_prediction.txt";

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
	//printf("%s\n", argv[1]);
	//printf("%s\n", argv[2]);
	
	strcpy(id, argv[1]);
	//printf("%s\n", id);
	strcpy(libSVM_d, argv[2]);
	
	//printf("%s\n", id);

	nfp = fopen(normalOuputFile, "wb+");
	if (nfp == NULL) {
		fprintf(stderr, "Can't create normal output file!\n");
		exit(1);
	}
	
	
	//printf("\nlibSVM Path: %s\n", libSVM_d);
	
	//printf("\nPrediction START...\n");
	
	
	fprintf(nfp, "%s\n", id);
	
	
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
	seqLength = strlen(fsequence);
	seqLength--;
	//printf("ID: %s --- Length: %d\n", id, seqLength);
			
	//----------------------------------------------------------------------------------------------------------------------------------
	// scale & predict (SL)
	// prepare svm-scale command (SL, W-21)
	svm_scale_command[0] = '\0';
	strcpy(svm_scale_command, libSVM_d);
	strcat(svm_scale_command, "/svm-scale -r ");
	strcat(svm_scale_command, "../Models/DisPredict_SVM/DisPredict_scale_range_SL477 ");
	strcat(svm_scale_command, "../Features/");
	strcat(svm_scale_command, id);
	strcat(svm_scale_command, "/");
	strcat(svm_scale_command, id);
	strcat(svm_scale_command, ".final.dispredict.input > ../Output/prediction/");
	strcat(svm_scale_command, id);
	strcat(svm_scale_command, "/");
	strcat(svm_scale_command, id);
	strcat(svm_scale_command, ".dispredict.sl477.scale");
	fprintf(nfp, "SVM SCALE Command: %s\n", svm_scale_command);
	system(svm_scale_command);
	//printf("ID: %s---SL477---SCALE----DONE\n", id);

	// prepare svm-predict command and execute (SL, W-21)
	svm_predict_command[0] = '\0';
	strcpy(svm_predict_command, libSVM_d);
	strcat(svm_predict_command, "/svm-predict -b 1 ");
	strcat(svm_predict_command, "../Output/prediction/");
	strcat(svm_predict_command, id);
	strcat(svm_predict_command, "/");
	strcat(svm_predict_command, id);
	strcat(svm_predict_command, ".dispredict.sl477.scale ");
	strcat(svm_predict_command, "../Models/DisPredict_SVM/DisPredict_SL477.model ");
	strcat(svm_predict_command, "../Output/prediction/");
	strcat(svm_predict_command, id);
	strcat(svm_predict_command, "/");
	strcat(svm_predict_command, id);
	strcat(svm_predict_command, ".dispredict.sl477.predict");
	strcat(svm_predict_command, " >> ../Output/log/log_");
	strcat(svm_predict_command, id);
	strcat(svm_predict_command, "_SL477_prediction");
	fprintf(nfp, "SVM PREDICT Command: %s\n", svm_predict_command);
	system(svm_predict_command);
		
	//printf("ID: %s---SL477---PREDICT----DONE\n", id);
	
	//----------------------------------------------------------------------------------------------------------------------------------
	// Process Output
	// open predicted output file (annotation + both probability)
	prediction_output_file_name[0] = '\0';
	strcpy(prediction_output_file_name, "../Output/prediction/");
	strcat(prediction_output_file_name, id);
	strcat(prediction_output_file_name, "/");
	strcat(prediction_output_file_name, id);
	strcat(prediction_output_file_name, ".dispredict.sl477.predict");
	
	pred = fopen(prediction_output_file_name, "r");
	if (prediction_output_file_name == NULL) {
		fprintf(stderr, "Can't open prediction output file!\n");
		exit(1);
	}
	fgets(output_line, sizeof output_line, pred);
		
	prediction_output_file_name_procss[0] = '\0';
	strcpy(prediction_output_file_name_procss, "../Output/prediction/");
	strcat(prediction_output_file_name_procss, id);
	strcat(prediction_output_file_name_procss, "/SL477/");
	strcat(prediction_output_file_name_procss, id);
	strcat(prediction_output_file_name_procss, ".drp");
	
	pred_procss = fopen(prediction_output_file_name_procss, "wb+");
	if (prediction_output_file_name_procss == NULL) {
		fprintf(stderr, "Can't open dpredict output file!\n");
		exit(1);
	}
	fprintf(pred_procss, "PFRMAT DR\n");
	fprintf(pred_procss, "TARGET %s\n", id);
	fprintf(pred_procss, "AUTHOR DisPredict: Sumaiya Iqbal and Md Tamjidul Hoque\n");
	fprintf(pred_procss, "REMARK SVM\n");
	fprintf(pred_procss, "METHOD SVM (SL477)\n");
	fprintf(pred_procss, "MODEL  1\n");
	//fprintf(pred_procss, "\n");

	//fprintf(pred_procss, ">%s\n", id);
	//fprintf(pred_procss, "AA O/D DisorderProbability\n");
	
	track = 0;
	while (fgets(output_line, sizeof output_line, pred) != NULL)
	{
			wline[0] = '\0';
			
			char* fresidue1 = substring(fsequence, track + 1, 1);		// residues	
			strcpy(wline, fresidue1);
			strcat(wline, "   ");
			track++;
			
			if (output_line[0] == '1')						// Write annotation
			{
				strcat(wline, "D");
			}
			else
			{
				strcat(wline, "O");
			}
			strcat(wline, "  ");

			char *tdp;
			tdp = strtok(output_line, " ");
			int token_track = 0;
			while (tdp != NULL)
			{
				if (token_track == 1)
				{
					strcat(wline, tdp);
					int sc = 0;
					while(sc < 20 - strlen(tdp))
					{
						strcat(wline, " ");
						sc++;
					}
												
				}
				if(token_track == 2)
				{
					//strcat(wline, tdp);
					break;
				}
				tdp = strtok(NULL, " ");
				token_track++;
			}
			fprintf(pred_procss, "%s\n", wline);				// write annotation and both probability
	}
	fprintf(pred_procss, "END\n");
	fclose(pred);
	fclose(pred_procss);
	//printf("ID: %s---PREDICTION OUTPUT PROCESSING----DONE\n", id);
	//----------------------------------------------------------------------------------------------------------------------------------
	//----------------------------------------------------------------------------------------------------------------------------------

	
	//----------------------------------------------------------------------------------------------------------------------------------
	rline[0] = '\0';
	
	//printf("Prediction DONE!!...\n");

	//printf("\n------------ DisPredict End------------------\n");
	
	fclose(nfp);
}