/*
Author: Sumaiya Iqbal
Part of balancedSSP software
---- predict secondary structure by binary model
---- extension of output (.binary.beta.predict/ .binary.coil.predict. / .binary.helix.predict)
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
char normalOuputFile[] = "../Output/log/binary_svm_prediction.txt";

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

	// BETA
	//==========================================================================================================================
	// scale & predict (beta)
	// prepare svm-scale command (beta, W-21)
	svm_scale_command[0] = '\0';
	strcpy(svm_scale_command, libSVM_d);
	strcat(svm_scale_command, "/svm-scale -r ");
	strcat(svm_scale_command, "../Models/binary/scale_range_f33_b1 ");
	strcat(svm_scale_command, "../Features/");
	strcat(svm_scale_command, id);
	strcat(svm_scale_command, "/");
	strcat(svm_scale_command, id);
	strcat(svm_scale_command, ".binary.input > ../Features/");
	strcat(svm_scale_command, id);
	strcat(svm_scale_command, "/");
	strcat(svm_scale_command, id);
	strcat(svm_scale_command, ".binary.beta.scale");
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
	strcat(svm_predict_command, ".binary.beta.scale ");
	strcat(svm_predict_command, "../Models/binary/train_b1_f33_W_15.scale.model ");
	strcat(svm_predict_command, "../Features/");
	strcat(svm_predict_command, id);
	strcat(svm_predict_command, "/");
	strcat(svm_predict_command, id);
	strcat(svm_predict_command, ".binary.beta.predict");
	strcat(svm_predict_command, " >> ../Output/log/log_");
	strcat(svm_predict_command, id);
	strcat(svm_predict_command, "_binary_beta_prediction");
	fprintf(nfp, "SVM PREDICT Command: %s\n", svm_predict_command);
	system(svm_predict_command);
	//==========================================================================================================================

	// COIL
	//==========================================================================================================================
	// scale & predict (coil)
	// prepare svm-scale command (coil, W-15)
	svm_scale_command[0] = '\0';
	strcpy(svm_scale_command, libSVM_d);
	strcat(svm_scale_command, "/svm-scale -r ");
	strcat(svm_scale_command, "../Models/binary/scale_range_f33_c1 ");
	strcat(svm_scale_command, "../Features/");
	strcat(svm_scale_command, id);
	strcat(svm_scale_command, "/");
	strcat(svm_scale_command, id);
	strcat(svm_scale_command, ".binary.input > ../Features/");
	strcat(svm_scale_command, id);
	strcat(svm_scale_command, "/");
	strcat(svm_scale_command, id);
	strcat(svm_scale_command, ".binary.coil.scale");
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
	strcat(svm_predict_command, ".binary.coil.scale ");
	strcat(svm_predict_command, "../Models/binary/train_c1_f33_W_15.scale.model ");
	strcat(svm_predict_command, "../Features/");
	strcat(svm_predict_command, id);
	strcat(svm_predict_command, "/");
	strcat(svm_predict_command, id);
	strcat(svm_predict_command, ".binary.coil.predict");
	strcat(svm_predict_command, " >> ../Output/log/log_");
	strcat(svm_predict_command, id);
	strcat(svm_predict_command, "_binary_coil_prediction");
	fprintf(nfp, "SVM PREDICT Command: %s\n", svm_predict_command);
	system(svm_predict_command);
	//==========================================================================================================================

	// HELIX
	//==========================================================================================================================
	// scale & predict (helix)
	// prepare svm-scale command (helix, W-15)
	svm_scale_command[0] = '\0';
	strcpy(svm_scale_command, libSVM_d);
	strcat(svm_scale_command, "/svm-scale -r ");
	strcat(svm_scale_command, "../Models/binary/scale_range_f33_h1 ");
	strcat(svm_scale_command, "../Features/");
	strcat(svm_scale_command, id);
	strcat(svm_scale_command, "/");
	strcat(svm_scale_command, id);
	strcat(svm_scale_command, ".binary.input > ../Features/");
	strcat(svm_scale_command, id);
	strcat(svm_scale_command, "/");
	strcat(svm_scale_command, id);
	strcat(svm_scale_command, ".binary.helix.scale");
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
	strcat(svm_predict_command, ".binary.helix.scale ");
	strcat(svm_predict_command, "../Models/binary/train_h1_f33_W_15.scale.model ");
	strcat(svm_predict_command, "../Features/");
	strcat(svm_predict_command, id);
	strcat(svm_predict_command, "/");
	strcat(svm_predict_command, id);
	strcat(svm_predict_command, ".binary.helix.predict");
	strcat(svm_predict_command, " >> ../Output/log/log_");
	strcat(svm_predict_command, id);
	strcat(svm_predict_command, "_binary_helix_prediction");
	fprintf(nfp, "SVM PREDICT Command: %s\n", svm_predict_command);
	system(svm_predict_command);
	//==========================================================================================================================

	fclose(nfp);
}