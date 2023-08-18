/*
Author: Sumaiya Iqbal
Part of balancedSSP software
Collect 55 features for initial ASA prediction by REGAd3p
*/

#define _CRT_SECURE_NO_DEPRECATE
#define LINE_SIZE 30000 //Maximum line size that can be read from the input file
#define BIGRAM_MONOGRAM_LOG_MEDIAN 6.0

#include<stdio.h>
#include<stdlib.h>
#include<string.h>
#include <malloc.h>
#include<ctype.h>
#include<math.h>

FILE *nfp;
FILE *fasta;
FILE *svm_predict;
FILE *feature2;
FILE *feature1;

char terminal[10];
char normalOuputFile[] = "../Output/log/log_feature_collection_for_ASA.txt";
char ppFile[] = "../AdditionalFiles/physiochemical_properties.txt";
char id[LINE_SIZE];
char svm_predict_file_name[LINE_SIZE];
char svm_prediction[LINE_SIZE];
char feature_set2_file_name[LINE_SIZE];
char feature_set1_file_name[LINE_SIZE];

char rline[LINE_SIZE];
char fseq[LINE_SIZE];
char all_feature[LINE_SIZE];
char feature_part_1[LINE_SIZE];
char all_feature_filename[300];
char fastaFileName[300];

int aa_counter = 0;
int pp_counter = 0;
int lineLength = 0;
int seqNumber = 0;
int seqLength = 0;
int seqNum = 0;
int local_counter = 0;

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
	// collect ID
	strcpy(id, argv[1]);

	nfp = fopen(normalOuputFile, "wb+");
	if (nfp == NULL) {
		fprintf(stderr, "Can't log output file!\n");
		exit(1);
	}

	//====================================================================================================================================
	// open svm predict file for ss probabilities
	strcpy(svm_predict_file_name, "../Features/");
	strcat(svm_predict_file_name, id);
	strcat(svm_predict_file_name, "/");
	strcat(svm_predict_file_name, id);
	strcat(svm_predict_file_name, ".svm.ss.predict");

	svm_predict = fopen(svm_predict_file_name, "r");					   // Open svm predict from step1 file for reading
	if (svm_predict == NULL) {
		fprintf(stderr, "Can't open svm predict file!\n");
		exit(1);
	}

	// skip header
	fgets(svm_prediction, sizeof svm_prediction, svm_predict);	// skipp header
	
	//====================================================================================================================================

	//====================================================================================================================================
	// prepare fasta file name, open and collect sequence length
	strcpy(fastaFileName, "../Features/");
	strcat(fastaFileName, id);
	strcat(fastaFileName, "/");
	strcat(fastaFileName, id);
	strcat(fastaFileName, ".fasta");


	//printf("Fasta File name: %s\n", fastaFileName);

	fasta = fopen(fastaFileName, "r");
	if (fasta == NULL) {
		fprintf(stderr, "Can't open fasta File!\n");
		exit(1);
	}

	fgets(fseq, sizeof fseq, fasta);	// skip header
	fgets(fseq, sizeof fseq, fasta);	// collect sequence
	fclose(fasta);

	int tr = 0;
	seqLength = 0;
	while ((fseq[tr] >= 'A') && (fseq[tr] <= 'Z'))
	{
		seqLength++;
		tr++;
	}
	//====================================================================================================================================

	//====================================================================================================================================
	// prepare feature set - step 2 file name and open to write
	strcpy(feature_set2_file_name, "../Features/");
	strcat(feature_set2_file_name, id);
	strcat(feature_set2_file_name, "/");
	strcat(feature_set2_file_name, id);
	strcat(feature_set2_file_name, ".ASA.features");
	fprintf(nfp, "feature_set2_file_name: %s\n", feature_set2_file_name);

	feature2 = fopen(feature_set2_file_name, "wb+");					// Open all feature2 file for writing
	if (feature2 == NULL) {
		fprintf(stderr, "Can't open rsa feature file!\n");
		exit(1);
	}
	//====================================================================================================================================

	//====================================================================================================================================
	// prepare seature set 1 file name and read on that file
	strcpy(feature_set1_file_name, "../Features/");
	strcat(feature_set1_file_name, id);
	strcat(feature_set1_file_name, "/");
	strcat(feature_set1_file_name, id);
	strcat(feature_set1_file_name, ".initialSS.features");
	

	feature1 = fopen(feature_set1_file_name, "r");			// Open feature step 1 file for reading
	if (feature1 == NULL) {
		fprintf(stderr, "Can't open feature1 file!\n");
		exit(1);
	}
	//====================================================================================================================================
	
	//====================================================================================================================================
	// Merger all features
	// Feature serial is: class, residue indication, pp, pssm, ss, sa, torsion angle, monogram, bigram, terminal
	int res_pos = 0;
	while (res_pos < seqLength)
	{
		//====================================================================================================================================
		// Take ASA annottaion --- dummy 0
		strcpy(all_feature, "0 ");
		//====================================================================================================================================

		//====================================================================================================================================
		// take feature set by step 1
		fgets(rline, sizeof rline, feature1);
	
		int tk_count = 0;
		char *tk;
		
		tk = strtok(rline, " ");
		while (tk != NULL)
		{
			tk_count++;
			if ((tk_count >= 2) && (tk_count <= 52))
			{
				strcat(all_feature, tk);						// take aa + pp + pssm + bi/mono + iu == 51 features
				strcat(all_feature, " ");
			}
			if (tk_count == 53)
			{
				strcpy(terminal, tk);					// save terminal
			}
			tk = strtok(NULL, " ");
		}
		//====================================================================================================================================
		
		//====================================================================================================================================

		// collection svm prediction probabilities -- SS0 information
		fgets(svm_prediction, sizeof svm_prediction, svm_predict);

		// collect SS0 probabilities

		char *ss_prob;
		int ss_count = 0;
		ss_prob = strtok(svm_prediction, " ");
		while (ss_prob != NULL)
		{
			ss_count++;
			if ((ss_count >= 2) && (ss_count <= 4))
			{
				strcat(all_feature, trim(ss_prob));						// take ss0 probabilities == 3 features
				strcat(all_feature, " ");
			}
			ss_prob = strtok(NULL, " ");
		}
		//====================================================================================================================================

		//====================================================================================================================================
		// take terminal information 
		strcat(all_feature, trim(terminal));
		//====================================================================================================================================

		//====================================================================================================================================
		fprintf(feature2, "%s\n", all_feature);
		//====================================================================================================================================

		res_pos++;
	}

	fclose(nfp);
	fclose(feature2);
	fclose(feature1);
	fclose(svm_predict);
}