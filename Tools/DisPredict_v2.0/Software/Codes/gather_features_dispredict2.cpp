/*
Author: Sumaiya Iqbal
Part of code for DisPredict 2 (Partial -- only with PSEE)
-- gather all feature per id into one file
*/
#define _CRT_SECURE_NO_DEPRECATE
#include<stdio.h>
#include<stdlib.h>
#include<string.h>
#include <malloc.h>
#include<ctype.h>
#include<math.h>
#define LINE_SIZE 30000 //Maximum line size that can be read from the input file
#define BIGRAM_MONOGRAM_LOG_MEDIAN 6.0
#define MAX_FEATURE 57
FILE *nfp;
FILE *allp;
FILE *alln;
FILE *spsee;
FILE *fasta;
char normalOuputFile[] = "../Output/log/log_feature_collection2.txt";
char fastaFileName[500];
char rline[LINE_SIZE];
char make_directory_command[LINE_SIZE];
char allFeature_file_prev[LINE_SIZE];
char allFeature_file_new[LINE_SIZE];
char sPSEE_filne_name[LINE_SIZE];
char fseq[LINE_SIZE];
char all_feature[LINE_SIZE];
char feature_part_1[LINE_SIZE];
char energy_line[LINE_SIZE];
char terminal[100];
char PSEE_val[500];
char id[200];
int lineLength = 0;
int seqNumber = 0;
int seqLength = 0;
int res_pos = 0;
int i = 0;
int l = 0;
int featureCount = 0;
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
	strcpy(id, argv[1]);
	
	nfp = fopen(normalOuputFile, "wb+");
	if (nfp == NULL) {
		fprintf(stderr, "Can't create normal output file!\n");
		exit(1);
	}
	
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
	//printf("Sequence No: %d ------------- ID: %s, Length: %d\n", seqNumber, id, seqLength);
	//====================================================================================================================================
	//====================================================================================================================================
	// prepare previous all feature file name and open to read
	strcpy(allFeature_file_prev, "../Features/");
	strcat(allFeature_file_prev, id);
	strcat(allFeature_file_prev, "/");
	strcat(allFeature_file_prev, id);
	strcat(allFeature_file_prev, ".all_dispredict_features");
	allp = fopen(allFeature_file_prev, "r");							// open annotated fasta file and read annotation
	if (allp == NULL) {
		fprintf(stderr, "Can't open previou all feature file!\n");
		exit(1);
	}
	//====================================================================================================================================
	//====================================================================================================================================
	// prepare new all feature file name and open to write
	strcpy(allFeature_file_new, "../Features/");
	strcat(allFeature_file_new, id);
	strcat(allFeature_file_new, "/");
	strcat(allFeature_file_new, id);
	strcat(allFeature_file_new, ".57pfeatures");
	alln = fopen(allFeature_file_new, "wb+");							// open annotated fasta file and read annotation
	if (alln == NULL) {
		fprintf(stderr, "Can't open new feature file file!\n");
		exit(1);
	}
	fprintf(alln, "O/D(1), AA(1), PP(7), PSSM(20), SS(3), ASA(1), dphi(1), dpsi(1), MG(1), BG(20), sPSEE(1), t(1)\n");
		
	//====================================================================================================================================
	//====================================================================================================================================
	// prepare sPSEE file name and open to read
	strcpy(sPSEE_filne_name, "../Features/");
	strcat(sPSEE_filne_name, id);
	strcat(sPSEE_filne_name, "/");
	strcat(sPSEE_filne_name, id);
	strcat(sPSEE_filne_name, ".PSEE");
	spsee = fopen(sPSEE_filne_name, "r");							// open annotated fasta file and read annotation
	if (spsee == NULL) {
		fprintf(nfp,"%s\n", id);
		fprintf(stderr, "Can't open energy file!\n");
		exit(1);
		//continue;
	}
	fgets(energy_line, sizeof energy_line, spsee);
	fgets(energy_line, sizeof energy_line, spsee);
	//====================================================================================================================================
		
	//====================================================================================================================================
	// Merger all features
	// Feature serial is: class, residue indication, pp, pssm, ss, sa, torsion angle, monogram, bigram, terminal
	res_pos = 0;
	featureCount = 0;
	i = 0;
	while (res_pos < seqLength)
	{
		//====================================================================================================================================
		//collect 56 features from previous feature file 
		fgets(feature_part_1, sizeof feature_part_1, allp);
		i = 0;
		char *tf;
		tf = strtok(feature_part_1, " ");
		while (tf != NULL)
		{
			i++;
			if (i == 1)
			{
				strcpy(all_feature, tf);
				strcat(all_feature, " ");
			}
			else if (i == 57)
			{
				strcpy(terminal, tf);
			}
			else if ((i>=2) && (i<=56))
			{
				strcat(all_feature, tf);
				strcat(all_feature, " ");
				featureCount++;
			}
				
			tf = strtok(NULL, " ");
		}
			
		//====================================================================================================================================
		//====================================================================================================================================
		//collect energy
		fgets(energy_line, sizeof energy_line, spsee);
		i = 0;
		char *te;
		te = strtok(energy_line, " ");
		while (te != NULL)
		{
			i++;
			if (i == 6)
			{
				strcpy(PSEE_val, te);
				l = strlen(PSEE_val);
				l--;
				PSEE_val[l] = '\0';
			}
				
			te = strtok(NULL, " ");
		}
		//====================================================================================================================================
		strcat(all_feature, PSEE_val);
		featureCount++;
		strcat(all_feature, " ");
		strcat(all_feature, terminal);
		strcat(all_feature, " ");
		featureCount++;
		//====================================================================================================================================
		fprintf(alln, "%s\n", all_feature);
		res_pos++;
	}
	if (res_pos < (seqLength - 1))
	{
		printf("Error!!\n");
		printf("Full Sequence Length NOT Covered!!\n");
		exit(1);
	}
	if (featureCount < MAX_FEATURE)
	{
		printf("Error!!\n");
		printf("57 features NOT Covered!!\n");
		exit(1);
	}
	rline[0] = '\0';
	//====================================================================================================================================
	fclose(allp);
	fclose(alln);
	fclose(spsee);
	seqLength = 0;
	fclose(nfp);
}