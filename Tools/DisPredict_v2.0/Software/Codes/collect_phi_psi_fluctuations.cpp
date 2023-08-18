/*
Author: Sumaiya Iqbal
Collect phi, psi fluctuation out from DAVAR output
Part of BalancedSSp software
*/

#define _CRT_SECURE_NO_DEPRECATE
#define LINE_SIZE 30000 //Maximum line size that can be read from the input file

#include<stdio.h>
#include<stdlib.h>
#include<string.h>
#include <malloc.h>
#include<ctype.h>
#include<math.h>

//FILE *i_listfp; // list file pointer
FILE *nfp;		// log file
FILE *dphi;
FILE *dpsi;
FILE *out1;
FILE *out2;
FILE *out3;
FILE *out4;
FILE *out5;
FILE *fasta;

//char inputFileList[] = "../Input/id_list.txt";
char normalOuputFile[] = "../Output/log/log_phi_psi_fluctuation_output_collection.txt";
char out1_file[500];
char out2_file[500];
char out3_file[500];
char out4_file[500];
char out5_file[500];
char dphi_file[] = "test.dphi";
char dpsi_file[] = "test.dpsi";
char rline[LINE_SIZE];
char wline[LINE_SIZE];
char r[5];
char out1line[LINE_SIZE];
char out2line[LINE_SIZE];
char out3line[LINE_SIZE];
char out4line[LINE_SIZE];
char out5line[LINE_SIZE];

int lineLength = 0;
int seqNumber = 0;
int seqLength = 0;

double outP = 0.0;
double outS = 0.0;
int run = 5;

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

int main()
{
	nfp = fopen(normalOuputFile, "wb+");
	if (nfp == NULL) {
		fprintf(stderr, "Can't create normal output file!\n");
		exit(1);
	}

	//printf("\nPhi/Psi Angle fluctuation prediction collection START.......\n");

	dphi = fopen(dphi_file, "wb+");							
	if (dphi == NULL) {
		fprintf(stderr, "Can't create test dphi file!\n");
		exit(1);
	}
	
	
	dpsi = fopen(dpsi_file, "wb+");
	if (dpsi == NULL) {
		fprintf(stderr, "Can't create test dpsi file!\n");
		exit(1);
	}
	
	
	//PHI
	//open 5 dphi out file to read
	strcpy(out1_file, "../AdditionalFiles/davar/dphi/round1/pred_test.out");
	out1 = fopen(out1_file, "r");							// open fasta file and read sequence
	if (out1 == NULL) {
		fprintf(stderr, "Can't open output1 file!\n");
		exit(1);
	}
	while(fgets(out1line, sizeof out1line, out1) != NULL)
	{
		seqLength++;
	}
	//printf("Length: %d\n", seqLength);
	fclose(out1);
	
	strcpy(out1_file, "../AdditionalFiles/davar/dphi/round1/pred_test.out");
	out1 = fopen(out1_file, "r");							// open fasta file and read sequence
	if (out1 == NULL) {
		fprintf(stderr, "Can't open output1 file!\n");
		exit(1);
	}
	
	strcpy(out2_file, "../AdditionalFiles/davar/dphi/round2/pred_test.out");
	out2 = fopen(out2_file, "r");							// open fasta file and read sequence
	if (out2 == NULL) {
		fprintf(stderr, "Can't open output2 file!\n");
		exit(1);
	}

	strcpy(out3_file, "../AdditionalFiles/davar/dphi/round3/pred_test.out");
	out3 = fopen(out3_file, "r");							// open fasta file and read sequence
	if (out3 == NULL) {
		fprintf(stderr, "Can't open output3 file!\n");
		exit(1);
	}

	strcpy(out4_file, "../AdditionalFiles/davar/dphi/round4/pred_test.out");
	out4 = fopen(out4_file, "r");							// open fasta file and read sequence
	if (out4 == NULL) {
		fprintf(stderr, "Can't open output1 file!\n");
		exit(1);
	}

	strcpy(out5_file, "../AdditionalFiles/davar/dphi/round5/pred_test.out");
	out5 = fopen(out5_file, "r");							// open fasta file and read sequence
	if (out5 == NULL) {
		fprintf(stderr, "Can't open output1 file!\n");
		exit(1);
	}
		
	int res_pos = 0;
	while (res_pos < seqLength)
	{
		outS = 0.0;
		fgets(out1line, sizeof out1line, out1);
		char * out1Sub = substring(out1line, 31, 12);
		outP = atof(out1Sub);
		outS = outS + outP;

		fgets(out2line, sizeof out2line, out2);
		char * out2Sub = substring(out2line, 31, 12);
		outP = atof(out2Sub);
		outS = outS + outP;

		fgets(out3line, sizeof out3line, out3);
		char * out3Sub = substring(out3line, 31, 12);
		outP = atof(out3Sub);
		outS = outS + outP;

		fgets(out4line, sizeof out4line, out4);
		char * out4Sub = substring(out4line, 31, 12);
		outP = atof(out4Sub);
		outS = outS + outP;

		fgets(out5line, sizeof out5line, out5);
		char * out5Sub = substring(out5line, 31, 12);
		outP = atof(out5Sub);
		outS = outS + outP;

		sprintf(wline, "%lf ", (double)outS/run);
		fprintf(dphi, "%s\n", wline);
		
		
		res_pos++;

		wline[0] = '\0';
		out1line[0] = '\0';
		out2line[0] = '\0';
		out3line[0] = '\0';
		out4line[0] = '\0';
		out5line[0] = '\0';
	}
	//=================================================================================

	fclose(out1);
	fclose(out2);
	fclose(out3);
	fclose(out4);
	fclose(out5);

	//PSI
	//open 5 dpsi out file to read
	strcpy(out1_file, "../AdditionalFiles/davar/dpsi/round1/pred_test.out");
	out1 = fopen(out1_file, "r");							// open fasta file and read sequence
	if (out1 == NULL) {
		fprintf(stderr, "Can't open output1 file!\n");
		exit(1);
	}

	strcpy(out2_file, "../AdditionalFiles/davar/dpsi/round2/pred_test.out");
	out2 = fopen(out2_file, "r");							// open fasta file and read sequence
	if (out2 == NULL) {
		fprintf(stderr, "Can't open output2 file!\n");
		exit(1);
	}

	strcpy(out3_file, "../AdditionalFiles/davar/dpsi/round3/pred_test.out");
	out3 = fopen(out3_file, "r");							// open fasta file and read sequence
	if (out3 == NULL) {
		fprintf(stderr, "Can't open output3 file!\n");
		exit(1);
	}

	strcpy(out4_file, "../AdditionalFiles/davar/dpsi/round4/pred_test.out");
	out4 = fopen(out4_file, "r");							// open fasta file and read sequence
	if (out4 == NULL) {
		fprintf(stderr, "Can't open output4 file!\n");
		exit(1);
	}

	strcpy(out5_file, "../AdditionalFiles/davar/dpsi/round5/pred_test.out");
	out5 = fopen(out5_file, "r");							// open fasta file and read sequence
	if (out5 == NULL) {
		fprintf(stderr, "Can't open output5 file!\n");
		exit(1);
	}

	outS = 0.0;
	res_pos = 0;
	while (res_pos < seqLength)
	{
		outS = 0.0;
		fgets(out1line, sizeof out1line, out1);
		char * out1Sub = substring(out1line, 31, 12);
		outP = atof(out1Sub);
		outS = outS + outP;

		fgets(out2line, sizeof out2line, out2);
		char * out2Sub = substring(out2line, 31, 12);
		outP = atof(out2Sub);
		outS = outS + outP;

		fgets(out3line, sizeof out3line, out3);
		char * out3Sub = substring(out3line, 31, 12);
		outP = atof(out3Sub);
		outS = outS + outP;

		fgets(out4line, sizeof out4line, out4);
		char * out4Sub = substring(out4line, 31, 12);
		outP = atof(out4Sub);
		outS = outS + outP;

		fgets(out5line, sizeof out5line, out5);
		char * out5Sub = substring(out5line, 31, 12);
		outP = atof(out5Sub);
		outS = outS + outP;

		sprintf(wline, "%lf ", (double)outS / run);
		fprintf(dpsi, "%s\n", wline);
		
		res_pos++;

		wline[0] = '\0';
		out1line[0] = '\0';
		out2line[0] = '\0';
		out3line[0] = '\0';
		out4line[0] = '\0';
		out5line[0] = '\0';
	}
	//=================================================================================

	fclose(out1);
	fclose(out2);
	fclose(out3);
	fclose(out4);
	fclose(out5);
	
	//printf("Phi/Psi Angle fluctuation prediction collection END.......\n");	
	
	fclose(nfp);
	fclose(dphi);
	fclose(dpsi);
}