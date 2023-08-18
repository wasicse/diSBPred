/*
Author: Sumaiya Iqbal
Part of balancedSSP sofware
--- Predict ASA by exact method with loading weights
--- Process formatted output
*/

#define _CRT_SECURE_NO_DEPRECATE
#define LINE_SIZE 40000					//Maximum line size that can be read from the input file
#define FEATURES 3466
#define KERNEL 3

#include<stdio.h>
#include<stdlib.h>
#include<string.h>
#include <malloc.h>
#include<ctype.h>

FILE *nfp;
FILE *X;
FILE *wgt;
FILE *fasta;
FILE *asap;

char normalOuputFile[] = "../Output/log/asa_prediction.txt";
char X_filename[200];
char weight_filename[] = "../Models/ASA_WEIGHT/weight.txt";
char fastaFileName[200];
char ASApFileName[200];

char id[10];
char fsequence[LINE_SIZE];
char featureLine[LINE_SIZE];
char wghtLine[100];
char ASAp_str[500];
char wline[500];

double featureVal = 0.0;
double wghtVal = 0.0;
int seqLength = 0;
int resPos = 0;
double ASAp_val = 0.0;
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
	// collect id 
	strcpy(id, argv[1]);
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
	// prepare feature file name, read
	strcpy(X_filename, "../Features/");
	strcat(X_filename, id);
	strcat(X_filename, "/");
	strcat(X_filename, id);
	strcat(X_filename, ".ASA.input");

	X = fopen(X_filename, "r");
	if (X == NULL) {
		fprintf(stderr, "Can't open ASA input File!\n");
		exit(1);
	}

	//==========================================================================================================================

	//==========================================================================================================================
	// prepare predicted ASA file name and write
	strcpy(ASApFileName, "../Output/prediction/");
	strcat(ASApFileName, id);
	strcat(ASApFileName, "/ASA/");
	strcat(ASApFileName, id);
	strcat(ASApFileName, ".ASAp");

	asap = fopen(ASApFileName, "wb+");
	if (asap == NULL) {
		fprintf(stderr, "Can't open feature File!\n");
		exit(1);
	}
	fprintf(asap, "Absolute Accessible Surface Area (ASA) Prediction Output by REGAd3p\n");
	fprintf(asap, "TARGET: %s, ", id);
	fprintf(asap, "Length: %d\n", seqLength);
	fprintf(asap, "\n");
	fprintf(asap, "SR#  AA  ASAp\n");
	//==========================================================================================================================

	//==========================================================================================================================
	// extend kernel to degree 3 polynomial with feature set internally and multiply with weight to get predicted ASA
	resPos = 0;
	while (resPos < seqLength)
	{
		ASAp_val = 0;
		fgets(featureLine, sizeof featureLine, X);						// collect feature
		
		wgt = fopen(weight_filename, "r");								// open weights for reading
		if (wgt == NULL) {
			fprintf(stderr, "Can't feature extended File!\n");
			exit(1);
		}

		char *tf;
		tf = strtok(featureLine, " ");
		int token_track = 1;
		while (tf != NULL)
		{
			if (token_track == 1)
			{
				featureVal = atof(tf);
				fgets(wghtLine, sizeof wghtLine, wgt);						// collect feature
				wghtVal = atof(wghtLine);

				ASAp_val = ASAp_val + (featureVal * wghtVal);
			}
			if ((token_track >= 2) && (token_track <= FEATURES))
			{
				featureVal = atof(tf);
			
				fgets(wghtLine, sizeof wghtLine, wgt);									// collect wight
				wghtVal = atof(wghtLine);
				ASAp_val = ASAp_val + (featureVal * wghtVal);							// 1st order	

				fgets(wghtLine, sizeof wghtLine, wgt);									// collect wight
				wghtVal = atof(wghtLine);
				ASAp_val = ASAp_val + (featureVal * featureVal * wghtVal);				// 2nd order				

				fgets(wghtLine, sizeof wghtLine, wgt);									// collect wight
				wghtVal = atof(wghtLine);
				ASAp_val = ASAp_val + (featureVal * featureVal * featureVal *wghtVal);	// 3rd order	
			}	
			
			tf = strtok(NULL, " ");
			token_track++;
		}
		
		fclose(wgt);

		//==========================================================================================================================
		// write serial number
		sprintf(wline, "%d", resPos + 1);
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
		char* fresidue1 = substring(fsequence, resPos + 1, 1);		// residues	
		strcat(wline, fresidue1);
		strcat(wline, "   ");
		resPos++;
		//==========================================================================================================================
		
		// write predicted ASA
		if (ASAp_val < 0.0)
		{
			ASAp_val = 0.0;
		}
		sprintf(ASAp_str, "%lf", ASAp_val);
		strcat(wline, ASAp_str);
		strcat(wline, "   ");

		fprintf(asap, "%s\n", trim(wline));				// write predicted asa
	}
	




	//==========================================================================================================================

	fclose(nfp);
	fclose(asap);
	fclose(X);
}








