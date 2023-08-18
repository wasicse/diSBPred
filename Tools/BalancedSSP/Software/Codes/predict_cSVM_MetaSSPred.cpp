/*
Author: Sumaiya Iqbal
Part of balancedSSP software
---- predict secondary structure by cSVM (conbining GA weight) + MetaSSPred
---- extension of output (.cSVM.ssp/ .MetaSSPred.ssp)
*/
#define _CRT_SECURE_NO_DEPRECATE
#define MAX_LINE_SIZE 30000 //Maximum line size that can be read from the input file
#define BETA 1															// window size
#define COIL 2
#define HELIX 3

#include<stdio.h>
#include<stdlib.h>
#include<string.h>
#include <malloc.h>
#include<ctype.h>


FILE *nfp;
FILE *fasta;
FILE *betap;
FILE *coilp;
FILE *helixp;
FILE *weightF;
FILE *spx;
FILE *cSVM;
FILE *meta;

char fsequence[MAX_LINE_SIZE];
char fastaFileName[MAX_LINE_SIZE];

char binary_beta_file_name[MAX_LINE_SIZE];
char binary_coil_file_name[MAX_LINE_SIZE];
char binary_helix_file_name[MAX_LINE_SIZE];
char cSVM_output_file_name[MAX_LINE_SIZE];
char meta_output_file_name[MAX_LINE_SIZE];
char spinex_output_file_name[MAX_LINE_SIZE];
char weightFileName[] = "../Models/cSVM/bestparameter_f33.471";
char output_line_beta[MAX_LINE_SIZE];
char output_line_coil[MAX_LINE_SIZE];
char output_line_helix[MAX_LINE_SIZE];
char output_line_csvm[MAX_LINE_SIZE];
char output_line_spx[MAX_LINE_SIZE];
char output_line_meta[MAX_LINE_SIZE];

char normalOuputFile[] = "../Output/log/cSVM_meta_svm_prediction.txt";

char rline[MAX_LINE_SIZE];
char wline[MAX_LINE_SIZE];
char wline1[MAX_LINE_SIZE];
char id[MAX_LINE_SIZE];

int lineLength = 0;
int seqNumber = 0;
int seqLength = 0;
int resPos = 0;
double wBeta = 0.0;
double wCoil = 0.0;
double wHelix = 0.0;
double sump = 0.0;

int cSVMClass = 0;
double cSVMProb[3];
int metaClass = 0;
double metaProb[3];

double pE = 0.0;
double pC = 0.0;
double pH = 0.0;
double wpE = 0.0;
double wpC = 0.0;
double wpH = 0.0;
int spc = 0;
int constl = 0;
char probability[50];

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

	fclose(fasta);
	//==========================================================================================================================


	//==========================================================================================================================
	// read wight file and store weights
	weightF = fopen(weightFileName, "r");
	if (weightF == NULL) {
		fprintf(stderr, "Can't open fasta File!\n");
		exit(1);
	}
	fgets(rline, sizeof rline, weightF);		// weight for beta
	wBeta = atof(trim(rline));
	
	fgets(rline, sizeof rline, weightF);		// weight for coil
	wCoil = atof(trim(rline));

	fgets(rline, sizeof rline, weightF);		// weight for helix
	wHelix = atof(trim(rline));
	fclose(weightF);
	//==========================================================================================================================

	//==========================================================================================================================
	// binary beta out file
	strcpy(binary_beta_file_name, "../Features/");
	strcat(binary_beta_file_name, id);
	strcat(binary_beta_file_name, "/");
	strcat(binary_beta_file_name, id);
	strcat(binary_beta_file_name, ".binary.beta.predict");

	betap = fopen(binary_beta_file_name, "r");
	if (betap == NULL) {
		fprintf(stderr, "Can't open beta output File!\n");
		exit(1);
	}
	fgets(output_line_beta, sizeof output_line_beta, betap);	// skip header

	//==========================================================================================================================

	//==========================================================================================================================
	// binary coil out file
	strcpy(binary_coil_file_name, "../Features/");
	strcat(binary_coil_file_name, id);
	strcat(binary_coil_file_name, "/");
	strcat(binary_coil_file_name, id);
	strcat(binary_coil_file_name, ".binary.coil.predict");

	coilp = fopen(binary_coil_file_name, "r");
	if (coilp == NULL) {
		fprintf(stderr, "Can't open coil output File!\n");
		exit(1);
	}
	fgets(output_line_coil, sizeof output_line_coil, coilp);	// skip header
	//==========================================================================================================================

	//==========================================================================================================================
	// binary helix out file
	strcpy(binary_helix_file_name, "../Features/");
	strcat(binary_helix_file_name, id);
	strcat(binary_helix_file_name, "/");
	strcat(binary_helix_file_name, id);
	strcat(binary_helix_file_name, ".binary.helix.predict");

	helixp = fopen(binary_helix_file_name, "r");
	if (helixp == NULL) {
		fprintf(stderr, "Can't open helix output File!\n");
		exit(1);
	}
	fgets(output_line_helix, sizeof output_line_helix, helixp);	// skip header
	//==========================================================================================================================

	//==========================================================================================================================
	// spinex file
	strcpy(spinex_output_file_name, "../Features/");
	strcat(spinex_output_file_name, id);
	strcat(spinex_output_file_name, "/");
	strcat(spinex_output_file_name, id);
	strcat(spinex_output_file_name, ".spXout");

	spx = fopen(spinex_output_file_name, "r");
	if (spx == NULL) {
		fprintf(stderr, "Can't open spinex output File!\n");
		exit(1);
	}
	fgets(output_line_spx, sizeof output_line_spx, spx);	// skip header
	//==========================================================================================================================

	//==========================================================================================================================
	// cSVM file
	strcpy(cSVM_output_file_name, "../Output/prediction/");
	strcat(cSVM_output_file_name, id);
	strcat(cSVM_output_file_name, "/cSVM/");
	strcat(cSVM_output_file_name, id);
	strcat(cSVM_output_file_name, ".cSVMpred");

	cSVM = fopen(cSVM_output_file_name, "wb+");
	if (cSVM == NULL) {
		fprintf(stderr, "Can't open csvm output File!\n");
		exit(1);
	}

	fprintf(cSVM, "Secondary Structure Prediction Output by cSVM (Part of BalancedSSP)\n");
	fprintf(cSVM, "TARGET: %s, ", id);
	fprintf(cSVM, "Length: %d\n", seqLength);
	fprintf(cSVM, "\n");
	fprintf(cSVM, "SR#  AA  cSVM  pBeta  pCoil  pHelix\n");
	//==========================================================================================================================

	//==========================================================================================================================
	// cSVM file
	strcpy(meta_output_file_name, "../Output/prediction/");
	strcat(meta_output_file_name, id);
	strcat(meta_output_file_name, "/MetaSSPred/");
	strcat(meta_output_file_name, id);
	strcat(meta_output_file_name, ".MetaSSpred");

	meta = fopen(meta_output_file_name, "wb+");
	if (meta == NULL) {
		fprintf(stderr, "Can't open metasspred output File!\n");
		exit(1);
	}

	fprintf(meta, "Secondary Structure Prediction Output by MetaSSPred (Part of BalancedSSP)\n");
	fprintf(meta, "TARGET: %s, ", id);
	fprintf(meta, "Length: %d\n", seqLength);
	fprintf(meta, "\n");
	fprintf(meta, "SR#  AA  MetaSSPred  pBeta  pCoil  pHelix\n");
	//==========================================================================================================================
	resPos = 0;
	while (resPos < seqLength)
	{
		
		fgets(output_line_beta, sizeof output_line_beta, betap);		// read prediction
		fgets(output_line_coil, sizeof output_line_coil, coilp);
		fgets(output_line_helix, sizeof output_line_helix, helixp);	
		
		
		char *tbp = strtok(output_line_beta, " ");		// collect beta probability		
		int token_track = 0;						// tokenize and cnvert to float
		while (tbp != NULL)
		{
			if (token_track == 1)
			{
				pE = atof(tbp);
			}
			
			tbp = strtok(NULL, " ");
			token_track++;
		}

		char *tcp = strtok(output_line_coil, " ");		// collect coil probability
		token_track = 0;						// tokenize and cnvert to float
		while (tcp != NULL)
		{
			if (token_track == 1)
			{
				pC = atof(tcp);
			}
		
			tcp = strtok(NULL, " ");			
			token_track++;
		}

		char *thp = strtok(output_line_helix, " ");		// collect helix probability
		token_track = 0;						// tokenize and cnvert to float
		while (thp != NULL)
		{
			if (token_track == 1)
			{
				pH = atof(thp);
			}
			thp = strtok(NULL, " ");
			token_track++;
		}

		wpE = pE + wBeta;							// assign weight
		wpC = pC + wCoil;
		wpH = pH + wHelix;
		sump = wpE + wpC + wpH;
		cSVMProb[0] = wpE / sump;
		cSVMProb[1] = wpC / sump;
		cSVMProb[2] = wpH / sump;

		//printf("done\n");

		if ((wpE >= wpC) && (wpE >= wpH))			// take maximum probability as result
		{
			cSVMClass = BETA;
		}
		else if ((wpC >= wpE) && (wpC >= wpH))
		{
			cSVMClass = COIL;
		}
		else if ((wpH >= wpE) && (wpH >= wpC))
		{
			cSVMClass = HELIX;
		}
		
		fgets(output_line_spx, sizeof output_line_spx, spx);
		char *spxout = substring(output_line_spx, 18, 1);
		
		//==========================================================================================================================
		// write serial number
		sprintf(wline, "%d", resPos + 1);
		sprintf(wline1, "%d", resPos + 1);
		spc = 0;
		constl = strlen(wline);
		while (spc < 5 - constl)
		{
			strcat(wline, " ");
			strcat(wline1, " ");
			spc++;
		}
		//==========================================================================================================================

		//==========================================================================================================================
		// write amino acid
		char* fresidue1 = substring(fsequence, resPos + 1, 1);		// residues	
		strcat(wline, fresidue1);
		strcat(wline1, fresidue1);
		strcat(wline, "   ");
		strcat(wline1, "   ");
		resPos++;
		//==========================================================================================================================

		//==========================================================================================================================
		if (cSVMClass == BETA)						// Write predicted annotation
		{
			strcat(wline, "E");
			strcat(wline1, "E");
		}
		else if (cSVMClass == COIL)
		{
			strcat(wline, "C");
			strcat(wline1, spxout);
		}
		else if (cSVMClass == HELIX)
		{
			strcat(wline, "H");
			strcat(wline1, spxout);
		}
		strcat(wline, "     ");
		strcat(wline1, "           ");
		//==========================================================================================================================

		//==========================================================================================================================
		sprintf(probability, "%1.3lf", cSVMProb[0]);	// cSVM probabilities
		strcat(wline, probability);
		strcat(wline, "  ");

		sprintf(probability, "%1.3lf", cSVMProb[1]);
		strcat(wline, probability);
		strcat(wline, "  ");

		sprintf(probability, "%1.3lf", cSVMProb[2]);
		strcat(wline, probability);
		strcat(wline, "  ");

		//==========================================================================================================================

		//==========================================================================================================================
		if (cSVMClass == BETA)
		{
			sprintf(probability, "%1.3lf", cSVMProb[0]);	// meta probabilities
			strcat(wline1, probability);
			strcat(wline1, "  ");

			sprintf(probability, "%1.3lf", cSVMProb[1]);
			strcat(wline1, probability);
			strcat(wline1, "  ");

			sprintf(probability, "%1.3lf", cSVMProb[2]);
			strcat(wline1, probability);
			strcat(wline1, "  ");
		}
		else
		{
			char * probE = substring(output_line_spx, 39, 5);
			strcat(wline1, probE);
			strcat(wline1, "  ");

			char * probC = substring(output_line_spx, 48, 5);
			strcat(wline1, probC);
			strcat(wline1, "  ");

			char * probH = substring(output_line_spx, 58, 5);
			strcat(wline1, probH);
			strcat(wline1, "  ");
		}

		//==========================================================================================================================
		fprintf(cSVM, "%s\n", trim(wline));				// write prediction of cSVM
		fprintf(meta, "%s\n", trim(wline1));			// write prediction of meta
	}

	fclose(betap);
	fclose(coilp);
	fclose(helixp);
	fclose(spx);
	fclose(cSVM);
	fclose(meta);
	fclose(nfp);
	
}