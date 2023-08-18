/*
Author: Sumaiya Iqbal
Part of balancedSSP software
--- computes bigram and monogram from PSSM
*/
#define _CRT_SECURE_NO_DEPRECATE
#define MAX_LINE_SIZE 15000 //Maximum line size that can be read from the input file

#include<stdio.h>
#include<stdlib.h>
#include<string.h>
#include <malloc.h>
#include<ctype.h>

FILE *mvalues;
FILE *bvalues;
FILE *pssm;
FILE *pssm1;
FILE *monogram;
FILE *bigram;
FILE *nfp;


char normalOuputFile[] = "../Output/log/log_generate_mg_bg.txt";
char pssmFileName[300];
char monogramFileName[300];
char bigramFileName[300];
char wline[MAX_LINE_SIZE];
char aa[] = "ARNDCQEGHILKMFPSTWYV";
char pssmLine[MAX_LINE_SIZE];
char pssmLine1[MAX_LINE_SIZE];
char cnt[200];
char id[MAX_LINE_SIZE];

int pssmLineCount = 0;
int lineLength = 0;
int monogramVal[20];
int bigramVal[20][20];
int seqLength = 0;


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
	//==============================================================================================
	// input protein id
	strcpy(id, argv[1]);
	//==============================================================================================

	//==============================================================================================
	// log file for monogram,bigram generation
	nfp = fopen(normalOuputFile, "wb+");
	if (nfp == NULL) {
		fprintf(stderr, "Can't create log output file!\n");
		exit(1);
	}
	//==============================================================================================

	
	//==============================================================================================
	// Monogram

	for (int i = 0; i < 20; i++)
	{
		monogramVal[i] = 0;
	}
	for (int i = 0; i < 20; i++)
	{
		for (int j = 0; j < 20; j++)
		{
			bigramVal[i][j] = 0;
		}

	}

	// specify PSSM file directory
	strcpy(pssmFileName, "../Features/");
	strcat(pssmFileName, id);
	strcat(pssmFileName, "/");
	strcat(pssmFileName, id);
	strcat(pssmFileName, ".pssm");
	fprintf(nfp, "%s\n", pssmFileName);
	//printf("%s\n", pssmFileName);

	pssm = fopen(pssmFileName, "r");
	if (pssm == NULL) {
		fprintf(stderr, "Can't open pssm file!\n");
		exit(1);
	}

	// Create Monogram file name
	strcpy(monogramFileName, "../Features/");
	strcat(monogramFileName, id);
	strcat(monogramFileName, "/");
	strcat(monogramFileName, id);
	strcat(monogramFileName, ".monogram");
	fprintf(nfp, "%s\n", monogramFileName);


	// Generate Monogram

	pssmLineCount = 0;
	while (!feof(pssm))
	{
		fgets(pssmLine, sizeof pssmLine, pssm);
		pssmLineCount++;
		if (pssmLineCount < 4)
		{
			continue;
		}
		if (strlen(pssmLine) > 60)
		{
			seqLength++;
			char *pssmSub = substring(pssmLine, 10, 60);
			//fprintf(nfp, "%s\n", pssmSub);

			for (int j = 0; j < 20; j++)
			{
				char* t;
				t = substring(pssmSub, (j * 3) + 1, 3);

				char* token;
				token = trim(t);
				int val = atoi(token);

				monogramVal[j] = monogramVal[j] + val;


				free(t);
			}
		}
		pssmLine[0] = '\0';
	}
	fclose(pssm);

	//Write on Monogram File
	monogram = fopen(monogramFileName, "wb+");
	fprintf(monogram, "A, R, N, D, C, Q, E, G, H, I, L, K, M, F, P, S, T, W, Y, V \n");
	wline[0] = '\0';
	char temp[10];
	for (int k = 0; k < 20; k++)
	{
		sprintf(temp, "%d", monogramVal[k]);
		strcat(wline, temp);
		strcat(wline, ",");
	}
	wline[strlen(wline) - 1] = '\0';
	fprintf(monogram, "%s\n", wline);
	fclose(monogram);

	//printf("Monogram Generated!!\n");
	//=======================================================================================================//


	//==============================================================================================
	// Bigram
	//Create Bigram File name
	strcpy(bigramFileName, "../Features/");
	strcat(bigramFileName, id);
	strcat(bigramFileName, "/");
	strcat(bigramFileName, id);
	strcat(bigramFileName, ".bigram");
	fprintf(nfp, "%s\n", bigramFileName);

	fprintf(nfp, "SEQUENCE LENGTH: %d\n", seqLength);
	// Generate Bigram
	for (int i = 0; i < 20; i++)										// Iterate for 20 possible amino acid									
	{
		for (int j = 0; j < 20; j++)									// Iterate for 20 possible amino acids (form pair)				
		{

			for (int k = 0; k < seqLength - 1; k++)						// Iterate for length of the sequence
			{
				pssm = fopen(pssmFileName, "r");						// Open PSSM file for fisrt scan
				if (pssm == NULL) {
					fprintf(stderr, "Can't open pssm file!\n");
					continue;
				}


				pssmLineCount = 0;										// Start counting line in PSSM file - 1st scan
				while (!feof(pssm))
				{
					pssmLineCount++;
					fgets(pssmLine, sizeof pssmLine, pssm);
					if (pssmLineCount == k + 4)							// Stop at 1st line ---- increment through loop
					{
						break;
					}
				}
				fclose(pssm);

				pssm1 = fopen(pssmFileName, "r");						// Open PSSM file for second scan
				pssmLineCount = 0;									// Start counting line in PSSM file - 1st scan	
				while (!feof(pssm1))
				{
					pssmLineCount++;
					fgets(pssmLine1, sizeof pssmLine1, pssm1);
					if (pssmLineCount == k + 5)						// Stop at 2nd line ---- increment through loop	
					{
						break;
					}
				}
				fclose(pssm1);

				char *pssmSub = substring(pssmLine, 10, 60);			// Parse PSSM information line
				char *pssmSub1 = substring(pssmLine1, 10, 60);

				char* t;
				t = substring(pssmSub, (i * 3) + 1, 3);					// Parse PSSM value
				char *token;
				token = trim(t);
				int val = atoi(token);

				char* t1;
				t1 = substring(pssmSub1, (j * 3) + 1, 3);				// Parse PSSM value
				char *token1;
				token1 = trim(t1);
				int val1 = atoi(token1);

				bigramVal[i][j] = bigramVal[i][j] + (val * val1);		// calculate bigram


				free(t);
				free(t1);

				pssmLine[0] = '\0';
				pssmLine1[0] = '\0';


			}
		}
	}


	//Write on Bigram File
	bigram = fopen(bigramFileName, "wb+");
	fprintf(bigram, "A, R, N, D, C, Q, E, G, H, I, L, K, M, F, P, S, T, W, Y, V \n");

	for (int m = 0; m < 20; m++)
	{
		for (int n = 0; n < 20; n++)
		{
			fprintf(bigram, "%d  ", bigramVal[m][n]);
		}
		fprintf(bigram, "\n");
	}

	fclose(bigram);

	//==============================================================================================

	fprintf(nfp, "Sequence Length: %d\n", seqLength);
	seqLength = 0;

	fclose(nfp);
}