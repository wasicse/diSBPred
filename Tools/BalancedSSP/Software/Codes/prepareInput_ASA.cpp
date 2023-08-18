/*
Author: Sumaiya Iqbal
Part of balancedSSP software
--- apply window of size 21 and prepare final input for ASA prediction by REGAd3p
--- This is input for exact method, so contains, bias + (55 * 21) columns
*/

#define _CRT_SECURE_NO_DEPRECATE
#define LINE_SIZE 30000 
#define WINDOW_SIZE 21
#define FEATURES 55
#define MAX_FEATURE_COUNT 2000

#include<stdio.h>
#include<stdlib.h>
#include<string.h>
#include <malloc.h>
#include<ctype.h>

FILE *nfp;
FILE *allf;
FILE *xfile;

char wline[LINE_SIZE];
char rline[LINE_SIZE];
char allFeature[LINE_SIZE];
char normalOuputFile[] = "../Output/log/log_windowing_input_ASA.txt";
char allFeatureFileName[400];
char x_FileName[400];
char directoryName[400];
char temp[100];
char id[100];

int seqNumber = 0;
int lineLength = 0;
int feature_count = 0;
int seqLength = 0;
int res_count = 0;

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
	
	//===============================================================================================================
	// collect id as parameter
	strcpy(id, argv[1]);
	//===============================================================================================================

	//===============================================================================================================
	nfp = fopen(normalOuputFile, "wb+");
	if (nfp == NULL) {
		fprintf(stderr, "Can't create normal output file!\n");
		exit(1);
	}
	//===============================================================================================================

	//===============================================================================================================
	// prepare input feature file name
	strcpy(x_FileName, "../Features/");
	strcat(x_FileName, id);
	strcat(x_FileName, "/");
	strcat(x_FileName, id);
	strcat(x_FileName, ".ASA.input");

	// open libsvm file
	xfile = fopen(x_FileName, "wb+");
	if (xfile == NULL) {
		fprintf(stderr, "Can't open libsvm output File!\n");
		exit(1);
	}
	//===============================================================================================================

	//===============================================================================================================
	//prepare all feature file names
	strcpy(allFeatureFileName, "../Features/");
	strcat(allFeatureFileName, id);
	strcat(allFeatureFileName, "/");
	strcat(allFeatureFileName, id);
	strcat(allFeatureFileName, ".ASA.features");

	// open all feature file
	allf = fopen(allFeatureFileName, "r");
	if (allf == NULL) {
		fprintf(stderr, "Can't open all feature file File!\n");
		exit(1);
	}
	
	// collect sequence length
	seqLength = 0;
	while (fgets(allFeature, sizeof allFeature, allf) != NULL)
	{
		seqLength++;
	}
	fclose(allf);
	//===============================================================================================================

	int win_split = WINDOW_SIZE / 2;
	int lineCount = 1;

	while (lineCount <= seqLength)
	{
		// take rsa
		allf = fopen(allFeatureFileName, "r");
		if (allf == NULL) {
			fprintf(stderr, "Can't open all feature file File!\n");
			exit(1);
		}
		int read_count = 1;
		while (fgets(allFeature, sizeof allFeature, allf) != NULL)
		{
			if (read_count == lineCount)
			{
				char *tf;
				tf = strtok(allFeature, " ");
				int token_track = 0;
				while (tf != NULL)
				{
					if (token_track == 0)
					{
						strcpy(wline, tf);
						strcat(wline, " ");
						break;
					}
					tf = strtok(NULL, " ");
					token_track++;
				}
			}
			read_count++;
		}

		fclose(allf);

		strcat(wline, "1 ");   // 0th coeeficient

		// take feature
		allf = fopen(allFeatureFileName, "r");
		if (allf == NULL) {
			fprintf(stderr, "Can't open all feature file File!\n");
			exit(1);
		}

		int start = lineCount - win_split;
		int end = lineCount + win_split;

		//fprintf(nfp, "Start: %d--------End : %d\n", start, end);

		feature_count = 1;
		int window_cover = 0;
		int line_scan = 0;
		while (window_cover < WINDOW_SIZE)
		{
			if (start < 1)
			{
				int feature_track = 1;
				while (feature_track <= FEATURES)
				{
					strcat(wline, "0 ");
					feature_track++;
					feature_count++;
				}

				window_cover++;
				start++;
			}
			else if (((start >= 1) || (end <= seqLength)) && (start <= seqLength))
			{
				fgets(allFeature, sizeof allFeature, allf);
				int l = strlen(allFeature);
				l--;
				allFeature[l] = '\0';
				line_scan++;

				if ((line_scan >= start) && (line_scan <= end))
				{
					char *t;
					t = strtok(allFeature, " ");
					int token_track = 0;
					while (t != NULL)
					{
						if (token_track > 0)
						{
							strcat(wline, t);
							strcat(wline, " ");
							feature_count++;
						}

						t = strtok(NULL, " ");
						token_track++;
					}
					//fprintf(nfp, "WLINE: %s\n", wline);
					window_cover++;
					start++;
				}
				else
				{
					//line_scan++;
					continue;
				}


			}
			else if (end > seqLength)
			{
				int feature_track = 1;
				while (feature_track <= FEATURES)
				{
					strcat(wline, "0 ");
					feature_track++;
					feature_count++;
				}

				window_cover++;
				end--;
				fprintf(nfp, "Overflow\n");
			}

		}
		fclose(allf);

		int spc = 0;
		while (1)
		{
			if (wline[spc] == ' ')
			{
				break;
			}
			spc++;
		}

		char *x_con = substring(wline, spc + 2, strlen(wline) - spc);
		fprintf(xfile, "%s\n", x_con);

		wline[0] = '\0';
		allFeature[0] = '\0';

		res_count++;
		lineCount++;
	}

	fclose(xfile);
	fclose(nfp);
}