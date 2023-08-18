/*
Author: Sumaiya Iqbal
Prepare final input for DisPredict by applying windowing
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

char normalOuputFile[] = "../Output/log/log_windowing_final_input.txt";
char fileIdList[] = "../Input/id_list.txt";
char rline[MAX_LINE_SIZE];
char wline[MAX_LINE_SIZE];
char allFeature[MAX_LINE_SIZE];
char allFeatureFileName[400];
char libSVMFileName[400];
char aa[] = "ARNDCQEGHILKMFPSTWYV";
char cnt[200];
char temp[100];
char id[MAX_LINE_SIZE];

int lineLength = 0;
int seqNumber = 0;
int feature_count = 0;
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
	strcpy(id, argv[1]);
	list = fopen(fileIdList, "r");
	if (list == NULL) {
		fprintf(stderr, "Can't open SL list File!\n");
		exit(1);
	}

	nfp = fopen(normalOuputFile, "wb+");
	if (nfp == NULL) {
		fprintf(stderr, "Can't create normal output file!\n");
		exit(1);
	}

	//printf("\nFinal Input (including windowing) Preparation START...\n");
	
	fprintf(nfp, "%s\n", id);
	
	//----------------------------------------------------------------------------------------------------------------------------------
	// generate 21 window size libsvm file for this id
	
	//prepare libSVM file names
	strcpy(libSVMFileName, "../Features/");
	strcat(libSVMFileName, id);
	strcat(libSVMFileName, "/");
	strcat(libSVMFileName, id);
	strcat(libSVMFileName, ".final.dispredict.input");

	//printf("DisPredict Input file name: %s\n", libSVMFileName);

	// open libsvm file
	svm = fopen(libSVMFileName, "wb+");
	if (svm == NULL) {
		fprintf(stderr, "Can't open all feature file File!\n");
		exit(1);
	}
	
	//prepare all feature file names
	strcpy(allFeatureFileName, "../Features/");
	strcat(allFeatureFileName, id);
	strcat(allFeatureFileName, "/");
	strcat(allFeatureFileName, id);
	strcat(allFeatureFileName, ".all_dispredict_features");
	//printf("All feature file name: %s\n", allFeatureFileName);

	

	// open all feature file
	allf = fopen(allFeatureFileName, "r");
	if (allf == NULL) {
		fprintf(stderr, "Can't open all feature file File!\n");
		exit(1);
	}

	seqLength = 0;
	while (fgets(allFeature, sizeof allFeature, allf) != NULL)
	{
		seqLength++;
	}
	//printf("Sequence Length : %d\n", seqLength);
	fclose(allf);

	int win_split = WINDOW_SIZE / 2;
	int lineCount = 1;

	while (lineCount <= seqLength)
	{
		// take class
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

		fprintf(nfp, "line : %d, wline (only class) : %s\n", lineCount, wline);
			
			
		// take feature
		allf = fopen(allFeatureFileName, "r");
		if (allf == NULL) {
			fprintf(stderr, "Can't open all feature file File!\n");
			exit(1);
		}

		int start = lineCount - win_split;
		int end = lineCount + win_split;

		fprintf(nfp, "Start: %d--------End : %d\n", start, end);
		
		feature_count = 1;
		int window_cover = 0;
		int line_scan = 0;
		while (window_cover < WINDOW_SIZE)
		{
			if (start < 1)
			{
				int feature_track = 1;
				while (feature_track <= NO_OF_FEATURES)
				{
					temp[0] = '\0';
					sprintf(temp, "%d", feature_count);
					strcat(wline, temp);
					strcat(wline, ":0 ");
					feature_track++;
					feature_count++;
					
				}
				
				window_cover++;
				start++;
				
				fprintf(nfp, "Underflow\n");
			}
			else if (((start >= 1) || (end <= seqLength)) && (start <= seqLength))
			{
				//printf("Within sequence\n");
				
				fgets(allFeature, sizeof allFeature, allf);
				int l = strlen(allFeature);
				l--;
				allFeature[l] = '\0';
				line_scan++;
				//fprintf(nfp, "LINE Scan : %d\n", line_scan);
				if ((line_scan >= start) && (line_scan <= end))
				{
					fprintf(nfp, "Within range\n");
					char *t;
					t = strtok(allFeature, " ");
					int token_track = 0;
					while (t != NULL)
					{
						if (token_track > 0)
						{
							temp[0] = '\0';
							sprintf(temp, "%d", feature_count);
							strcat(wline, temp);
							strcat(wline, ":");
							strcat(wline, t);
							strcat(wline, " ");
							feature_count++;
						}
						
						t = strtok(NULL, " ");
						token_track++;
					}
					fprintf(nfp, "WLINE: %s\n", wline);
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
				while (feature_track <= NO_OF_FEATURES)
				{
					temp[0] = '\0';
					sprintf(temp, "%d", feature_count);
					strcat(wline, temp);
					strcat(wline, ":0 ");
					feature_track++;
					feature_count++;

				}

				window_cover++;
				end--;
				fprintf(nfp, "Overflow\n");
			}
			
		}
		fclose(allf);
		fprintf(svm, "%s\n", wline);
		wline[0] = '\0';
		allFeature[0] = '\0';

		
		lineCount++;
	}
	
	rline[0] = '\0';
	fclose(svm);
	//printf("%s--------DONE\n", rline);

	//printf("\nFinal Input (including windowing) Preparation END...\n");	
	fclose(list);
	fclose(nfp);
}