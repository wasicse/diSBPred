/*
Component of balancedSSP
Author: Sumaiya Iqbal
Process Input:
--- Create directory for each protein ID
--- Copy fasta into that directory
*/

#define _CRT_SECURE_NO_DEPRECATE

#include<stdio.h>
#include<stdlib.h>
#include<string.h>
#include<ctype.h>

#define MAX_LINE_SIZE 30000

FILE *nfp;
FILE *fasta;

char normalOuputFile[] = "../Output/log/log_processInput.txt";
char make_directory_command[MAX_LINE_SIZE];
char cp_fasta_file_command[MAX_LINE_SIZE];
char id[MAX_LINE_SIZE];
char fastaFileName[MAX_LINE_SIZE];
char fsequence[MAX_LINE_SIZE];

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
	//===================================================================================================
	// open log file 
	nfp = fopen(normalOuputFile, "wb+");
	if (nfp == NULL) {
		fprintf(stderr, "Can't create normal output file!\n");
		exit(1);
	}
	//===================================================================================================

	//===================================================================================================
	// take input protein id
	strcpy(id, argv[1]);
	//===================================================================================================

	//===================================================================================================
	// make directory for this protein id
	
	strcpy(make_directory_command, "mkdir ../Features/");
	strcat(make_directory_command, id);
	fprintf(nfp, "Make Directory Command: %s\n", make_directory_command);
	system(make_directory_command);

	//===================================================================================================
	
	//===================================================================================================
	// copy fasta file into the directory
	strcpy(cp_fasta_file_command, "cp ../Input/FASTA/");
	strcat(cp_fasta_file_command, id);
	strcat(cp_fasta_file_command, ".fasta");
	strcat(cp_fasta_file_command, " ../Features/");
	strcat(cp_fasta_file_command, id);
	fprintf(nfp, "Copy fasta file command: %s\n", cp_fasta_file_command);
	system(cp_fasta_file_command);
	
	//===================================================================================================

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

	printf("Length: %d\n", seqLength);

	fclose(nfp);
}