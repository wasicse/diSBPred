/*
Author: Sumaiya Iqbal
Part of balancedSSp software
Collect 52 features for initial SS prediction by REGAd3p
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
FILE *all;
FILE *dr;
FILE *bi;
FILE *mono;
FILE *fasta;
FILE *pssm;
FILE *pp;		// physical properties file
FILE *ius;
FILE *iul;


char normalOuputFile[] = "../Output/log/log_feature_collection_for_SS.txt";
char ppFile[] = "../AdditionalFiles/physiochemical_properties.txt";
char id[LINE_SIZE];
char make_directory_command[LINE_SIZE];
char gather_fasta_seq_command[LINE_SIZE];
char gather_annotated_fasta_seq_command[LINE_SIZE];
char gather_pssm_command[LINE_SIZE];
char gather_33feature_command[LINE_SIZE];
char gather_monogram_command[LINE_SIZE];
char gather_bigram_command[LINE_SIZE];
char iupredSFile[LINE_SIZE];
char iupredLFile[LINE_SIZE];

char iusline[LINE_SIZE];
char iulline[LINE_SIZE];
char fseq[LINE_SIZE];
char all_feature[LINE_SIZE];
char feature_part_1[LINE_SIZE];
char ppline[LINE_SIZE];
char m_value[30];
char b_value[30];
char residue[4];
char ssc[4];
char cnt[200];
char m_list[LINE_SIZE];
char b_list[LINE_SIZE];
double monogram[20];
double pssm_values_20[20];
double bigram[20][20];
char aa[] = "ARNDCQEGHILKMFPSTWYV";
double pp_val[20][7];
double pssm_value = 0.0;
char pp_value[30];

char pssm_val_str[30];
char all_feature_filename[300];
char fastaFileName[300];
char annotated_fasta_filename[300];
char dr_feature_filename[300];
char monogram_filename[300];
char bigram_filename[300];
char pssmFile[300];
char pssmline[1000];

int aa_counter = 0;
int pp_counter = 0;
int lineLength = 0;
int seqNumber = 0;
int seqLength = 0;
int seqNum = 0;
int local_counter = 0;
int pssm_counter_new = 0;

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
	// prepare pp file name, open to read, collect pp

	pp = fopen(ppFile, "r");													// open pp file
	if (pp == NULL) {
		fprintf(stderr, "Can't open pp file!\n");
		exit(1);
	}

	aa_counter = 0;
	while (aa_counter < 20)														// read for 20 amino acids
	{
		fgets(ppline, sizeof ppline, pp);
		pp_counter = 0;															// collect 7 properties
		char *t_pp;
		t_pp = strtok(ppline, " ");
		while (t_pp != NULL)
		{
			if (pp_counter > 0)													// skip first column (AA)
			{
				pp_val[aa_counter][pp_counter - 1] = (double)atof(t_pp);		// Collect each pp
			}
			pp_counter++;
			t_pp = strtok(NULL, " ");
		}
		aa_counter++;
	}

	ppline[0] = '\0';

	fclose(pp);
	//====================================================================================================================================

	//====================================================================================================================================
	// prepare fasta file name, open and collect sequence length
	fastaFileName[0] = '\0';
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
	// prepare pssm file name, open to read
	strcpy(pssmFile, "../Features/");
	strcat(pssmFile, id);
	strcat(pssmFile, "/");
	strcat(pssmFile, id);
	strcat(pssmFile, ".pssm");

	pssm = fopen(pssmFile, "r");							// open fasta file and read sequence
	if (pssm == NULL) {
		fprintf(stderr, "Can't open pssm file!\n");
		exit(1);
	}

	fgets(pssmline, sizeof pssmline, pssm);				// skip header
	fgets(pssmline, sizeof pssmline, pssm);				// skip header
	fgets(pssmline, sizeof pssmline, pssm);				// skip header
	//====================================================================================================================================

	//====================================================================================================================================
	// prepare monogram fasta file name
	strcpy(monogram_filename, "../Features/");
	strcat(monogram_filename, id);
	strcat(monogram_filename, "/");
	strcat(monogram_filename, id);
	strcat(monogram_filename, ".monogram");


	mono = fopen(monogram_filename, "r");								// open monogram file
	if (mono == NULL) {
		fprintf(stderr, "Can't open monogram file!\n");
		exit(1);
	}

	local_counter = 0;
	while (!feof(mono))
	{
		local_counter++;
		fgets(m_list, sizeof m_list, mono);								// collect monogram list
		if (local_counter == 2)
		{
			break;
		}
	}

	int i = 0;
	char *t;
	t = strtok(m_list, ",");
	while (t != NULL)
	{
		monogram[i] = (double)atoi(t) / exp(BIGRAM_MONOGRAM_LOG_MEDIAN);											// Collect each monogram
		i++;
		t = strtok(NULL, ",");
	}

	m_list[0] = '\0';

	fclose(mono);
	//====================================================================================================================================


	//====================================================================================================================================
	// prepare bigram file name
	strcpy(bigram_filename, "../Features/");
	strcat(bigram_filename, id);
	strcat(bigram_filename, "/");
	strcat(bigram_filename, id);
	strcat(bigram_filename, ".bigram");

	bi = fopen(bigram_filename, "r");											// open bigram file
	if (bi == NULL) {
		fprintf(stderr, "Can't open bigram file!\n");
		exit(1);
	}

	local_counter = 0;
	while (!feof(bi))
	{
		local_counter++;
		fgets(b_list, sizeof b_list, bi);										// collect bigram list
		if ((local_counter > 1) && (local_counter <= 21))
		{
			int i = 0;
			char *tb;
			tb = strtok(b_list, " ");
			while (tb != NULL)
			{
				bigram[local_counter - 2][i] = (double)atoi(tb) / exp(BIGRAM_MONOGRAM_LOG_MEDIAN);						// Collect each bigram
				i++;
				tb = strtok(NULL, " ");
			}
		}
		b_list[0] = '\0';
	}

	fclose(bi);
	
	//====================================================================================================================================
	// prepare iupred short output file, open to read
	strcpy(iupredSFile, "../Features/");
	strcat(iupredSFile, id);
	strcat(iupredSFile, "/");
	strcat(iupredSFile, id);
	strcat(iupredSFile, ".iupredS");

	ius = fopen(iupredSFile, "r");							// open iupreds output file and read sequence
	if (ius == NULL) {
		fprintf(stderr, "Can't open iupred short file!\n");
		exit(1);
	}

	fgets(iusline, sizeof iusline, ius);				// skip header
	fgets(iusline, sizeof iusline, ius);
	fgets(iusline, sizeof iusline, ius);
	fgets(iusline, sizeof iusline, ius);
	fgets(iusline, sizeof iusline, ius);
	fgets(iusline, sizeof iusline, ius);
	fgets(iusline, sizeof iusline, ius);
	fgets(iusline, sizeof iusline, ius);
	fgets(iusline, sizeof iusline, ius);
	//====================================================================================================================================

	//====================================================================================================================================
	// prepare iupred long output file, open to read
	strcpy(iupredLFile, "../Features/");
	strcat(iupredLFile, id);
	strcat(iupredLFile, "/");
	strcat(iupredLFile, id);
	strcat(iupredLFile, ".iupredL");

	iul = fopen(iupredLFile, "r");							// open iupreds output file and read sequence
	if (iul == NULL) {
		fprintf(stderr, "Can't open iupred long file!\n");
		exit(1);
	}

	fgets(iulline, sizeof iulline, iul);				// skip header
	fgets(iulline, sizeof iulline, iul);
	fgets(iulline, sizeof iulline, iul);
	fgets(iulline, sizeof iulline, iul);
	fgets(iulline, sizeof iulline, iul);
	fgets(iulline, sizeof iulline, iul);
	fgets(iulline, sizeof iulline, iul);
	fgets(iulline, sizeof iulline, iul);
	fgets(iulline, sizeof iulline, iul);
	//====================================================================================================================================

	//====================================================================================================================================
	// prepare all feature file name and write on that file
	strcpy(all_feature_filename, "../Features/");
	strcat(all_feature_filename, id);
	strcat(all_feature_filename, "/");
	strcat(all_feature_filename, id);
	strcat(all_feature_filename, ".initialSS.features");

	all = fopen(all_feature_filename, "wb+");					// Open all feature file for writing
	if (all == NULL) {
		fprintf(stderr, "Can't open all feature file!\n");
		exit(1);
	}

	//====================================================================================================================================
	// Merger all features
	// Feature serial is: class, residue indication, pp, pssm, ss, sa, torsion angle, monogram, bigram, terminal
	int res_pos = 0;
	while (res_pos < seqLength)
	{

		all_feature[0] = '\0';
		//====================================================================================================================================
		// Take class (annotation) --- dummy 1
		strcpy(all_feature, "1 ");
		//====================================================================================================================================

		//====================================================================================================================================

		char res = fseq[res_pos];
		int j = 0;
		for (; j < 20; j++)
		{
			if (res == aa[j])
			{
				break;
			}
		}
		sprintf(residue, "%d", j + 1);
		strcat(all_feature, residue);
		strcat(all_feature, " ");
		//====================================================================================================================================

		//====================================================================================================================================

		// Take pp value
		int k = 0;
		while (k < 7)
		{
			sprintf(pp_value, "%.2lf", pp_val[j][k]);
			strcat(all_feature, pp_value);
			strcat(all_feature, " ");
			k++;
			pp_value[0] = '\0';
		}
		//====================================================================================================================================

		//====================================================================================================================================
		// Take pssm value
		fgets(pssmline, sizeof pssmline, pssm);
		
		char *pssm_val = substring(pssmline, 10, 60);
		pssm_counter_new = 0;
		char *tpssm;
		tpssm = strtok(pssm_val, " ");
		while (tpssm != NULL)
		{
			pssm_value = (double)atoi(tpssm) / 9.0;						// Collect each pssm (normalized)
			pssm_values_20[pssm_counter_new] = pssm_value;

			tpssm = strtok(NULL, " ");
			pssm_counter_new++;
		}

		int y = 0;
		if (pssm_counter_new < 20)
		{
			pssm_val = substring(pssmline, 10, 60);
			printf("No 20 PSSMs by tokenize!!! --- don't worry, managed\n");

			for (y = 0; y < 20; y++)
			{
				char *t;
				t = substring(pssm_val, (y * 3) + 1, 3);
				char *token;
				token = trim(t);
				pssm_value = (double)atoi(token) / 9.0;
				pssm_values_20[y] = pssm_value;
				free(t);
			}
		}

		for (y = 0; y < 20; y++)
		{
			sprintf(pssm_val_str, "%.3lf", pssm_values_20[y]);
			strcat(all_feature, pssm_val_str);
			strcat(all_feature, " ");
		}
		//====================================================================================================================================

		//====================================================================================================================================
		sprintf(m_value, "%lf", monogram[j]);					// Take monogram value
		//strcat(all_feature, " ");
		strcat(all_feature, m_value);
		strcat(all_feature, " ");
		m_value[0] = '\0';
		//====================================================================================================================================


		//====================================================================================================================================
		k = 0;
		while (k < 20)
		{
			sprintf(b_value, "%lf", bigram[j][k]);				// Take bigram value
			//strcat(all_feature, " ");
			strcat(all_feature, b_value);
			strcat(all_feature, " ");
			k++;
			b_value[0] = '\0';
		}
		//====================================================================================================================================

		//====================================================================================================================================
		// Take iupred short prediction value

		fgets(iusline, sizeof iusline, ius);

		// take short probability
		char *p_ius = substring(iusline, 13, 6);
		strcat(all_feature, p_ius);
		strcat(all_feature, " ");


		// Take iupred long prediction value

		fgets(iulline, sizeof iulline, iul);

		// take long probability
		char *p_iul = substring(iulline, 13, 6);
		strcat(all_feature, p_iul);
		strcat(all_feature, " ");
		//====================================================================================================================================

		//====================================================================================================================================
		// take terminal information 
		if (res_pos < 0) // left out of window
		{
			strcat(all_feature, "0"); // 3 character for terminal indication (-1 for start (N-terminal))
		}
		else if (res_pos == 0) // first
		{
			strcat(all_feature, "-1"); // 3 character for terminal indication (-1 for start (N-terminal))
		}
		else if (res_pos == 0 + 1) //second
		{
			strcat(all_feature, "-0.8");
		}
		else if (res_pos == 0 + 2) //third
		{
			strcat(all_feature, "-0.6");
		}
		else if (res_pos == 0 + 3) //fourth
		{

			strcat(all_feature, "-0.4");
		}
		else if (res_pos == 0 + 4) //fifth
		{
			strcat(all_feature, "-0.2");
		}
		else if (res_pos == (seqLength - 1)) //last
		{
			strcat(all_feature, "+1"); // 3 character for terminal indication (+1 for end (C-terminal))
		}
		else if (res_pos == (seqLength - 2)) // second last
		{
			strcat(all_feature, "+0.8"); // 3 character for terminal indication (+1 for end (C-terminal))
		}
		else if (res_pos == (seqLength - 3)) // third last
		{
			strcat(all_feature, "+0.6"); // 3 character for terminal indication (+1 for end (C-terminal))
		}
		else if (res_pos == (seqLength - 4)) // fourth last
		{
			strcat(all_feature, "+0.4"); // 3 character for terminal indication (+1 for end (C-terminal))
		}
		else if (res_pos == (seqLength - 5)) // fifth last
		{
			strcat(all_feature, "+0.2"); // 3 character for terminal indication (+1 for end (C-terminal))
		}
		else if (res_pos >= seqLength) // right out of window
		{
			strcat(all_feature, "0"); // 3 character for terminal indication (-1 for start (N-terminal))
		}
		else if ((res_pos > 4) && (res_pos < (seqLength - 5)))
		{
			strcat(all_feature, "0"); // 2 character for terminal indication (0 for Non terminal)
		}


		//====================================================================================================================================
		fprintf(nfp, "%s\n", all_feature);
		fprintf(all, "%s\n", all_feature);
		res_pos++;
	}

	fclose(all);
	fclose(pssm);
	fclose(ius);
	fclose(iul);
	fclose(nfp);
}