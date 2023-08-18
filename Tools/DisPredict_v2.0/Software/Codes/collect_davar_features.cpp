#define _CRT_SECURE_NO_DEPRECATE
#define LINE_SIZE 30000 //Maximum line size that can be read from the input file
#define BIGRAM_MONOGRAM_LOG_MEDIAN 6.0

#include<stdio.h>
#include<stdlib.h>
#include<string.h>
#include <malloc.h>
#include<ctype.h>
#include<math.h>

FILE *nfp;		// log file
FILE *all;		// all feature file
FILE *pp;		// physical properties file
FILE *fasta;	// fasta file
FILE *pssm;		// PSSM file
FILE *spx;		// SPINE-X file
FILE *asa_ex;		// ASA exetended conformation information
FILE *ius;

char normalOuputFile[] = "../Output/log/log_davar_feature_collection.txt";
char ppFile[] = "../AdditionalFiles/physiochemical_properties.txt";
char rsa_conformationFile[] = "../AdditionalFiles/RSA_extended_conformation.txt";
char fastaFile[LINE_SIZE];
char pssmFile[LINE_SIZE];
char spxoutFile[LINE_SIZE];
char iupredSFile[LINE_SIZE];
char all_feature_filename[300];
char id[LINE_SIZE];

char rline[LINE_SIZE];
char all_feature[LINE_SIZE];
char ppline[LINE_SIZE];
char pssmline[LINE_SIZE];
char spxline[LINE_SIZE];
char iusline[LINE_SIZE];
char rsa_exntdline[100];
char m_value[30];
char b_value[30];
char pp_value[30];
char pssm_val_str[30];
char p_ius_val_str[30];
char p_e_val_str[30];
char p_c_val_str[30];
char p_h_val_str[30];
char asa_val_str[30];
char phi_val_str[30];
char psi_val_str[30];
double pssm_value = 0.0;
char residue[4];
char cnt[200];
char m_list[LINE_SIZE];
char b_list[LINE_SIZE];
char fseq[LINE_SIZE];
double monogram[20];
double bigram[20][20];
double pp_val[20][7];
double rsa_extnd_val[20];
char aa[] = "ARNDCQEGHILKMFPSTWYV";
double p_ius_value = 0.0;
double p_e_value = 0.0;
double p_c_value = 0.0;
double p_h_value = 0.0;
double asa_value = 0.0;
double phi_value = 0.0;
double psi_value = 0.0;
double pssm_values_20[20];

int aa_counter = 0;
int pp_counter = 0;
int local_counter = 0;
int lineLength = 0;
int seqNumber = 0;
int seqLength = 0;
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
	strcpy(id, argv[1]);
	
	nfp = fopen(normalOuputFile, "wb+");
	if (nfp == NULL) {
		fprintf(stderr, "Can't create normal output file!\n");
		exit(1);
	}

	//printf("\nPhi/Psi Angle fluctuation prediction Feature Collection START.......\n");
	//====================================================================================================================================
	// prepare pp file name, open to read, collect pp

	pp = fopen(ppFile, "r");								// open pp file
	if (pp == NULL) {
		fprintf(stderr, "Can't open pp file!\n");
		exit(1);
	}

	aa_counter = 0;
	while (aa_counter < 20)									// read for 20 amino acids
	{
		fgets(ppline, sizeof ppline, pp);					
		pp_counter = 0;									// collect 7 properties
		char *t_pp;
		t_pp = strtok(ppline, " ");
		while (t_pp != NULL)
		{
			if (pp_counter > 0)								// skip first column (AA)
			{
				pp_val[aa_counter][pp_counter - 1] = (double)atof(t_pp);	// Collect each pp
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
	// prepare rsa extended conformational value file name, open to read, collect rsa_entnd values

	asa_ex = fopen(rsa_conformationFile, "r");								// open pp file
	if (asa_ex == NULL) {
		fprintf(stderr, "Can't open rsa extended conformation values file!\n");
		exit(1);
	}

	aa_counter = 0;
	while (aa_counter < 20)									// read for 20 amino acids
	{
		fgets(rsa_exntdline, sizeof rsa_exntdline, asa_ex);
		pp_counter = 0;									// collect 7 properties
		char *t_sax;
		t_sax = strtok(rsa_exntdline, ":");
		while (t_sax != NULL)
		{
			if (pp_counter == 1)								// skip first column (AA)
			{
				rsa_extnd_val[aa_counter] = (double)atof(t_sax);	// Collect each pp
			}
			pp_counter++;
			t_sax = strtok(NULL, ":");
		}
		aa_counter++;
	}

	rsa_exntdline[0] = '\0';

	fclose(asa_ex);
	//====================================================================================================================================

	// collect ID
	
	fprintf(nfp, "ID found: %s\n", id);
	//printf("ID: %s\n", id);

	//====================================================================================================================================
	// prepare fasta file name, open, read and collect length
	strcpy(fastaFile, "../Features/");
	strcat(fastaFile, id);
	strcat(fastaFile, "/");
	strcat(fastaFile, id);
	strcat(fastaFile, ".fasta");

	fasta = fopen(fastaFile, "r");							// open fasta file and read sequence
	if (fasta == NULL) {
		fprintf(stderr, "Can't open annotated fasta file!\n");
		exit(1);
	}

	fgets(fseq, sizeof fseq, fasta);				// skip header
	fgets(fseq, sizeof fseq, fasta);				// collect fasta sequence

	seqLength = strlen(fseq);						// collect length
	seqLength--;
	fseq[seqLength] = '\0';

	fclose(fasta);
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
	// prepare spx output file, open to read
	strcpy(spxoutFile, "../Features/");
	strcat(spxoutFile, id);
	strcat(spxoutFile, "/");
	strcat(spxoutFile, id);
	strcat(spxoutFile, ".spXout");

	spx = fopen(spxoutFile, "r");							// open fasta file and read sequence
	if (spx == NULL) {
		fprintf(stderr, "Can't open pssm file!\n");
		exit(1);
	}

	fgets(spxline, sizeof spxline, spx);				// skip header
	
	//====================================================================================================================================

	//====================================================================================================================================
	// prepare iupred short output file, open to read
	strcpy(iupredSFile, "../Features/");
	strcat(iupredSFile, id);
	strcat(iupredSFile, "/");
	strcat(iupredSFile, id);
	strcat(iupredSFile, ".iupredS");

	ius = fopen(iupredSFile, "r");							// open iupreds output file and read sequence
	if (ius == NULL) {
		fprintf(stderr, "Can't open pssm file!\n");
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
	
	// prepare all feature file name and open to write on that file
	strcpy(all_feature_filename, "../Features//");
	strcat(all_feature_filename, id);
	strcat(all_feature_filename, "/");
	strcat(all_feature_filename, id);
	strcat(all_feature_filename, ".davar_features");
	fprintf(nfp, "All feature file name: %s\n", all_feature_filename);

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
		
		//====================================================================================================================================
		// Take class (annotation) -- dummy (+1 == disorder) for test protein
		//strcpy(all_feature, "+1 ");
		
		//====================================================================================================================================

		//====================================================================================================================================
		// Take residue information
		
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
		strcpy(all_feature, residue);
		strcat(all_feature, " ");
			
		//====================================================================================================================================

		//====================================================================================================================================
		// Take pp value
		int k = 0;
		while (k < 7)
		{
			sprintf(pp_value, "%.2lf", pp_val[j][k]);				
			if(k == 0)
			{
				strcpy(all_feature, pp_value);
				strcat(all_feature, " ");
				k++;	
			}
			else
			{
				strcat(all_feature, pp_value);
				strcat(all_feature, " ");
				k++;
			}
			
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
		// Take spinex output
		fgets(spxline, sizeof spxline, spx);
		
		// take strand probability
		char *p_e = substring(spxline, 39, 6);
		p_e_value = (double)(2.0 * atof(p_e)) - 1.0;			// scale within (-1, 1)
		sprintf(p_e_val_str, "%.4lf", p_e_value);
		strcat(all_feature, p_e_val_str);
		strcat(all_feature, " ");

		// take coil probability
		char *p_c = substring(spxline, 48, 6);
		p_c_value = (double)(2.0 * atof(p_c)) - 1.0;			// scale within (-1, 1)
		sprintf(p_c_val_str, "%.4lf", p_c_value);
		strcat(all_feature, p_c_val_str);
		strcat(all_feature, " ");

		// take strand probability
		char *p_h = substring(spxline, 58, 6);
		p_h_value = (double)(2.0 * atof(p_h)) - 1.0;			// scale within (-1, 1)
		sprintf(p_h_val_str, "%.4lf", p_h_value);
		strcat(all_feature, p_h_val_str);
		strcat(all_feature, " ");

		// take phi
		char *p_phi = substring(spxline, 21, 6);
		phi_value = (double)atof(p_phi) / 180.0;			// normalize by 180.0 
		sprintf(phi_val_str, "%.4lf", phi_value);
		strcat(all_feature, phi_val_str);
		strcat(all_feature, " ");

		// take psi
		char *p_psi = substring(spxline, 29, 6);
		psi_value = (double)atof(p_psi) / 180.0;			// normalize by 180.0 
		sprintf(psi_val_str, "%.4lf", psi_value);
		strcat(all_feature, psi_val_str);
		strcat(all_feature, " ");

		// take asa
		char *p_asa = substring(spxline, 88, 6);
		asa_value = (double)(2.0 * (atof(p_asa) / rsa_extnd_val[j])) - 1.0;			// normalize by extended conformation and scale within (-1, 1)
		sprintf(asa_val_str, "%.4lf", asa_value);
		strcat(all_feature, asa_val_str);
		strcat(all_feature, " ");
		
		p_e_val_str[0] = '\0';
		p_c_val_str[0] = '\0';
		p_h_val_str[0] = '\0';
		asa_val_str[0] = '\0';
		phi_val_str[0] = '\0';
		psi_val_str[0] = '\0';
		//====================================================================================================================================

		//====================================================================================================================================
		// Take iupred short prediction value
		
		fgets(iusline, sizeof iusline, ius);
		
		// take strand probability
		char *p_ius = substring(iusline, 13, 6);
		p_ius_value = (double)(2.0 * atof(p_ius)) - 1.0;			// scale within (-1, 1)
		sprintf(p_ius_val_str, "%.4lf", p_ius_value);
		strcat(all_feature, p_ius_val_str);
		strcat(all_feature, " ");
		
		strcat(all_feature, "0.0 0.0");									// fill up gap
		//====================================================================================================================================

		
		//====================================================================================================================================
		fprintf(nfp, "FEATURE\n");
		fprintf(nfp, "%s\n", all_feature);
		fprintf(all, "%s\n", all_feature);
		res_pos++;
	}


	rline[0] = '\0';
	//====================================================================================================================================

	fclose(all);
	fclose(pssm);
	fclose(spx);
	fclose(ius);
	fprintf(nfp, "\n");

	seqLength = 0;

	rline[0] = '\0';
	all_feature_filename[0] = '\0';
	
	fclose(nfp);

	//printf("Phi/Psi Angle fluctuation prediction Feature Collection DONE!!.......\n");
}