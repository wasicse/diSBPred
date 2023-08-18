#define _CRT_SECURE_NO_DEPRECATE
#define LINE_SIZE 15000 //Maximum line size that can be read from the input file
#define BIGRAM_MONOGRAM_LOG_MEDIAN 6.0

#include<stdio.h>
#include<stdlib.h>
#include<string.h>
#include <malloc.h>
#include<ctype.h>
#include<math.h>

FILE *i_listfp; // list file pointer
FILE *nfp;		// log file
FILE *all;		// all feature file
FILE *pp;		// physical properties file
FILE *fasta;	// fasta file
FILE *pssm;		// PSSM file
FILE *bi;		// bigram file	
FILE *mono;		// monogram file
FILE *spx;		// SPINE-X file
FILE *dphi;
FILE *dpsi;
FILE *asa_ex;		// ASA exetended conformation information

char inputFileList[] = "../Input/id_list.txt";
char normalOuputFile[] = "../Output/log/log_feature_collection.txt";
char ppFile[] = "../AdditionalFiles/physiochemical_properties.txt";
char rsa_conformationFile[] = "../AdditionalFiles/RSA_extended_conformation.txt";
char fastaFile[LINE_SIZE];
char pssmFile[LINE_SIZE];
char spxoutFile[LINE_SIZE];
char dphiFile[300];
char dpsiFile[300];
char monogram_filename[300];
char bigram_filename[300];
char all_feature_filename[300];
char id[LINE_SIZE];

int pssm_counter_new = 0;
char rline[LINE_SIZE];
char all_feature[LINE_SIZE];
char ppline[LINE_SIZE];
char pssmline[LINE_SIZE];
char spxline[LINE_SIZE];
char rsa_exntdline[100];
char philine[100];
char psiline[100];
char m_value[30];
char b_value[30];
char pp_value[30];
char pssm_val_str[30];
char p_e_val_str[30];
char p_c_val_str[30];
char p_h_val_str[30];
char asa_val_str[30];
double pssm_value = 0.0;
double pssm_values_20[20];
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
double p_e_value = 0.0;
double p_c_value = 0.0;
double p_h_value = 0.0;
double asa_value = 0.0;

int aa_counter = 0;
int pp_counter = 0;
int local_counter = 0;
int lineLength = 0;
int seqNumber = 0;
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
	i_listfp = fopen(inputFileList, "r");
	if (i_listfp == NULL) {
		fprintf(stderr, "Can't open id list file!\n");
		exit(1);
	}

	nfp = fopen(normalOuputFile, "wb+");
	if (nfp == NULL) {
		fprintf(stderr, "Can't create normal output file!\n");
		exit(1);
	}

	//printf("\nDisPredict Feature Collection START.......\n");
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

	
	fprintf(nfp, "ID found: %s\n", id);
	

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
	// prepare dphi output file, open to read
	strcpy(dphiFile, "../Features/");
	strcat(dphiFile, id);
	strcat(dphiFile, "/");
	strcat(dphiFile, id);
	strcat(dphiFile, ".dphi");

	dphi = fopen(dphiFile, "r");							// open fasta file and read sequence
	if (dphi == NULL) {
		fprintf(stderr, "Can't open pssm file!\n");
		exit(1);
	}
		
	//====================================================================================================================================
	
	//====================================================================================================================================
	// prepare dpsi output file, open to read
	strcpy(dpsiFile, "../Features/");
	strcat(dpsiFile, id);
	strcat(dpsiFile, "/");
	strcat(dpsiFile, id);
	strcat(dpsiFile, ".dpsi");

	dpsi = fopen(dpsiFile, "r");							// open fasta file and read sequence
	if (dpsi == NULL) {
		fprintf(stderr, "Can't open pssm file!\n");
		exit(1);
	}
	
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

	//fprintf(nfp, "MLIST FOUND: %s\n", m_list);
	int i = 0;
	char *t;
	t = strtok(m_list, ",");
	while (t != NULL)
	{
		monogram[i] = (double)atoi(t) / exp(BIGRAM_MONOGRAM_LOG_MEDIAN);											// Collect each monogram
		//fprintf(nfp, "MLIST TOKEN FOUND: %d\n", monogram[i]);
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
	// prepare all feature file name and open to write on that file
	strcpy(all_feature_filename, "../Features/");
	strcat(all_feature_filename, id);
	strcat(all_feature_filename, "/");
	strcat(all_feature_filename, id);
	strcat(all_feature_filename, ".all_dispredict_features");
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
		strcpy(all_feature, "+1 ");
		
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
		//====================================================================================================================================
		
		//====================================================================================================================================
		// Take phi fluctuation value
		fgets(philine, sizeof philine, dphi);
		int phi_l = strlen(philine);
		phi_l = phi_l - 2;
		char* phi_sub = substring(philine, 1, phi_l);
		strcat(all_feature, phi_sub);
		strcat(all_feature, " ");
		//====================================================================================================================================
		
		//====================================================================================================================================
		// Take psi fluctuation value
		fgets(psiline, sizeof psiline, dpsi);
		int psi_l = strlen(psiline);
		psi_l = psi_l - 2;
		char* psi_sub = substring(psiline, 1, psi_l);
		strcat(all_feature, psi_sub);
		strcat(all_feature, " ");
		
		//====================================================================================================================================
		
		//====================================================================================================================================
		// Take monogram value
		sprintf(m_value, "%.4lf", monogram[j]);					
		strcat(all_feature, m_value);
		strcat(all_feature, " ");
		m_value[0] = '\0';
		//====================================================================================================================================

		//====================================================================================================================================
		// Take bigram value
		k = 0;
		while (k < 20)
		{
			sprintf(b_value, "%.4lf", bigram[j][k]);				
			strcat(all_feature, b_value);
			strcat(all_feature, " ");
			k++;
			b_value[0] = '\0';
		}
		//====================================================================================================================================

		//====================================================================================================================================
		// take terminal information 
		if (res_pos < 0) // left out of window
		{
			strcat(all_feature, "0 "); // 3 character for terminal indication (-1 for start (N-terminal))
		}
		else if (res_pos == 0) // first
		{
			strcat(all_feature, "-1 "); // 3 character for terminal indication (-1 for start (N-terminal))
		}
		else if (res_pos == 0 + 1) //second
		{
			strcat(all_feature, "-0.8 ");
		}
		else if (res_pos == 0 + 2) //third
		{
			strcat(all_feature, "-0.6 ");
		}
		else if (res_pos == 0 + 3) //fourth
		{

			strcat(all_feature, "-0.4 ");
		}
		else if (res_pos == 0 + 4) //fifth
		{
			strcat(all_feature, "-0.2 ");
		}
		else if (res_pos == (seqLength - 1)) //last
		{
			strcat(all_feature, "+1 "); // 3 character for terminal indication (+1 for end (C-terminal))
		}
		else if (res_pos == (seqLength - 2)) // second last
		{
			strcat(all_feature, "+0.8 "); // 3 character for terminal indication (+1 for end (C-terminal))
		}
		else if (res_pos == (seqLength - 3)) // third last
		{
			strcat(all_feature, "+0.6 "); // 3 character for terminal indication (+1 for end (C-terminal))
		}
		else if (res_pos == (seqLength - 4)) // fourth last
		{
			strcat(all_feature, "+0.4 "); // 3 character for terminal indication (+1 for end (C-terminal))
		}
		else if (res_pos == (seqLength - 5)) // fifth last
		{
			strcat(all_feature, "+0.2 "); // 3 character for terminal indication (+1 for end (C-terminal))
		}
		else if (res_pos >= seqLength) // right out of window
		{
			strcat(all_feature, "0 "); // 3 character for terminal indication (-1 for start (N-terminal))
		}
		else if ((res_pos > 4) && (res_pos < (seqLength - 5)))
		{
			strcat(all_feature, "0 "); // 2 character for terminal indication (0 for Non terminal)
		}


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
	fclose(dphi);
	fclose(dpsi);
	fprintf(nfp, "\n");

	seqLength = 0;

	rline[0] = '\0';
	all_feature_filename[0] = '\0';
	

	fclose(i_listfp);
	fclose(nfp);

	//printf("DisPredict Feature Collection DONE!!.......\n");
}