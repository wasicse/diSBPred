/*
Author: Sumaiya Iqbal
Part of balancedSSp software
Collect 33 features for final balanced SS prediction by REGAd3p
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

FILE *i_listfp; // disprot database file pointer
FILE *nfp; // normal output
FILE *all; // contains all feature
FILE *dr; //
FILE *annf;
FILE *bi;
FILE *mono;
FILE *fasta;
FILE *pssm;
FILE *pp;		// physical properties file
FILE *dispred;
FILE *asa; // accessible surface area
FILE *dphi;
FILE *dpsi;
FILE *probSeq;

char normalOuputFile[] = "../Output/log/log_feature_collection_for_balanced_SS.txt"; // may change
char ppFile[] = "../AdditionalFiles/physiochemical_properties.txt";
char rline[LINE_SIZE];
char make_directory_command[LINE_SIZE];
char gather_fasta_seq_command[LINE_SIZE];
char gather_annotated_fasta_seq_command[LINE_SIZE];
char gather_pssm_command[LINE_SIZE];
char gather_29feature_command[LINE_SIZE];
char gather_monogram_command[LINE_SIZE];
char gather_bigram_command[LINE_SIZE];
char id[100];

char fseq[LINE_SIZE]; // fasta
char all_feature[LINE_SIZE]; // for a residue
char feature_part_1[LINE_SIZE]; // for spineD,  not required for me
char ppline[LINE_SIZE];
char dispred_output_line[LINE_SIZE];
char m_value[30];
char b_value[30];
char residue[4];
char cnt[200];
char m_list[LINE_SIZE];
char b_list[LINE_SIZE];
char annotation[LINE_SIZE];
double monogram[20];
double bigram[20][20];
double pssm_values_20[20];
char aa[] = "ARNDCQEGHILKMFPSTWYV";
double pp_val[20][7];
double pssm_value = 0.0;
char pp_value[30];
//double asa_val = 0.0;
char asa_line[LINE_SIZE];
char dphi_line[LINE_SIZE];
char dpsi_line[LINE_SIZE];
char problem[] = "nan"; 
int pssm_counter_new = 0;

char pssm_val_str[30];
char all_feature_filename[300];
char fastaFileName[300];
char annotated_fasta_filename[300];
char dr_feature_filename[300];
char monogram_filename[300];
char bigram_filename[300];
char pssmFile[300];
char pssmline[1000];
char dispred_output_filename[LINE_SIZE];
char asa_fileName[300];
char dphi_fileName[300];
char dpsi_fileName[300];


int aa_counter = 0;
int pp_counter = 0;
int lineLength = 0;
int seqNumber = 0;
int seqLength = 0;
int seqNum = 0;
int local_counter = 0;

/*Start of SUB-ROUTINE: substring - It returns a pointer to the substring */
char* substring(char *string, int position, int length)
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

//===============================end subroutines=========================================================//


int main(int argc, char *argv[])
{
	// collect ID
	strcpy(id, argv[1]);

	nfp = fopen(normalOuputFile, "wb+");
	if (nfp == NULL) {
		fprintf(stderr, "Can't create normal output file!\n");
		exit(1);
	}
	
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
		fgets(ppline, sizeof ppline, pp);					// read line and store it to ppline
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
	} // end of pp reading loop

	ppline[0] = '\0';

	fclose(pp);
	//====================================================================================================================================

	//====================================================================================================================================
	// prepare fasta file name, open and collect sequence length
	fastaFileName[0] = '\0';
	strcpy(fastaFileName, "../Features/");  // may need change
	strcat(fastaFileName, id);
	strcat(fastaFileName, "/");
	strcat(fastaFileName, id);
	strcat(fastaFileName, ".fasta");
	
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
	strcpy(pssmFile, "../Features/"); // may need change here 
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
		
	// prepare dispredict output file name, open to read
	strcpy(dispred_output_filename, "../Features/"); // may need change
	strcat(dispred_output_filename, id);
	strcat(dispred_output_filename, "/");
	strcat(dispred_output_filename, id);
	strcat(dispred_output_filename, ".drp");

	dispred = fopen(dispred_output_filename, "r");							// open fasta file and read sequence
	if (dispred == NULL) {
		fprintf(stderr, "Can't open dispredict output file!\n");
		exit(1);
	}

	fgets(dispred_output_line, sizeof dispred_output_line, dispred);				// skip header
	fgets(dispred_output_line, sizeof dispred_output_line, dispred);				// skip header
	fgets(dispred_output_line, sizeof dispred_output_line, dispred);				// skip header
	fgets(dispred_output_line, sizeof dispred_output_line, dispred);				// skip header
	fgets(dispred_output_line, sizeof dispred_output_line, dispred);				// skip header
	fgets(dispred_output_line, sizeof dispred_output_line, dispred);				// skip header
	
	
	//====================================================================================================================================
				
	//====================================================================================================================================
	
	// prepare predicted RSA filename and open======================================================
	strcpy(asa_fileName, "../Features/"); //change
	strcat(asa_fileName, id);
	strcat(asa_fileName, "/");
	strcat(asa_fileName, id);
	strcat(asa_fileName, ".ASAp");

	asa = fopen(asa_fileName, "r");							// open prsa file and read value
	if (asa == NULL) {
		fprintf(stderr, "Can't open predicted rsa file!\n");
		exit(1);
	}
	
	//skip header
	fgets(asa_line, sizeof asa_line, asa);				// skip header
	fgets(asa_line, sizeof asa_line, asa);				// skip header
	fgets(asa_line, sizeof asa_line, asa);				// skip header
	fgets(asa_line, sizeof asa_line, asa);				// skip header
	//====================================================================================================================================
		
	// prepare predicted dphi filename and open======================================================
	strcpy(dphi_fileName, "../Features/"); //change
	strcat(dphi_fileName, id);
	strcat(dphi_fileName, "/");
	strcat(dphi_fileName, id);
	strcat(dphi_fileName, ".dphi");

	dphi = fopen(dphi_fileName, "r");							// open dphi file to read value
	if (dphi == NULL) {
		fprintf(stderr, "Can't open spineX predicted dphi file!\n");
		exit(1);
	}
	
	//====================================================================================================================================
		
	// prepare predicted dpsi filename and open======================================================
	strcpy(dpsi_fileName, "../Features/"); //change
	strcat(dpsi_fileName, id);
	strcat(dpsi_fileName, "/");
	strcat(dpsi_fileName, id);
	strcat(dpsi_fileName, ".dpsi");

	dpsi = fopen(dpsi_fileName, "r");							// open dpsi file to read value
	if (dpsi == NULL) {
		fprintf(stderr, "Can't open spineX predicted dpsi file!\n");
		exit(1);
	}
	
	//====================================================================================================================================
		
		
	// prepare all feature file name and write on that file
	strcpy(all_feature_filename, "../Features/"); 
	strcat(all_feature_filename, id);
	strcat(all_feature_filename, "/");
	strcat(all_feature_filename, id);
	strcat(all_feature_filename, ".BalancedSSP.f33"); // change here, ultimate feature file for seq, VVI 
	
	all = fopen(all_feature_filename, "wb+");					// Open all feature file for writing
	if (all == NULL) {
		fprintf(stderr, "Can't open all feature file!\n");
		exit(1);
	}
	//====================================================================================================================================
		
		
	// Merger all features, WRITE here
	// Feature serial is: class, residue indication, pp, pssm, ss, sa, torsion angle, monogram, bigram, dispred, RSA, dphi, dpsi, terminal 
	int res_pos = 0;
	while (res_pos < seqLength) // start traversing along the sequence ******************************************************************************************************
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
			if (res == aa[j]) // match residue with 20 aa
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
			
			
		//add disorder probability====================================================================================================================================
		fgets(dispred_output_line, sizeof dispred_output_line, dispred);
		char* drp_val = trim(substring(dispred_output_line, 8, strlen(dispred_output_line) - 2)); // to exclude endline character
		strcat(all_feature, drp_val);
		strcat(all_feature, " ");
		
		//====================================================================================================================================
		
		
		//add rsa copied===================================================================================================================================
		
		//fgets(asa_line, sizeof asa_line, asa);		//skip header
		//fgets(asa_line, sizeof asa_line, asa);		//skip header
		fgets(asa_line, sizeof asa_line, asa);		
		char* asa_val = trim(substring(asa_line, 10, strlen(asa_line) - 10)); // to exclude endline character change
		strcat(all_feature, asa_val);
		strcat(all_feature, " ");
		
		//====================================================================================================================================
			
			
		//add dphi===================================================================================================================================
		
		fgets(dphi_line, sizeof dphi_line, dphi);		
		char* dphi_val = trim(substring(dphi_line, 1, strlen(dphi_line) - 1)); // to exclude endline character change
		if(strcmp(problem,dphi_val) == 0)
		{
			printf("NaN in phi fluctuation input!!!");
			exit(1);
		}
		strcat(all_feature, dphi_val);
		strcat(all_feature, " ");
		
		//====================================================================================================================================
		
		//add dpsi===================================================================================================================================
		fgets(dpsi_line, sizeof dpsi_line, dpsi);		
		char* dpsi_val = trim(substring(dpsi_line, 1, strlen(dpsi_line) - 1)); // to exclude endline character change
		if(strcmp(problem,dpsi_val) == 0)
		{
			printf("NaN in psi fluctuation input!!!");
			exit(1);
		}
		strcat(all_feature, dpsi_val);
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
	} // end traversing along the sequence 


	rline[0] = '\0';
	//====================================================================================================================================

	fclose(all);
	fclose(pssm);
	fclose(dispred);
	fclose(asa);
	fclose(dphi);
	fclose(dpsi);
	

	fprintf(nfp, "\n");

	seqLength = 0;

	rline[0] = '\0';
	
	gather_fasta_seq_command[0] = '\0';
	gather_pssm_command[0] = '\0';
	gather_29feature_command[0] = '\0';
	
	all_feature_filename[0] = '\0';


	
	fclose(nfp);

	
} // end of main