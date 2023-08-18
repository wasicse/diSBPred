/*
Author: Sumaiya iqbal, CSE, UNO # 2450707
Calculate sequence composition based Position Specific Estimated Energy (sPSEE)
Neighborhood size can be varried by specifying the start and end the region around both side of the target residue
*/


#define _CRT_SECURE_NO_DEPRECATE
#define LINE_SIZE 30000

#include<stdio.h>
#include<stdlib.h>
#include<string.h>
#include<malloc.h>
#include<ctype.h>
#include<math.h>

FILE *nfp;				// normal output file
FILE *fasta;			// indivudual fatsa file
FILE *asa;				// individual predicted ASA file
FILE *asa_ex;			// ASA exetended conformation information
FILE *energyF;			// energy matrix file
FILE *sPSEIE_output;	// output

char normalOuputFile[] = "../Output/log/log_generate_psee.txt";
char rsa_conformationFile[] = "../AdditionalFiles/EASA_new.txt";
char energyFile[] = "../AdditionalFiles/DT_predicted_energy.txt";
char asa_fileName[300];
char sPSEIE_output_fileName[300];
char commnd[LINE_SIZE];

char aa[] = "ARNDCQEGHILKMFPSTWYV";
char rline[LINE_SIZE];
char fseq[LINE_SIZE];
char sPSEIE_output_line[100];
char asa_line[300];
char energy_line[500];
char fastaFileName[LINE_SIZE]; 
char rsa_exntdline[100];
double rsa_extnd_val[20];
double energy_val[20][20];

double res_asap = 0.0;
double res_exposure = 0.0;
double Nres_exposure = 0.0;
double res_burial = 0.0;
double Nres_burial = 0.0;
double res_res_energy = 0.0;
int lineLength = 0;
int seqNumber = 0;
int seqLength = 0;
int aa_counter = 0;
int pp_counter = 0;
int region_length = 0;

double sPSEIE = 0.0;						// energy output
double sPSEIE_normalized = 0.0;						// energy output
double lNE = 0.0;							// energy from left neighborhood
double rNE = 0.0;							// energy from right neighborhood
double lNE_norm = 0.0;
double rNE_norm = 0.0;
int radius = 0;

int sRange = 0;
int eRange = 0;
int lNstart = 0;							// start position of left neighborhood region
int	lNend = 0;								// end position of left nighborhood region

int rNstart = 0;							// start position of right neighborhood region
int	rNend = 0;								// end position of right nighborhood region

char sN[100];
char eN[100];
char lE_string[100];
char rE_string[100];
char tE_string[100];
char total_E_string[100];
char normalized_E_string[100];
char id[200];

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

//===============================end subroutines=========================================================//

int main(int argc, char *argv[])
{
	if (argc < 4)
	{
		printf("You have to give two inputs --- start of range and end of range!!!\n");
	}

	//====================================================================================================================================
	strcpy(id, argv[1]);
	strcpy(sN, argv[2]);			// start of neighborhood region
	sRange = atoi(sN);
	strcpy(eN, argv[3]);			// end of neighborhood region
	eRange = atoi(eN);


	region_length = 2 * (eRange - sRange + 1);	// total length of the region
	//====================================================================================================================================

	//====================================================================================================================================
	nfp = fopen(normalOuputFile, "wb+");
	if (nfp == NULL) {
		fprintf(stderr, "Can't create normal output file!\n");
		exit(1);
	}
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

	//====================================================================================================================================
	// read energy matrix and store

	energyF = fopen(energyFile, "r");								// open energy matrix file
	if (energyF == NULL) {
		fprintf(stderr, "Can't open energy matrix file!\n");
		exit(1);
	}
	fgets(energy_line, sizeof energy_line, energyF);				// skip header

	aa_counter = 0;
	while (aa_counter < 20)											// read for 20 amino acids
	{
		fgets(energy_line, sizeof energy_line, energyF);
		pp_counter = 0;												// counter
		char *t_energy;
		t_energy = strtok(energy_line, " ");
		while (t_energy != NULL)
		{
			if (pp_counter > 0)								// skip first column (AA)
			{
				energy_val[aa_counter][pp_counter - 1] = (double)atof(t_energy);	// Collect each pp
			}
			pp_counter++;
			t_energy = strtok(NULL, " ");
		}
		aa_counter++;
	}
	//printf("%c and %c---------------%lf\n", aa[2], aa[5],energy_val[2][5]);
	energy_line[0] = '\0';

	fclose(energyF);
	//====================================================================================================================================

	//====================================================================================================================================
	// prepare fasta file name, open and collect sequence length
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
	//printf("Sequence No: %d ------------- ID: %s, Length: %d\n", seqNumber, id, seqLength);
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

	double* asap_storage;
	asap_storage = (double *)malloc((seqLength + 1)*sizeof(double));
	if (asap_storage == NULL)
	{
		printf("Unable to allocate memory.\n");
		exit(EXIT_FAILURE);
	}
	int t = 0;
	while (t < seqLength)
	{
		fgets(asa_line, sizeof asa_line, asa);
		char* asa_val = trim(substring(asa_line, 10, strlen(asa_line) - 10)); // to exclude endline character change
		res_asap = (double)atof(asa_val);
		*(asap_storage + t) = res_asap;
		t++;
	}

	fclose(asa);
	//====================================================================================================================================

	//====================================================================================================================================
	// prepare output file name and open to write======================================================
	strcpy(sPSEIE_output_fileName, "../Features/"); //change
	strcat(sPSEIE_output_fileName, id);
	strcat(sPSEIE_output_fileName, "/");
	strcat(sPSEIE_output_fileName, id);
	strcat(sPSEIE_output_fileName, ".PSEE");
		
	sPSEIE_output = fopen(sPSEIE_output_fileName, "wb+");							// open prsa file and read value
	if (sPSEIE_output == NULL) {
		fprintf(stderr, "Can't open output file!\n");
		exit(1);
	}
	fprintf(sPSEIE_output, ">%s\n", id);
	fprintf(sPSEIE_output, "#SR AA left_neighborhood_energy right_neighborhood_energy PSEE PSEE_per_residue\n");		// header
	//printf("cheking!!!\n");
	//====================================================================================================================================

	// Start of residue specific energy calculation
	int res = 1;
	while (res <= seqLength)
	{
			
		//====================================================================================================================================
		char residue = fseq[res-1];
		int j = 0;
		for (; j < 20; j++)
		{
			if (residue == aa[j])
			{
				break;
			}
		}
		sprintf(sPSEIE_output_line, "%d", res);
		strcat(sPSEIE_output_line, " ");
		char* string_res = substring(fseq, res, 1);
		strcat(sPSEIE_output_line, string_res);
		strcat(sPSEIE_output_line, " ");
		//====================================================================================================================================
			
		//====================================================================================================================================

		sPSEIE = 0;
		lNE = 0.0;							// energy from left neighborhood
		rNE = 0.0;							// energy from right neighborhood

		lNstart = 0;						// start position of left neighborhood region
		lNend = 0;							// end position of left nighborhood region

		rNstart = 0;						// start position of right neighborhood region
		rNend = 0;							// end position of right nighborhood region

		//====================================================================================================================================
		// Compute start and end of left neighborhood subsequence
		if ((res - sRange) <= 0)				// no neighborhood region at the left of target
		{
			lNE = 0;
			lNstart = 0;
			lNend = 0;
		}
		else
		{
			if ((res - eRange) > 0)				// maximum neighborhood with 100 residues present
			{
				lNstart = res - eRange;
				lNend = res - sRange; 
			}
			else
			{
				lNstart = 1;
				lNend = res - sRange;
			}
		}

		// Now compute energy contribution by this left neighborhood region
		// calculate only exists
		if ((lNstart > 0) && (lNend > 0))
		{
			int eT = 0;
			for (eT = lNstart; eT <= lNend; eT++)
			{
				// find current neighborhood residue (AA) type
				char Nresidue = fseq[eT - 1];
				int k = 0;
				for (; k < 20; k++)
				{
					if (Nresidue == aa[k])
					{
						break;
					}
				}
				radius = lNend - lNstart + 1;
				// calculate solvent exposure of the neighborhood residue
				Nres_exposure = (*(asap_storage + (eT-1))) / rsa_extnd_val[k];

				// calculate solvent exposure of the neighborhood residue
				Nres_burial = 1 - Nres_exposure;

				// extract interaction energy with this neighborhood residue
				res_res_energy = energy_val[j][k];

				// calculte energy contribution from the left part of neighborhood
				lNE = lNE + (res_res_energy * Nres_burial);
				lNE_norm = (double)(lNE / radius);
			}
		}
		else
		{
			lNE = 0;
		}
		//====================================================================================================================================

		//====================================================================================================================================
		// Compute start and end of right neighborhood subsequence
		if ((seqLength - res) < sRange)				// no neighborhood region at the right of target
		{
			rNE = 0;
			rNstart = 0;
			rNend = 0;
		}
		else
		{
			if ((seqLength - res) >= eRange)				// maximum neighborhood with 100 residues present
			{
				rNstart = res + sRange;
				rNend = res + eRange;
			}
			else
			{
				rNstart = res + sRange;
				rNend = seqLength;
			}
		}

		// Now compute energy contribution by this left neighborhood region
		// calculate only exists
		if ((rNstart > 0) && (rNend > 0))
		{
			int eT = 0;
			for (eT = rNstart; eT <= rNend; eT++)
			{
				// find current neighborhood residue (AA) type
				char Nresidue = fseq[eT - 1];
				int k = 0;
				for (; k < 20; k++)
				{
					if (Nresidue == aa[k])
					{
						break;
					}
				}
					
				radius = rNend - rNstart + 1;
				// calculate solvent exposure of the neighborhood residue
				Nres_exposure = (*(asap_storage + (eT - 1))) / rsa_extnd_val[k];

				// calculate solvent exposure of the neighborhood residue
				Nres_burial = 1 - Nres_exposure;

				// extract interaction energy with this neighborhood residue
				res_res_energy = energy_val[j][k];

				// calculte energy contribution from the left part of neighborhood
				rNE = rNE + (res_res_energy * Nres_burial);
				rNE_norm = (double)(rNE / radius);
			}
		}
		else
		{
			rNE = 0;
		}
		//====================================================================================================================================
		// calculate sPSEIE
		// calculate solvent exposure of the target residue
		res_exposure = (*(asap_storage + (res - 1))) / rsa_extnd_val[j];

		// calculate solvent exposure of the target residue
		res_burial = 1 - res_exposure;

		// calculte sPSEIE
		sPSEIE = res_burial * (lNE + rNE);
		sPSEIE_normalized = res_burial * (lNE_norm + rNE_norm);
		//printf("sPSEIE: %lf\n", sPSEIE);


		// setup writing
		sprintf(lE_string, "%lf", lNE);
		sprintf(rE_string, "%lf", rNE);
		sprintf(total_E_string, "%lf", sPSEIE);
		sprintf(normalized_E_string, "%lf", sPSEIE_normalized);
			
		strcat(sPSEIE_output_line, lE_string);
		strcat(sPSEIE_output_line, " ");
		strcat(sPSEIE_output_line, rE_string);
		strcat(sPSEIE_output_line, " ");
		strcat(sPSEIE_output_line, total_E_string);
		strcat(sPSEIE_output_line, " ");
		strcat(sPSEIE_output_line, normalized_E_string);

		fprintf(sPSEIE_output, "%s\n", sPSEIE_output_line);
		res++;

		//exit(1);
	}
		//fclose(asa);
	fclose(sPSEIE_output);

	fclose(nfp);
}