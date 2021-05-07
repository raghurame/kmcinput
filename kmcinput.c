#include <stdio.h>
#include <string.h>
#include <math.h>
#include <unistd.h>
#include <stdlib.h>
#include "kmclib.h"

#include <time.h>

int main(int argc, char const *argv[])
{
	if (argc == 1)
	{
		printf("\nINSTRUCTIONS:\n~~~~~~~~~~~~~\n\nargv[0] = program\nargv[1] = input dump file name\nargv[2] = input dihedral file name\nargv[3] = max timeframes to run\n\n");
		exit(1);
	}

	FILE *inputDump, *inputDihedral, *outputJumps, *configJumps, *configJumps2;
	inputDump = fopen (argv[1], "r");
	inputDihedral = fopen (argv[2], "r");
	outputJumps = fopen ("possibleJumps.csv", "w");
	configJumps = fopen ("configJumps.csv", "w");
	configJumps2 = fopen ("configJumps2.csv", "w");


	int nDihedrals = countDihedrals (inputDihedral), nAtoms = countAtoms (inputDump), isEOF_dump = 0, isEOF_dihedral = 0, init = 0, nBackboneDihedrals, nTimeframes = 0, nLattice = 0, maxtf = atoi (argv[3]);
	printf("%d atoms detected...\n", nAtoms);
	printf("%d dihedrals detected...\n", nDihedrals);

	int transition_CR, transition_CL, transition_CC, transition_CM;
	int transition_RR, transition_RL, transition_RC, transition_RM;
	int transition_LR, transition_LL, transition_LC, transition_LM;
	int transition_MR, transition_ML, transition_MC, transition_MM;

	while ((!isEOF_dihedral) && (!isEOF_dump))
	{

		DUMPINFO *dump;
		dump = (DUMPINFO *) malloc (nAtoms * sizeof (DUMPINFO));

		DIHEDRALINFO *dihedral;
		dihedral = (DIHEDRALINFO *) malloc (nDihedrals * sizeof (DIHEDRALINFO));

		LATTICE_PROBABILITIES conigurations;

		KMC_LATTICE *lattice, *lattice_tminus1;

		KMC_JUMPS possibleJumps;

		// TRANSITIONS *data1, *data2;
		// data1 = (TRANSITIONS *) malloc (maxtf * sizeof(TRANSITIONS));
		// data2 = (TRANSITIONS *) malloc (maxtf * sizeof(TRANSITIONS));

		dump = readDump (inputDump, nAtoms, &isEOF_dump);
		dihedral = readDihedral (inputDihedral, nDihedrals, &isEOF_dihedral, 1);

		if (nTimeframes == 0)
		{
			nBackboneDihedrals = countBackboneDihedrals (dump, nAtoms, dihedral, nDihedrals, 3);
			lattice = (KMC_LATTICE *) malloc (nBackboneDihedrals * sizeof (KMC_LATTICE));
			lattice_tminus1 = (KMC_LATTICE *) malloc (nBackboneDihedrals * sizeof (KMC_LATTICE));

			fprintf(configJumps, "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t\n", 
				"coil to right",
				"coil to left",
				"coil to coil",
				"coil to meso",
				"right to right",
				"right to left",
				"right to coil",
				"right to meso",
				"left to right",
				"left to left",
				"left to coil",
				"left to meso"
				"meso to right",
				"meso to left",
				"meso to coil",
				"meso to meso");

			printf("reading dump file...\n");
			printf("reading dihedral file...\n");
			printf("printing configuration jumps...\n");
		}

		if (nTimeframes > 0)
			lattice_tminus1 = lattice;

		lattice = computeKMCLattice (dump, nAtoms, dihedral, nDihedrals, nBackboneDihedrals, 3, &nLattice);

		if (nTimeframes > 0)
		{
			// Algorithm 1
			possibleJumps = computeKCMJumps (lattice, lattice_tminus1, nLattice);

			fprintf(configJumps, "%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\n", 
				possibleJumps.transition_CR,
				possibleJumps.transition_CL,
				possibleJumps.transition_CC,
				possibleJumps.transition_CM,
				possibleJumps.transition_RR,
				possibleJumps.transition_RL,
				possibleJumps.transition_RC,
				possibleJumps.transition_RM,
				possibleJumps.transition_LR,
				possibleJumps.transition_LL,
				possibleJumps.transition_LC,
				possibleJumps.transition_LM,
				possibleJumps.transition_MR,
				possibleJumps.transition_ML,
				possibleJumps.transition_MC,
				possibleJumps.transition_MM);

			// data1[nTimeframes-1].CR = possibleJumps.transition_CR; data1[nTimeframes-1].CL = possibleJumps.transition_CL; data1[nTimeframes-1].CC = possibleJumps.transition_CC; data1[nTimeframes-1].CM = possibleJumps.transition_CM; data1[nTimeframes-1].RR = possibleJumps.transition_RR; data1[nTimeframes-1].RL = possibleJumps.transition_RL; data1[nTimeframes-1].RC = possibleJumps.transition_RC; data1[nTimeframes-1].RM = possibleJumps.transition_RM; data1[nTimeframes-1].LR = possibleJumps.transition_LR; data1[nTimeframes-1].LL = possibleJumps.transition_LL; data1[nTimeframes-1].LC = possibleJumps.transition_LC; data1[nTimeframes-1].LM = possibleJumps.transition_LM; data1[nTimeframes-1].MR = possibleJumps.transition_MR; data1[nTimeframes-1].ML = possibleJumps.transition_ML; data1[nTimeframes-1].MC = possibleJumps.transition_MC; data1[nTimeframes-1].MM = possibleJumps.transition_MM; 

			// Algorithm 2
			possibleJumps = computeKCMJumps2 (lattice, lattice_tminus1, nLattice);

			fprintf(configJumps2, "%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\n", 
				possibleJumps.transition_CR,
				possibleJumps.transition_CL,
				possibleJumps.transition_CC,
				possibleJumps.transition_CM,
				possibleJumps.transition_RR,
				possibleJumps.transition_RL,
				possibleJumps.transition_RC,
				possibleJumps.transition_RM,
				possibleJumps.transition_LR,
				possibleJumps.transition_LL,
				possibleJumps.transition_LC,
				possibleJumps.transition_LM,
				possibleJumps.transition_MR,
				possibleJumps.transition_ML,
				possibleJumps.transition_MC,
				possibleJumps.transition_MM);

			// data2[nTimeframes-1].CR = possibleJumps.transition_CR; data2[nTimeframes-1].CL = possibleJumps.transition_CL; data2[nTimeframes-1].CC = possibleJumps.transition_CC; data2[nTimeframes-1].CM = possibleJumps.transition_CM; data2[nTimeframes-1].RR = possibleJumps.transition_RR; data2[nTimeframes-1].RL = possibleJumps.transition_RL; data2[nTimeframes-1].RC = possibleJumps.transition_RC; data2[nTimeframes-1].RM = possibleJumps.transition_RM; data2[nTimeframes-1].LR = possibleJumps.transition_LR; data2[nTimeframes-1].LL = possibleJumps.transition_LL; data2[nTimeframes-1].LC = possibleJumps.transition_LC; data2[nTimeframes-1].LM = possibleJumps.transition_LM; data2[nTimeframes-1].MR = possibleJumps.transition_MR; data2[nTimeframes-1].ML = possibleJumps.transition_ML; data2[nTimeframes-1].MC = possibleJumps.transition_MC; data2[nTimeframes-1].MM = possibleJumps.transition_MM; 

			// To do: Print running average, reduced numbers
			// NOTE: Do the above seperately
		}
		
		nTimeframes++;
		if (!(nTimeframes % 100))
		{
			printf("nTimeframes scanned so far... %d\n", nTimeframes);
		}

		if (nTimeframes == maxtf)
		{
			exit(1);
		}

		free (dump);
		free (dihedral);
	}

	fclose (inputDump);
	fclose (inputDihedral);
	fclose (outputJumps);
	fclose (configJumps);
	fclose (configJumps2);

	return 0;
}