#include <stdio.h>
#include <string.h>
#include <math.h>
#include <unistd.h>
#include <stdlib.h>
#include <assert.h>
#include <time.h>
#include "kmclib.h"

int countDihedrals (FILE *inputDihedral)
{
	char waste[5000];
	int nDihedrals = 0;

	for (int i = 0; i < 4; ++i) {
		fgets (waste, 3000, inputDihedral); }

	sscanf (waste, "%d\n", &nDihedrals);

	return nDihedrals;
}

int countAtoms (FILE *inputDump)
{
	char waste[5000];
	int nAtoms = 0;

	for (int i = 0; i < 4; ++i)
		fgets (waste, 3000, inputDump);

	sscanf (waste, "%d\n", &nAtoms);

	return nAtoms;
}

DUMPINFO *readDump (FILE *inputDump, int nAtoms, int *isEOF_dump)
{
	DUMPINFO *dump;
	dump = (DUMPINFO *) malloc (nAtoms * sizeof (DUMPINFO));

	char lineString[5000];

	if (fgets (lineString, 5000, inputDump) == NULL)
	{

		(*isEOF_dump) = 1;
		return dump;
	}
	else
	{
		for (int i = 0; i < 8; ++i)
			fgets (lineString, 5000, inputDump);

		for (int i = 0; i < nAtoms; ++i)
		{
			fgets (lineString, 5000, inputDump);
			sscanf (lineString, "%d %d %d %f %f %f %*f %*f %*f %d %d %d\n", &dump[i].sino, &dump[i].atomType, &dump[i].molType, &dump[i].x, &dump[i].y, &dump[i].z, &dump[i].ix, &dump[i].iy, &dump[i].iz);
			// printf("%d %d %d %f %f %f %d %d %d\n", dump[i].sino, dump[i].atomType, dump[i].molType, dump[i].x, dump[i].y, dump[i].z, dump[i].ix, dump[i].iy, dump[i].iz);
			// sleep(1);
		}

		(*isEOF_dump) = 0;
		return dump;
	}
}

DIHEDRALINFO *readDihedral (FILE *inputDihedral, int nDihedrals, int *isEOF_dihedral, int sort)
{
	DIHEDRALINFO *dihedral, tempDihedral;
	dihedral = (DIHEDRALINFO *) malloc (nDihedrals * sizeof (DIHEDRALINFO));

	char lineString[5000];

	if (fgets (lineString, 5000, inputDihedral) == NULL)
	{
		(*isEOF_dihedral) = 1;
		return dihedral;
	}
	else
	{
		for (int i = 0; i < 8; ++i)
			fgets (lineString, 5000, inputDihedral);

		for (int i = 0; i < nDihedrals; ++i)
		{
			fgets (lineString, 5000, inputDihedral);
			sscanf (lineString, "%d %f %d %d %d %d\n", &dihedral[i].sino, &dihedral[i].dihedralAngle, &dihedral[i].atom1, &dihedral[i].atom2, &dihedral[i].atom3, &dihedral[i].atom4);
			// printf("%d %f %d %d %d %d\n", dihedral[i].sino, dihedral[i].dihedralAngle, dihedral[i].atom1, dihedral[i].atom2, dihedral[i].atom3, dihedral[i].atom4);
			// sleep(1);
		}

		(*isEOF_dihedral) = 0;

		if (sort == 1)
		{
			for (int i = 0; i < nDihedrals; ++i)
			{
				for (int j = i + 1; j < nDihedrals; ++j)
				{
					if ((dihedral[i].atom1 + dihedral[i].atom2 + dihedral[i].atom3 + dihedral[i].atom4) > (dihedral[j].atom1 + dihedral[j].atom2 + dihedral[j].atom3 + dihedral[j].atom4))
					{
						tempDihedral = dihedral[j];
						dihedral[j] = dihedral[i];
						dihedral[i] = tempDihedral;
					}
				}
			}

			for (int i = 0; i < nDihedrals; ++i)
				dihedral[i].sino = i + 1;
		}

		return dihedral;
	}
}

int *checkPendant (DUMPINFO *dump, int nAtoms, int pendantMolID)
{
	int *isPendant;
	isPendant = (int *) malloc ((nAtoms + 1) * sizeof (int));

	for (int i = 0; i < nAtoms; ++i)
	{
		if (dump[i].molType == pendantMolID)
			isPendant [i + 1] = 1;
		else
			isPendant [i + 1] = 0;
	}

	return isPendant;
}

int countBackboneDihedrals (DUMPINFO *dump, int nAtoms, DIHEDRALINFO *dihedral, int nDihedrals, int pendantMolID)
{
	int *isPendant, nBackboneDihedrals = 0;
	isPendant = (int *) malloc ((nAtoms + 1) * sizeof (int));

	isPendant = checkPendant (dump, nAtoms, pendantMolID);

	for (int i = 0; i < nDihedrals; ++i)
	{
		if ((isPendant [dihedral[i].atom1] == 0) && (isPendant [dihedral[i].atom2] == 0) && (isPendant [dihedral[i].atom3] == 0) && (isPendant [dihedral[i].atom4] == 0))
			nBackboneDihedrals++;
	}

	free (isPendant);

	printf("nBackboneDihedrals: %d\n", nBackboneDihedrals);

	return nBackboneDihedrals;
}

KMC_LATTICE *computeKMCLattice (DUMPINFO *dump, int nAtoms, DIHEDRALINFO *dihedral, int nDihedrals, int nBackboneDihedrals, int pendantMolID, int *nLattice)
{
	
	int *isPendant, nLattice_local = 0 /*, nLattice = 0*/;
	isPendant = (int *) malloc ((nAtoms + 1) * sizeof (int));
	isPendant = checkPendant (dump, nAtoms, pendantMolID);	

	assert (nLattice);

	KMC_LATTICE *lattice;
	lattice = (KMC_LATTICE *) malloc (nBackboneDihedrals * sizeof (KMC_LATTICE));

	for (int i = 0; i < nDihedrals; ++i)
	{
		if ((isPendant [dihedral[i].atom1] == 0) && (isPendant [dihedral[i].atom2] == 0) && (isPendant [dihedral[i].atom3] == 0) && (isPendant [dihedral[i].atom4] == 0))
		{
			lattice[nLattice_local].sino1 = dihedral[i].atom2;
			lattice[nLattice_local].sino2 = dihedral[i].atom3;
			lattice[nLattice_local].atomType1 = dump[(dihedral[i].atom2 - 1)].atomType;
			lattice[nLattice_local].atomType2 = dump[(dihedral[i].atom3 - 1)].atomType;
			lattice[nLattice_local].molType1 = dump[(dihedral[i].atom2 - 1)].molType;
			lattice[nLattice_local].molType2 = dump[(dihedral[i].atom3 - 1)].molType;
			lattice[nLattice_local].ix1 = dump[(dihedral[i].atom2 - 1)].ix;
			lattice[nLattice_local].ix2 = dump[(dihedral[i].atom3 - 1)].ix;
			lattice[nLattice_local].iy1 = dump[(dihedral[i].atom2 - 1)].iy;
			lattice[nLattice_local].iy2 = dump[(dihedral[i].atom3 - 1)].iy;
			lattice[nLattice_local].iz1 = dump[(dihedral[i].atom2 - 1)].iz;
			lattice[nLattice_local].iz2 = dump[(dihedral[i].atom3 - 1)].iz;
			lattice[nLattice_local].x1 = dump[(dihedral[i].atom2 - 1)].x;
			lattice[nLattice_local].x2 = dump[(dihedral[i].atom3 - 1)].x;
			lattice[nLattice_local].y1 = dump[(dihedral[i].atom2 - 1)].y;
			lattice[nLattice_local].y2 = dump[(dihedral[i].atom3 - 1)].y;
			lattice[nLattice_local].z1 = dump[(dihedral[i].atom2 - 1)].z;
			lattice[nLattice_local].z2 = dump[(dihedral[i].atom3 - 1)].z;
			lattice[nLattice_local].dihedralAngle = dihedral[i].dihedralAngle;

			nLattice_local++;
		}
	}

	printf("nLattice_local: %d\n", nLattice_local);

	*nLattice = nLattice_local;

	// printf("nLattice: %d\n", (*nLattice));
	free (isPendant);
	return lattice;
}

LATTICE_PROBABILITIES computeProbabilities (KMC_LATTICE *lattice, int nBackboneDihedrals)
{
	LATTICE_PROBABILITIES probabilities;

	probabilities.GmTGp = 0;
	probabilities.GmTT = 0;
	probabilities.GmTGm = 0;
	probabilities.GpTGm = 0;
	probabilities.GpTT = 0;
	probabilities.GpTGp = 0;
	probabilities.TTGp = 0;
	probabilities.TTGm = 0;
	probabilities.TTT = 0;
	probabilities.TGpGm = 0;
	probabilities.TGpGp = 0;
	probabilities.TGpT = 0;
	probabilities.GmGpT = 0;
	probabilities.GmGpGp = 0;
	probabilities.GmGpGm = 0;
	probabilities.GpGpT = 0;
	probabilities.GpGpGp = 0;
	probabilities.GpGpGm = 0;
	probabilities.TGmT = 0;
	probabilities.TGmGp = 0;
	probabilities.TGmGm = 0;
	probabilities.GmGmT = 0;
	probabilities.GmGmGp = 0;
	probabilities.GmGmGm = 0;
	probabilities.GpGmT = 0;
	probabilities.GpGmGp = 0;
	probabilities.GpGmGm = 0;

	for (int i = 0; i < nBackboneDihedrals; ++i)
	{
		if ((i == 0) || (i == nBackboneDihedrals - 1))
		{
			
		}
		else
		{
			// GmTGp
			if (
				(lattice[i-1].dihedralAngle < -38) && (lattice[i-1].dihedralAngle > -116) && 

				(((lattice[i].dihedralAngle < -117) && (lattice[i].dihedralAngle > -180)) || 
					((lattice[i].dihedralAngle < 180) && (lattice[i].dihedralAngle > 116))) && 

				(lattice[i+1].dihedralAngle < 115) && (lattice[i+1].dihedralAngle > 34)
				)
			{
				probabilities.GmTGp++;
			}
			// GmTT
			if (
				(lattice[i-1].dihedralAngle < -38) && (lattice[i-1].dihedralAngle > -116) &&

				(((lattice[i].dihedralAngle < -117) && (lattice[i].dihedralAngle > -180)) || 
					((lattice[i].dihedralAngle < 180) && (lattice[i].dihedralAngle > 116))) && 

				(((lattice[i].dihedralAngle < -117) && (lattice[i].dihedralAngle > -180)) || 
					((lattice[i].dihedralAngle < 180) && (lattice[i].dihedralAngle > 116)))
				)
			{
				probabilities.GmTT++;
			}
			// GmTGm
			if (
				(lattice[i-1].dihedralAngle < -38) && (lattice[i-1].dihedralAngle > -116) &&

				(((lattice[i].dihedralAngle < -117) && (lattice[i].dihedralAngle > -180)) || 
					((lattice[i].dihedralAngle < 180) && (lattice[i].dihedralAngle > 116))) &&

				(lattice[i-1].dihedralAngle < -38) && (lattice[i-1].dihedralAngle > -116)
				)
			{
				probabilities.GmTGm++;
			}

			// GpTGm
			if (
				(lattice[i+1].dihedralAngle < 115) && (lattice[i+1].dihedralAngle > 34) &&

				(((lattice[i].dihedralAngle < -117) && (lattice[i].dihedralAngle > -180)) || 
					((lattice[i].dihedralAngle < 180) && (lattice[i].dihedralAngle > 116))) &&

				(lattice[i-1].dihedralAngle < -38) && (lattice[i-1].dihedralAngle > -116)
				)
			{
				probabilities.GpTGm++;
			}
			// GpTT
			if (
				(lattice[i+1].dihedralAngle < 115) && (lattice[i+1].dihedralAngle > 34) &&

				(((lattice[i].dihedralAngle < -117) && (lattice[i].dihedralAngle > -180)) || 
					((lattice[i].dihedralAngle < 180) && (lattice[i].dihedralAngle > 116))) &&

				(((lattice[i].dihedralAngle < -117) && (lattice[i].dihedralAngle > -180)) || 
					((lattice[i].dihedralAngle < 180) && (lattice[i].dihedralAngle > 116)))
				)
			{
				probabilities.GpTT++;
			}
			// GpTGp
			if (
				(lattice[i+1].dihedralAngle < 115) && (lattice[i+1].dihedralAngle > 34) &&

				(((lattice[i].dihedralAngle < -117) && (lattice[i].dihedralAngle > -180)) || 
					((lattice[i].dihedralAngle < 180) && (lattice[i].dihedralAngle > 116))) &&

				(lattice[i+1].dihedralAngle < 115) && (lattice[i+1].dihedralAngle > 34)
				)
			{
				probabilities.GpTGp++;
			}

			// TTGp
			if (
				(((lattice[i].dihedralAngle < -117) && (lattice[i].dihedralAngle > -180)) || 
					((lattice[i].dihedralAngle < 180) && (lattice[i].dihedralAngle > 116))) &&

				(((lattice[i].dihedralAngle < -117) && (lattice[i].dihedralAngle > -180)) || 
					((lattice[i].dihedralAngle < 180) && (lattice[i].dihedralAngle > 116))) &&

				(lattice[i+1].dihedralAngle < 115) && (lattice[i+1].dihedralAngle > 34)
				)
			{
				probabilities.TTGp++;
			}
			// TTGm
			if (
				(((lattice[i].dihedralAngle < -117) && (lattice[i].dihedralAngle > -180)) || 
					((lattice[i].dihedralAngle < 180) && (lattice[i].dihedralAngle > 116))) &&

				(((lattice[i].dihedralAngle < -117) && (lattice[i].dihedralAngle > -180)) || 
					((lattice[i].dihedralAngle < 180) && (lattice[i].dihedralAngle > 116))) &&

				(lattice[i-1].dihedralAngle < -38) && (lattice[i-1].dihedralAngle > -116)
				)
			{
				probabilities.TTGm++;
			}
			// TTT
			if (
				(((lattice[i].dihedralAngle < -117) && (lattice[i].dihedralAngle > -180)) || 
					((lattice[i].dihedralAngle < 180) && (lattice[i].dihedralAngle > 116))) &&

				(((lattice[i].dihedralAngle < -117) && (lattice[i].dihedralAngle > -180)) || 
					((lattice[i].dihedralAngle < 180) && (lattice[i].dihedralAngle > 116))) &&

				(((lattice[i].dihedralAngle < -117) && (lattice[i].dihedralAngle > -180)) || 
					((lattice[i].dihedralAngle < 180) && (lattice[i].dihedralAngle > 116)))
				)
			{
				probabilities.TTT++;
			}

			// TGpGm
			if (
				(((lattice[i].dihedralAngle < -117) && (lattice[i].dihedralAngle > -180)) || 
					((lattice[i].dihedralAngle < 180) && (lattice[i].dihedralAngle > 116))) &&

				(lattice[i+1].dihedralAngle < 115) && (lattice[i+1].dihedralAngle > 34) &&

				(lattice[i-1].dihedralAngle < -38) && (lattice[i-1].dihedralAngle > -116)
				)
			{
				probabilities.TGpGm++;
			}
			// TGpGp
			if (
				(((lattice[i].dihedralAngle < -117) && (lattice[i].dihedralAngle > -180)) || 
					((lattice[i].dihedralAngle < 180) && (lattice[i].dihedralAngle > 116))) &&

				(lattice[i+1].dihedralAngle < 115) && (lattice[i+1].dihedralAngle > 34) &&

				(lattice[i+1].dihedralAngle < 115) && (lattice[i+1].dihedralAngle > 34)
				)
			{
				probabilities.TGpGp++;
			}
			// TGpT
			if (
				(((lattice[i].dihedralAngle < -117) && (lattice[i].dihedralAngle > -180)) || 
					((lattice[i].dihedralAngle < 180) && (lattice[i].dihedralAngle > 116))) &&

				(lattice[i+1].dihedralAngle < 115) && (lattice[i+1].dihedralAngle > 34) &&

				(((lattice[i].dihedralAngle < -117) && (lattice[i].dihedralAngle > -180)) || 
					((lattice[i].dihedralAngle < 180) && (lattice[i].dihedralAngle > 116)))
				)
			{
				probabilities.TGpT++;
			}

			// GmGpT
			if (
				(lattice[i-1].dihedralAngle < -38) && (lattice[i-1].dihedralAngle > -116) &&

				(lattice[i+1].dihedralAngle < 115) && (lattice[i+1].dihedralAngle > 34) &&

				(((lattice[i].dihedralAngle < -117) && (lattice[i].dihedralAngle > -180)) || 
					((lattice[i].dihedralAngle < 180) && (lattice[i].dihedralAngle > 116)))
				)
			{
				probabilities.GmGpT++;
			}
			// GmGpGp
			if (
				(lattice[i-1].dihedralAngle < -38) && (lattice[i-1].dihedralAngle > -116) &&

				(lattice[i+1].dihedralAngle < 115) && (lattice[i+1].dihedralAngle > 34) &&

				(lattice[i+1].dihedralAngle < 115) && (lattice[i+1].dihedralAngle > 34)
				)
			{
				probabilities.GmGpGp++;
			}
			// GmGpGm
			if (
				(lattice[i-1].dihedralAngle < -38) && (lattice[i-1].dihedralAngle > -116) &&

				(lattice[i+1].dihedralAngle < 115) && (lattice[i+1].dihedralAngle > 34) &&

				(lattice[i-1].dihedralAngle < -38) && (lattice[i-1].dihedralAngle > -116)
				)
			{
				probabilities.GmGpGm++;
			}

			// GpGpT
			if (
				(lattice[i+1].dihedralAngle < 115) && (lattice[i+1].dihedralAngle > 34) &&

				(lattice[i+1].dihedralAngle < 115) && (lattice[i+1].dihedralAngle > 34) &&

				(((lattice[i].dihedralAngle < -117) && (lattice[i].dihedralAngle > -180)) || 
					((lattice[i].dihedralAngle < 180) && (lattice[i].dihedralAngle > 116)))
				)
			{
				probabilities.GpGpT++;
			}
			// GpGpGp
			if (
				(lattice[i+1].dihedralAngle < 115) && (lattice[i+1].dihedralAngle > 34) &&

				(lattice[i+1].dihedralAngle < 115) && (lattice[i+1].dihedralAngle > 34) &&

				(lattice[i+1].dihedralAngle < 115) && (lattice[i+1].dihedralAngle > 34)
				)
			{
				probabilities.GpGpGp++;
			}
			// GpGpGm
			if (
				(lattice[i+1].dihedralAngle < 115) && (lattice[i+1].dihedralAngle > 34) &&

				(lattice[i+1].dihedralAngle < 115) && (lattice[i+1].dihedralAngle > 34) &&

				(lattice[i-1].dihedralAngle < -38) && (lattice[i-1].dihedralAngle > -116)
				)
			{
				probabilities.GpGpGm++;
			}

			// TGmT
			if (
				(((lattice[i].dihedralAngle < -117) && (lattice[i].dihedralAngle > -180)) || 
					((lattice[i].dihedralAngle < 180) && (lattice[i].dihedralAngle > 116))) &&

				(lattice[i-1].dihedralAngle < -38) && (lattice[i-1].dihedralAngle > -116) &&

				(((lattice[i].dihedralAngle < -117) && (lattice[i].dihedralAngle > -180)) || 
					((lattice[i].dihedralAngle < 180) && (lattice[i].dihedralAngle > 116)))
				)
			{
				probabilities.TGmT++;
			}
			// TGmGp
			if (
				(((lattice[i].dihedralAngle < -117) && (lattice[i].dihedralAngle > -180)) || 
					((lattice[i].dihedralAngle < 180) && (lattice[i].dihedralAngle > 116))) &&

				(lattice[i-1].dihedralAngle < -38) && (lattice[i-1].dihedralAngle > -116) &&

				(lattice[i+1].dihedralAngle < 115) && (lattice[i+1].dihedralAngle > 34)
				)
			{
				probabilities.TGmGp++;
			}
			// TGmGm
			if (
				(((lattice[i].dihedralAngle < -117) && (lattice[i].dihedralAngle > -180)) || 
					((lattice[i].dihedralAngle < 180) && (lattice[i].dihedralAngle > 116))) &&

				(lattice[i-1].dihedralAngle < -38) && (lattice[i-1].dihedralAngle > -116) &&

				(lattice[i-1].dihedralAngle < -38) && (lattice[i-1].dihedralAngle > -116)
				)
			{
				probabilities.TGmGm++;
			}

			// GmGmGp
			if (
				(lattice[i-1].dihedralAngle < -38) && (lattice[i-1].dihedralAngle > -116) &&

				(lattice[i-1].dihedralAngle < -38) && (lattice[i-1].dihedralAngle > -116) &&

				(lattice[i+1].dihedralAngle < 115) && (lattice[i+1].dihedralAngle > 34)
				)
			{
				probabilities.GmGmGp++;
			}
			// GmGmGm
			if (
				(lattice[i-1].dihedralAngle < -38) && (lattice[i-1].dihedralAngle > -116) &&

				(lattice[i-1].dihedralAngle < -38) && (lattice[i-1].dihedralAngle > -116) &&

				(lattice[i-1].dihedralAngle < -38) && (lattice[i-1].dihedralAngle > -116) 
				)
			{
				probabilities.GmGmGm++;
			}
			// GmGmT
			if (
				(lattice[i-1].dihedralAngle < -38) && (lattice[i-1].dihedralAngle > -116) &&

				(lattice[i-1].dihedralAngle < -38) && (lattice[i-1].dihedralAngle > -116) &&

				(((lattice[i].dihedralAngle < -117) && (lattice[i].dihedralAngle > -180)) || 
					((lattice[i].dihedralAngle < 180) && (lattice[i].dihedralAngle > 116)))
				)
			{
				probabilities.GmGmT++;
			}

			// GpGmT
			if (
				(lattice[i+1].dihedralAngle < 115) && (lattice[i+1].dihedralAngle > 34) &&

				(lattice[i-1].dihedralAngle < -38) && (lattice[i-1].dihedralAngle > -116) &&

				(((lattice[i].dihedralAngle < -117) && (lattice[i].dihedralAngle > -180)) || 
					((lattice[i].dihedralAngle < 180) && (lattice[i].dihedralAngle > 116)))
				)
			{
				probabilities.GpGmT++;
			}
			// GpGmGm
			if (
				(lattice[i+1].dihedralAngle < 115) && (lattice[i+1].dihedralAngle > 34) &&

				(lattice[i-1].dihedralAngle < -38) && (lattice[i-1].dihedralAngle > -116) &&

				(lattice[i-1].dihedralAngle < -38) && (lattice[i-1].dihedralAngle > -116)
				)
			{
				probabilities.GpGmGm++;
			}
			// GpGmGp
			if (
				(lattice[i+1].dihedralAngle < 115) && (lattice[i+1].dihedralAngle > 34) &&

				(lattice[i-1].dihedralAngle < -38) && (lattice[i-1].dihedralAngle > -116) &&

				(lattice[i+1].dihedralAngle < 115) && (lattice[i+1].dihedralAngle > 34)
				)
			{
				probabilities.GpGmGp++;
			}

		}
	}

	return probabilities;
}

int chiralChecking (int dihedralAngle)
{
	int chirality;

	if (dihedralAngle > -116 && dihedralAngle < -38)
	{
		chirality = -1;
	}
	else if (dihedralAngle > 34 && dihedralAngle < 115)
	{
		chirality = 1;
	}
	else if ((dihedralAngle > -180 && dihedralAngle < -117) || (dihedralAngle > 116 && dihedralAngle < 180))
	{
		chirality = 0;
	}
	else
	{
		chirality = 5;
	}

	return chirality;
}

KMC_JUMPS computeKCMJumps (KMC_LATTICE *lattice, KMC_LATTICE *lattice_tminus1, int nLattice)
{
	KMC_JUMPS possibleJumps;
	// fprintf(stdout, "%s\n", "Entered KMC_JUMP function");

	possibleJumps.t1Gp = 0;
	possibleJumps.t1Gm = 0;
	possibleJumps.t1T = 0;
	possibleJumps.t2Gp = 0;
	possibleJumps.t2Gm = 0;
	possibleJumps.t2T = 0;
	possibleJumps.t3Gp = 0;
	possibleJumps.t3Gm = 0;
	possibleJumps.t3T = 0;
	possibleJumps.t4Gp = 0;
	possibleJumps.t4Gm = 0;
	possibleJumps.t4T = 0;
	possibleJumps.t5Gp = 0;
	possibleJumps.t5Gm = 0;
	possibleJumps.t5T = 0;
	possibleJumps.t6Gp = 0;
	possibleJumps.t6Gm = 0;
	possibleJumps.t6T = 0;
	possibleJumps.t7Gp = 0;
	possibleJumps.t7Gm = 0;
	possibleJumps.t7T = 0;
	possibleJumps.t8Gp = 0;
	possibleJumps.t8Gm = 0;
	possibleJumps.t8T = 0;
	possibleJumps.t9Gp = 0;
	possibleJumps.t9Gm = 0;
	possibleJumps.t9T = 0;

	possibleJumps.Gp1Gp = 0;
	possibleJumps.Gp1Gm = 0;
	possibleJumps.Gp1T = 0;
	possibleJumps.Gp2Gp = 0;
	possibleJumps.Gp2Gm = 0;
	possibleJumps.Gp2T = 0;
	possibleJumps.Gp3Gp = 0;
	possibleJumps.Gp3Gm = 0;
	possibleJumps.Gp3T = 0;
	possibleJumps.Gp4Gp = 0;
	possibleJumps.Gp4Gm = 0;
	possibleJumps.Gp4T = 0;
	possibleJumps.Gp5Gp = 0;
	possibleJumps.Gp5Gm = 0;
	possibleJumps.Gp5T = 0;
	possibleJumps.Gp6Gp = 0;
	possibleJumps.Gp6Gm = 0;
	possibleJumps.Gp6T = 0;
	possibleJumps.Gp7Gp = 0;
	possibleJumps.Gp7Gm = 0;
	possibleJumps.Gp7T = 0;
	possibleJumps.Gp8Gp = 0;
	possibleJumps.Gp8Gm = 0;
	possibleJumps.Gp8T = 0;
	possibleJumps.Gp9Gp = 0;
	possibleJumps.Gp9Gm = 0;
	possibleJumps.Gp9T = 0;

	possibleJumps.Gm1Gp = 0;
	possibleJumps.Gm1Gm = 0;
	possibleJumps.Gm1T = 0;
	possibleJumps.Gm2Gp = 0;
	possibleJumps.Gm2Gm = 0;
	possibleJumps.Gm2T = 0;
	possibleJumps.Gm3Gp = 0;
	possibleJumps.Gm3Gm = 0;
	possibleJumps.Gm3T = 0;
	possibleJumps.Gm4Gp = 0;
	possibleJumps.Gm4Gm = 0;
	possibleJumps.Gm4T = 0;
	possibleJumps.Gm5Gp = 0;
	possibleJumps.Gm5Gm = 0;
	possibleJumps.Gm5T = 0;
	possibleJumps.Gm6Gp = 0;
	possibleJumps.Gm6Gm = 0;
	possibleJumps.Gm6T = 0;
	possibleJumps.Gm7Gp = 0;
	possibleJumps.Gm7Gm = 0;
	possibleJumps.Gm7T = 0;
	possibleJumps.Gm8Gp = 0;
	possibleJumps.Gm8Gm = 0;
	possibleJumps.Gm8T = 0;
	possibleJumps.Gm9Gp = 0;
	possibleJumps.Gm9Gm = 0;
	possibleJumps.Gm9T = 0;

	for (int i = 1; i < nLattice - 1; ++i)
	{
		if (
			(chiralChecking (lattice_tminus1[i-1].dihedralAngle) == -1) &&
			(chiralChecking (lattice_tminus1[i].dihedralAngle) == 0) &&
			(chiralChecking (lattice_tminus1[i+1].dihedralAngle) == -1) &&
			(chiralChecking (lattice[i-1].dihedralAngle) == -1) &&
			(chiralChecking (lattice[i].dihedralAngle) == 1) &&
			(chiralChecking (lattice[i+1].dihedralAngle) == -1)
			)
		{
			possibleJumps.t1Gp++;
		}
		else if (
			(chiralChecking (lattice_tminus1[i-1].dihedralAngle) == -1) &&
			(chiralChecking (lattice_tminus1[i].dihedralAngle) == 0) &&
			(chiralChecking (lattice_tminus1[i+1].dihedralAngle) == -1) &&
			(chiralChecking (lattice[i-1].dihedralAngle) == -1) &&
			(chiralChecking (lattice[i].dihedralAngle) == -1) &&
			(chiralChecking (lattice[i+1].dihedralAngle) == -1)
			)
		{
			possibleJumps.t1Gm++;
		}
		else if (
			(chiralChecking (lattice_tminus1[i-1].dihedralAngle) == -1) &&
			(chiralChecking (lattice_tminus1[i].dihedralAngle) == 0) &&
			(chiralChecking (lattice_tminus1[i+1].dihedralAngle) == -1) &&
			(chiralChecking (lattice[i-1].dihedralAngle) == -1) &&
			(chiralChecking (lattice[i].dihedralAngle) == 0) &&
			(chiralChecking (lattice[i+1].dihedralAngle) == -1)
			)
		{
			possibleJumps.t1T++;
		}
		else if (
			(chiralChecking (lattice_tminus1[i-1].dihedralAngle) == -1) &&
			(chiralChecking (lattice_tminus1[i].dihedralAngle) == 0) &&
			(chiralChecking (lattice_tminus1[i+1].dihedralAngle) == 1) &&
			(chiralChecking (lattice[i-1].dihedralAngle) == -1) &&
			(chiralChecking (lattice[i].dihedralAngle) == 1) &&
			(chiralChecking (lattice[i+1].dihedralAngle) == 1)
			)
		{
			possibleJumps.t2Gp++;
		}
		else if (
			(chiralChecking (lattice_tminus1[i-1].dihedralAngle) == -1) &&
			(chiralChecking (lattice_tminus1[i].dihedralAngle) == 0) &&
			(chiralChecking (lattice_tminus1[i+1].dihedralAngle) == 1) &&
			(chiralChecking (lattice[i-1].dihedralAngle) == -1) &&
			(chiralChecking (lattice[i].dihedralAngle) == -1) &&
			(chiralChecking (lattice[i+1].dihedralAngle) == 1)
			)
		{
			possibleJumps.t2Gm++;
		}
		else if (
			(chiralChecking (lattice_tminus1[i-1].dihedralAngle) == -1) &&
			(chiralChecking (lattice_tminus1[i].dihedralAngle) == 0) &&
			(chiralChecking (lattice_tminus1[i+1].dihedralAngle) == 1) &&
			(chiralChecking (lattice[i-1].dihedralAngle) == -1) &&
			(chiralChecking (lattice[i].dihedralAngle) == 0) &&
			(chiralChecking (lattice[i+1].dihedralAngle) == 1)
			)
		{
			possibleJumps.t2T++;
		}
		else if (
			(chiralChecking (lattice_tminus1[i-1].dihedralAngle) == -1) &&
			(chiralChecking (lattice_tminus1[i].dihedralAngle) == 0) &&
			(chiralChecking (lattice_tminus1[i+1].dihedralAngle) == 0) &&
			(chiralChecking (lattice[i-1].dihedralAngle) == -1) &&
			(chiralChecking (lattice[i].dihedralAngle) == 1) &&
			(chiralChecking (lattice[i+1].dihedralAngle) == 0)
			)
		{
			possibleJumps.t3Gp++;
		}
		else if (
			(chiralChecking (lattice_tminus1[i-1].dihedralAngle) == -1) &&
			(chiralChecking (lattice_tminus1[i].dihedralAngle) == 0) &&
			(chiralChecking (lattice_tminus1[i+1].dihedralAngle) == 0) &&
			(chiralChecking (lattice[i-1].dihedralAngle) == -1) &&
			(chiralChecking (lattice[i].dihedralAngle) == -1) &&
			(chiralChecking (lattice[i+1].dihedralAngle) == 0)
			)
		{
			possibleJumps.t3Gm++;
		}
		else if (
			(chiralChecking (lattice_tminus1[i-1].dihedralAngle) == -1) &&
			(chiralChecking (lattice_tminus1[i].dihedralAngle) == 0) &&
			(chiralChecking (lattice_tminus1[i+1].dihedralAngle) == 0) &&
			(chiralChecking (lattice[i-1].dihedralAngle) == -1) &&
			(chiralChecking (lattice[i].dihedralAngle) == 0) &&
			(chiralChecking (lattice[i+1].dihedralAngle) == 0)
			)
		{
			possibleJumps.t3T++;
		}
		else if (
			(chiralChecking (lattice_tminus1[i-1].dihedralAngle) == 1) &&
			(chiralChecking (lattice_tminus1[i].dihedralAngle) == 0) &&
			(chiralChecking (lattice_tminus1[i+1].dihedralAngle) == -1) &&
			(chiralChecking (lattice[i-1].dihedralAngle) == 1) &&
			(chiralChecking (lattice[i].dihedralAngle) == 1) &&
			(chiralChecking (lattice[i+1].dihedralAngle) == -1)
			)
		{
			possibleJumps.t4Gp++;
		}
		else if (
			(chiralChecking (lattice_tminus1[i-1].dihedralAngle) == 1) &&
			(chiralChecking (lattice_tminus1[i].dihedralAngle) == 0) &&
			(chiralChecking (lattice_tminus1[i+1].dihedralAngle) == -1) &&
			(chiralChecking (lattice[i-1].dihedralAngle) == 1) &&
			(chiralChecking (lattice[i].dihedralAngle) == -1) &&
			(chiralChecking (lattice[i+1].dihedralAngle) == -1)
			)
		{
			possibleJumps.t4Gm++;
		}
		else if (
			(chiralChecking (lattice_tminus1[i-1].dihedralAngle) == 1) &&
			(chiralChecking (lattice_tminus1[i].dihedralAngle) == 0) &&
			(chiralChecking (lattice_tminus1[i+1].dihedralAngle) == -1) &&
			(chiralChecking (lattice[i-1].dihedralAngle) == 1) &&
			(chiralChecking (lattice[i].dihedralAngle) == 0) &&
			(chiralChecking (lattice[i+1].dihedralAngle) == -1)
			)
		{
			possibleJumps.t4T++;
		}
		else if (
			(chiralChecking (lattice_tminus1[i-1].dihedralAngle) == 1) &&
			(chiralChecking (lattice_tminus1[i].dihedralAngle) == 0) &&
			(chiralChecking (lattice_tminus1[i+1].dihedralAngle) == 1) &&
			(chiralChecking (lattice[i-1].dihedralAngle) == 1) &&
			(chiralChecking (lattice[i].dihedralAngle) == 1) &&
			(chiralChecking (lattice[i+1].dihedralAngle) == 1)
			)
		{
			possibleJumps.t5Gp++;
		}
		else if (
			(chiralChecking (lattice_tminus1[i-1].dihedralAngle) == 1) &&
			(chiralChecking (lattice_tminus1[i].dihedralAngle) == 0) &&
			(chiralChecking (lattice_tminus1[i+1].dihedralAngle) == 1) &&
			(chiralChecking (lattice[i-1].dihedralAngle) == 1) &&
			(chiralChecking (lattice[i].dihedralAngle) == -1) &&
			(chiralChecking (lattice[i+1].dihedralAngle) == 1)
			)
		{
			possibleJumps.t5Gm++;
		}
		else if (
			(chiralChecking (lattice_tminus1[i-1].dihedralAngle) == 1) &&
			(chiralChecking (lattice_tminus1[i].dihedralAngle) == 0) &&
			(chiralChecking (lattice_tminus1[i+1].dihedralAngle) == 1) &&
			(chiralChecking (lattice[i-1].dihedralAngle) == 1) &&
			(chiralChecking (lattice[i].dihedralAngle) == 0) &&
			(chiralChecking (lattice[i+1].dihedralAngle) == 1)
			)
		{
			possibleJumps.t5T++;
		}
		else if (
			(chiralChecking (lattice_tminus1[i-1].dihedralAngle) == 1) &&
			(chiralChecking (lattice_tminus1[i].dihedralAngle) == 0) &&
			(chiralChecking (lattice_tminus1[i+1].dihedralAngle) == 0) &&
			(chiralChecking (lattice[i-1].dihedralAngle) == 1) &&
			(chiralChecking (lattice[i].dihedralAngle) == 1) &&
			(chiralChecking (lattice[i+1].dihedralAngle) == 0)
			)
		{
			possibleJumps.t6Gp++;
		}
		else if (
			(chiralChecking (lattice_tminus1[i-1].dihedralAngle) == 1) &&
			(chiralChecking (lattice_tminus1[i].dihedralAngle) == 0) &&
			(chiralChecking (lattice_tminus1[i+1].dihedralAngle) == 0) &&
			(chiralChecking (lattice[i-1].dihedralAngle) == 1) &&
			(chiralChecking (lattice[i].dihedralAngle) == -1) &&
			(chiralChecking (lattice[i+1].dihedralAngle) == 0)
			)
		{
			possibleJumps.t6Gm++;
		}
		else if (
			(chiralChecking (lattice_tminus1[i-1].dihedralAngle) == 1) &&
			(chiralChecking (lattice_tminus1[i].dihedralAngle) == 0) &&
			(chiralChecking (lattice_tminus1[i+1].dihedralAngle) == 0) &&
			(chiralChecking (lattice[i-1].dihedralAngle) == 1) &&
			(chiralChecking (lattice[i].dihedralAngle) == 0) &&
			(chiralChecking (lattice[i+1].dihedralAngle) == 0)
			)
		{
			possibleJumps.t6T++;
		}
		else if (
			(chiralChecking (lattice_tminus1[i-1].dihedralAngle) == 0) &&
			(chiralChecking (lattice_tminus1[i].dihedralAngle) == 0) &&
			(chiralChecking (lattice_tminus1[i+1].dihedralAngle) == -1) &&
			(chiralChecking (lattice[i-1].dihedralAngle) == 0) &&
			(chiralChecking (lattice[i].dihedralAngle) == 1) &&
			(chiralChecking (lattice[i+1].dihedralAngle) == -1)
			)
		{
			possibleJumps.t7Gp++;
		}
		else if (
			(chiralChecking (lattice_tminus1[i-1].dihedralAngle) == 0) &&
			(chiralChecking (lattice_tminus1[i].dihedralAngle) == 0) &&
			(chiralChecking (lattice_tminus1[i+1].dihedralAngle) == -1) &&
			(chiralChecking (lattice[i-1].dihedralAngle) == 0) &&
			(chiralChecking (lattice[i].dihedralAngle) == -1) &&
			(chiralChecking (lattice[i+1].dihedralAngle) == -1)
			)
		{
			possibleJumps.t7Gm++;
		}
		else if (
			(chiralChecking (lattice_tminus1[i-1].dihedralAngle) == 0) &&
			(chiralChecking (lattice_tminus1[i].dihedralAngle) == 0) &&
			(chiralChecking (lattice_tminus1[i+1].dihedralAngle) == -1) &&
			(chiralChecking (lattice[i-1].dihedralAngle) == 0) &&
			(chiralChecking (lattice[i].dihedralAngle) == 0) &&
			(chiralChecking (lattice[i+1].dihedralAngle) == -1)
			)
		{
			possibleJumps.t7T++;
		}
		else if (
			(chiralChecking (lattice_tminus1[i-1].dihedralAngle) == 0) &&
			(chiralChecking (lattice_tminus1[i].dihedralAngle) == 0) &&
			(chiralChecking (lattice_tminus1[i+1].dihedralAngle) == 1) &&
			(chiralChecking (lattice[i-1].dihedralAngle) == 0) &&
			(chiralChecking (lattice[i].dihedralAngle) == 1) &&
			(chiralChecking (lattice[i+1].dihedralAngle) == 1)
			)
		{
			possibleJumps.t8Gp++;
		}
		else if (
			(chiralChecking (lattice_tminus1[i-1].dihedralAngle) == 0) &&
			(chiralChecking (lattice_tminus1[i].dihedralAngle) == 0) &&
			(chiralChecking (lattice_tminus1[i+1].dihedralAngle) == 1) &&
			(chiralChecking (lattice[i-1].dihedralAngle) == 0) &&
			(chiralChecking (lattice[i].dihedralAngle) == -1) &&
			(chiralChecking (lattice[i+1].dihedralAngle) == 1)
			)
		{
			possibleJumps.t8Gm++;
		}
		else if (
			(chiralChecking (lattice_tminus1[i-1].dihedralAngle) == 0) &&
			(chiralChecking (lattice_tminus1[i].dihedralAngle) == 0) &&
			(chiralChecking (lattice_tminus1[i+1].dihedralAngle) == 1) &&
			(chiralChecking (lattice[i-1].dihedralAngle) == 0) &&
			(chiralChecking (lattice[i].dihedralAngle) == 0) &&
			(chiralChecking (lattice[i+1].dihedralAngle) == 1)
			)
		{
			possibleJumps.t8T++;
		}
		else if (
			(chiralChecking (lattice_tminus1[i-1].dihedralAngle) == 0) &&
			(chiralChecking (lattice_tminus1[i].dihedralAngle) == 0) &&
			(chiralChecking (lattice_tminus1[i+1].dihedralAngle) == 0) &&
			(chiralChecking (lattice[i-1].dihedralAngle) == 0) &&
			(chiralChecking (lattice[i].dihedralAngle) == 1) &&
			(chiralChecking (lattice[i+1].dihedralAngle) == 0)
			)
		{
			possibleJumps.t9Gp++;
		}
		else if (
			(chiralChecking (lattice_tminus1[i-1].dihedralAngle) == 0) &&
			(chiralChecking (lattice_tminus1[i].dihedralAngle) == 0) &&
			(chiralChecking (lattice_tminus1[i+1].dihedralAngle) == 0) &&
			(chiralChecking (lattice[i-1].dihedralAngle) == 0) &&
			(chiralChecking (lattice[i].dihedralAngle) == -1) &&
			(chiralChecking (lattice[i+1].dihedralAngle) == 0)
			)
		{
			possibleJumps.t9Gm++;
		}
		else if (
			(chiralChecking (lattice_tminus1[i-1].dihedralAngle) == 0) &&
			(chiralChecking (lattice_tminus1[i].dihedralAngle) == 0) &&
			(chiralChecking (lattice_tminus1[i+1].dihedralAngle) == 0) &&
			(chiralChecking (lattice[i-1].dihedralAngle) == 0) &&
			(chiralChecking (lattice[i].dihedralAngle) == 0) &&
			(chiralChecking (lattice[i+1].dihedralAngle) == 0)
			)
		{
			possibleJumps.t9T++;
		}
		else if (
			(chiralChecking (lattice_tminus1[i-1].dihedralAngle) == -1) &&
			(chiralChecking (lattice_tminus1[i].dihedralAngle) == -1) &&
			(chiralChecking (lattice_tminus1[i+1].dihedralAngle) == -1) &&
			(chiralChecking (lattice[i-1].dihedralAngle) == -1) &&
			(chiralChecking (lattice[i].dihedralAngle) == 1) &&
			(chiralChecking (lattice[i+1].dihedralAngle) == -1)
			)
		{
			possibleJumps.Gm1Gp++;
		}
		else if (
			(chiralChecking (lattice_tminus1[i-1].dihedralAngle) == -1) &&
			(chiralChecking (lattice_tminus1[i].dihedralAngle) == -1) &&
			(chiralChecking (lattice_tminus1[i+1].dihedralAngle) == -1) &&
			(chiralChecking (lattice[i-1].dihedralAngle) == -1) &&
			(chiralChecking (lattice[i].dihedralAngle) == -1) &&
			(chiralChecking (lattice[i+1].dihedralAngle) == -1)
			)
		{
			possibleJumps.Gm1Gm++;
		}
		else if (
			(chiralChecking (lattice_tminus1[i-1].dihedralAngle) == -1) &&
			(chiralChecking (lattice_tminus1[i].dihedralAngle) == -1) &&
			(chiralChecking (lattice_tminus1[i+1].dihedralAngle) == -1) &&
			(chiralChecking (lattice[i-1].dihedralAngle) == -1) &&
			(chiralChecking (lattice[i].dihedralAngle) == 0) &&
			(chiralChecking (lattice[i+1].dihedralAngle) == -1)
			)
		{
			possibleJumps.Gm1T++;
		}
		else if (
			(chiralChecking (lattice_tminus1[i-1].dihedralAngle) == -1) &&
			(chiralChecking (lattice_tminus1[i].dihedralAngle) == -1) &&
			(chiralChecking (lattice_tminus1[i+1].dihedralAngle) == 1) &&
			(chiralChecking (lattice[i-1].dihedralAngle) == -1) &&
			(chiralChecking (lattice[i].dihedralAngle) == 1) &&
			(chiralChecking (lattice[i+1].dihedralAngle) == 1)
			)
		{
			possibleJumps.Gm2Gp++;
		}
		else if (
			(chiralChecking (lattice_tminus1[i-1].dihedralAngle) == -1) &&
			(chiralChecking (lattice_tminus1[i].dihedralAngle) == -1) &&
			(chiralChecking (lattice_tminus1[i+1].dihedralAngle) == 1) &&
			(chiralChecking (lattice[i-1].dihedralAngle) == -1) &&
			(chiralChecking (lattice[i].dihedralAngle) == -1) &&
			(chiralChecking (lattice[i+1].dihedralAngle) == 1)
			)
		{
			possibleJumps.Gm2Gm++;
		}
		else if (
			(chiralChecking (lattice_tminus1[i-1].dihedralAngle) == -1) &&
			(chiralChecking (lattice_tminus1[i].dihedralAngle) == -1) &&
			(chiralChecking (lattice_tminus1[i+1].dihedralAngle) == 1) &&
			(chiralChecking (lattice[i-1].dihedralAngle) == -1) &&
			(chiralChecking (lattice[i].dihedralAngle) == 0) &&
			(chiralChecking (lattice[i+1].dihedralAngle) == 1)
			)
		{
			possibleJumps.Gm2T++;
		}
		else if (
			(chiralChecking (lattice_tminus1[i-1].dihedralAngle) == -1) &&
			(chiralChecking (lattice_tminus1[i].dihedralAngle) == -1) &&
			(chiralChecking (lattice_tminus1[i+1].dihedralAngle) == 0) &&
			(chiralChecking (lattice[i-1].dihedralAngle) == -1) &&
			(chiralChecking (lattice[i].dihedralAngle) == 1) &&
			(chiralChecking (lattice[i+1].dihedralAngle) == 0)
			)
		{
			possibleJumps.Gm3Gp++;
		}
		else if (
			(chiralChecking (lattice_tminus1[i-1].dihedralAngle) == -1) &&
			(chiralChecking (lattice_tminus1[i].dihedralAngle) == -1) &&
			(chiralChecking (lattice_tminus1[i+1].dihedralAngle) == 0) &&
			(chiralChecking (lattice[i-1].dihedralAngle) == -1) &&
			(chiralChecking (lattice[i].dihedralAngle) == -1) &&
			(chiralChecking (lattice[i+1].dihedralAngle) == 0)
			)
		{
			possibleJumps.Gm3Gm++;
		}
		else if (
			(chiralChecking (lattice_tminus1[i-1].dihedralAngle) == -1) &&
			(chiralChecking (lattice_tminus1[i].dihedralAngle) == -1) &&
			(chiralChecking (lattice_tminus1[i+1].dihedralAngle) == 0) &&
			(chiralChecking (lattice[i-1].dihedralAngle) == -1) &&
			(chiralChecking (lattice[i].dihedralAngle) == 0) &&
			(chiralChecking (lattice[i+1].dihedralAngle) == 0)
			)
		{
			possibleJumps.Gm3T++;
		}
		else if (
			(chiralChecking (lattice_tminus1[i-1].dihedralAngle) == 1) &&
			(chiralChecking (lattice_tminus1[i].dihedralAngle) == -1) &&
			(chiralChecking (lattice_tminus1[i+1].dihedralAngle) == -1) &&
			(chiralChecking (lattice[i-1].dihedralAngle) == 1) &&
			(chiralChecking (lattice[i].dihedralAngle) == 1) &&
			(chiralChecking (lattice[i+1].dihedralAngle) == -1)
			)
		{
			possibleJumps.Gm4Gp++;
		}
		else if (
			(chiralChecking (lattice_tminus1[i-1].dihedralAngle) == 1) &&
			(chiralChecking (lattice_tminus1[i].dihedralAngle) == -1) &&
			(chiralChecking (lattice_tminus1[i+1].dihedralAngle) == -1) &&
			(chiralChecking (lattice[i-1].dihedralAngle) == 1) &&
			(chiralChecking (lattice[i].dihedralAngle) == -1) &&
			(chiralChecking (lattice[i+1].dihedralAngle) == -1)
			)
		{
			possibleJumps.Gm4Gm++;
		}
		else if (
			(chiralChecking (lattice_tminus1[i-1].dihedralAngle) == 1) &&
			(chiralChecking (lattice_tminus1[i].dihedralAngle) == -1) &&
			(chiralChecking (lattice_tminus1[i+1].dihedralAngle) == -1) &&
			(chiralChecking (lattice[i-1].dihedralAngle) == 1) &&
			(chiralChecking (lattice[i].dihedralAngle) == 0) &&
			(chiralChecking (lattice[i+1].dihedralAngle) == -1)
			)
		{
			possibleJumps.Gm4T++;
		}
		else if (
			(chiralChecking (lattice_tminus1[i-1].dihedralAngle) == 1) &&
			(chiralChecking (lattice_tminus1[i].dihedralAngle) == -1) &&
			(chiralChecking (lattice_tminus1[i+1].dihedralAngle) == 1) &&
			(chiralChecking (lattice[i-1].dihedralAngle) == 1) &&
			(chiralChecking (lattice[i].dihedralAngle) == 1) &&
			(chiralChecking (lattice[i+1].dihedralAngle) == 1)
			)
		{
			possibleJumps.Gm5Gp++;
		}
		else if (
			(chiralChecking (lattice_tminus1[i-1].dihedralAngle) == 1) &&
			(chiralChecking (lattice_tminus1[i].dihedralAngle) == -1) &&
			(chiralChecking (lattice_tminus1[i+1].dihedralAngle) == 1) &&
			(chiralChecking (lattice[i-1].dihedralAngle) == 1) &&
			(chiralChecking (lattice[i].dihedralAngle) == -1) &&
			(chiralChecking (lattice[i+1].dihedralAngle) == 1)
			)
		{
			possibleJumps.Gm5Gm++;
		}
		else if (
			(chiralChecking (lattice_tminus1[i-1].dihedralAngle) == 1) &&
			(chiralChecking (lattice_tminus1[i].dihedralAngle) == -1) &&
			(chiralChecking (lattice_tminus1[i+1].dihedralAngle) == 1) &&
			(chiralChecking (lattice[i-1].dihedralAngle) == 1) &&
			(chiralChecking (lattice[i].dihedralAngle) == 0) &&
			(chiralChecking (lattice[i+1].dihedralAngle) == 1)
			)
		{
			possibleJumps.Gm5T++;
		}
		else if (
			(chiralChecking (lattice_tminus1[i-1].dihedralAngle) == 1) &&
			(chiralChecking (lattice_tminus1[i].dihedralAngle) == -1) &&
			(chiralChecking (lattice_tminus1[i+1].dihedralAngle) == 0) &&
			(chiralChecking (lattice[i-1].dihedralAngle) == 1) &&
			(chiralChecking (lattice[i].dihedralAngle) == 1) &&
			(chiralChecking (lattice[i+1].dihedralAngle) == 0)
			)
		{
			possibleJumps.Gm6Gp++;
		}
		else if (
			(chiralChecking (lattice_tminus1[i-1].dihedralAngle) == 1) &&
			(chiralChecking (lattice_tminus1[i].dihedralAngle) == -1) &&
			(chiralChecking (lattice_tminus1[i+1].dihedralAngle) == 0) &&
			(chiralChecking (lattice[i-1].dihedralAngle) == 1) &&
			(chiralChecking (lattice[i].dihedralAngle) == -1) &&
			(chiralChecking (lattice[i+1].dihedralAngle) == 0)
			)
		{
			possibleJumps.Gm6Gm++;
		}
		else if (
			(chiralChecking (lattice_tminus1[i-1].dihedralAngle) == 1) &&
			(chiralChecking (lattice_tminus1[i].dihedralAngle) == -1) &&
			(chiralChecking (lattice_tminus1[i+1].dihedralAngle) == 0) &&
			(chiralChecking (lattice[i-1].dihedralAngle) == 1) &&
			(chiralChecking (lattice[i].dihedralAngle) == 0) &&
			(chiralChecking (lattice[i+1].dihedralAngle) == 0)
			)
		{
			possibleJumps.Gm6T++;
		}
		else if (
			(chiralChecking (lattice_tminus1[i-1].dihedralAngle) == 0) &&
			(chiralChecking (lattice_tminus1[i].dihedralAngle) == -1) &&
			(chiralChecking (lattice_tminus1[i+1].dihedralAngle) == -1) &&
			(chiralChecking (lattice[i-1].dihedralAngle) == 0) &&
			(chiralChecking (lattice[i].dihedralAngle) == 1) &&
			(chiralChecking (lattice[i+1].dihedralAngle) == -1)
			)
		{
			possibleJumps.Gm7Gp++;
		}
		else if (
			(chiralChecking (lattice_tminus1[i-1].dihedralAngle) == 0) &&
			(chiralChecking (lattice_tminus1[i].dihedralAngle) == -1) &&
			(chiralChecking (lattice_tminus1[i+1].dihedralAngle) == -1) &&
			(chiralChecking (lattice[i-1].dihedralAngle) == 0) &&
			(chiralChecking (lattice[i].dihedralAngle) == -1) &&
			(chiralChecking (lattice[i+1].dihedralAngle) == -1)
			)
		{
			possibleJumps.Gm7Gm++;
		}
		else if (
			(chiralChecking (lattice_tminus1[i-1].dihedralAngle) == 0) &&
			(chiralChecking (lattice_tminus1[i].dihedralAngle) == -1) &&
			(chiralChecking (lattice_tminus1[i+1].dihedralAngle) == -1) &&
			(chiralChecking (lattice[i-1].dihedralAngle) == 0) &&
			(chiralChecking (lattice[i].dihedralAngle) == 0) &&
			(chiralChecking (lattice[i+1].dihedralAngle) == -1)
			)
		{
			possibleJumps.Gm7T++;
		}
		else if (
			(chiralChecking (lattice_tminus1[i-1].dihedralAngle) == 0) &&
			(chiralChecking (lattice_tminus1[i].dihedralAngle) == -1) &&
			(chiralChecking (lattice_tminus1[i+1].dihedralAngle) == 1) &&
			(chiralChecking (lattice[i-1].dihedralAngle) == 0) &&
			(chiralChecking (lattice[i].dihedralAngle) == 1) &&
			(chiralChecking (lattice[i+1].dihedralAngle) == 1)
			)
		{
			possibleJumps.Gm8Gp++;
		}
		else if (
			(chiralChecking (lattice_tminus1[i-1].dihedralAngle) == 0) &&
			(chiralChecking (lattice_tminus1[i].dihedralAngle) == -1) &&
			(chiralChecking (lattice_tminus1[i+1].dihedralAngle) == 1) &&
			(chiralChecking (lattice[i-1].dihedralAngle) == 0) &&
			(chiralChecking (lattice[i].dihedralAngle) == -1) &&
			(chiralChecking (lattice[i+1].dihedralAngle) == 1)
			)
		{
			possibleJumps.Gm8Gm++;
		}
		else if (
			(chiralChecking (lattice_tminus1[i-1].dihedralAngle) == 0) &&
			(chiralChecking (lattice_tminus1[i].dihedralAngle) == -1) &&
			(chiralChecking (lattice_tminus1[i+1].dihedralAngle) == 1) &&
			(chiralChecking (lattice[i-1].dihedralAngle) == 0) &&
			(chiralChecking (lattice[i].dihedralAngle) == 0) &&
			(chiralChecking (lattice[i+1].dihedralAngle) == 1)
			)
		{
			possibleJumps.Gm8T++;
		}
		else if (
			(chiralChecking (lattice_tminus1[i-1].dihedralAngle) == 0) &&
			(chiralChecking (lattice_tminus1[i].dihedralAngle) == -1) &&
			(chiralChecking (lattice_tminus1[i+1].dihedralAngle) == 0) &&
			(chiralChecking (lattice[i-1].dihedralAngle) == 0) &&
			(chiralChecking (lattice[i].dihedralAngle) == 1) &&
			(chiralChecking (lattice[i+1].dihedralAngle) == 0)
			)
		{
			possibleJumps.Gm9Gp++;
		}
		else if (
			(chiralChecking (lattice_tminus1[i-1].dihedralAngle) == 0) &&
			(chiralChecking (lattice_tminus1[i].dihedralAngle) == -1) &&
			(chiralChecking (lattice_tminus1[i+1].dihedralAngle) == 0) &&
			(chiralChecking (lattice[i-1].dihedralAngle) == 0) &&
			(chiralChecking (lattice[i].dihedralAngle) == -1) &&
			(chiralChecking (lattice[i+1].dihedralAngle) == 0)
			)
		{
			possibleJumps.Gm9Gm++;
		}
		else if (
			(chiralChecking (lattice_tminus1[i-1].dihedralAngle) == 0) &&
			(chiralChecking (lattice_tminus1[i].dihedralAngle) == -1) &&
			(chiralChecking (lattice_tminus1[i+1].dihedralAngle) == 0) &&
			(chiralChecking (lattice[i-1].dihedralAngle) == 0) &&
			(chiralChecking (lattice[i].dihedralAngle) == 0) &&
			(chiralChecking (lattice[i+1].dihedralAngle) == 0)
			)
		{
			possibleJumps.Gm9T++;
		}
		else if (
			(chiralChecking (lattice_tminus1[i-1].dihedralAngle) == -1) &&
			(chiralChecking (lattice_tminus1[i].dihedralAngle) == 1) &&
			(chiralChecking (lattice_tminus1[i+1].dihedralAngle) == -1) &&
			(chiralChecking (lattice[i-1].dihedralAngle) == -1) &&
			(chiralChecking (lattice[i].dihedralAngle) == 1) &&
			(chiralChecking (lattice[i+1].dihedralAngle) == -1)
			)
		{
			possibleJumps.Gp1Gp++;
		}
		else if (
			(chiralChecking (lattice_tminus1[i-1].dihedralAngle) == -1) &&
			(chiralChecking (lattice_tminus1[i].dihedralAngle) == 1) &&
			(chiralChecking (lattice_tminus1[i+1].dihedralAngle) == -1) &&
			(chiralChecking (lattice[i-1].dihedralAngle) == -1) &&
			(chiralChecking (lattice[i].dihedralAngle) == -1) &&
			(chiralChecking (lattice[i+1].dihedralAngle) == -1)
			)
		{
			possibleJumps.Gp1Gm++;
		}
		else if (
			(chiralChecking (lattice_tminus1[i-1].dihedralAngle) == -1) &&
			(chiralChecking (lattice_tminus1[i].dihedralAngle) == 1) &&
			(chiralChecking (lattice_tminus1[i+1].dihedralAngle) == -1) &&
			(chiralChecking (lattice[i-1].dihedralAngle) == -1) &&
			(chiralChecking (lattice[i].dihedralAngle) == 0) &&
			(chiralChecking (lattice[i+1].dihedralAngle) == -1)
			)
		{
			possibleJumps.Gp1T++;
		}
		else if (
			(chiralChecking (lattice_tminus1[i-1].dihedralAngle) == -1) &&
			(chiralChecking (lattice_tminus1[i].dihedralAngle) == 1) &&
			(chiralChecking (lattice_tminus1[i+1].dihedralAngle) == 1) &&
			(chiralChecking (lattice[i-1].dihedralAngle) == -1) &&
			(chiralChecking (lattice[i].dihedralAngle) == 1) &&
			(chiralChecking (lattice[i+1].dihedralAngle) == 1)
			)
		{
			possibleJumps.Gp2Gp++;
		}
		else if (
			(chiralChecking (lattice_tminus1[i-1].dihedralAngle) == -1) &&
			(chiralChecking (lattice_tminus1[i].dihedralAngle) == 1) &&
			(chiralChecking (lattice_tminus1[i+1].dihedralAngle) == 1) &&
			(chiralChecking (lattice[i-1].dihedralAngle) == -1) &&
			(chiralChecking (lattice[i].dihedralAngle) == -1) &&
			(chiralChecking (lattice[i+1].dihedralAngle) == 1)
			)
		{
			possibleJumps.Gp2Gm++;
		}
		else if (
			(chiralChecking (lattice_tminus1[i-1].dihedralAngle) == -1) &&
			(chiralChecking (lattice_tminus1[i].dihedralAngle) == 1) &&
			(chiralChecking (lattice_tminus1[i+1].dihedralAngle) == 1) &&
			(chiralChecking (lattice[i-1].dihedralAngle) == -1) &&
			(chiralChecking (lattice[i].dihedralAngle) == 0) &&
			(chiralChecking (lattice[i+1].dihedralAngle) == 1)
			)
		{
			possibleJumps.Gp2T++;
		}
		else if (
			(chiralChecking (lattice_tminus1[i-1].dihedralAngle) == -1) &&
			(chiralChecking (lattice_tminus1[i].dihedralAngle) == 1) &&
			(chiralChecking (lattice_tminus1[i+1].dihedralAngle) == 0) &&
			(chiralChecking (lattice[i-1].dihedralAngle) == -1) &&
			(chiralChecking (lattice[i].dihedralAngle) == 1) &&
			(chiralChecking (lattice[i+1].dihedralAngle) == 0)
			)
		{
			possibleJumps.Gp3Gp++;
		}
		else if (
			(chiralChecking (lattice_tminus1[i-1].dihedralAngle) == -1) &&
			(chiralChecking (lattice_tminus1[i].dihedralAngle) == 1) &&
			(chiralChecking (lattice_tminus1[i+1].dihedralAngle) == 0) &&
			(chiralChecking (lattice[i-1].dihedralAngle) == -1) &&
			(chiralChecking (lattice[i].dihedralAngle) == -1) &&
			(chiralChecking (lattice[i+1].dihedralAngle) == 0)
			)
		{
			possibleJumps.Gp3Gm++;
		}
		else if (
			(chiralChecking (lattice_tminus1[i-1].dihedralAngle) == -1) &&
			(chiralChecking (lattice_tminus1[i].dihedralAngle) == 1) &&
			(chiralChecking (lattice_tminus1[i+1].dihedralAngle) == 0) &&
			(chiralChecking (lattice[i-1].dihedralAngle) == -1) &&
			(chiralChecking (lattice[i].dihedralAngle) == 0) &&
			(chiralChecking (lattice[i+1].dihedralAngle) == 0)
			)
		{
			possibleJumps.Gp3T++;
		}
		else if (
			(chiralChecking (lattice_tminus1[i-1].dihedralAngle) == 1) &&
			(chiralChecking (lattice_tminus1[i].dihedralAngle) == 1) &&
			(chiralChecking (lattice_tminus1[i+1].dihedralAngle) == -1) &&
			(chiralChecking (lattice[i-1].dihedralAngle) == 1) &&
			(chiralChecking (lattice[i].dihedralAngle) == 1) &&
			(chiralChecking (lattice[i+1].dihedralAngle) == -1)
			)
		{
			possibleJumps.Gp4Gp++;
		}
		else if (
			(chiralChecking (lattice_tminus1[i-1].dihedralAngle) == 1) &&
			(chiralChecking (lattice_tminus1[i].dihedralAngle) == 1) &&
			(chiralChecking (lattice_tminus1[i+1].dihedralAngle) == -1) &&
			(chiralChecking (lattice[i-1].dihedralAngle) == 1) &&
			(chiralChecking (lattice[i].dihedralAngle) == -1) &&
			(chiralChecking (lattice[i+1].dihedralAngle) == -1)
			)
		{
			possibleJumps.Gp4Gm++;
		}
		else if (
			(chiralChecking (lattice_tminus1[i-1].dihedralAngle) == 1) &&
			(chiralChecking (lattice_tminus1[i].dihedralAngle) == 1) &&
			(chiralChecking (lattice_tminus1[i+1].dihedralAngle) == -1) &&
			(chiralChecking (lattice[i-1].dihedralAngle) == 1) &&
			(chiralChecking (lattice[i].dihedralAngle) == 0) &&
			(chiralChecking (lattice[i+1].dihedralAngle) == -1)
			)
		{
			possibleJumps.Gp4T++;
		}
		else if (
			(chiralChecking (lattice_tminus1[i-1].dihedralAngle) == 1) &&
			(chiralChecking (lattice_tminus1[i].dihedralAngle) == 1) &&
			(chiralChecking (lattice_tminus1[i+1].dihedralAngle) == 1) &&
			(chiralChecking (lattice[i-1].dihedralAngle) == 1) &&
			(chiralChecking (lattice[i].dihedralAngle) == 1) &&
			(chiralChecking (lattice[i+1].dihedralAngle) == 1)
			)
		{
			possibleJumps.Gp5Gp++;
		}
		else if (
			(chiralChecking (lattice_tminus1[i-1].dihedralAngle) == 1) &&
			(chiralChecking (lattice_tminus1[i].dihedralAngle) == 1) &&
			(chiralChecking (lattice_tminus1[i+1].dihedralAngle) == 1) &&
			(chiralChecking (lattice[i-1].dihedralAngle) == 1) &&
			(chiralChecking (lattice[i].dihedralAngle) == -1) &&
			(chiralChecking (lattice[i+1].dihedralAngle) == 1)
			)
		{
			possibleJumps.Gp5Gm++;
		}
		else if (
			(chiralChecking (lattice_tminus1[i-1].dihedralAngle) == 1) &&
			(chiralChecking (lattice_tminus1[i].dihedralAngle) == 1) &&
			(chiralChecking (lattice_tminus1[i+1].dihedralAngle) == 1) &&
			(chiralChecking (lattice[i-1].dihedralAngle) == 1) &&
			(chiralChecking (lattice[i].dihedralAngle) == 0) &&
			(chiralChecking (lattice[i+1].dihedralAngle) == 1)
			)
		{
			possibleJumps.Gp5T++;
		}
		else if (
			(chiralChecking (lattice_tminus1[i-1].dihedralAngle) == 1) &&
			(chiralChecking (lattice_tminus1[i].dihedralAngle) == 1) &&
			(chiralChecking (lattice_tminus1[i+1].dihedralAngle) == 0) &&
			(chiralChecking (lattice[i-1].dihedralAngle) == 1) &&
			(chiralChecking (lattice[i].dihedralAngle) == 1) &&
			(chiralChecking (lattice[i+1].dihedralAngle) == 0)
			)
		{
			possibleJumps.Gp6Gp++;
		}
		else if (
			(chiralChecking (lattice_tminus1[i-1].dihedralAngle) == 1) &&
			(chiralChecking (lattice_tminus1[i].dihedralAngle) == 1) &&
			(chiralChecking (lattice_tminus1[i+1].dihedralAngle) == 0) &&
			(chiralChecking (lattice[i-1].dihedralAngle) == 1) &&
			(chiralChecking (lattice[i].dihedralAngle) == -1) &&
			(chiralChecking (lattice[i+1].dihedralAngle) == 0)
			)
		{
			possibleJumps.Gp6Gm++;
		}
		else if (
			(chiralChecking (lattice_tminus1[i-1].dihedralAngle) == 1) &&
			(chiralChecking (lattice_tminus1[i].dihedralAngle) == 1) &&
			(chiralChecking (lattice_tminus1[i+1].dihedralAngle) == 0) &&
			(chiralChecking (lattice[i-1].dihedralAngle) == 1) &&
			(chiralChecking (lattice[i].dihedralAngle) == 0) &&
			(chiralChecking (lattice[i+1].dihedralAngle) == 0)
			)
		{
			possibleJumps.Gp6T++;
		}
		else if (
			(chiralChecking (lattice_tminus1[i-1].dihedralAngle) == 0) &&
			(chiralChecking (lattice_tminus1[i].dihedralAngle) == 1) &&
			(chiralChecking (lattice_tminus1[i+1].dihedralAngle) == -1) &&
			(chiralChecking (lattice[i-1].dihedralAngle) == 0) &&
			(chiralChecking (lattice[i].dihedralAngle) == 1) &&
			(chiralChecking (lattice[i+1].dihedralAngle) == -1)
			)
		{
			possibleJumps.Gp7Gp++;
		}
		else if (
			(chiralChecking (lattice_tminus1[i-1].dihedralAngle) == 0) &&
			(chiralChecking (lattice_tminus1[i].dihedralAngle) == 1) &&
			(chiralChecking (lattice_tminus1[i+1].dihedralAngle) == -1) &&
			(chiralChecking (lattice[i-1].dihedralAngle) == 0) &&
			(chiralChecking (lattice[i].dihedralAngle) == -1) &&
			(chiralChecking (lattice[i+1].dihedralAngle) == -1)
			)
		{
			possibleJumps.Gp7Gm++;
		}
		else if (
			(chiralChecking (lattice_tminus1[i-1].dihedralAngle) == 0) &&
			(chiralChecking (lattice_tminus1[i].dihedralAngle) == 1) &&
			(chiralChecking (lattice_tminus1[i+1].dihedralAngle) == -1) &&
			(chiralChecking (lattice[i-1].dihedralAngle) == 0) &&
			(chiralChecking (lattice[i].dihedralAngle) == 0) &&
			(chiralChecking (lattice[i+1].dihedralAngle) == -1)
			)
		{
			possibleJumps.Gp7T++;
		}
		else if (
			(chiralChecking (lattice_tminus1[i-1].dihedralAngle) == 0) &&
			(chiralChecking (lattice_tminus1[i].dihedralAngle) == 1) &&
			(chiralChecking (lattice_tminus1[i+1].dihedralAngle) == 1) &&
			(chiralChecking (lattice[i-1].dihedralAngle) == 0) &&
			(chiralChecking (lattice[i].dihedralAngle) == 1) &&
			(chiralChecking (lattice[i+1].dihedralAngle) == 1)
			)
		{
			possibleJumps.Gp8Gp++;
		}
		else if (
			(chiralChecking (lattice_tminus1[i-1].dihedralAngle) == 0) &&
			(chiralChecking (lattice_tminus1[i].dihedralAngle) == 1) &&
			(chiralChecking (lattice_tminus1[i+1].dihedralAngle) == 1) &&
			(chiralChecking (lattice[i-1].dihedralAngle) == 0) &&
			(chiralChecking (lattice[i].dihedralAngle) == -1) &&
			(chiralChecking (lattice[i+1].dihedralAngle) == 1)
			)
		{
			possibleJumps.Gp8Gm++;
		}
		else if (
			(chiralChecking (lattice_tminus1[i-1].dihedralAngle) == 0) &&
			(chiralChecking (lattice_tminus1[i].dihedralAngle) == 1) &&
			(chiralChecking (lattice_tminus1[i+1].dihedralAngle) == 1) &&
			(chiralChecking (lattice[i-1].dihedralAngle) == 0) &&
			(chiralChecking (lattice[i].dihedralAngle) == 0) &&
			(chiralChecking (lattice[i+1].dihedralAngle) == 1)
			)
		{
			possibleJumps.Gp8T++;
		}
		else if (
			(chiralChecking (lattice_tminus1[i-1].dihedralAngle) == 0) &&
			(chiralChecking (lattice_tminus1[i].dihedralAngle) == 1) &&
			(chiralChecking (lattice_tminus1[i+1].dihedralAngle) == 0) &&
			(chiralChecking (lattice[i-1].dihedralAngle) == 0) &&
			(chiralChecking (lattice[i].dihedralAngle) == 1) &&
			(chiralChecking (lattice[i+1].dihedralAngle) == 0)
			)
		{
			possibleJumps.Gp9Gp++;
		}
		else if (
			(chiralChecking (lattice_tminus1[i-1].dihedralAngle) == 0) &&
			(chiralChecking (lattice_tminus1[i].dihedralAngle) == 1) &&
			(chiralChecking (lattice_tminus1[i+1].dihedralAngle) == 0) &&
			(chiralChecking (lattice[i-1].dihedralAngle) == 0) &&
			(chiralChecking (lattice[i].dihedralAngle) == -1) &&
			(chiralChecking (lattice[i+1].dihedralAngle) == 0)
			)
		{
			possibleJumps.Gp9Gm++;
		}
		else if (
			(chiralChecking (lattice_tminus1[i-1].dihedralAngle) == 0) &&
			(chiralChecking (lattice_tminus1[i].dihedralAngle) == 1) &&
			(chiralChecking (lattice_tminus1[i+1].dihedralAngle) == 0) &&
			(chiralChecking (lattice[i-1].dihedralAngle) == 0) &&
			(chiralChecking (lattice[i].dihedralAngle) == 0) &&
			(chiralChecking (lattice[i+1].dihedralAngle) == 0)
			)
		{
			possibleJumps.Gp9T++;
		}
		// if ((lattice[i-1].dihedralAngle < -38) && (lattice[i-1].dihedralAngle > -116) && (((lattice[i].dihedralAngle < -117) && (lattice[i].dihedralAngle > -180)) || ((lattice[i].dihedralAngle < 180) && (lattice[i].dihedralAngle > 116))) && (lattice[i+1].dihedralAngle < 115) && (lattice[i+1].dihedralAngle > 34) && (lattice_tminus1[i-1].dihedralAngle < -38) && (lattice_tminus1[i-1].dihedralAngle > -116) && (lattice_tminus1[i].dihedralAngle < 115) && (lattice_tminus1[i].dihedralAngle > 34) &&	(lattice_tminus1[i+1].dihedralAngle < 115) && (lattice_tminus1[i+1].dihedralAngle > 34))
		// {
		// 	possibleJumps.t1Gp++;
		// }
		// else if ((lattice[i-1].dihedralAngle < -38) && (lattice[i-1].dihedralAngle > -116) && (((lattice[i].dihedralAngle < -117) && (lattice[i].dihedralAngle > -180)) || ((lattice[i].dihedralAngle < 180) && (lattice[i].dihedralAngle > 116))) && (lattice[i+1].dihedralAngle < 115) && (lattice[i+1].dihedralAngle > 34) && (lattice_tminus1[i-1].dihedralAngle < -38) && (lattice_tminus1[i-1].dihedralAngle > -116) && (lattice_tminus1[i].dihedralAngle < -38) && (lattice_tminus1[i].dihedralAngle > -116) && (lattice_tminus1[i+1].dihedralAngle < 115) && (lattice_tminus1[i+1].dihedralAngle > 34))
		// {
		// 	possibleJumps.t1Gm++;
		// }
		// else if ((lattice[i-1].dihedralAngle < -38) && (lattice[i-1].dihedralAngle > -116) && (((lattice[i].dihedralAngle < -117) && (lattice[i].dihedralAngle > -180)) || ((lattice[i].dihedralAngle < 180) && (lattice[i].dihedralAngle > 116))) && (lattice[i+1].dihedralAngle < 115) && (lattice[i+1].dihedralAngle > 34) && (lattice_tminus1[i-1].dihedralAngle < -38) && (lattice_tminus1[i-1].dihedralAngle > -116) && (((lattice_tminus1[i].dihedralAngle < -117) && (lattice_tminus1[i].dihedralAngle > -180)) || ((lattice_tminus1[i].dihedralAngle < 180) && (lattice_tminus1[i].dihedralAngle > 116))) && (lattice_tminus1[i+1].dihedralAngle < 115) && (lattice_tminus1[i+1].dihedralAngle > 34))
		// {
		// 	possibleJumps.t1T++;
		// }
		// else if ((lattice[i-1].dihedralAngle < -38) && (lattice[i-1].dihedralAngle > -116) && (((lattice[i].dihedralAngle < -117) && (lattice[i].dihedralAngle > -180)) || ((lattice[i].dihedralAngle < 180) && (lattice[i].dihedralAngle > 116))) && (((lattice[i+1].dihedralAngle < -117) && (lattice[i+1].dihedralAngle > -180)) || ((lattice[i+1].dihedralAngle < 180) && (lattice[i+1].dihedralAngle > 116))) && (lattice_tminus1[i-1].dihedralAngle < -38) && (lattice_tminus1[i-1].dihedralAngle > -116) && (lattice_tminus1[i].dihedralAngle < 115) && (lattice_tminus1[i].dihedralAngle > 34) && (((lattice_tminus1[i+1].dihedralAngle < -117) && (lattice_tminus1[i+1].dihedralAngle > -180)) || ((lattice_tminus1[i+1].dihedralAngle < 180) && (lattice_tminus1[i+1].dihedralAngle > 116))))
		// {
		// 	possibleJumps.t2Gp++;
		// }
		// else if ((lattice[i-1].dihedralAngle < -38) && (lattice[i-1].dihedralAngle > -116) && (((lattice[i].dihedralAngle < -117) && (lattice[i].dihedralAngle > -180)) || ((lattice[i].dihedralAngle < 180) && (lattice[i].dihedralAngle > 116))) && (((lattice[i+1].dihedralAngle < -117) && (lattice[i+1].dihedralAngle > -180)) || ((lattice[i+1].dihedralAngle < 180) && (lattice[i+1].dihedralAngle > 116))) && (lattice_tminus1[i-1].dihedralAngle < -38) && (lattice_tminus1[i-1].dihedralAngle > -116) && (lattice_tminus1[i].dihedralAngle < -38) && (lattice_tminus1[i].dihedralAngle > -116) &&	(((lattice_tminus1[i+1].dihedralAngle < -117) && (lattice_tminus1[i+1].dihedralAngle > -180)) || ((lattice_tminus1[i+1].dihedralAngle < 180) && (lattice_tminus1[i+1].dihedralAngle > 116))))
		// {
		// 	possibleJumps.t2Gm++;
		// }
		// else if ((lattice[i-1].dihedralAngle < -38) && (lattice[i-1].dihedralAngle > -116) && (((lattice[i].dihedralAngle < -117) && (lattice[i].dihedralAngle > -180)) || ((lattice[i].dihedralAngle < 180) && (lattice[i].dihedralAngle > 116))) && (((lattice[i+1].dihedralAngle < -117) && (lattice[i+1].dihedralAngle > -180)) || ((lattice[i+1].dihedralAngle < 180) && (lattice[i+1].dihedralAngle > 116))) && (lattice_tminus1[i-1].dihedralAngle < -38) && (lattice_tminus1[i-1].dihedralAngle > -116) && (((lattice_tminus1[i].dihedralAngle < -117) && (lattice_tminus1[i].dihedralAngle > -180)) || ((lattice_tminus1[i].dihedralAngle < 180) && (lattice_tminus1[i].dihedralAngle > 116))) && (((lattice_tminus1[i+1].dihedralAngle < -117) && (lattice_tminus1[i+1].dihedralAngle > -180)) || ((lattice_tminus1[i+1].dihedralAngle < 180) && (lattice_tminus1[i+1].dihedralAngle > 116))))
		// {
		// 	possibleJumps.t2T++;
		// }
		// else if ((lattice[i-1].dihedralAngle < -38) && (lattice[i-1].dihedralAngle > -116) && (((lattice[i].dihedralAngle < -117) && (lattice[i].dihedralAngle > -180)) || ((lattice[i].dihedralAngle < 180) && (lattice[i].dihedralAngle > 116))) && (lattice[i+1].dihedralAngle < -38) && (lattice[i+1].dihedralAngle > -116) && (lattice_tminus1[i-1].dihedralAngle < -38) && (lattice_tminus1[i-1].dihedralAngle > -116) && (lattice_tminus1[i].dihedralAngle < 115) && (lattice_tminus1[i].dihedralAngle > 34) && (lattice_tminus1[i+1].dihedralAngle < -38) && (lattice_tminus1[i+1].dihedralAngle > -116))
		// {
		// 	possibleJumps.t3Gp++;
		// }
		// else if ((lattice[i-1].dihedralAngle < -38) && (lattice[i-1].dihedralAngle > -116) && (((lattice[i].dihedralAngle < -117) && (lattice[i].dihedralAngle > -180)) || ((lattice[i].dihedralAngle < 180) && (lattice[i].dihedralAngle > 116))) && (lattice[i+1].dihedralAngle < -38) && (lattice[i+1].dihedralAngle > -116) && (lattice_tminus1[i-1].dihedralAngle < -38) && (lattice_tminus1[i-1].dihedralAngle > -116) && (lattice_tminus1[i].dihedralAngle < -38) && (lattice_tminus1[i].dihedralAngle > -116) && (lattice_tminus1[i+1].dihedralAngle < -38) && (lattice_tminus1[i+1].dihedralAngle > -116))
		// {
		// 	possibleJumps.t3Gm++;
		// }
		// else if ((lattice[i-1].dihedralAngle < -38) && (lattice[i-1].dihedralAngle > -116) && (((lattice[i].dihedralAngle < -117) && (lattice[i].dihedralAngle > -180)) || ((lattice[i].dihedralAngle < 180) && (lattice[i].dihedralAngle > 116))) && (lattice[i+1].dihedralAngle < -38) && (lattice[i+1].dihedralAngle > -116) && (lattice_tminus1[i-1].dihedralAngle < -38) && (lattice_tminus1[i-1].dihedralAngle > -116) &&	(((lattice_tminus1[i].dihedralAngle < -117) && (lattice_tminus1[i].dihedralAngle > -180)) || ((lattice_tminus1[i].dihedralAngle < 180) && (lattice_tminus1[i].dihedralAngle > 116))) &&	(lattice_tminus1[i+1].dihedralAngle < -38) && (lattice_tminus1[i+1].dihedralAngle > -116))
		// {
		// 	possibleJumps.t3T++;
		// }
		// else if ((lattice[i-1].dihedralAngle < 115) && (lattice[i-1].dihedralAngle > 34) &&	(((lattice[i].dihedralAngle < -117) && (lattice[i].dihedralAngle > -180)) || ((lattice[i].dihedralAngle < 180) && (lattice[i].dihedralAngle > 116))) && (lattice[i+1].dihedralAngle < -38) && (lattice[i+1].dihedralAngle > -116) && (lattice_tminus1[i-1].dihedralAngle < 115) && (lattice_tminus1[i-1].dihedralAngle > 34) && (lattice_tminus1[i].dihedralAngle < 115) && (lattice_tminus1[i].dihedralAngle > 34) && (lattice_tminus1[i+1].dihedralAngle < -38) && (lattice_tminus1[i+1].dihedralAngle > -116))
		// {
		// 	possibleJumps.t4Gp++;
		// }
		// else if ((lattice[i-1].dihedralAngle < 115) && (lattice[i-1].dihedralAngle > 34) && (((lattice[i].dihedralAngle < -117) && (lattice[i].dihedralAngle > -180)) || ((lattice[i].dihedralAngle < 180) && (lattice[i].dihedralAngle > 116))) && (lattice[i+1].dihedralAngle < -38) && (lattice[i+1].dihedralAngle > -116) && (lattice_tminus1[i-1].dihedralAngle < 115) && (lattice_tminus1[i-1].dihedralAngle > 34) && (lattice_tminus1[i].dihedralAngle < -38) && (lattice_tminus1[i].dihedralAngle > -116) && (lattice_tminus1[i+1].dihedralAngle < -38) && (lattice_tminus1[i+1].dihedralAngle > -116))
		// {
		// 	possibleJumps.t4Gm++;
		// }
		// else if ((lattice[i-1].dihedralAngle < 115) && (lattice[i-1].dihedralAngle > 34) &&	(((lattice[i].dihedralAngle < -117) && (lattice[i].dihedralAngle > -180)) || ((lattice[i].dihedralAngle < 180) && (lattice[i].dihedralAngle > 116))) && (lattice[i+1].dihedralAngle < -38) && (lattice[i+1].dihedralAngle > -116) && (lattice_tminus1[i-1].dihedralAngle < 115) && (lattice_tminus1[i-1].dihedralAngle > 34) && (((lattice_tminus1[i].dihedralAngle < -117) && (lattice_tminus1[i].dihedralAngle > -180)) || ((lattice_tminus1[i].dihedralAngle < 180) && (lattice_tminus1[i].dihedralAngle > 116))) && (lattice_tminus1[i+1].dihedralAngle < -38) && (lattice_tminus1[i+1].dihedralAngle > -116))
		// {
		// 	possibleJumps.t4T++;
		// }
		// else if ((lattice[i-1].dihedralAngle < 115) && (lattice[i-1].dihedralAngle > 34) &&	(((lattice[i].dihedralAngle < -117) && (lattice[i].dihedralAngle > -180)) || ((lattice[i].dihedralAngle < 180) && (lattice[i].dihedralAngle > 116))) &&	(((lattice[i+1].dihedralAngle < -117) && (lattice[i+1].dihedralAngle > -180)) || ((lattice[i+1].dihedralAngle < 180) && (lattice[i+1].dihedralAngle > 116))) && (lattice_tminus1[i-1].dihedralAngle < 115) && (lattice_tminus1[i-1].dihedralAngle > 34) && (lattice_tminus1[i].dihedralAngle < 115) && (lattice_tminus1[i].dihedralAngle > 34) && (((lattice_tminus1[i+1].dihedralAngle < -117) && (lattice_tminus1[i+1].dihedralAngle > -180)) || ((lattice_tminus1[i+1].dihedralAngle < 180) && (lattice_tminus1[i+1].dihedralAngle > 116))))
		// {
		// 	possibleJumps.t5Gp++;
		// }
		// else if ((lattice[i-1].dihedralAngle < 115) && (lattice[i-1].dihedralAngle > 34) &&	(((lattice[i].dihedralAngle < -117) && (lattice[i].dihedralAngle > -180)) || ((lattice[i].dihedralAngle < 180) && (lattice[i].dihedralAngle > 116))) && (((lattice[i+1].dihedralAngle < -117) && (lattice[i+1].dihedralAngle > -180)) || ((lattice[i+1].dihedralAngle < 180) && (lattice[i+1].dihedralAngle > 116))) &&	(lattice_tminus1[i-1].dihedralAngle < 115) && (lattice_tminus1[i-1].dihedralAngle > 34) && (lattice_tminus1[i].dihedralAngle < -38) && (lattice_tminus1[i].dihedralAngle > -116) &&	(((lattice_tminus1[i+1].dihedralAngle < -117) && (lattice_tminus1[i+1].dihedralAngle > -180)) || ((lattice_tminus1[i+1].dihedralAngle < 180) && (lattice_tminus1[i+1].dihedralAngle > 116))))
		// {
		// 	possibleJumps.t5Gm++;
		// }
		// else if ((lattice[i-1].dihedralAngle < 115) && (lattice[i-1].dihedralAngle > 34) &&	(((lattice[i].dihedralAngle < -117) && (lattice[i].dihedralAngle > -180)) || ((lattice[i].dihedralAngle < 180) && (lattice[i].dihedralAngle > 116))) &&	(((lattice[i+1].dihedralAngle < -117) && (lattice[i+1].dihedralAngle > -180)) || ((lattice[i+1].dihedralAngle < 180) && (lattice[i+1].dihedralAngle > 116))) &&	(lattice_tminus1[i-1].dihedralAngle < 115) && (lattice_tminus1[i-1].dihedralAngle > 34) && (((lattice_tminus1[i].dihedralAngle < -117) && (lattice_tminus1[i].dihedralAngle > -180)) || ((lattice_tminus1[i].dihedralAngle < 180) && (lattice_tminus1[i].dihedralAngle > 116))) && (((lattice_tminus1[i+1].dihedralAngle < -117) && (lattice_tminus1[i+1].dihedralAngle > -180)) || ((lattice_tminus1[i+1].dihedralAngle < 180) && (lattice_tminus1[i+1].dihedralAngle > 116))))
		// {
		// 	possibleJumps.t5T++;
		// }
		// else if ((lattice[i-1].dihedralAngle < 115) && (lattice[i-1].dihedralAngle > 34) &&	(((lattice[i].dihedralAngle < -117) && (lattice[i].dihedralAngle > -180)) || ((lattice[i].dihedralAngle < 180) && (lattice[i].dihedralAngle > 116))) &&	(lattice[i+1].dihedralAngle < 115) && (lattice[i+1].dihedralAngle > 34) && (lattice_tminus1[i-1].dihedralAngle < 115) && (lattice_tminus1[i-1].dihedralAngle > 34) && (lattice_tminus1[i].dihedralAngle < 115) && (lattice_tminus1[i].dihedralAngle > 34) && (lattice_tminus1[i+1].dihedralAngle < 115) && (lattice_tminus1[i+1].dihedralAngle > 34))
		// {
		// 	possibleJumps.t6Gp++;
		// }
		// else if ((lattice[i-1].dihedralAngle < 115) && (lattice[i-1].dihedralAngle > 34) &&	(((lattice[i].dihedralAngle < -117) && (lattice[i].dihedralAngle > -180)) || ((lattice[i].dihedralAngle < 180) && (lattice[i].dihedralAngle > 116))) && (lattice[i+1].dihedralAngle < 115) && (lattice[i+1].dihedralAngle > 34) && (lattice_tminus1[i-1].dihedralAngle < 115) && (lattice_tminus1[i-1].dihedralAngle > 34) && (((lattice_tminus1[i].dihedralAngle < -117) && (lattice_tminus1[i].dihedralAngle > -180)) || ((lattice_tminus1[i].dihedralAngle < 180) && (lattice_tminus1[i].dihedralAngle > 116))) && (lattice_tminus1[i+1].dihedralAngle < 115) && (lattice_tminus1[i+1].dihedralAngle > 34))
		// {
		// 	possibleJumps.t6T++;
		// }
		// else if ((lattice[i-1].dihedralAngle < 115) && (lattice[i-1].dihedralAngle > 34) &&	(((lattice[i].dihedralAngle < -117) && (lattice[i].dihedralAngle > -180)) || ((lattice[i].dihedralAngle < 180) && (lattice[i].dihedralAngle > 116))) &&	(lattice[i+1].dihedralAngle < 115) && (lattice[i+1].dihedralAngle > 34) && (lattice_tminus1[i-1].dihedralAngle < 115) && (lattice_tminus1[i-1].dihedralAngle > 34) && (lattice_tminus1[i].dihedralAngle < -38) && (lattice_tminus1[i].dihedralAngle > -116) && (lattice_tminus1[i+1].dihedralAngle < 115) && (lattice_tminus1[i+1].dihedralAngle > 34))
		// {
		// 	possibleJumps.t6Gm++;
		// }
		// else if ((((lattice[i-1].dihedralAngle < -117) && (lattice[i-1].dihedralAngle > -180)) || ((lattice[i-1].dihedralAngle < 180) && (lattice[i-1].dihedralAngle > 116))) && (((lattice[i].dihedralAngle < -117) && (lattice[i].dihedralAngle > -180)) || ((lattice[i].dihedralAngle < 180) && (lattice[i].dihedralAngle > 116))) && (lattice[i+1].dihedralAngle < 115) && (lattice[i+1].dihedralAngle > 34) && (((lattice_tminus1[i-1].dihedralAngle < -117) && (lattice_tminus1[i-1].dihedralAngle > -180)) || ((lattice_tminus1[i-1].dihedralAngle < 180) && (lattice_tminus1[i-1].dihedralAngle > 116))) &&	(lattice_tminus1[i+1].dihedralAngle < 115) && (lattice_tminus1[i+1].dihedralAngle > 34) && (lattice_tminus1[i+1].dihedralAngle < 115) && (lattice_tminus1[i+1].dihedralAngle > 34))
		// {
		// 	possibleJumps.t7Gp++;
		// }
		// else if ((((lattice[i-1].dihedralAngle < -117) && (lattice[i-1].dihedralAngle > -180)) || ((lattice[i-1].dihedralAngle < 180) && (lattice[i-1].dihedralAngle > 116))) && (((lattice[i].dihedralAngle < -117) && (lattice[i].dihedralAngle > -180)) || ((lattice[i].dihedralAngle < 180) && (lattice[i].dihedralAngle > 116))) && (lattice[i+1].dihedralAngle < 115) && (lattice[i+1].dihedralAngle > 34) &&	(((lattice_tminus1[i-1].dihedralAngle < -117) && (lattice_tminus1[i-1].dihedralAngle > -180)) || ((lattice_tminus1[i-1].dihedralAngle < 180) && (lattice_tminus1[i-1].dihedralAngle > 116))) && (((lattice_tminus1[i].dihedralAngle < -117) && (lattice_tminus1[i].dihedralAngle > -180)) || ((lattice_tminus1[i].dihedralAngle < 180) && (lattice_tminus1[i].dihedralAngle > 116))) && (lattice_tminus1[i+1].dihedralAngle < 115) && (lattice_tminus1[i+1].dihedralAngle > 34))
		// {
		// 	possibleJumps.t7T++;
		// }
		// else if ((((lattice[i-1].dihedralAngle < -117) && (lattice[i-1].dihedralAngle > -180)) || ((lattice[i-1].dihedralAngle < 180) && (lattice[i-1].dihedralAngle > 116))) && (((lattice[i].dihedralAngle < -117) && (lattice[i].dihedralAngle > -180)) || ((lattice[i].dihedralAngle < 180) && (lattice[i].dihedralAngle > 116))) && (lattice[i+1].dihedralAngle < 115) && (lattice[i+1].dihedralAngle > 34) &&	(((lattice_tminus1[i-1].dihedralAngle < -117) && (lattice_tminus1[i-1].dihedralAngle > -180)) || ((lattice_tminus1[i-1].dihedralAngle < 180) && (lattice_tminus1[i-1].dihedralAngle > 116))) &&	(lattice_tminus1[i].dihedralAngle < -38) && (lattice_tminus1[i].dihedralAngle > -116) && (lattice_tminus1[i+1].dihedralAngle < 115) && (lattice_tminus1[i+1].dihedralAngle > 34))
		// {
		// 	possibleJumps.t7Gm++;
		// }
		// else if ((((lattice[i-1].dihedralAngle < -117) && (lattice[i-1].dihedralAngle > -180)) || ((lattice[i-1].dihedralAngle < 180) && (lattice[i-1].dihedralAngle > 116))) && (((lattice[i].dihedralAngle < -117) && (lattice[i].dihedralAngle > -180)) || ((lattice[i].dihedralAngle < 180) && (lattice[i].dihedralAngle > 116))) && (lattice[i+1].dihedralAngle < -38) && (lattice[i+1].dihedralAngle > -116) && (((lattice_tminus1[i-1].dihedralAngle < -117) && (lattice_tminus1[i-1].dihedralAngle > -180)) || ((lattice_tminus1[i-1].dihedralAngle < 180) && (lattice_tminus1[i-1].dihedralAngle > 116))) && (lattice_tminus1[i].dihedralAngle < 115) && (lattice_tminus1[i].dihedralAngle > 34) && (lattice_tminus1[i+1].dihedralAngle < -38) && (lattice_tminus1[i+1].dihedralAngle > -116))
		// {
		// 	possibleJumps.t8Gp++;
		// }
		// else if ((((lattice[i-1].dihedralAngle < -117) && (lattice[i-1].dihedralAngle > -180)) || ((lattice[i-1].dihedralAngle < 180) && (lattice[i-1].dihedralAngle > 116))) && (((lattice[i].dihedralAngle < -117) && (lattice[i].dihedralAngle > -180)) || ((lattice[i].dihedralAngle < 180) && (lattice[i].dihedralAngle > 116))) && (lattice[i+1].dihedralAngle < -38) && (lattice[i+1].dihedralAngle > -116) && (((lattice_tminus1[i-1].dihedralAngle < -117) && (lattice_tminus1[i-1].dihedralAngle > -180)) || ((lattice_tminus1[i-1].dihedralAngle < 180) && (lattice_tminus1[i-1].dihedralAngle > 116))) && (((lattice_tminus1[i].dihedralAngle < -117) && (lattice_tminus1[i].dihedralAngle > -180)) || ((lattice_tminus1[i].dihedralAngle < 180) && (lattice_tminus1[i].dihedralAngle > 116))) && (lattice_tminus1[i+1].dihedralAngle < -38) && (lattice_tminus1[i+1].dihedralAngle > -116))
		// {
		// 	possibleJumps.t8T++;
		// }
		// else if ((((lattice[i-1].dihedralAngle < -117) && (lattice[i-1].dihedralAngle > -180)) || ((lattice[i-1].dihedralAngle < 180) && (lattice[i-1].dihedralAngle > 116))) && (((lattice[i].dihedralAngle < -117) && (lattice[i].dihedralAngle > -180)) || ((lattice[i].dihedralAngle < 180) && (lattice[i].dihedralAngle > 116))) && (lattice[i+1].dihedralAngle < -38) && (lattice[i+1].dihedralAngle > -116) && (((lattice_tminus1[i-1].dihedralAngle < -117) && (lattice_tminus1[i-1].dihedralAngle > -180)) || ((lattice_tminus1[i-1].dihedralAngle < 180) && (lattice_tminus1[i-1].dihedralAngle > 116))) && (lattice_tminus1[i].dihedralAngle < -38) && (lattice_tminus1[i].dihedralAngle > -116) && (lattice_tminus1[i+1].dihedralAngle < -38) && (lattice_tminus1[i+1].dihedralAngle > -116))
		// {
		// 	possibleJumps.t8Gm;
		// }
		// else if ((((lattice[i-1].dihedralAngle < -117) && (lattice[i-1].dihedralAngle > -180)) || ((lattice[i-1].dihedralAngle < 180) && (lattice[i-1].dihedralAngle > 116))) && (((lattice[i].dihedralAngle < -117) && (lattice[i].dihedralAngle > -180)) || ((lattice[i].dihedralAngle < 180) && (lattice[i].dihedralAngle > 116))) && (((lattice[i+1].dihedralAngle < -117) && (lattice[i+1].dihedralAngle > -180)) || ((lattice[i+1].dihedralAngle < 180) && (lattice[i+1].dihedralAngle > 116))) && (((lattice_tminus1[i-1].dihedralAngle < -117) && (lattice_tminus1[i-1].dihedralAngle > -180)) || ((lattice_tminus1[i-1].dihedralAngle < 180) && (lattice_tminus1[i-1].dihedralAngle > 116))) && (lattice_tminus1[i+1].dihedralAngle < 115) && (lattice_tminus1[i+1].dihedralAngle > 34) &&	(((lattice_tminus1[i+1].dihedralAngle < -117) && (lattice_tminus1[i+1].dihedralAngle > -180)) || ((lattice_tminus1[i+1].dihedralAngle < 180) && (lattice_tminus1[i+1].dihedralAngle > 116))))
		// {
		// 	possibleJumps.t9Gp++;
		// }
		// else if ((((lattice[i-1].dihedralAngle < -117) && (lattice[i-1].dihedralAngle > -180)) || ((lattice[i-1].dihedralAngle < 180) && (lattice[i-1].dihedralAngle > 116))) && (((lattice[i].dihedralAngle < -117) && (lattice[i].dihedralAngle > -180)) || ((lattice[i].dihedralAngle < 180) && (lattice[i].dihedralAngle > 116))) && (((lattice[i+1].dihedralAngle < -117) && (lattice[i+1].dihedralAngle > -180)) || ((lattice[i+1].dihedralAngle < 180) && (lattice[i+1].dihedralAngle > 116))) && (((lattice_tminus1[i-1].dihedralAngle < -117) && (lattice_tminus1[i-1].dihedralAngle > -180)) || ((lattice_tminus1[i-1].dihedralAngle < 180) && (lattice_tminus1[i-1].dihedralAngle > 116))) && (((lattice_tminus1[i].dihedralAngle < -117) && (lattice_tminus1[i].dihedralAngle > -180)) || ((lattice_tminus1[i].dihedralAngle < 180) && (lattice_tminus1[i].dihedralAngle > 116))) && (((lattice_tminus1[i+1].dihedralAngle < -117) && (lattice_tminus1[i+1].dihedralAngle > -180)) || ((lattice_tminus1[i+1].dihedralAngle < 180) && (lattice_tminus1[i+1].dihedralAngle > 116))))
		// {
		// 	possibleJumps.t9T++;
		// }
		// else if ((((lattice[i-1].dihedralAngle < -117) && (lattice[i-1].dihedralAngle > -180)) || ((lattice[i-1].dihedralAngle < 180) && (lattice[i-1].dihedralAngle > 116))) && (((lattice[i].dihedralAngle < -117) && (lattice[i].dihedralAngle > -180)) || ((lattice[i].dihedralAngle < 180) && (lattice[i].dihedralAngle > 116))) && (((lattice[i+1].dihedralAngle < -117) && (lattice[i+1].dihedralAngle > -180)) || ((lattice[i+1].dihedralAngle < 180) && (lattice[i+1].dihedralAngle > 116))) && (((lattice_tminus1[i-1].dihedralAngle < -117) && (lattice_tminus1[i-1].dihedralAngle > -180)) || ((lattice_tminus1[i-1].dihedralAngle < 180) && (lattice_tminus1[i-1].dihedralAngle > 116))) && (lattice_tminus1[i-1].dihedralAngle < -38) && (lattice_tminus1[i-1].dihedralAngle > -116) && (((lattice_tminus1[i+1].dihedralAngle < -117) && (lattice_tminus1[i+1].dihedralAngle > -180)) || ((lattice_tminus1[i+1].dihedralAngle < 180) && (lattice_tminus1[i+1].dihedralAngle > 116))))
		// {
		// 	possibleJumps.t9Gm++;
		// }
		// else if ((((lattice[i-1].dihedralAngle < -117) && (lattice[i-1].dihedralAngle > -180)) || ((lattice[i-1].dihedralAngle < 180) && (lattice[i-1].dihedralAngle > 116))) && (lattice[i].dihedralAngle < 115) && (lattice[i].dihedralAngle > 34) &&	(lattice[i+1].dihedralAngle < -38) && (lattice[i+1].dihedralAngle > -116) && (((lattice_tminus1[i-1].dihedralAngle < -117) && (lattice_tminus1[i-1].dihedralAngle > -180)) || ((lattice_tminus1[i-1].dihedralAngle < 180) && (lattice_tminus1[i-1].dihedralAngle > 116))) && (lattice_tminus1[i].dihedralAngle < 115) && (lattice_tminus1[i].dihedralAngle > 34) &&	(lattice_tminus1[i+1].dihedralAngle < -38) && (lattice_tminus1[i+1].dihedralAngle > -116))
		// {
		// 	possibleJumps.Gp1Gp++;
		// }
		// else if ((((lattice[i-1].dihedralAngle < -117) && (lattice[i-1].dihedralAngle > -180)) || ((lattice[i-1].dihedralAngle < 180) && (lattice[i-1].dihedralAngle > 116))) && (lattice[i].dihedralAngle < 115) && (lattice[i].dihedralAngle > 34) &&	(lattice[i+1].dihedralAngle < -38) && (lattice[i+1].dihedralAngle > -116) && (((lattice_tminus1[i-1].dihedralAngle < -117) && (lattice_tminus1[i-1].dihedralAngle > -180)) || ((lattice_tminus1[i-1].dihedralAngle < 180) && (lattice_tminus1[i-1].dihedralAngle > 116))) && (lattice_tminus1[i].dihedralAngle < -38) && (lattice_tminus1[i].dihedralAngle > -116) && (lattice_tminus1[i+1].dihedralAngle < -38) && (lattice_tminus1[i+1].dihedralAngle > -116))
		// {
		// 	possibleJumps.Gp1Gm++;
		// }
		// else if ((((lattice[i-1].dihedralAngle < -117) && (lattice[i-1].dihedralAngle > -180)) || ((lattice[i-1].dihedralAngle < 180) && (lattice[i-1].dihedralAngle > 116))) && (lattice[i].dihedralAngle < 115) && (lattice[i].dihedralAngle > 34) && (lattice[i+1].dihedralAngle < -38) && (lattice[i+1].dihedralAngle > -116) && (((lattice_tminus1[i-1].dihedralAngle < -117) && (lattice_tminus1[i-1].dihedralAngle > -180)) || ((lattice_tminus1[i-1].dihedralAngle < 180) && (lattice_tminus1[i-1].dihedralAngle > 116))) && (((lattice_tminus1[i].dihedralAngle < -117) && (lattice_tminus1[i].dihedralAngle > -180)) || ((lattice_tminus1[i].dihedralAngle < 180) && (lattice_tminus1[i].dihedralAngle > 116))) && (lattice_tminus1[i+1].dihedralAngle < -38) && (lattice_tminus1[i+1].dihedralAngle > -116))
		// {
		// 	possibleJumps.Gp1T++;
		// }
		// else if ((((lattice[i-1].dihedralAngle < -117) && (lattice[i-1].dihedralAngle > -180)) || ((lattice[i-1].dihedralAngle < 180) && (lattice[i-1].dihedralAngle > 116))) && (lattice[i].dihedralAngle < 115) && (lattice[i].dihedralAngle > 34) && (lattice[i+1].dihedralAngle < 115) && (lattice[i+1].dihedralAngle > 34) && (((lattice_tminus1[i-1].dihedralAngle < -117) && (lattice_tminus1[i-1].dihedralAngle > -180)) || ((lattice_tminus1[i-1].dihedralAngle < 180) && (lattice_tminus1[i-1].dihedralAngle > 116))) && (lattice_tminus1[i].dihedralAngle < 115) && (lattice_tminus1[i].dihedralAngle > 34) && (lattice_tminus1[i+1].dihedralAngle < 115) && (lattice_tminus1[i+1].dihedralAngle > 34))
		// {
		// 	possibleJumps.Gp2Gp++;
		// }
		// else if ((((lattice[i-1].dihedralAngle < -117) && (lattice[i-1].dihedralAngle > -180)) || ((lattice[i-1].dihedralAngle < 180) && (lattice[i-1].dihedralAngle > 116))) && (lattice[i].dihedralAngle < 115) && (lattice[i].dihedralAngle > 34) && (lattice[i+1].dihedralAngle < 115) && (lattice[i+1].dihedralAngle > 34) && (((lattice_tminus1[i-1].dihedralAngle < -117) && (lattice_tminus1[i-1].dihedralAngle > -180)) || ((lattice_tminus1[i-1].dihedralAngle < 180) && (lattice_tminus1[i-1].dihedralAngle > 116))) && (((lattice_tminus1[i].dihedralAngle < -117) && (lattice_tminus1[i].dihedralAngle > -180)) || ((lattice_tminus1[i].dihedralAngle < 180) && (lattice_tminus1[i].dihedralAngle > 116))) && (lattice_tminus1[i+1].dihedralAngle < 115) && (lattice_tminus1[i+1].dihedralAngle > 34))
		// {
		// 	possibleJumps.Gp2T++;
		// }
		// else if ((((lattice[i-1].dihedralAngle < -117) && (lattice[i-1].dihedralAngle > -180)) || ((lattice[i-1].dihedralAngle < 180) && (lattice[i-1].dihedralAngle > 116))) && (lattice[i].dihedralAngle < 115) && (lattice[i].dihedralAngle > 34) &&	(lattice[i+1].dihedralAngle < 115) && (lattice[i+1].dihedralAngle > 34) && (((lattice_tminus1[i-1].dihedralAngle < -117) && (lattice_tminus1[i-1].dihedralAngle > -180)) || ((lattice_tminus1[i-1].dihedralAngle < 180) && (lattice_tminus1[i-1].dihedralAngle > 116))) && (lattice_tminus1[i].dihedralAngle < -38) && (lattice_tminus1[i].dihedralAngle > -116) &&	(lattice_tminus1[i+1].dihedralAngle < 115) && (lattice_tminus1[i+1].dihedralAngle > 34))
		// {
		// 	possibleJumps.Gp2Gm++;
		// }
		// else if ((((lattice[i-1].dihedralAngle < -117) && (lattice[i-1].dihedralAngle > -180)) || ((lattice[i-1].dihedralAngle < 180) && (lattice[i-1].dihedralAngle > 116))) && (lattice[i].dihedralAngle < 115) && (lattice[i].dihedralAngle > 34) &&	(((lattice[i+1].dihedralAngle < -117) && (lattice[i+1].dihedralAngle > -180)) || ((lattice[i+1].dihedralAngle < 180) && (lattice[i+1].dihedralAngle > 116))) &&	(((lattice_tminus1[i-1].dihedralAngle < -117) && (lattice_tminus1[i-1].dihedralAngle > -180)) || ((lattice_tminus1[i-1].dihedralAngle < 180) && (lattice_tminus1[i-1].dihedralAngle > 116))) &&	(((lattice_tminus1[i].dihedralAngle < -117) && (lattice_tminus1[i].dihedralAngle > -180)) || ((lattice_tminus1[i].dihedralAngle < 180) && (lattice_tminus1[i].dihedralAngle > 116))) && (((lattice_tminus1[i+1].dihedralAngle < -117) && (lattice_tminus1[i+1].dihedralAngle > -180)) || ((lattice_tminus1[i+1].dihedralAngle < 180) && (lattice_tminus1[i+1].dihedralAngle > 116))))
		// {
		// 	possibleJumps.Gp3T++;
		// }
		// else if ((((lattice[i-1].dihedralAngle < -117) && (lattice[i-1].dihedralAngle > -180)) || ((lattice[i-1].dihedralAngle < 180) && (lattice[i-1].dihedralAngle > 116))) && (lattice[i].dihedralAngle < 115) && (lattice[i].dihedralAngle > 34) &&	(((lattice[i+1].dihedralAngle < -117) && (lattice[i+1].dihedralAngle > -180)) || ((lattice[i+1].dihedralAngle < 180) && (lattice[i+1].dihedralAngle > 116))) &&	(((lattice_tminus1[i-1].dihedralAngle < -117) && (lattice_tminus1[i-1].dihedralAngle > -180)) || ((lattice_tminus1[i-1].dihedralAngle < 180) && (lattice_tminus1[i-1].dihedralAngle > 116))) &&	(lattice_tminus1[i].dihedralAngle < 115) && (lattice_tminus1[i].dihedralAngle > 34) && (((lattice_tminus1[i+1].dihedralAngle < -117) && (lattice_tminus1[i+1].dihedralAngle > -180)) || ((lattice_tminus1[i+1].dihedralAngle < 180) && (lattice_tminus1[i+1].dihedralAngle > 116))))
		// {
		// 	possibleJumps.Gp3Gp++;
		// }
		// else if ((((lattice[i-1].dihedralAngle < -117) && (lattice[i-1].dihedralAngle > -180)) || ((lattice[i-1].dihedralAngle < 180) && (lattice[i-1].dihedralAngle > 116))) && (lattice[i].dihedralAngle < 115) && (lattice[i].dihedralAngle > 34) &&	(((lattice[i+1].dihedralAngle < -117) && (lattice[i+1].dihedralAngle > -180)) || ((lattice[i+1].dihedralAngle < 180) && (lattice[i+1].dihedralAngle > 116))) &&	(((lattice_tminus1[i-1].dihedralAngle < -117) && (lattice_tminus1[i-1].dihedralAngle > -180)) || ((lattice_tminus1[i-1].dihedralAngle < 180) && (lattice_tminus1[i-1].dihedralAngle > 116))) &&	(lattice_tminus1[i].dihedralAngle < -38) && (lattice_tminus1[i].dihedralAngle > -116) && (((lattice_tminus1[i+1].dihedralAngle < -117) && (lattice_tminus1[i+1].dihedralAngle > -180)) || ((lattice_tminus1[i+1].dihedralAngle < 180) && (lattice_tminus1[i+1].dihedralAngle > 116))))
		// {
		// 	possibleJumps.Gp3Gm++;
		// }
		// else if ((lattice[i-1].dihedralAngle < -38) && (lattice[i-1].dihedralAngle > -116) && (lattice[i].dihedralAngle < 115) && (lattice[i].dihedralAngle > 34) && (((lattice[i+1].dihedralAngle < -117) && (lattice[i+1].dihedralAngle > -180)) || ((lattice[i+1].dihedralAngle < 180) && (lattice[i+1].dihedralAngle > 116))) && (lattice_tminus1[i-1].dihedralAngle < -38) && (lattice_tminus1[i-1].dihedralAngle > -116) && (((lattice_tminus1[i].dihedralAngle < -117) && (lattice_tminus1[i].dihedralAngle > -180)) || ((lattice_tminus1[i].dihedralAngle < 180) && (lattice_tminus1[i].dihedralAngle > 116))) && (((lattice_tminus1[i+1].dihedralAngle < -117) && (lattice_tminus1[i+1].dihedralAngle > -180)) || ((lattice_tminus1[i+1].dihedralAngle < 180) && (lattice_tminus1[i+1].dihedralAngle > 116))))
		// {
		// 	possibleJumps.Gp4T++;
		// }
		// else if ((lattice[i-1].dihedralAngle < -38) && (lattice[i-1].dihedralAngle > -116) && (lattice[i].dihedralAngle < 115) && (lattice[i].dihedralAngle > 34) && (((lattice[i+1].dihedralAngle < -117) && (lattice[i+1].dihedralAngle > -180)) || ((lattice[i+1].dihedralAngle < 180) && (lattice[i+1].dihedralAngle > 116))) && (lattice_tminus1[i-1].dihedralAngle < -38) && (lattice_tminus1[i-1].dihedralAngle > -116) && (lattice_tminus1[i].dihedralAngle < -38) && (lattice_tminus1[i].dihedralAngle > -116) && (((lattice_tminus1[i+1].dihedralAngle < -117) && (lattice_tminus1[i+1].dihedralAngle > -180)) || ((lattice_tminus1[i+1].dihedralAngle < 180) && (lattice_tminus1[i+1].dihedralAngle > 116))))
		// {
		// 	possibleJumps.Gp4Gm++;
		// }
		// else if ((lattice[i-1].dihedralAngle < -38) && (lattice[i-1].dihedralAngle > -116) && (lattice[i].dihedralAngle < 115) && (lattice[i].dihedralAngle > 34) && (((lattice[i+1].dihedralAngle < -117) && (lattice[i+1].dihedralAngle > -180)) || ((lattice[i+1].dihedralAngle < 180) && (lattice[i+1].dihedralAngle > 116))) && (lattice_tminus1[i-1].dihedralAngle < -38) && (lattice_tminus1[i-1].dihedralAngle > -116) && (lattice_tminus1[i].dihedralAngle < 115) && (lattice_tminus1[i].dihedralAngle > 34) && (((lattice_tminus1[i+1].dihedralAngle < -117) && (lattice_tminus1[i+1].dihedralAngle > -180)) || ((lattice_tminus1[i+1].dihedralAngle < 180) && (lattice_tminus1[i+1].dihedralAngle > 116))))
		// {
		// 	possibleJumps.Gp4Gp++;
		// }
		// else if ((lattice[i-1].dihedralAngle < -38) && (lattice[i-1].dihedralAngle > -116) && (lattice[i].dihedralAngle < 115) && (lattice[i].dihedralAngle > 34) && (lattice[i+1].dihedralAngle < 115) && (lattice[i+1].dihedralAngle > 34) &&	(lattice_tminus1[i-1].dihedralAngle < -38) && (lattice_tminus1[i-1].dihedralAngle > -116) && (((lattice_tminus1[i].dihedralAngle < -117) && (lattice_tminus1[i].dihedralAngle > -180)) || ((lattice_tminus1[i].dihedralAngle < 180) && (lattice_tminus1[i].dihedralAngle > 116))) && (lattice_tminus1[i+1].dihedralAngle < 115) && (lattice_tminus1[i+1].dihedralAngle > 34))
		// {
		// 	possibleJumps.Gp5T++;
		// }
		// else if ((lattice[i-1].dihedralAngle < -38) && (lattice[i-1].dihedralAngle > -116) && (lattice[i].dihedralAngle < 115) && (lattice[i].dihedralAngle > 34) && (lattice[i+1].dihedralAngle < 115) && (lattice[i+1].dihedralAngle > 34) &&	(lattice_tminus1[i-1].dihedralAngle < -38) && (lattice_tminus1[i-1].dihedralAngle > -116) && (lattice_tminus1[i].dihedralAngle < -38) && (lattice_tminus1[i].dihedralAngle > -116) && (lattice_tminus1[i+1].dihedralAngle < 115) && (lattice_tminus1[i+1].dihedralAngle > 34))
		// {
		// 	possibleJumps.Gp5Gm++;
		// }
		// else if ((lattice[i-1].dihedralAngle < -38) && (lattice[i-1].dihedralAngle > -116) && (lattice[i].dihedralAngle < 115) && (lattice[i].dihedralAngle > 34) && (lattice[i+1].dihedralAngle < 115) && (lattice[i+1].dihedralAngle > 34) &&	(lattice_tminus1[i-1].dihedralAngle < -38) && (lattice_tminus1[i-1].dihedralAngle > -116) && (lattice_tminus1[i].dihedralAngle < 115) && (lattice_tminus1[i].dihedralAngle > 34) &&	(lattice_tminus1[i+1].dihedralAngle < 115) && (lattice_tminus1[i+1].dihedralAngle > 34))
		// {
		// 	possibleJumps.Gp5Gp++;
		// }
		// else if ((lattice[i-1].dihedralAngle < -38) && (lattice[i-1].dihedralAngle > -116) && (lattice[i].dihedralAngle < 115) && (lattice[i].dihedralAngle > 34) && (lattice[i-1].dihedralAngle < -38) && (lattice[i-1].dihedralAngle > -116) && (lattice_tminus1[i-1].dihedralAngle < -38) && (lattice_tminus1[i-1].dihedralAngle > -116) && (((lattice_tminus1[i].dihedralAngle < -117) && (lattice_tminus1[i].dihedralAngle > -180)) || ((lattice_tminus1[i].dihedralAngle < 180) && (lattice_tminus1[i].dihedralAngle > 116))) && (lattice_tminus1[i-1].dihedralAngle < -38) && (lattice_tminus1[i-1].dihedralAngle > -116))
		// {
		// 	possibleJumps.Gp6T++;
		// }
		// else if ((lattice[i-1].dihedralAngle < -38) && (lattice[i-1].dihedralAngle > -116) && (lattice[i].dihedralAngle < 115) && (lattice[i].dihedralAngle > 34) && (lattice[i-1].dihedralAngle < -38) && (lattice[i-1].dihedralAngle > -116) && (lattice_tminus1[i-1].dihedralAngle < -38) && (lattice_tminus1[i-1].dihedralAngle > -116) && (lattice_tminus1[i].dihedralAngle < -38) && (lattice_tminus1[i].dihedralAngle > -116) &&	(lattice_tminus1[i-1].dihedralAngle < -38) && (lattice_tminus1[i-1].dihedralAngle > -116))
		// {
		// 	possibleJumps.Gp6Gm++;
		// }
		// else if ((lattice[i-1].dihedralAngle < -38) && (lattice[i-1].dihedralAngle > -116) && (lattice[i].dihedralAngle < 115) && (lattice[i].dihedralAngle > 34) && (lattice[i-1].dihedralAngle < -38) && (lattice[i-1].dihedralAngle > -116) && (lattice_tminus1[i-1].dihedralAngle < -38) && (lattice_tminus1[i-1].dihedralAngle > -116) && (lattice_tminus1[i].dihedralAngle < 115) && (lattice_tminus1[i].dihedralAngle > 34) && (lattice_tminus1[i-1].dihedralAngle < -38) && (lattice_tminus1[i-1].dihedralAngle > -116))
		// {
		// 	possibleJumps.Gp6Gp++;
		// }
		// else if ((lattice[i-1].dihedralAngle < 115) && (lattice[i-1].dihedralAngle > 34) &&	(lattice[i].dihedralAngle < 115) && (lattice[i].dihedralAngle > 34) && (((lattice[i+1].dihedralAngle < -117) && (lattice[i+1].dihedralAngle > -180)) || ((lattice[i+1].dihedralAngle < 180) && (lattice[i+1].dihedralAngle > 116))) && (lattice_tminus1[i-1].dihedralAngle < 115) && (lattice_tminus1[i-1].dihedralAngle > 34) && (lattice_tminus1[i].dihedralAngle < 115) && (lattice_tminus1[i].dihedralAngle > 34) && (((lattice_tminus1[i+1].dihedralAngle < -117) && (lattice_tminus1[i+1].dihedralAngle > -180)) || ((lattice_tminus1[i+1].dihedralAngle < 180) && (lattice_tminus1[i+1].dihedralAngle > 116))))
		// {
		// 	possibleJumps.Gp7Gp++;
		// }
		// else if ((lattice[i-1].dihedralAngle < 115) && (lattice[i-1].dihedralAngle > 34) &&	(lattice[i].dihedralAngle < 115) && (lattice[i].dihedralAngle > 34) && (((lattice[i+1].dihedralAngle < -117) && (lattice[i+1].dihedralAngle > -180)) || ((lattice[i+1].dihedralAngle < 180) && (lattice[i+1].dihedralAngle > 116))) && (lattice_tminus1[i-1].dihedralAngle < 115) && (lattice_tminus1[i-1].dihedralAngle > 34) && (lattice_tminus1[i].dihedralAngle < -38) && (lattice_tminus1[i].dihedralAngle > -116) && (((lattice_tminus1[i+1].dihedralAngle < -117) && (lattice_tminus1[i+1].dihedralAngle > -180)) || ((lattice_tminus1[i+1].dihedralAngle < 180) && (lattice_tminus1[i+1].dihedralAngle > 116))))
		// {
		// 	possibleJumps.Gp7Gm++;
		// }
		// else if ((lattice[i-1].dihedralAngle < 115) && (lattice[i-1].dihedralAngle > 34) && (lattice[i].dihedralAngle < 115) && (lattice[i].dihedralAngle > 34) && (((lattice[i+1].dihedralAngle < -117) && (lattice[i+1].dihedralAngle > -180)) || ((lattice[i+1].dihedralAngle < 180) && (lattice[i+1].dihedralAngle > 116))) && (lattice_tminus1[i-1].dihedralAngle < 115) && (lattice_tminus1[i-1].dihedralAngle > 34) && (((lattice_tminus1[i].dihedralAngle < -117) && (lattice_tminus1[i].dihedralAngle > -180)) || ((lattice_tminus1[i].dihedralAngle < 180) && (lattice_tminus1[i].dihedralAngle > 116))) && (((lattice_tminus1[i+1].dihedralAngle < -117) && (lattice_tminus1[i+1].dihedralAngle > -180)) || ((lattice_tminus1[i+1].dihedralAngle < 180) && (lattice_tminus1[i+1].dihedralAngle > 116))))
		// {
		// 	possibleJumps.Gp7T++;
		// }
		// else if ((lattice[i-1].dihedralAngle < 115) && (lattice[i-1].dihedralAngle > 34) &&	(lattice[i].dihedralAngle < 115) && (lattice[i].dihedralAngle > 34) && (lattice[i+1].dihedralAngle < 115) && (lattice[i+1].dihedralAngle > 34) && (lattice_tminus1[i-1].dihedralAngle < 115) && (lattice_tminus1[i-1].dihedralAngle > 34) && (((lattice_tminus1[i].dihedralAngle < -117) && (lattice_tminus1[i].dihedralAngle > -180)) || ((lattice_tminus1[i].dihedralAngle < 180) && (lattice_tminus1[i].dihedralAngle > 116))) && (lattice_tminus1[i+1].dihedralAngle < 115) && (lattice_tminus1[i+1].dihedralAngle > 34))
		// {
		// 	possibleJumps.Gp8T++;
		// }
		// else if ((lattice[i-1].dihedralAngle < 115) && (lattice[i-1].dihedralAngle > 34) &&	(lattice[i].dihedralAngle < 115) && (lattice[i].dihedralAngle > 34) && (lattice[i+1].dihedralAngle < 115) && (lattice[i+1].dihedralAngle > 34) && (lattice_tminus1[i-1].dihedralAngle < 115) && (lattice_tminus1[i-1].dihedralAngle > 34) && (lattice_tminus1[i].dihedralAngle < 115) && (lattice_tminus1[i].dihedralAngle > 34) &&	(lattice_tminus1[i+1].dihedralAngle < 115) && (lattice_tminus1[i+1].dihedralAngle > 34))
		// {
		// 	possibleJumps.Gp8Gp++;
		// }
		// else if ((lattice[i-1].dihedralAngle < 115) && (lattice[i-1].dihedralAngle > 34) &&	(lattice[i].dihedralAngle < 115) && (lattice[i].dihedralAngle > 34) && (lattice[i+1].dihedralAngle < 115) && (lattice[i+1].dihedralAngle > 34) && (lattice_tminus1[i-1].dihedralAngle < 115) && (lattice_tminus1[i-1].dihedralAngle > 34) && (lattice_tminus1[i].dihedralAngle < -38) && (lattice_tminus1[i].dihedralAngle > -116) && (lattice_tminus1[i+1].dihedralAngle < 115) && (lattice_tminus1[i+1].dihedralAngle > 34))
		// {
		// 	possibleJumps.Gp8Gm++;
		// }
		// else if ((lattice[i-1].dihedralAngle < 115) && (lattice[i-1].dihedralAngle > 34) &&	(lattice[i].dihedralAngle < 115) && (lattice[i].dihedralAngle > 34) && (lattice[i+1].dihedralAngle < -38) && (lattice[i+1].dihedralAngle > -116) &&	(lattice_tminus1[i-1].dihedralAngle < 115) && (lattice_tminus1[i-1].dihedralAngle > 34) && (lattice_tminus1[i].dihedralAngle < -38) && (lattice_tminus1[i].dihedralAngle > -116) && (lattice_tminus1[i+1].dihedralAngle < -38) && (lattice_tminus1[i+1].dihedralAngle > -116))
		// {
		// 	possibleJumps.Gp9Gm++;
		// }
		// else if ((lattice[i-1].dihedralAngle < 115) && (lattice[i-1].dihedralAngle > 34) &&	(lattice[i].dihedralAngle < 115) && (lattice[i].dihedralAngle > 34) && (lattice[i+1].dihedralAngle < -38) && (lattice[i+1].dihedralAngle > -116) &&	(lattice_tminus1[i-1].dihedralAngle < 115) && (lattice_tminus1[i-1].dihedralAngle > 34) && (((lattice_tminus1[i].dihedralAngle < -117) && (lattice_tminus1[i].dihedralAngle > -180)) || ((lattice_tminus1[i].dihedralAngle < 180) && (lattice_tminus1[i].dihedralAngle > 116))) && (lattice_tminus1[i+1].dihedralAngle < -38) && (lattice_tminus1[i+1].dihedralAngle > -116))
		// {
		// 	possibleJumps.Gp9T++;
		// }
		// else if ((lattice[i-1].dihedralAngle < 115) && (lattice[i-1].dihedralAngle > 34) &&	(lattice[i].dihedralAngle < 115) && (lattice[i].dihedralAngle > 34) && (lattice[i+1].dihedralAngle < -38) && (lattice[i+1].dihedralAngle > -116) && (lattice_tminus1[i-1].dihedralAngle < 115) && (lattice_tminus1[i-1].dihedralAngle > 34) && (lattice_tminus1[i].dihedralAngle < 115) && (lattice_tminus1[i].dihedralAngle > 34) && (lattice_tminus1[i+1].dihedralAngle < -38) && (lattice_tminus1[i+1].dihedralAngle > -116))
		// {
		// 	possibleJumps.Gp9Gp++;
		// }
		// else if ((((lattice[i-1].dihedralAngle < -117) && (lattice[i-1].dihedralAngle > -180)) || ((lattice[i-1].dihedralAngle < 180) && (lattice[i-1].dihedralAngle > 116))) && (lattice[i].dihedralAngle < -38) && (lattice[i].dihedralAngle > -116) && (((lattice[i+1].dihedralAngle < -117) && (lattice[i+1].dihedralAngle > -180)) || ((lattice[i+1].dihedralAngle < 180) && (lattice[i+1].dihedralAngle > 116))) && (((lattice_tminus1[i-1].dihedralAngle < -117) && (lattice_tminus1[i-1].dihedralAngle > -180)) || ((lattice_tminus1[i-1].dihedralAngle < 180) && (lattice_tminus1[i-1].dihedralAngle > 116))) && (lattice_tminus1[i].dihedralAngle < 115) && (lattice_tminus1[i].dihedralAngle > 34) && (((lattice_tminus1[i+1].dihedralAngle < -117) && (lattice_tminus1[i+1].dihedralAngle > -180)) || ((lattice_tminus1[i+1].dihedralAngle < 180) && (lattice_tminus1[i+1].dihedralAngle > 116))))
		// {
		// 	possibleJumps.Gm1Gp++;
		// }
		// else if ((((lattice[i-1].dihedralAngle < -117) && (lattice[i-1].dihedralAngle > -180)) || ((lattice[i-1].dihedralAngle < 180) && (lattice[i-1].dihedralAngle > 116))) && (lattice[i].dihedralAngle < -38) && (lattice[i].dihedralAngle > -116) && (((lattice[i+1].dihedralAngle < -117) && (lattice[i+1].dihedralAngle > -180)) || ((lattice[i+1].dihedralAngle < 180) && (lattice[i+1].dihedralAngle > 116))) && (((lattice_tminus1[i-1].dihedralAngle < -117) && (lattice_tminus1[i-1].dihedralAngle > -180)) || ((lattice_tminus1[i-1].dihedralAngle < 180) && (lattice_tminus1[i-1].dihedralAngle > 116))) && (((lattice_tminus1[i].dihedralAngle < -117) && (lattice_tminus1[i].dihedralAngle > -180)) || ((lattice_tminus1[i].dihedralAngle < 180) && (lattice_tminus1[i].dihedralAngle > 116))) && (((lattice_tminus1[i+1].dihedralAngle < -117) && (lattice_tminus1[i+1].dihedralAngle > -180)) || ((lattice_tminus1[i+1].dihedralAngle < 180) && (lattice_tminus1[i+1].dihedralAngle > 116))))
		// {
		// 	possibleJumps.Gm1T++;
		// }
		// else if ((((lattice[i-1].dihedralAngle < -117) && (lattice[i-1].dihedralAngle > -180)) || ((lattice[i-1].dihedralAngle < 180) && (lattice[i-1].dihedralAngle > 116))) && (lattice[i].dihedralAngle < -38) && (lattice[i].dihedralAngle > -116) && (((lattice[i+1].dihedralAngle < -117) && (lattice[i+1].dihedralAngle > -180)) || ((lattice[i+1].dihedralAngle < 180) && (lattice[i+1].dihedralAngle > 116))) && (((lattice_tminus1[i-1].dihedralAngle < -117) && (lattice_tminus1[i-1].dihedralAngle > -180)) || ((lattice_tminus1[i-1].dihedralAngle < 180) && (lattice_tminus1[i-1].dihedralAngle > 116))) && (lattice_tminus1[i].dihedralAngle < -38) && (lattice_tminus1[i].dihedralAngle > -116) && (((lattice_tminus1[i+1].dihedralAngle < -117) && (lattice_tminus1[i+1].dihedralAngle > -180)) || ((lattice_tminus1[i+1].dihedralAngle < 180) && (lattice_tminus1[i+1].dihedralAngle > 116))))
		// {
		// 	possibleJumps.Gm1Gm++;
		// }
		// else if ((((lattice[i-1].dihedralAngle < -117) && (lattice[i-1].dihedralAngle > -180)) || ((lattice[i-1].dihedralAngle < 180) && (lattice[i-1].dihedralAngle > 116))) && (lattice[i].dihedralAngle < -38) && (lattice[i].dihedralAngle > -116) && (lattice[i+1].dihedralAngle < 115) && (lattice[i+1].dihedralAngle > 34) && (((lattice_tminus1[i-1].dihedralAngle < -117) && (lattice_tminus1[i-1].dihedralAngle > -180)) || ((lattice_tminus1[i-1].dihedralAngle < 180) && (lattice_tminus1[i-1].dihedralAngle > 116))) && (((lattice_tminus1[i].dihedralAngle < -117) && (lattice_tminus1[i].dihedralAngle > -180)) || ((lattice_tminus1[i].dihedralAngle < 180) && (lattice_tminus1[i].dihedralAngle > 116))) && (lattice_tminus1[i+1].dihedralAngle < 115) && (lattice_tminus1[i+1].dihedralAngle > 34))
		// {
		// 	possibleJumps.Gm2T++;
		// }
		// else if ((((lattice[i-1].dihedralAngle < -117) && (lattice[i-1].dihedralAngle > -180)) || ((lattice[i-1].dihedralAngle < 180) && (lattice[i-1].dihedralAngle > 116))) && (lattice[i].dihedralAngle < -38) && (lattice[i].dihedralAngle > -116) && (lattice[i+1].dihedralAngle < 115) && (lattice[i+1].dihedralAngle > 34) && (((lattice_tminus1[i-1].dihedralAngle < -117) && (lattice_tminus1[i-1].dihedralAngle > -180)) || ((lattice_tminus1[i-1].dihedralAngle < 180) && (lattice_tminus1[i-1].dihedralAngle > 116))) && (lattice_tminus1[i].dihedralAngle < 115) && (lattice_tminus1[i].dihedralAngle > 34) &&	(lattice_tminus1[i+1].dihedralAngle < 115) && (lattice_tminus1[i+1].dihedralAngle > 34))
		// {
		// 	possibleJumps.Gm2Gp++;
		// }
		// else if ((((lattice[i-1].dihedralAngle < -117) && (lattice[i-1].dihedralAngle > -180)) || ((lattice[i-1].dihedralAngle < 180) && (lattice[i-1].dihedralAngle > 116))) && (lattice[i].dihedralAngle < -38) && (lattice[i].dihedralAngle > -116) && (lattice[i+1].dihedralAngle < 115) && (lattice[i+1].dihedralAngle > 34) && (((lattice_tminus1[i-1].dihedralAngle < -117) && (lattice_tminus1[i-1].dihedralAngle > -180)) || ((lattice_tminus1[i-1].dihedralAngle < 180) && (lattice_tminus1[i-1].dihedralAngle > 116))) && (lattice_tminus1[i].dihedralAngle < -38) && (lattice_tminus1[i].dihedralAngle > -116) && (lattice_tminus1[i+1].dihedralAngle < 115) && (lattice_tminus1[i+1].dihedralAngle > 34))
		// {
		// 	possibleJumps.Gm2Gm++;
		// }
		// else if ((((lattice[i-1].dihedralAngle < -117) && (lattice[i-1].dihedralAngle > -180)) || ((lattice[i-1].dihedralAngle < 180) && (lattice[i-1].dihedralAngle > 116))) && (lattice[i].dihedralAngle < -38) && (lattice[i].dihedralAngle > -116) && (lattice[i+1].dihedralAngle < -38) && (lattice[i+1].dihedralAngle > -116) && (((lattice_tminus1[i-1].dihedralAngle < -117) && (lattice_tminus1[i-1].dihedralAngle > -180)) || ((lattice_tminus1[i-1].dihedralAngle < 180) && (lattice_tminus1[i-1].dihedralAngle > 116))) && (lattice_tminus1[i].dihedralAngle < 115) && (lattice_tminus1[i].dihedralAngle > 34) && (lattice_tminus1[i+1].dihedralAngle < -38) && (lattice_tminus1[i+1].dihedralAngle > -116))
		// {
		// 	possibleJumps.Gm3Gp++;
		// }
		// else if ((((lattice[i-1].dihedralAngle < -117) && (lattice[i-1].dihedralAngle > -180)) || ((lattice[i-1].dihedralAngle < 180) && (lattice[i-1].dihedralAngle > 116))) && (lattice[i].dihedralAngle < -38) && (lattice[i].dihedralAngle > -116) && (lattice[i+1].dihedralAngle < -38) && (lattice[i+1].dihedralAngle > -116) && (((lattice_tminus1[i-1].dihedralAngle < -117) && (lattice_tminus1[i-1].dihedralAngle > -180)) || ((lattice_tminus1[i-1].dihedralAngle < 180) && (lattice_tminus1[i-1].dihedralAngle > 116))) && (lattice_tminus1[i].dihedralAngle < -38) && (lattice_tminus1[i].dihedralAngle > -116) && (lattice_tminus1[i+1].dihedralAngle < -38) && (lattice_tminus1[i+1].dihedralAngle > -116))
		// {
		// 	possibleJumps.Gm3Gm++;
		// }
		// else if ((((lattice[i-1].dihedralAngle < -117) && (lattice[i-1].dihedralAngle > -180)) || ((lattice[i-1].dihedralAngle < 180) && (lattice[i-1].dihedralAngle > 116))) && (lattice[i].dihedralAngle < -38) && (lattice[i].dihedralAngle > -116) && (lattice[i+1].dihedralAngle < -38) && (lattice[i+1].dihedralAngle > -116) && (((lattice_tminus1[i-1].dihedralAngle < -117) && (lattice_tminus1[i-1].dihedralAngle > -180)) || ((lattice_tminus1[i-1].dihedralAngle < 180) && (lattice_tminus1[i-1].dihedralAngle > 116))) && (((lattice_tminus1[i].dihedralAngle < -117) && (lattice_tminus1[i].dihedralAngle > -180)) || ((lattice_tminus1[i].dihedralAngle < 180) && (lattice_tminus1[i].dihedralAngle > 116))) && (lattice_tminus1[i+1].dihedralAngle < -38) && (lattice_tminus1[i+1].dihedralAngle > -116))
		// {
		// 	possibleJumps.Gm3T++;
		// }
		// else if ((lattice[i-1].dihedralAngle < -38) && (lattice[i-1].dihedralAngle > -116) && (lattice[i].dihedralAngle < -38) && (lattice[i].dihedralAngle > -116) && (((lattice[i+1].dihedralAngle < -117) && (lattice[i+1].dihedralAngle > -180)) || ((lattice[i+1].dihedralAngle < 180) && (lattice[i+1].dihedralAngle > 116))) && (lattice_tminus1[i-1].dihedralAngle < -38) && (lattice_tminus1[i-1].dihedralAngle > -116) && (lattice_tminus1[i].dihedralAngle < -38) && (lattice_tminus1[i].dihedralAngle > -116) && (((lattice_tminus1[i+1].dihedralAngle < -117) && (lattice_tminus1[i+1].dihedralAngle > -180)) || ((lattice_tminus1[i+1].dihedralAngle < 180) && (lattice_tminus1[i+1].dihedralAngle > 116))))
		// {
		// 	possibleJumps.Gm4Gm++;
		// }
		// else if ((lattice[i-1].dihedralAngle < -38) && (lattice[i-1].dihedralAngle > -116) && (lattice[i].dihedralAngle < -38) && (lattice[i].dihedralAngle > -116) && (((lattice[i+1].dihedralAngle < -117) && (lattice[i+1].dihedralAngle > -180)) || ((lattice[i+1].dihedralAngle < 180) && (lattice[i+1].dihedralAngle > 116))) && (lattice_tminus1[i-1].dihedralAngle < -38) && (lattice_tminus1[i-1].dihedralAngle > -116) && (lattice_tminus1[i].dihedralAngle < 115) && (lattice_tminus1[i].dihedralAngle > 34) && (((lattice_tminus1[i+1].dihedralAngle < -117) && (lattice_tminus1[i+1].dihedralAngle > -180)) || ((lattice_tminus1[i+1].dihedralAngle < 180) && (lattice_tminus1[i+1].dihedralAngle > 116))))
		// {
		// 	possibleJumps.Gm4Gp++;
		// }
		// else if ((lattice[i-1].dihedralAngle < -38) && (lattice[i-1].dihedralAngle > -116) && (lattice[i].dihedralAngle < -38) && (lattice[i].dihedralAngle > -116) && (((lattice[i+1].dihedralAngle < -117) && (lattice[i+1].dihedralAngle > -180)) || ((lattice[i+1].dihedralAngle < 180) && (lattice[i+1].dihedralAngle > 116))) && (lattice_tminus1[i-1].dihedralAngle < -38) && (lattice_tminus1[i-1].dihedralAngle > -116) && (((lattice_tminus1[i].dihedralAngle < -117) && (lattice_tminus1[i].dihedralAngle > -180)) || ((lattice_tminus1[i].dihedralAngle < 180) && (lattice_tminus1[i].dihedralAngle > 116))) && (((lattice_tminus1[i+1].dihedralAngle < -117) && (lattice_tminus1[i+1].dihedralAngle > -180)) || ((lattice_tminus1[i+1].dihedralAngle < 180) && (lattice_tminus1[i+1].dihedralAngle > 116))))
		// {
		// 	possibleJumps.Gm4T++;
		// }
		// else if ((lattice[i-1].dihedralAngle < -38) && (lattice[i-1].dihedralAngle > -116) && (lattice[i].dihedralAngle < -38) && (lattice[i].dihedralAngle > -116) && (lattice[i+1].dihedralAngle < 115) && (lattice[i+1].dihedralAngle > 34) && (lattice_tminus1[i-1].dihedralAngle < -38) && (lattice_tminus1[i-1].dihedralAngle > -116) && (lattice_tminus1[i].dihedralAngle < 115) && (lattice_tminus1[i].dihedralAngle > 34) && (lattice_tminus1[i+1].dihedralAngle < 115) && (lattice_tminus1[i+1].dihedralAngle > 34))
		// {
		// 	possibleJumps.Gm5Gp++;
		// }
		// else if ((lattice[i-1].dihedralAngle < -38) && (lattice[i-1].dihedralAngle > -116) && (lattice[i].dihedralAngle < -38) && (lattice[i].dihedralAngle > -116) && (lattice[i+1].dihedralAngle < 115) && (lattice[i+1].dihedralAngle > 34) && (lattice_tminus1[i-1].dihedralAngle < -38) && (lattice_tminus1[i-1].dihedralAngle > -116) && (lattice_tminus1[i].dihedralAngle < -38) && (lattice_tminus1[i].dihedralAngle > -116) &&	(lattice_tminus1[i+1].dihedralAngle < 115) && (lattice_tminus1[i+1].dihedralAngle > 34))
		// {
		// 	possibleJumps.Gm5Gm++;
		// }
		// else if ((lattice[i-1].dihedralAngle < -38) && (lattice[i-1].dihedralAngle > -116) && (lattice[i].dihedralAngle < -38) && (lattice[i].dihedralAngle > -116) && (lattice[i+1].dihedralAngle < 115) && (lattice[i+1].dihedralAngle > 34) && (lattice_tminus1[i-1].dihedralAngle < -38) && (lattice_tminus1[i-1].dihedralAngle > -116) && (((lattice_tminus1[i].dihedralAngle < -117) && (lattice_tminus1[i].dihedralAngle > -180)) || ((lattice_tminus1[i].dihedralAngle < 180) && (lattice_tminus1[i].dihedralAngle > 116))) && (lattice_tminus1[i+1].dihedralAngle < 115) && (lattice_tminus1[i+1].dihedralAngle > 34))
		// {
		// 	possibleJumps.Gm5T++;
		// }
		// else if ((lattice[i-1].dihedralAngle < -38) && (lattice[i-1].dihedralAngle > -116) && (lattice[i].dihedralAngle < -38) && (lattice[i].dihedralAngle > -116) && (lattice[i+1].dihedralAngle < -38) && (lattice[i+1].dihedralAngle > -116) &&	(lattice_tminus1[i-1].dihedralAngle < -38) && (lattice_tminus1[i-1].dihedralAngle > -116) && (lattice_tminus1[i].dihedralAngle < 115) && (lattice_tminus1[i].dihedralAngle > 34) && (lattice_tminus1[i+1].dihedralAngle < -38) && (lattice_tminus1[i+1].dihedralAngle > -116))
		// {
		// 	possibleJumps.Gm6Gp++;
		// }
		// else if ((lattice[i-1].dihedralAngle < -38) && (lattice[i-1].dihedralAngle > -116) && (lattice[i].dihedralAngle < -38) && (lattice[i].dihedralAngle > -116) && (lattice[i+1].dihedralAngle < -38) && (lattice[i+1].dihedralAngle > -116) &&	(lattice_tminus1[i-1].dihedralAngle < -38) && (lattice_tminus1[i-1].dihedralAngle > -116) && (lattice_tminus1[i].dihedralAngle < -38) && (lattice_tminus1[i].dihedralAngle > -116) && (lattice_tminus1[i+1].dihedralAngle < -38) && (lattice_tminus1[i+1].dihedralAngle > -116))
		// {
		// 	possibleJumps.Gm6Gm++;
		// }
		// else if ((lattice[i-1].dihedralAngle < -38) && (lattice[i-1].dihedralAngle > -116) && (lattice[i].dihedralAngle < -38) && (lattice[i].dihedralAngle > -116) && (lattice[i+1].dihedralAngle < -38) && (lattice[i+1].dihedralAngle > -116) &&	(lattice_tminus1[i-1].dihedralAngle < -38) && (lattice_tminus1[i-1].dihedralAngle > -116) && (((lattice_tminus1[i].dihedralAngle < -117) && (lattice_tminus1[i].dihedralAngle > -180)) || ((lattice_tminus1[i].dihedralAngle < 180) && (lattice_tminus1[i].dihedralAngle > 116))) && (lattice_tminus1[i+1].dihedralAngle < -38) && (lattice_tminus1[i+1].dihedralAngle > -116))
		// {
		// 	possibleJumps.Gm6T++;
		// }
		// else if ((lattice[i-1].dihedralAngle < 115) && (lattice[i-1].dihedralAngle > 34) &&	(lattice[i].dihedralAngle < -38) && (lattice[i].dihedralAngle > -116) && (((lattice[i+1].dihedralAngle < -117) && (lattice[i+1].dihedralAngle > -180)) || ((lattice[i+1].dihedralAngle < 180) && (lattice[i+1].dihedralAngle > 116))) && (lattice_tminus1[i-1].dihedralAngle < 115) && (lattice_tminus1[i-1].dihedralAngle > 34) && (((lattice_tminus1[i].dihedralAngle < -117) && (lattice_tminus1[i].dihedralAngle > -180)) || ((lattice_tminus1[i].dihedralAngle < 180) && (lattice_tminus1[i].dihedralAngle > 116))) &&	(((lattice_tminus1[i+1].dihedralAngle < -117) && (lattice_tminus1[i+1].dihedralAngle > -180)) || ((lattice_tminus1[i+1].dihedralAngle < 180) && (lattice_tminus1[i+1].dihedralAngle > 116))))
		// {
		// 	possibleJumps.Gm7T++;
		// }
		// else if ((lattice[i-1].dihedralAngle < 115) && (lattice[i-1].dihedralAngle > 34) &&	(lattice[i].dihedralAngle < -38) && (lattice[i].dihedralAngle > -116) && (((lattice[i+1].dihedralAngle < -117) && (lattice[i+1].dihedralAngle > -180)) || ((lattice[i+1].dihedralAngle < 180) && (lattice[i+1].dihedralAngle > 116))) && (lattice_tminus1[i-1].dihedralAngle < 115) && (lattice_tminus1[i-1].dihedralAngle > 34) && (lattice_tminus1[i].dihedralAngle < 115) && (lattice_tminus1[i].dihedralAngle > 34) && (((lattice_tminus1[i+1].dihedralAngle < -117) && (lattice_tminus1[i+1].dihedralAngle > -180)) || ((lattice_tminus1[i+1].dihedralAngle < 180) && (lattice_tminus1[i+1].dihedralAngle > 116))))
		// {
		// 	possibleJumps.Gm7Gp++;
		// }
		// else if ((lattice[i-1].dihedralAngle < 115) && (lattice[i-1].dihedralAngle > 34) &&	(lattice[i].dihedralAngle < -38) && (lattice[i].dihedralAngle > -116) && (((lattice[i+1].dihedralAngle < -117) && (lattice[i+1].dihedralAngle > -180)) || ((lattice[i+1].dihedralAngle < 180) && (lattice[i+1].dihedralAngle > 116))) && (lattice_tminus1[i-1].dihedralAngle < 115) && (lattice_tminus1[i-1].dihedralAngle > 34) && (lattice_tminus1[i].dihedralAngle < -38) && (lattice_tminus1[i].dihedralAngle > -116) && (((lattice_tminus1[i+1].dihedralAngle < -117) && (lattice_tminus1[i+1].dihedralAngle > -180)) || ((lattice_tminus1[i+1].dihedralAngle < 180) && (lattice_tminus1[i+1].dihedralAngle > 116))))
		// {
		// 	possibleJumps.Gm7Gm;
		// }
		// else if ((lattice[i-1].dihedralAngle < 115) && (lattice[i-1].dihedralAngle > 34) && (lattice[i].dihedralAngle < -38) && (lattice[i].dihedralAngle > -116) && (lattice[i+1].dihedralAngle < 115) && (lattice[i+1].dihedralAngle > 34) && (lattice_tminus1[i-1].dihedralAngle < 115) && (lattice_tminus1[i-1].dihedralAngle > 34) && (((lattice_tminus1[i].dihedralAngle < -117) && (lattice_tminus1[i].dihedralAngle > -180)) || ((lattice_tminus1[i].dihedralAngle < 180) && (lattice_tminus1[i].dihedralAngle > 116))) && (lattice_tminus1[i+1].dihedralAngle < 115) && (lattice_tminus1[i+1].dihedralAngle > 34))
		// {
		// 	possibleJumps.Gm8T++;
		// }
		// else if ((lattice[i-1].dihedralAngle < 115) && (lattice[i-1].dihedralAngle > 34) &&	(lattice[i].dihedralAngle < -38) && (lattice[i].dihedralAngle > -116) && (lattice[i+1].dihedralAngle < 115) && (lattice[i+1].dihedralAngle > 34) &&	(lattice_tminus1[i-1].dihedralAngle < 115) && (lattice_tminus1[i-1].dihedralAngle > 34) && (lattice_tminus1[i].dihedralAngle < 115) && (lattice_tminus1[i].dihedralAngle > 34) && (lattice_tminus1[i+1].dihedralAngle < 115) && (lattice_tminus1[i+1].dihedralAngle > 34))
		// {
		// 	possibleJumps.Gm8Gp++;
		// }
		// else if ((lattice[i-1].dihedralAngle < 115) && (lattice[i-1].dihedralAngle > 34) && (lattice[i].dihedralAngle < -38) && (lattice[i].dihedralAngle > -116) && (lattice[i+1].dihedralAngle < 115) && (lattice[i+1].dihedralAngle > 34) &&	(lattice_tminus1[i-1].dihedralAngle < 115) && (lattice_tminus1[i-1].dihedralAngle > 34) && (lattice_tminus1[i].dihedralAngle < -38) && (lattice_tminus1[i].dihedralAngle > -116) && (lattice_tminus1[i+1].dihedralAngle < 115) && (lattice_tminus1[i+1].dihedralAngle > 34))
		// {
		// 	possibleJumps.Gm8Gm++;
		// }
		// else if ((lattice[i-1].dihedralAngle < 115) && (lattice[i-1].dihedralAngle > 34) &&	(lattice[i].dihedralAngle < -38) && (lattice[i].dihedralAngle > -116) && (lattice[i+1].dihedralAngle < -38) && (lattice[i+1].dihedralAngle > -116) && (lattice_tminus1[i-1].dihedralAngle < 115) && (lattice_tminus1[i-1].dihedralAngle > 34) && (((lattice_tminus1[i].dihedralAngle < -117) && (lattice_tminus1[i].dihedralAngle > -180)) || ((lattice_tminus1[i].dihedralAngle < 180) && (lattice_tminus1[i].dihedralAngle > 116))) && (lattice_tminus1[i+1].dihedralAngle < -38) && (lattice_tminus1[i+1].dihedralAngle > -116))
		// {
		// 	possibleJumps.Gm9T++;
		// }
		// else if ((lattice[i-1].dihedralAngle < 115) && (lattice[i-1].dihedralAngle > 34) &&	(lattice[i].dihedralAngle < -38) && (lattice[i].dihedralAngle > -116) && (lattice[i+1].dihedralAngle < -38) && (lattice[i+1].dihedralAngle > -116) && (lattice_tminus1[i-1].dihedralAngle < 115) && (lattice_tminus1[i-1].dihedralAngle > 34) && (lattice_tminus1[i].dihedralAngle < 115) && (lattice_tminus1[i].dihedralAngle > 34) && (lattice_tminus1[i+1].dihedralAngle < -38) && (lattice_tminus1[i+1].dihedralAngle > -116))
		// {
		// 	possibleJumps.Gm9Gp++;
		// }
		// else if ((lattice[i-1].dihedralAngle < 115) && (lattice[i-1].dihedralAngle > 34) &&	(lattice[i].dihedralAngle < -38) && (lattice[i].dihedralAngle > -116) && (lattice[i+1].dihedralAngle < -38) && (lattice[i+1].dihedralAngle > -116) && (lattice_tminus1[i-1].dihedralAngle < 115) && (lattice_tminus1[i-1].dihedralAngle > 34) && (lattice_tminus1[i].dihedralAngle < -38) && (lattice_tminus1[i].dihedralAngle > -116) && (lattice_tminus1[i+1].dihedralAngle < -38) && (lattice_tminus1[i+1].dihedralAngle > -116))
		// {
		// 	possibleJumps.Gm9Gm++;
		// }
	}

	possibleJumps.transition_CR = possibleJumps.Gp2T + possibleJumps.Gp4T + possibleJumps.Gp5T + possibleJumps.Gp5T + possibleJumps.Gp6T + possibleJumps.Gp8T + possibleJumps.Gm2T + possibleJumps.Gm4T + possibleJumps.Gm5T + possibleJumps.Gm5T + possibleJumps.Gm6T + possibleJumps.Gm8T;

	possibleJumps.transition_CL = possibleJumps.Gm1T + possibleJumps.Gm1T + possibleJumps.Gm2T + possibleJumps.Gm3T + possibleJumps.Gm4T + possibleJumps.Gm7T + possibleJumps.Gp1T + possibleJumps.Gp1T + possibleJumps.Gp2T + possibleJumps.Gp3T + possibleJumps.Gp4T + possibleJumps.Gp7T;

	possibleJumps.transition_CC = possibleJumps.Gp1Gp + possibleJumps.Gp1Gp + possibleJumps.Gp1Gm + possibleJumps.Gp1Gm + possibleJumps.Gp2Gp + possibleJumps.Gp2Gp + possibleJumps.Gp2Gm + possibleJumps.Gp2Gm + possibleJumps.Gp3Gp + possibleJumps.Gp3Gm + possibleJumps.Gp4Gm + possibleJumps.Gp4Gm + possibleJumps.Gp4Gp + possibleJumps.Gp4Gp + possibleJumps.Gp5Gp + possibleJumps.Gp5Gp + possibleJumps.Gp5Gm + possibleJumps.Gp5Gm + possibleJumps.Gp6Gp + possibleJumps.Gp6Gm + possibleJumps.Gp7Gp + possibleJumps.Gp7Gm + possibleJumps.Gp8Gp + possibleJumps.Gp8Gm + possibleJumps.Gm1Gm + possibleJumps.Gm1Gm + possibleJumps.Gm1Gp + possibleJumps.Gm1Gp + possibleJumps.Gm2Gp + possibleJumps.Gm2Gp + possibleJumps.Gm2Gm + possibleJumps.Gm2Gm + possibleJumps.Gm3Gm + possibleJumps.Gm3Gp + possibleJumps.Gm4Gm + possibleJumps.Gm4Gm + possibleJumps.Gm4Gp + possibleJumps.Gm4Gp + possibleJumps.Gm5Gp + possibleJumps.Gm5Gp + possibleJumps.Gm5Gm + possibleJumps.Gm5Gm + possibleJumps.Gm6Gm + possibleJumps.Gm6Gp + possibleJumps.Gm7Gm + possibleJumps.Gm7Gp + possibleJumps.Gm8Gp + possibleJumps.Gm8Gm;

	possibleJumps.transition_CM = 0;

	possibleJumps.transition_RR = possibleJumps.t2T + possibleJumps.t4T + possibleJumps.t5T + possibleJumps.t5T + possibleJumps.t6T + possibleJumps.t8T + possibleJumps.Gp3Gp + possibleJumps.Gp6Gp + possibleJumps.Gp7Gp + possibleJumps.Gp8Gp + possibleJumps.Gp9Gp + possibleJumps.Gp9Gp;

	possibleJumps.transition_RL = possibleJumps.Gp3Gm + possibleJumps.Gp6Gm + possibleJumps.Gp7Gm + possibleJumps.Gp8Gm + possibleJumps.Gp9Gm + possibleJumps.Gp9Gm;

	possibleJumps.transition_RC = possibleJumps.t2Gp + possibleJumps.t2Gm + possibleJumps.t4Gm + possibleJumps.t4Gp + possibleJumps.t5Gp + possibleJumps.t5Gp + possibleJumps.t5Gm + possibleJumps.t5Gm + possibleJumps.t6Gp + possibleJumps.t6Gm + possibleJumps.t8Gp + possibleJumps.t8Gm;

	possibleJumps.transition_RM = possibleJumps.Gp3T + possibleJumps.Gp6T + possibleJumps.Gp7T + possibleJumps.Gp8T + possibleJumps.Gp9T + possibleJumps.Gp9T;

	possibleJumps.transition_LR = possibleJumps.Gm3Gp + possibleJumps.Gm6Gp + possibleJumps.Gm7Gp + possibleJumps.Gm8Gp + possibleJumps.Gm9Gp + possibleJumps.Gm9Gp;

	possibleJumps.transition_LL = possibleJumps.t1T + possibleJumps.t1T + possibleJumps.t2T + possibleJumps.t3T + possibleJumps.t4T + possibleJumps.t7T + possibleJumps.Gm3Gm + possibleJumps.Gm6Gm + possibleJumps.Gm7Gm + possibleJumps.Gm8Gm + possibleJumps.Gm9Gm + possibleJumps.Gm9Gm;

	possibleJumps.transition_LC = possibleJumps.t1Gm + possibleJumps.t1Gp + possibleJumps.t1Gm + possibleJumps.t1Gp + possibleJumps.t2Gp + possibleJumps.t2Gm + possibleJumps.t3Gp + possibleJumps.t3Gm + possibleJumps.t4Gp + possibleJumps.t4Gm + possibleJumps.t7Gp + possibleJumps.t7Gm;

	possibleJumps.transition_LM = possibleJumps.Gm3T + possibleJumps.Gm6T + possibleJumps.Gm7T + possibleJumps.Gm8T + possibleJumps.Gm9T + possibleJumps.Gm9T;

	possibleJumps.transition_MR = possibleJumps.t3Gp + possibleJumps.t6Gp + possibleJumps.t7Gp + possibleJumps.t8Gp + possibleJumps.t9Gp + possibleJumps.t9Gp;

	possibleJumps.transition_ML = possibleJumps.t3Gm + possibleJumps.t6Gm + possibleJumps.t7Gm + possibleJumps.t8Gm + possibleJumps.t9Gm + possibleJumps.t9Gm;

	possibleJumps.transition_MC = 0;

	possibleJumps.transition_MM = possibleJumps.t3T + possibleJumps.t6T + possibleJumps.t7T + possibleJumps.t8T + possibleJumps.t9T + possibleJumps.t9T;
	
	return possibleJumps;
}

KMC_JUMPS computeKCMJumps2 (KMC_LATTICE *lattice, KMC_LATTICE *lattice_tminus1, int nLattice)
{
	KMC_JUMPS possibleJumps;
	int sequenceIndex_t = 0, sequenceIndex_tminus1 = 0; // 1 to 27, for each dihedral trios
	// int configurationIndex = 0, configurationIndex_tminus1 = 0; // 1 to 4, for coil, R, L and mesophase configuration
	int configurationIndexL = 0, configurationIndexR = 0, configurationIndexM = 0, configurationIndexC = 0;
	int configurationIndexL_tminus1 = 0, configurationIndexR_tminus1 = 0, configurationIndexM_tminus1 = 0, configurationIndexC_tminus1 = 0;

	// Checking for lattice_tminus1
	for (int i = 1; i < nLattice - 1; ++i)
	{
		// Re-initialization
		sequenceIndex_t = 0; sequenceIndex_tminus1 = 0;
		configurationIndexL = 0; configurationIndexR = 0; configurationIndexM = 0; configurationIndexC = 0; configurationIndexL_tminus1 = 0; configurationIndexR_tminus1 = 0; configurationIndexM_tminus1 = 0; configurationIndexC_tminus1 = 0;

		if (
			(chiralChecking (lattice_tminus1[i-1].dihedralAngle) == -1) &&
			(chiralChecking (lattice_tminus1[i].dihedralAngle) == 0) &&
			(chiralChecking (lattice_tminus1[i+1].dihedralAngle) == -1)
			)
		{
			sequenceIndex_tminus1 = 1;
		}
		else if (
			(chiralChecking (lattice_tminus1[i-1].dihedralAngle) == -1) &&
			(chiralChecking (lattice_tminus1[i].dihedralAngle) == 0) &&
			(chiralChecking (lattice_tminus1[i+1].dihedralAngle) == 1)
			)
		{
			sequenceIndex_tminus1 = 2;
		}
		else if (
			(chiralChecking (lattice_tminus1[i-1].dihedralAngle) == -1) &&
			(chiralChecking (lattice_tminus1[i].dihedralAngle) == 0) &&
			(chiralChecking (lattice_tminus1[i+1].dihedralAngle) == 0)
			)
		{
			sequenceIndex_tminus1 = 3;
		}
		else if (
			(chiralChecking (lattice_tminus1[i-1].dihedralAngle) == 1) &&
			(chiralChecking (lattice_tminus1[i].dihedralAngle) == 0) &&
			(chiralChecking (lattice_tminus1[i+1].dihedralAngle) == -1)
			)
		{
			sequenceIndex_tminus1 = 4;
		}
		else if (
			(chiralChecking (lattice_tminus1[i-1].dihedralAngle) == 1) &&
			(chiralChecking (lattice_tminus1[i].dihedralAngle) == 0) &&
			(chiralChecking (lattice_tminus1[i+1].dihedralAngle) == 1)
			)
		{
			sequenceIndex_tminus1 = 5;
		}
		else if (
			(chiralChecking (lattice_tminus1[i-1].dihedralAngle) == 1) &&
			(chiralChecking (lattice_tminus1[i].dihedralAngle) == 0) &&
			(chiralChecking (lattice_tminus1[i+1].dihedralAngle) == 0)
			)
		{
			sequenceIndex_tminus1 = 6;
		}
		else if (
			(chiralChecking (lattice_tminus1[i-1].dihedralAngle) == 0) &&
			(chiralChecking (lattice_tminus1[i].dihedralAngle) == 0) &&
			(chiralChecking (lattice_tminus1[i+1].dihedralAngle) == -1)
			)
		{
			sequenceIndex_tminus1 = 7;
		}
		else if (
			(chiralChecking (lattice_tminus1[i-1].dihedralAngle) == 0) &&
			(chiralChecking (lattice_tminus1[i].dihedralAngle) == 0) &&
			(chiralChecking (lattice_tminus1[i+1].dihedralAngle) == 1)
			)
		{
			sequenceIndex_tminus1 = 8;
		}
		else if (
			(chiralChecking (lattice_tminus1[i-1].dihedralAngle) == 0) &&
			(chiralChecking (lattice_tminus1[i].dihedralAngle) == 0) &&
			(chiralChecking (lattice_tminus1[i+1].dihedralAngle) == 0)
			)
		{
			sequenceIndex_tminus1 = 9;
		}
		else if (
			(chiralChecking (lattice_tminus1[i-1].dihedralAngle) == -1) &&
			(chiralChecking (lattice_tminus1[i].dihedralAngle) == -1) &&
			(chiralChecking (lattice_tminus1[i+1].dihedralAngle) == -1)
			)
		{
			sequenceIndex_tminus1 = 10;
		}
		else if (
			(chiralChecking (lattice_tminus1[i-1].dihedralAngle) == -1) &&
			(chiralChecking (lattice_tminus1[i].dihedralAngle) == -1) &&
			(chiralChecking (lattice_tminus1[i+1].dihedralAngle) == 1)
			)
		{
			sequenceIndex_tminus1 = 11;
		}
		else if (
			(chiralChecking (lattice_tminus1[i-1].dihedralAngle) == -1) &&
			(chiralChecking (lattice_tminus1[i].dihedralAngle) == -1) &&
			(chiralChecking (lattice_tminus1[i+1].dihedralAngle) == 0)
			)
		{
			sequenceIndex_tminus1 = 12;
		}
		else if (
			(chiralChecking (lattice_tminus1[i-1].dihedralAngle) == 1) &&
			(chiralChecking (lattice_tminus1[i].dihedralAngle) == -1) &&
			(chiralChecking (lattice_tminus1[i+1].dihedralAngle) == -1)
			)
		{
			sequenceIndex_tminus1 = 13;
		}
		else if (
			(chiralChecking (lattice_tminus1[i-1].dihedralAngle) == 1) &&
			(chiralChecking (lattice_tminus1[i].dihedralAngle) == -1) &&
			(chiralChecking (lattice_tminus1[i+1].dihedralAngle) == 1)
			)
		{
			sequenceIndex_tminus1 = 14;
		}
		else if (
			(chiralChecking (lattice_tminus1[i-1].dihedralAngle) == 1) &&
			(chiralChecking (lattice_tminus1[i].dihedralAngle) == -1) &&
			(chiralChecking (lattice_tminus1[i+1].dihedralAngle) == 0)
			)
		{
			sequenceIndex_tminus1 = 15;
		}
		else if (
			(chiralChecking (lattice_tminus1[i-1].dihedralAngle) == 0) &&
			(chiralChecking (lattice_tminus1[i].dihedralAngle) == -1) &&
			(chiralChecking (lattice_tminus1[i+1].dihedralAngle) == -1)
			)
		{
			sequenceIndex_tminus1 = 16;
		}
		else if (
			(chiralChecking (lattice_tminus1[i-1].dihedralAngle) == 0) &&
			(chiralChecking (lattice_tminus1[i].dihedralAngle) == -1) &&
			(chiralChecking (lattice_tminus1[i+1].dihedralAngle) == 1)
			)
		{
			sequenceIndex_tminus1 = 17;
		}
		else if (
			(chiralChecking (lattice_tminus1[i-1].dihedralAngle) == 0) &&
			(chiralChecking (lattice_tminus1[i].dihedralAngle) == -1) &&
			(chiralChecking (lattice_tminus1[i+1].dihedralAngle) == 0)
			)
		{
			sequenceIndex_tminus1 = 18;
		}
		else if (
			(chiralChecking (lattice_tminus1[i-1].dihedralAngle) == -1) &&
			(chiralChecking (lattice_tminus1[i].dihedralAngle) == 1) &&
			(chiralChecking (lattice_tminus1[i+1].dihedralAngle) == -1)
			)
		{
			sequenceIndex_tminus1 = 19;
		}
		else if (
			(chiralChecking (lattice_tminus1[i-1].dihedralAngle) == -1) &&
			(chiralChecking (lattice_tminus1[i].dihedralAngle) == 1) &&
			(chiralChecking (lattice_tminus1[i+1].dihedralAngle) == 1)
			)
		{
			sequenceIndex_tminus1 = 20;
		}
		else if (
			(chiralChecking (lattice_tminus1[i-1].dihedralAngle) == -1) &&
			(chiralChecking (lattice_tminus1[i].dihedralAngle) == 1) &&
			(chiralChecking (lattice_tminus1[i+1].dihedralAngle) == 0)
			)
		{
			sequenceIndex_tminus1 = 21;
		}
		else if (
			(chiralChecking (lattice_tminus1[i-1].dihedralAngle) == 1) &&
			(chiralChecking (lattice_tminus1[i].dihedralAngle) == 1) &&
			(chiralChecking (lattice_tminus1[i+1].dihedralAngle) == -1)
			)
		{
			sequenceIndex_tminus1 = 22;
		}
		else if (
			(chiralChecking (lattice_tminus1[i-1].dihedralAngle) == 1) &&
			(chiralChecking (lattice_tminus1[i].dihedralAngle) == 1) &&
			(chiralChecking (lattice_tminus1[i+1].dihedralAngle) == 1)
			)
		{
			sequenceIndex_tminus1 = 23;
		}
		else if (
			(chiralChecking (lattice_tminus1[i-1].dihedralAngle) == 1) &&
			(chiralChecking (lattice_tminus1[i].dihedralAngle) == 1) &&
			(chiralChecking (lattice_tminus1[i+1].dihedralAngle) == 0)
			)
		{
			sequenceIndex_tminus1 = 24;
		}
		else if (
			(chiralChecking (lattice_tminus1[i-1].dihedralAngle) == 0) &&
			(chiralChecking (lattice_tminus1[i].dihedralAngle) == 1) &&
			(chiralChecking (lattice_tminus1[i+1].dihedralAngle) == -1)
			)
		{
			sequenceIndex_tminus1 = 25;
		}
		else if (
			(chiralChecking (lattice_tminus1[i-1].dihedralAngle) == 0) &&
			(chiralChecking (lattice_tminus1[i].dihedralAngle) == 1) &&
			(chiralChecking (lattice_tminus1[i+1].dihedralAngle) == 1)
			)
		{
			sequenceIndex_tminus1 = 26;
		}
		else if (
			(chiralChecking (lattice_tminus1[i-1].dihedralAngle) == 0) &&
			(chiralChecking (lattice_tminus1[i].dihedralAngle) == 1) &&
			(chiralChecking (lattice_tminus1[i+1].dihedralAngle) == 0)
			)
		{
			sequenceIndex_tminus1 = 27;
		}

		// Checking for lattice
		if (
			(chiralChecking (lattice[i-1].dihedralAngle) == -1) &&
			(chiralChecking (lattice[i].dihedralAngle) == 0) &&
			(chiralChecking (lattice[i+1].dihedralAngle) == -1)
			)
		{
			sequenceIndex_t = 1;
		}
		else if (
			(chiralChecking (lattice[i-1].dihedralAngle) == -1) &&
			(chiralChecking (lattice[i].dihedralAngle) == 0) &&
			(chiralChecking (lattice[i+1].dihedralAngle) == 1)
			)
		{
			sequenceIndex_t = 2;
		}
		else if (
			(chiralChecking (lattice[i-1].dihedralAngle) == -1) &&
			(chiralChecking (lattice[i].dihedralAngle) == 0) &&
			(chiralChecking (lattice[i+1].dihedralAngle) == 0)
			)
		{
			sequenceIndex_t = 3;
		}
		else if (
			(chiralChecking (lattice[i-1].dihedralAngle) == 1) &&
			(chiralChecking (lattice[i].dihedralAngle) == 0) &&
			(chiralChecking (lattice[i+1].dihedralAngle) == -1)
			)
		{
			sequenceIndex_t = 4;
		}
		else if (
			(chiralChecking (lattice[i-1].dihedralAngle) == 1) &&
			(chiralChecking (lattice[i].dihedralAngle) == 0) &&
			(chiralChecking (lattice[i+1].dihedralAngle) == 1)
			)
		{
			sequenceIndex_t = 5;
		}
		else if (
			(chiralChecking (lattice[i-1].dihedralAngle) == 1) &&
			(chiralChecking (lattice[i].dihedralAngle) == 0) &&
			(chiralChecking (lattice[i+1].dihedralAngle) == 0)
			)
		{
			sequenceIndex_t = 6;
		}
		else if (
			(chiralChecking (lattice[i-1].dihedralAngle) == 0) &&
			(chiralChecking (lattice[i].dihedralAngle) == 0) &&
			(chiralChecking (lattice[i+1].dihedralAngle) == -1)
			)
		{
			sequenceIndex_t = 7;
		}
		else if (
			(chiralChecking (lattice[i-1].dihedralAngle) == 0) &&
			(chiralChecking (lattice[i].dihedralAngle) == 0) &&
			(chiralChecking (lattice[i+1].dihedralAngle) == 1)
			)
		{
			sequenceIndex_t = 8;
		}
		else if (
			(chiralChecking (lattice[i-1].dihedralAngle) == 0) &&
			(chiralChecking (lattice[i].dihedralAngle) == 0) &&
			(chiralChecking (lattice[i+1].dihedralAngle) == 0)
			)
		{
			sequenceIndex_t = 9;
		}
		else if (
			(chiralChecking (lattice[i-1].dihedralAngle) == -1) &&
			(chiralChecking (lattice[i].dihedralAngle) == -1) &&
			(chiralChecking (lattice[i+1].dihedralAngle) == -1)
			)
		{
			sequenceIndex_t = 10;
		}
		else if (
			(chiralChecking (lattice[i-1].dihedralAngle) == -1) &&
			(chiralChecking (lattice[i].dihedralAngle) == -1) &&
			(chiralChecking (lattice[i+1].dihedralAngle) == 1)
			)
		{
			sequenceIndex_t = 11;
		}
		else if (
			(chiralChecking (lattice[i-1].dihedralAngle) == -1) &&
			(chiralChecking (lattice[i].dihedralAngle) == -1) &&
			(chiralChecking (lattice[i+1].dihedralAngle) == 0)
			)
		{
			sequenceIndex_t = 12;
		}
		else if (
			(chiralChecking (lattice[i-1].dihedralAngle) == 1) &&
			(chiralChecking (lattice[i].dihedralAngle) == -1) &&
			(chiralChecking (lattice[i+1].dihedralAngle) == -1)
			)
		{
			sequenceIndex_t = 13;
		}
		else if (
			(chiralChecking (lattice[i-1].dihedralAngle) == 1) &&
			(chiralChecking (lattice[i].dihedralAngle) == -1) &&
			(chiralChecking (lattice[i+1].dihedralAngle) == 1)
			)
		{
			sequenceIndex_t = 14;
		}
		else if (
			(chiralChecking (lattice[i-1].dihedralAngle) == 1) &&
			(chiralChecking (lattice[i].dihedralAngle) == -1) &&
			(chiralChecking (lattice[i+1].dihedralAngle) == 0)
			)
		{
			sequenceIndex_t = 15;
		}
		else if (
			(chiralChecking (lattice[i-1].dihedralAngle) == 0) &&
			(chiralChecking (lattice[i].dihedralAngle) == -1) &&
			(chiralChecking (lattice[i+1].dihedralAngle) == -1)
			)
		{
			sequenceIndex_t = 16;
		}
		else if (
			(chiralChecking (lattice[i-1].dihedralAngle) == 0) &&
			(chiralChecking (lattice[i].dihedralAngle) == -1) &&
			(chiralChecking (lattice[i+1].dihedralAngle) == 1)
			)
		{
			sequenceIndex_t = 17;
		}
		else if (
			(chiralChecking (lattice[i-1].dihedralAngle) == 0) &&
			(chiralChecking (lattice[i].dihedralAngle) == -1) &&
			(chiralChecking (lattice[i+1].dihedralAngle) == 0)
			)
		{
			sequenceIndex_t = 18;
		}
		else if (
			(chiralChecking (lattice[i-1].dihedralAngle) == -1) &&
			(chiralChecking (lattice[i].dihedralAngle) == 1) &&
			(chiralChecking (lattice[i+1].dihedralAngle) == -1)
			)
		{
			sequenceIndex_t = 19;
		}
		else if (
			(chiralChecking (lattice[i-1].dihedralAngle) == -1) &&
			(chiralChecking (lattice[i].dihedralAngle) == 1) &&
			(chiralChecking (lattice[i+1].dihedralAngle) == 1)
			)
		{
			sequenceIndex_t = 20;
		}
		else if (
			(chiralChecking (lattice[i-1].dihedralAngle) == -1) &&
			(chiralChecking (lattice[i].dihedralAngle) == 1) &&
			(chiralChecking (lattice[i+1].dihedralAngle) == 0)
			)
		{
			sequenceIndex_t = 21;
		}
		else if (
			(chiralChecking (lattice[i-1].dihedralAngle) == 1) &&
			(chiralChecking (lattice[i].dihedralAngle) == 1) &&
			(chiralChecking (lattice[i+1].dihedralAngle) == -1)
			)
		{
			sequenceIndex_t = 22;
		}
		else if (
			(chiralChecking (lattice[i-1].dihedralAngle) == 1) &&
			(chiralChecking (lattice[i].dihedralAngle) == 1) &&
			(chiralChecking (lattice[i+1].dihedralAngle) == 1)
			)
		{
			sequenceIndex_t = 23;
		}
		else if (
			(chiralChecking (lattice[i-1].dihedralAngle) == 1) &&
			(chiralChecking (lattice[i].dihedralAngle) == 1) &&
			(chiralChecking (lattice[i+1].dihedralAngle) == 0)
			)
		{
			sequenceIndex_t = 24;
		}
		else if (
			(chiralChecking (lattice[i-1].dihedralAngle) == 0) &&
			(chiralChecking (lattice[i].dihedralAngle) == 1) &&
			(chiralChecking (lattice[i+1].dihedralAngle) == -1)
			)
		{
			sequenceIndex_t = 25;
		}
		else if (
			(chiralChecking (lattice[i-1].dihedralAngle) == 0) &&
			(chiralChecking (lattice[i].dihedralAngle) == 1) &&
			(chiralChecking (lattice[i+1].dihedralAngle) == 1)
			)
		{
			sequenceIndex_t = 26;
		}
		else if (
			(chiralChecking (lattice[i-1].dihedralAngle) == 0) &&
			(chiralChecking (lattice[i].dihedralAngle) == 1) &&
			(chiralChecking (lattice[i+1].dihedralAngle) == 0)
			)
		{
			sequenceIndex_t = 27;
		}

		// Finding configuration index for tminus1
		if (sequenceIndex_tminus1 == 11 ||
			sequenceIndex_tminus1 == 12 ||
			sequenceIndex_tminus1 == 13 ||
			sequenceIndex_tminus1 == 14 ||
			sequenceIndex_tminus1 == 15 ||
			sequenceIndex_tminus1 == 16 ||
			sequenceIndex_tminus1 == 17 ||
			sequenceIndex_tminus1 == 20 ||
			sequenceIndex_tminus1 == 21 ||
			sequenceIndex_tminus1 == 22 ||
			sequenceIndex_tminus1 == 23 ||
			sequenceIndex_tminus1 == 24 ||
			sequenceIndex_tminus1 == 25 ||
			sequenceIndex_tminus1 == 26
			)
		{
			configurationIndexC_tminus1 = 1;
			// configurationIndex_tminus1 = 1;
		}
		if (sequenceIndex_tminus1 == 2 ||
			sequenceIndex_tminus1 == 4 ||
			sequenceIndex_tminus1 == 6 ||
			sequenceIndex_tminus1 == 8 ||
			sequenceIndex_tminus1 == 21 ||
			sequenceIndex_tminus1 == 24 ||
			sequenceIndex_tminus1 == 25 ||
			sequenceIndex_tminus1 == 26
			)
		{
			configurationIndexR_tminus1 = 1;
			// configurationIndex_tminus1 = 2;
		}
		if (sequenceIndex_tminus1 == 2 ||
			sequenceIndex_tminus1 == 3 ||
			sequenceIndex_tminus1 == 4 ||
			sequenceIndex_tminus1 == 7 ||
			sequenceIndex_tminus1 == 12 ||
			sequenceIndex_tminus1 == 15 ||
			sequenceIndex_tminus1 == 16 ||
			sequenceIndex_tminus1 == 17
			)
		{
			configurationIndexL_tminus1 = 1;
			// configurationIndex_tminus1 = 3;
		}
		if (sequenceIndex_tminus1 == 3 ||
			sequenceIndex_tminus1 == 6 ||
			sequenceIndex_tminus1 == 7 ||
			sequenceIndex_tminus1 == 8
			)
		{
			configurationIndexM_tminus1 = 1;
			// configurationIndex_tminus1 = 4;
		}

		// Configuration index for t
		if (sequenceIndex_t == 11 ||
			sequenceIndex_t == 12 ||
			sequenceIndex_t == 13 ||
			sequenceIndex_t == 14 ||
			sequenceIndex_t == 15 ||
			sequenceIndex_t == 16 ||
			sequenceIndex_t == 17 ||
			sequenceIndex_t == 20 ||
			sequenceIndex_t == 21 ||
			sequenceIndex_t == 22 ||
			sequenceIndex_t == 23 ||
			sequenceIndex_t == 24 ||
			sequenceIndex_t == 25 ||
			sequenceIndex_t == 26
			)
		{
			configurationIndexC = 1;
			// configurationIndex = 1;
		}
		if (sequenceIndex_t == 2 ||
			sequenceIndex_t == 4 ||
			sequenceIndex_t == 6 ||
			sequenceIndex_t == 8 ||
			sequenceIndex_t == 21 ||
			sequenceIndex_t == 24 ||
			sequenceIndex_t == 25 ||
			sequenceIndex_t == 26
			)
		{
			configurationIndexR = 1;
			// configurationIndex = 2;
		}
		if (sequenceIndex_t == 2 ||
			sequenceIndex_t == 3 ||
			sequenceIndex_t == 4 ||
			sequenceIndex_t == 7 ||
			sequenceIndex_t == 12 ||
			sequenceIndex_t == 15 ||
			sequenceIndex_t == 16 ||
			sequenceIndex_t == 17
			)
		{
			configurationIndexL = 1;
			// configurationIndex = 3;
		}
		if (sequenceIndex_t == 3 ||
			sequenceIndex_t == 6 ||
			sequenceIndex_t == 7 ||
			sequenceIndex_t == 8
			)
		{
			configurationIndexM = 1;
			// configurationIndex = 4;
		}

		// Counting the transitions
		if (configurationIndexC_tminus1 == 1 &&
			configurationIndexR == 1
			)
		{
			possibleJumps.transition_CR++;
		}
		if (configurationIndexC_tminus1 == 1 &&
			configurationIndexL == 1
			)
		{
			possibleJumps.transition_CL++;
		}
		if (configurationIndexC_tminus1 == 1 &&
			configurationIndexC == 1
			)
		{
			possibleJumps.transition_CC++;
		}
		if (configurationIndexC_tminus1 == 1 &&
			configurationIndexM == 1
			)
		{
			possibleJumps.transition_CM++;
		}
		if (configurationIndexR_tminus1 == 1 &&
			configurationIndexR == 1
			)
		{
			possibleJumps.transition_RR++;
		}
		if (configurationIndexR_tminus1 == 1 &&
			configurationIndexL == 1
			)
		{
			possibleJumps.transition_RL++;
		}
		if (configurationIndexR_tminus1 == 1 &&
			configurationIndexC == 1
			)
		{
			possibleJumps.transition_RC++;
		}
		if (configurationIndexR_tminus1 == 1 &&
			configurationIndexM == 1
			)
		{
			possibleJumps.transition_RM++;
		}
		if (configurationIndexL_tminus1 == 1 &&
			configurationIndexR == 1
			)
		{
			possibleJumps.transition_LR++;
		}
		if (configurationIndexL_tminus1 == 1 &&
			configurationIndexL == 1
			)
		{
			possibleJumps.transition_LL++;
		}
		if (configurationIndexL_tminus1 == 1 &&
			configurationIndexC == 1
			)
		{
			possibleJumps.transition_LC++;
		}
		if (configurationIndexL_tminus1 == 1 &&
			configurationIndexM == 1
			)
		{
			possibleJumps.transition_LM++;
		}
		if (configurationIndexM_tminus1 == 1 &&
			configurationIndexR == 1
			)
		{
			possibleJumps.transition_MR++;
		}
		if (configurationIndexM_tminus1 == 1 &&
			configurationIndexL == 1
			)
		{
			possibleJumps.transition_ML++;
		}
		if (configurationIndexM_tminus1 == 4 &&
			configurationIndexC == 1
			)
		{
			possibleJumps.transition_MC++;
		}
		if (configurationIndexM_tminus1 == 1 &&
			configurationIndexM == 1
			)
		{
			possibleJumps.transition_MM++;
		}

		// Counting for double config sequences, such as G- T G-, which contains 2xL helix
		// Checking for double 'L'
		if (
			(
				sequenceIndex_tminus1 == 1 || sequenceIndex_tminus1 == 18
				) &&
			(
				sequenceIndex_t == 2 || sequenceIndex_t == 3 || sequenceIndex_t == 4 || sequenceIndex_t == 7 || sequenceIndex_t == 12 || sequenceIndex_t == 15 || sequenceIndex_t == 16 || sequenceIndex_t == 17
				)
			)
		{
			possibleJumps.transition_LL++;
		}
		if (
			(
				sequenceIndex_tminus1 == 1 || sequenceIndex_tminus1 == 18
				) &&
			(
				sequenceIndex_t == 1 || sequenceIndex_t == 18
				)
			)
		{
			possibleJumps.transition_LL += 2;
		}
		if (
			(
				sequenceIndex_tminus1 == 1 || sequenceIndex_tminus1 == 18
				) &&
			(
				sequenceIndex_t == 2 || sequenceIndex_t == 4 || sequenceIndex_t == 6 || sequenceIndex_t == 8 || sequenceIndex_t == 21  || sequenceIndex_t == 24 || sequenceIndex_t == 25 || sequenceIndex_t == 26
				)
			)
		{
			possibleJumps.transition_LR++;
		}
		if (
			(
				sequenceIndex_tminus1 == 1 || sequenceIndex_tminus1 == 18
				) &&
			(
				sequenceIndex_t == 27 || sequenceIndex_t == 5
				)
			)
		{
			possibleJumps.transition_LR += 2;
		}
		if (
			(
				sequenceIndex_tminus1 == 1 || sequenceIndex_tminus1 == 18
				) &&
			(
				sequenceIndex_t == 12 || sequenceIndex_t == 15 || sequenceIndex_t == 16 || sequenceIndex_t == 17 || sequenceIndex_t == 21 || sequenceIndex_t == 24 || sequenceIndex_t == 25 || sequenceIndex_t == 26
				)
			)
		{
			possibleJumps.transition_LC++;
		}
		if (
			(
				sequenceIndex_tminus1 == 1 || sequenceIndex_tminus1 == 18
				) &&
			(
				sequenceIndex_t == 10 || sequenceIndex_t == 11 || sequenceIndex_t == 13 || sequenceIndex_t == 14 || sequenceIndex_t == 19 || sequenceIndex_t == 20 || sequenceIndex_t == 22 || sequenceIndex_t == 23
				)
			)
		{
			possibleJumps.transition_LC += 2;
		}
		if (
			(
				sequenceIndex_tminus1 == 1 || sequenceIndex_tminus1 == 18
				) &&
			(
				sequenceIndex_t == 3 || sequenceIndex_t == 6 || sequenceIndex_t == 7 || sequenceIndex_t == 8
				)
			)
		{
			possibleJumps.transition_LM++;
		}
		if (
			(
				sequenceIndex_tminus1 == 1 || sequenceIndex_tminus1 == 18
				) &&
			(
				sequenceIndex_t == 9
				)
			)
		{
			possibleJumps.transition_LM += 2;
		}

		// Checking for double 'R'
		if (
			(
				sequenceIndex_tminus1 == 5 || sequenceIndex_tminus1 == 27
				) &&
			(
				sequenceIndex_t == 2 || sequenceIndex_t == 3 || sequenceIndex_t == 4 || sequenceIndex_t == 7 || sequenceIndex_t == 12 || sequenceIndex_t == 15 || sequenceIndex_t == 16 || sequenceIndex_t == 17
				)
			)
		{
			possibleJumps.transition_RL++;
		}
		if (
			(
				sequenceIndex_tminus1 == 5 || sequenceIndex_tminus1 == 27
				) &&
			(
				sequenceIndex_t == 1 || sequenceIndex_t == 18
				)
			)
		{
			possibleJumps.transition_RL += 2;
		}
		if (
			(
				sequenceIndex_tminus1 == 5 || sequenceIndex_tminus1 == 27
				) &&
			(
				sequenceIndex_t == 2 || sequenceIndex_t == 4 || sequenceIndex_t == 6 || sequenceIndex_t == 8 || sequenceIndex_t == 21 || sequenceIndex_t == 24 || sequenceIndex_t == 25 || sequenceIndex_t == 26
				)
			)
		{
			possibleJumps.transition_RR++;
		}
		if (
			(
				sequenceIndex_tminus1 == 5 || sequenceIndex_tminus1 == 27
				) &&
			(
				sequenceIndex_t == 5 || sequenceIndex_t == 27
				)
			)
		{
			possibleJumps.transition_RR += 2;
		}
		if (
			(
				sequenceIndex_tminus1 == 5 || sequenceIndex_tminus1 == 27
				) &&
			(
				sequenceIndex_t == 12 || sequenceIndex_t == 15 || sequenceIndex_t == 16 || sequenceIndex_t == 17 || sequenceIndex_t == 21 || sequenceIndex_t == 24 || sequenceIndex_t == 25 || sequenceIndex_t == 26
				)
			)
		{
			possibleJumps.transition_RC++;
		}
		if (
			(
				sequenceIndex_tminus1 == 5 || sequenceIndex_tminus1 == 27
				) &&
			(
				sequenceIndex_t == 10 || sequenceIndex_t == 11 || sequenceIndex_t == 13 || sequenceIndex_t == 14 || sequenceIndex_t == 19 || sequenceIndex_t == 20 || sequenceIndex_t == 22 || sequenceIndex_t == 23
				)
			)
		{
			possibleJumps.transition_RC += 2;
		}
		if (
			(
				sequenceIndex_tminus1 == 5 || sequenceIndex_tminus1 == 27
				) &&
			(
				sequenceIndex_t == 3 || sequenceIndex_t == 6 || sequenceIndex_t == 7 || sequenceIndex_t == 8
				)
			)
		{
			possibleJumps.transition_RM++;
		}
		if (
			(
				sequenceIndex_tminus1 == 5 || sequenceIndex_tminus1 == 27
				) &&
			(
				sequenceIndex_t == 9
				)
			)
		{
			possibleJumps.transition_RM += 2;
		}

		// Checking for double 'C' (coil)
		if (
			(
				sequenceIndex_tminus1 == 10 || sequenceIndex_tminus1 == 11 || sequenceIndex_tminus1 == 13 || sequenceIndex_tminus1 == 14 || sequenceIndex_tminus1 == 19 || sequenceIndex_tminus1 == 20 || sequenceIndex_tminus1 == 22 || sequenceIndex_tminus1 == 23
				) &&
			(
				sequenceIndex_t == 2 || sequenceIndex_t == 3 || sequenceIndex_t == 4 || sequenceIndex_t == 7 || sequenceIndex_t == 12 || sequenceIndex_t == 15 || sequenceIndex_t == 16 || sequenceIndex_t == 17
				)
			)
		{
			possibleJumps.transition_CL++;
		}
		if (
			(
				sequenceIndex_tminus1 == 10 || sequenceIndex_tminus1 == 11 || sequenceIndex_tminus1 == 13 || sequenceIndex_tminus1 == 14 || sequenceIndex_tminus1 == 19 || sequenceIndex_tminus1 == 20 || sequenceIndex_tminus1 == 22 || sequenceIndex_tminus1 == 23
				) &&
			(
				sequenceIndex_t == 1 || sequenceIndex_t == 18
				)
			)
		{
			possibleJumps.transition_CL += 2;
		}
		if (
			(
				sequenceIndex_tminus1 == 10 || sequenceIndex_tminus1 == 11 || sequenceIndex_tminus1 == 13 || sequenceIndex_tminus1 == 14 || sequenceIndex_tminus1 == 19 || sequenceIndex_tminus1 == 20 || sequenceIndex_tminus1 == 22 || sequenceIndex_tminus1 == 23
				) &&
			(
				sequenceIndex_t == 2 || sequenceIndex_t == 4 || sequenceIndex_t == 6 || sequenceIndex_t == 8 || sequenceIndex_t == 21 || sequenceIndex_t == 24 || sequenceIndex_t == 25 || sequenceIndex_t == 26
				)
			)
		{
			possibleJumps.transition_CR++;
		}
		if (
			(
				sequenceIndex_tminus1 == 10 || sequenceIndex_tminus1 == 11 || sequenceIndex_tminus1 == 13 || sequenceIndex_tminus1 == 14 || sequenceIndex_tminus1 == 19 || sequenceIndex_tminus1 == 20 || sequenceIndex_tminus1 == 22 || sequenceIndex_tminus1 == 23
				) &&
			(
				sequenceIndex_t == 5 || sequenceIndex_t == 27
				)
			)
		{
			possibleJumps.transition_CR += 2;
		}
		if (
			(
				sequenceIndex_tminus1 == 10 || sequenceIndex_tminus1 == 11 || sequenceIndex_tminus1 == 13 || sequenceIndex_tminus1 == 14 || sequenceIndex_tminus1 == 19 || sequenceIndex_tminus1 == 20 || sequenceIndex_tminus1 == 22 || sequenceIndex_tminus1 == 23
				) &&
			(
				sequenceIndex_t == 12 || sequenceIndex_t == 15 || sequenceIndex_t == 16 || sequenceIndex_t == 17 || sequenceIndex_t == 21 || sequenceIndex_t == 24 || sequenceIndex_t == 25 || sequenceIndex_t == 26
				)
			)
		{
			possibleJumps.transition_CC++;
		}
		if (
			(
				sequenceIndex_tminus1 == 10 || sequenceIndex_tminus1 == 11 || sequenceIndex_tminus1 == 13 || sequenceIndex_tminus1 == 14 || sequenceIndex_tminus1 == 19 || sequenceIndex_tminus1 == 20 || sequenceIndex_tminus1 == 22 || sequenceIndex_tminus1 == 23
				) &&
			(
				sequenceIndex_t == 10 || sequenceIndex_t == 11 || sequenceIndex_t == 13 || sequenceIndex_t == 14 || sequenceIndex_t == 19 || sequenceIndex_t == 20 || sequenceIndex_t == 22 || sequenceIndex_t == 23
				)
			)
		{
			possibleJumps.transition_CC += 2;
		}
		if (
			(
				sequenceIndex_tminus1 == 10 || sequenceIndex_tminus1 == 11 || sequenceIndex_tminus1 == 13 || sequenceIndex_tminus1 == 14 || sequenceIndex_tminus1 == 19 || sequenceIndex_tminus1 == 20 || sequenceIndex_tminus1 == 22 || sequenceIndex_tminus1 == 23
				) &&
			(
				sequenceIndex_t == 3 || sequenceIndex_t == 6 || sequenceIndex_t == 7 || sequenceIndex_t == 8
				)
			)
		{
			possibleJumps.transition_CM++;
		}
		if (
			(
				sequenceIndex_tminus1 == 10 || sequenceIndex_tminus1 == 11 || sequenceIndex_tminus1 == 13 || sequenceIndex_tminus1 == 14 || sequenceIndex_tminus1 == 19 || sequenceIndex_tminus1 == 20 || sequenceIndex_tminus1 == 22 || sequenceIndex_tminus1 == 23
				) &&
			(
				sequenceIndex_t == 9
				)
			)
		{
			possibleJumps.transition_CM += 2;
		}

		// Checking for double 'M' (mesophase)
		if (
			(
				sequenceIndex_tminus1 == 9
				) &&
			(
				sequenceIndex_t == 2 || sequenceIndex_t == 3 || sequenceIndex_t == 4 || sequenceIndex_t == 7 || sequenceIndex_t == 12 || sequenceIndex_t == 15 || sequenceIndex_t == 16 || sequenceIndex_t == 17
				)
			)
		{
			possibleJumps.transition_ML++;
		}
		if (
			(
				sequenceIndex_tminus1 == 9
				) &&
			(
				sequenceIndex_t == 1 || sequenceIndex_t == 18
				)
			)
		{
			possibleJumps.transition_ML += 2;
		}
		if (
			(
				sequenceIndex_tminus1 == 9
				) &&
			(
				sequenceIndex_t == 2 || sequenceIndex_t == 4 || sequenceIndex_t == 6 || sequenceIndex_t == 8 || sequenceIndex_t == 21 || sequenceIndex_t == 24 || sequenceIndex_t == 25 || sequenceIndex_t == 26
				)
			)
		{
			possibleJumps.transition_MR++;
		}
		if (
			(
				sequenceIndex_tminus1 == 9
				) &&
			(
				sequenceIndex_t == 5 || sequenceIndex_t == 27
				)
			)
		{
			possibleJumps.transition_MR += 2;
		}
		if (
			(
				sequenceIndex_tminus1 == 9
				) &&
			(
				sequenceIndex_t == 12 || sequenceIndex_t == 15 || sequenceIndex_t == 16 || sequenceIndex_t == 17 || sequenceIndex_t == 21 || sequenceIndex_t == 24 || sequenceIndex_t == 25 || sequenceIndex_t == 26
				)
			)
		{
			possibleJumps.transition_MC++;
		}
		if (
			(
				sequenceIndex_tminus1 == 9
				) &&
			(
				sequenceIndex_t == 10 || sequenceIndex_t == 11 || sequenceIndex_t == 13 || sequenceIndex_t == 14 || sequenceIndex_t == 19 || sequenceIndex_t == 20 || sequenceIndex_t == 22 || sequenceIndex_t == 23
				)
			)
		{
			possibleJumps.transition_MC += 2;
		}
		if (
			(
				sequenceIndex_tminus1 == 9
				) &&
			(
				sequenceIndex_t == 3 || sequenceIndex_t == 6 || sequenceIndex_t == 7 || sequenceIndex_t == 8
				)
			)
		{
			possibleJumps.transition_MM++;
		}
		if (
			(
				sequenceIndex_tminus1 == 9
				) &&
			(
				sequenceIndex_t == 9
				)
			)
		{
			possibleJumps.transition_MM += 2;
		}
	}
	
	return possibleJumps;
}
