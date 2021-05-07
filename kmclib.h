#ifndef POSTPROCESS_H
#define POSTPROCESS_H

typedef struct transitions
{
	int CR, CL, CC, CM, RR, RL, RC, RM, LR, LL, LC, LM, MR, ML, MC, MM;
} TRANSITIONS;

typedef struct kmc_jumps
{
	int t1Gp, t1Gm, t1T, t2Gp, t2Gm, t2T, t3Gp, t3Gm, t3T, t4Gp, t4Gm, t4T, t5Gp, t5Gm, t5T, t6Gp, t6Gm, t6T, t7Gp, t7Gm, t7T, t8Gp, t8Gm, t8T, t9Gp, t9Gm, t9T, Gp1Gp, Gp1Gm, Gp1T, Gp2Gp, Gp2Gm, Gp2T, Gp3Gp, Gp3Gm, Gp3T, Gp4Gp, Gp4Gm, Gp4T, Gp5Gp, Gp5Gm, Gp5T, Gp6Gp, Gp6Gm, Gp6T, Gp7Gp, Gp7Gm, Gp7T, Gp8Gp, Gp8Gm, Gp8T, Gp9Gp, Gp9Gm, Gp9T, Gm1Gp, Gm1Gm, Gm1T, Gm2Gp, Gm2Gm, Gm2T, Gm3Gp, Gm3Gm, Gm3T, Gm4Gp, Gm4Gm, Gm4T, Gm5Gp, Gm5Gm, Gm5T, Gm6Gp, Gm6Gm, Gm6T, Gm7Gp, Gm7Gm, Gm7T, Gm8Gp, Gm8Gm, Gm8T, Gm9Gp, Gm9Gm, Gm9T;
	int transition_CR, transition_CL, transition_CC, transition_CM;
	int transition_RR, transition_RL, transition_RC, transition_RM;
	int transition_LR, transition_LL, transition_LC, transition_LM;
	int transition_MR, transition_ML, transition_MC, transition_MM;
} KMC_JUMPS;

typedef struct latticeProbabilities
{
	int GmTGp, GmTT, GmTGm, GpTGm, GpTT, GpTGp, TTGp, TTGm, TTT;
	int TGpGm, TGpGp, TGpT, GmGpT, GmGpGp, GmGpGm, GpGpT, GpGpGp, GpGpGm;
	int TGmT, TGmGp, TGmGm, GmGmT, GmGmGp, GmGmGm, GpGmT, GpGmGp, GpGmGm;
} LATTICE_PROBABILITIES;

typedef struct dumpInfo
{
	int sino, atomType, molType, ix, iy, iz;
	float x, y, z;
} DUMPINFO;

typedef struct diheralInfo
{
	int sino, atom1, atom2, atom3, atom4;
	float dihedralAngle;
} DIHEDRALINFO;

typedef struct kmc_lattice
{
	int sino1, sino2, atomType1, atomType2, molType1, molType2, ix1, ix2, iy1, iy2, iz1, iz2;
	float x1, x2, y1, y2, z1, z2, dihedralAngle;
} KMC_LATTICE;

int countDihedrals (FILE *inputDihedral);
int countAtoms (FILE *inputDump);
DUMPINFO *readDump (FILE *inputDump, int nAtoms, int *isEOF_dump);
DIHEDRALINFO *readDihedral (FILE *inputDihedral, int nDihedrals, int *isEOF_dihedral, int sort);
int *checkPendant (DUMPINFO *dump, int nAtoms, int pendantMolID);
int countBackboneDihedrals (DUMPINFO *dump, int nAtoms, DIHEDRALINFO *dihedral, int nDihedrals, int pendantMolID);
KMC_LATTICE *computeKMCLattice (DUMPINFO *dump, int nAtoms, DIHEDRALINFO *dihedral, int nDihedrals, int nBackboneDihedrals, int pendantMolID, int *nLattice);
LATTICE_PROBABILITIES computeProbabilities (KMC_LATTICE *lattice, int nBackboneDihedrals);
KMC_JUMPS computeKCMJumps (KMC_LATTICE *lattice, KMC_LATTICE *lattice_tminus1, int nLattice);
KMC_JUMPS computeKCMJumps2 (KMC_LATTICE *lattice, KMC_LATTICE *lattice_tminus1, int nLattice);

#endif