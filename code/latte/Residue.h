// This is a -*- C++ -*- header file.

#ifndef RESIDUE_H
#define RESIDUE_H

ZZ Residue(listCone* cones, int numOfVars);

// Returns a -1 if dot product in denominator is 0, 1 otherwise
int Residue_Single_Cone (listCone *cones, int numOfVars,
			 ZZ *Random_Lambda, ZZ *Total_Lattice_Points,
			 ZZ *Ten_Power);

#endif
