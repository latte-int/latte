void printVector(vec_ZZ, int);
void printListVector(listVector*, int);
listVector* Changeproj(listVector* basis, listVector* equ, int numOfVars);
void printRationalVector(rationalVector*, int);
void printCone(listCone*, int);
void printListCone(listCone*, int);
void printVectorToFile(ofstream &, vec_ZZ, int);
void printListVectorToFile(ofstream &, listVector*, int);
void printVectorToFileWithoutBrackets(ofstream &, vec_ZZ, int);
void printListVectorToFileWithoutBrackets(ofstream &, listVector*, int);
void printRationalVectorToFile(ofstream &, rationalVector*, int);
void printRationalVectorToFileWithoutBrackets(ofstream &, rationalVector*, 
					      int);
void printConeToFile(ofstream &,listCone*, int);
void printListConeToFile(char*, listCone*, int); 
void printResidueFile(char*, listCone*, int);
ZZ Residue(listCone* cones, int numOfVars);

// Returns a -1 if dot product in denominator is 0, 1 otherwise
int Residue_Single_Cone (listCone *cones, int numOfVars, ZZ *Random_Lambda, ZZ *Total_Lattice_Points, ZZ *Ten_Power);

