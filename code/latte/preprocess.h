#ifndef PREPROCESS__H
#define PREPROCESS__H
int ihermite(vec_ZZ *S, vec_ZZ *U, vec_ZZ* rhs, int m, int n);
listVector* preprocessProblem(listVector*, listVector*, vec_ZZ**, int*, vec_ZZ&, mat_ZZ &, char*, int);
 listVector* TransformToDualCone(listVector* matrix, int& numOfVars);
void dilateListVector(listVector* basis, int numOfVars, int dil);
vec_ZZ transpose(vec_ZZ mat, int numOfVars, int numOfRows);
vec_ZZ ProjectingUp(mat_ZZ ProjU, vec_ZZ cost, int numOfVars);
vec_RR ProjectingUpRR(mat_RR ProjU, vec_RR cost, int numOfVars);
void Interior(listVector* basis);
#endif
