#ifndef PREPROCESS__H
#define PREPROCESS__H
int ihermite(vector *S, vector *U, vector* rhs, int m, int n);
listVector* preprocessProblem(listVector*, listVector*, vector**, int*, vector&, mat_ZZ &, char*, int);
 listVector* TransformToDualCone(listVector* matrix, int& numOfVars);
void dilateListVector(listVector* basis, int numOfVars, int dil);
vector transpose(vector mat, int numOfVars, int numOfRows);
vector ProjectingUp(mat_ZZ ProjU, vector cost, int numOfVars);
vec_RR ProjectingUpRR(mat_RR ProjU, vec_RR cost, int numOfVars);
void Interior(listVector* basis);
#endif
