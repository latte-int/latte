#include "../flags.h"

listCone* readListCone(rationalVector*, int);
vector* createConeDecMatrix(listCone*, int, int);
listCone* decomposeCones(listCone*, int, unsigned int Flags, char *File_Name);

// Added by Peter/David
// 
void decomposeCones_Single (listCone *, int, int degree, unsigned int flags, char *File_Name);




