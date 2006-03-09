#include <stdlib.h>
#include "myheader.h"
#include "ramon.h"
/* ----------------------------------------------------------------- */
listCone* createListCone() {
  listCone* z;

  z = new listCone;
//    z = (listCone*)malloc(sizeof(listCone));
  if (z==0) exit(0);

  z->coefficient=1;
  z->vertex=0;
  z->rays=0;
  z->facets=0;
  z->determinant = 0;
  z->latticePoints=0;
  z->rest=0;

  return (z);
}
/* ----------------------------------------------------------------- */
int lengthListCone(listCone* LIST) {
  int len=0;

  while (LIST) {len++; LIST = LIST->rest;}
  return (len);
}
/* ----------------------------------------------------------------- */
void freeCone(listCone *cone)
{
  delete cone->vertex;
  freeListVector(cone->rays);
  freeListVector(cone->facets);
  freeListVector(cone->latticePoints);
  delete cone;
}

void freeListCone(listCone *list)
{
  while (list != NULL) {
    listCone *rest = list->rest;
    freeCone(list);
    list = rest;
  }
}
