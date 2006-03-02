#ifndef CONE_H
#define CONE_H

listCone* createListCone();
int lengthListCone(listCone*);

/* Free the whole list of cones. */
void freeListCone(listCone *list);

#endif
