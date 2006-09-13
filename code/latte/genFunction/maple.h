rationalVector* addRationalVectorsWithUpperBoundOne(rationalVector*, 
						    rationalVector*, int);
rationalVector* subRationalVector(rationalVector*, rationalVector*, int);
listVector* readListVector(char*);
vec_ZZ movePoint(vec_ZZ, rationalVector*, rationalVector*, vec_ZZ*, int, int);
listVector* pointsInParallelepiped(rationalVector*, listVector*, int);
void writeTermToFile(FILE*, vec_ZZ, int);
void writeTermOfGeneratingFunctionToFile(FILE*, listCone*, int);
void createGeneratingFunctionAsMapleInput(char*, listCone*, int);
void createGeneratingFunctionAsMapleInputGrob(listCone* cones, 
					      int numOfVars, ofstream & out);
