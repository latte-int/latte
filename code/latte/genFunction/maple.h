rationalVector* addRationalVectorsWithUpperBoundOne(rationalVector*, 
						    rationalVector*, int);
rationalVector* subRationalVector(rationalVector*, rationalVector*, int);
listVector* readListVector(char*);
vector* subtractRowFromRow(vector*, int, int, int, vector*, int);
rationalVector* solveLinearSystem(vector*, vector, int, int);
vector movePoint(vector, rationalVector*, rationalVector*, vector*, int, int);
listVector* pointsInParallelepiped(rationalVector*, listVector*, int);
listVector* pointsInParallelepipedOfUnimodularCone(rationalVector*, 
						   listVector*, int);
void writeTermToFile(FILE*, vector, int);
void writeTermOfGeneratingFunctionToFile(FILE*, listCone*, int);
void createGeneratingFunctionAsMapleInput(char*, listCone*, int);
void createGeneratingFunctionAsMapleInputGrob(listCone* cones, 
					      int numOfVars, ofstream & out);
