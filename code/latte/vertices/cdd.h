void createCddIneFile(listVector*, int);
void createCddExtFile(listVector*, int);
void createCddIneLPFile(listVector* matrix, int numOfVars, vector & cost);
listVector* createListOfInequalities(listVector*, int);
listCone* readCddExtFile();
rationalVector* ReadLpsFile(int numOfVars);
listCone* readCddEadFile(listCone*, int);
listCone* readCddEadFileFromVrep(listCone* cones, int numOfVars);
listCone* computeVertexCones(char*, listVector*, int);
rationalVector* LP(listVector* matrix, vector & cost, int numOfVars);
listCone* CopyListCones(listCone* RudyCones, int numOfVars, 
			rationalVector* Opt_vertex);
listCone* CopyListCones(listCone* RudyCones, int numOfVars); 
listCone* computeVertexConesFromVrep(char* fileName, listVector* matrix, 
				     int numOfVars);
listCone* computeVertexConesViaLrs(char* fileName, listVector* matrix, 
				   int numOfVars);
void CreatExtEadFile();
