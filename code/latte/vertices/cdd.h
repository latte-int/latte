// This is a -*- C++ -*- header file.

#ifndef VERTICES_CDD__H
#define VERTICES_CDD__H

void createCddIneFile(listVector*, int);
void createCddExtFile(listVector*, int);
void createCddIneLPFile(listVector* matrix, int numOfVars, vec_ZZ & cost);
listVector* createListOfInequalities(listVector*, int);
listCone* readCddExtFile();
rationalVector* ReadLpsFile(int numOfVars, bool verbose = true);
listCone* readCddEadFile(listCone*, int);
listCone* readCddEadFileFromVrep(listCone* cones, int numOfVars);
listCone* computeVertexCones(char*, listVector*, int);
rationalVector* LP(listVector* matrix, vec_ZZ& cost, int numOfVars,
		   bool verbose = true);
listCone* CopyListCones(listCone* RudyCones, int numOfVars, 
			rationalVector* Opt_vertex);
listCone* CopyListCones(listCone* RudyCones, int numOfVars); 
listCone* computeVertexConesFromVrep(char* fileName, int numOfVars);
listCone* computeVertexConesViaLrs(char* fileName, listVector* matrix, 
				   int numOfVars);
void CreatExtEadFile();

#endif
