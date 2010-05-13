#define COEFF_MAX 10000

#include "PolyRep.h"
#include "newIntegration.h"
#include "../timing.h"
#include <iostream>
#include <fstream>
#include <NTL/ZZ.h>

using namespace std;

int main(int argc, char *argv[])
{
	if (argc < 3) { cout << "Usage: " << argv[0] << " fileIn fileOut IntegrateMode (polynomial/linear)" << endl; return 1; }
	_linFormSum _forms;
	_monomialSum _monomials;
        simplexZZ mySimplex;
	ifstream myStream (argv[1], ios::in);
	ofstream outStream(argv[2], ios::out);
	if (!myStream.is_open() || !outStream.is_open()) { cout << "Error opening file " << argv[1] << ", please make sure it is spelled correctly." << endl; return 1; }
	int polyCount = 0;
	ZZ a,b;
	float loadTime, newTime, lastTime;
	Timer myTimer("Integration timer");
	string line,line2;
        loadTime=0;
	while (!myStream.eof())
	{
		polyCount++;
		getline(myStream,line);delSpace(line);
		getline(myStream,line2);delSpace(line2);
 		if ((line.length()==0)||(line2.length()==0)) break;
		lastTime = myTimer.get_seconds();
		myTimer.start();
		if (!strcmp(argv[3],"linear"))
   		{_loadLinForms(_forms, line);
		convertToSimplex(mySimplex,line2);
		_integrateLinFormSum(a, b, _forms, mySimplex);
		}
		else if (!strcmp(argv[3],"polynomial"))
		{_loadMonomials(_monomials, line);
		convertToSimplex(mySimplex, line2);
		//_integrateMonomialSum(a,b,_monomials, mySimplex);
		};
                outStream<<a<<"/"<<b<<endl;
		myTimer.stop();	
		newTime = myTimer.get_seconds();
		FILE* markFile = fopen("integration/mark.txt","w");
		if ( ((!strcmp(argv[3],"linear"))&&(newTime-lastTime>10)) || ((!strcmp(argv[3],"polynomial"))&&(newTime-lastTime>600)) )
		{
			FILE* myFile = fopen("integration/benchmarksh.txt","a");
			fprintf(myFile, "%10.4f", newTime-lastTime);
			fclose(myFile);
			fprintf(markFile, "1");
			fclose(markFile);
			exit(0);
		}
		else loadTime += (newTime - lastTime);	
        };
        myStream.close();
        outStream.close();
	loadTime = loadTime/polyCount;
	FILE* myFile = fopen("integration/benchmarksh.txt","a");
        if (loadTime<0) loadTime=-loadTime;
	fprintf(myFile, "%10.4f", loadTime);
	fclose(myFile);
        FILE* markFile=fopen("integration/mark.txt","w");
	fprintf(markFile, "0");
	fclose(markFile);
	return 0;
}
