//This program reads a simplex and a linear form a specified file and computes the integral over that simplex
#include "myIntegration.h"

#include <iostream>
#include <fstream>
#include <NTL/vec_ZZ.h>
#include <NTL/vec_vec_ZZ.h>
#include <NTL/ZZ.h>

NTL_CLIENT;

using namespace std;

int main(int argc, char *argv[])
{
  int d,m,i;
  vec_ZZ lForm; //for now, suppose there is only one linear form with power m.
  vec_vec_ZZ s;
        if (argc < 2) { cout << "Usage: ./integrate filename" << endl; return 1; }
	ifstream myStream (argv[1]);
	if (!myStream.is_open()) { cout << "Error opening file " << argv[1] << ", please make sure it is spelled correctly." << endl; return 1; }
	myStream >> d;	                                  	  //d is the dimension
	lForm.SetLength(d);
	s.SetLength(d+1);
	for (i=0;i<=d;i++) s[i].SetLength(d);
	myStream >> m;                                            //raised to mth power; 
	for (i=0;i<d;i++) myStream >>lForm[i];                    //linear form
	for (i=0;i<=d;i++)
	  {for (int j=0;j<d;j++)
	    {
	      myStream>>s[i][j];                                  //d+1 vectors representing vertices
	    };
	  };
	print_Integrate(d, m, lForm, s);
	myStream.close();
	return 0; 
}
