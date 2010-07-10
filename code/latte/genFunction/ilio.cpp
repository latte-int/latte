#include <iostream>
#include <NTL/mat_ZZ.h>
#include <NTL/vec_vec_ZZ.h>
#include <fstream>

#include "matrix_ops.h"
#include "config.h"

#include "NTL_to_LiDIA.h"
#include "lidia-include.h"

using namespace LiDIA;
using namespace std;


NTL_CLIENT

NTL_io_matrix_decl(ZZ,vec_ZZ,vec_vec_ZZ,mat_ZZ)

int main(int argc, char *argv[])
{
    int rdim, cdim, loops, mySeed;
    int val;
    
    double nTime, lTime;
    double oldTime;

    cout << "Testing Ilio Smith normal form algorithm..." << endl;
    
    nTime = lTime = 0;
    if (argc < 3)
    {
     //cout << "Usage: ilio numRuns dim [seed]" << endl;
     //return 1;
     //assume testing
	loops = 1;
	rdim = cdim = 11;
    }
    else
    {
	loops = atoi(argv[1]);
	rdim = cdim = atoi(argv[2]);
    }
    if (argc == 4) { mySeed = atoi(argv[3]); } else { mySeed = time(NULL); }
    mat_ZZ M, L, R, X;
    vec_ZZ snf;
    snf.SetLength(rdim);
    srand ( mySeed );
    for (int iter = 0; iter < loops; iter++)
    {
	M = ident_mat_ZZ(rdim);
	L = ident_mat_ZZ(rdim);
	R = ident_mat_ZZ(rdim);
	X = ident_mat_ZZ(rdim);
	snf[0] = to_ZZ(1);
	for (int i = 0; i < rdim; i++)
	{
	    for (int j = 0; j < i; j++)
	    {
		R[j][i] = to_ZZ(rand() % 1000);
		L[i][j] = to_ZZ(rand() % 1000);
	    }
	    if (i > 0) { snf[i] = snf[i-1] * to_ZZ(rand() % 9  + 1); M[i][i] = snf[i]; }
	}

	M = L*M;
	M = M*R;

	ident(L, rdim);
	ident(R, rdim);
	oldTime = GetTime();
	SmithNormalFormLidia(M, L, R); 
	lTime += (GetTime() - oldTime);
	ident(L, rdim);
	ident(R, rdim);
	oldTime = GetTime();
	X = SmithNormalFormIlio(M, L, R); 
	nTime += (GetTime() - oldTime);

	ZZ g;
	  
	for (int i = 0; i < rdim; ++ i) {
				      
	    for (int j = i + 1; j < rdim; ++ j) {
					      
		if (IsOne(snf[i]))  break;
			      
		else if (IsZero(snf[j])) continue;
		      
		else if (IsZero(snf[i])) {
		    std::swap (snf[i], snf[j]);
		}
		      
		else {
      
		    ZZ x, y, myNeg;
		    vec_ZZ tmp1;
			      
		    tmp1.SetLength(rdim);
									      
		    XGCD (g, y, x, snf[j], snf[i]);
			      
		    snf[j] = snf[j] / g;
		      
		    snf[j] = snf[j] * snf[i];
      
		    snf[i] = g;
		}
	    }
	}
      
	M = L * M;
	M = M * R;
	for (int i = 0; i < rdim; i++)
	{
	    if (M[i][i] != X[i][i]) { cerr << "Wrong diagonal values - " << M[i][i] << "=/=" << X[i][i] << endl; return 0; }
	}
	if (M != X)
	{
	    cerr << "Incorrect transform matrices." << endl;
	    cerr << "M is: " << endl;
	    for (int i = 0; i < rdim; i++)
	    {
		for (int j = 0; j < rdim; j++)
		{
		    cerr << M[i][j] << " ";
		}
		cerr << endl;
	    }
	    cerr << "X is: " << endl;
	    for (int i = 0; i < rdim; i++)
	    {
		for (int j = 0; j < rdim; j++)
		{
		    cerr << X[i][j] << " ";
		}
		cerr << endl;
	    }
	    M.kill(); L.kill(); R.kill(); X.kill(); snf.kill();
	    return 1;
	}
	if (abs(determinant(L)) != 1 || abs(determinant(R)) != 1)
	{
	    cerr << "Wrong determinant, " << determinant(L) << " or " << determinant(R) << endl;
	    M.kill(); L.kill(); R.kill(); X.kill(); snf.kill();
	    return 1;
	}
    }
    cout << "LiDIA total of " << lTime << "s. compared to NTL total of " << nTime << "s." << endl;
    
    M.kill(); L.kill(); R.kill(); X.kill(); snf.kill();
    return 0; 
}

