 /**
  * The goal is to test the automatic building of
  * polymake files, dual polymake files, and v- and h-reps.
  */
 
 #include "../../buildPolytopes/BuildPolytope.h"
 #include <iostream>
 
 using namespace std;
 
 int main(void)
 {
 
 	cout << "goig to test building polymake and latte files" << endl;
 	BuildPolytope p;
 	
 	p.setIntegerPoints(true);
 	p.forDebugging();
 	p.setBaseFileName("pie");

 	cout << "making polymake" << endl;	
 	p.buildPolymakeFile();
 	cout << "making latte H-rep" << endl;
	p.buildLatteHRepFile();
	cout << "making latte V-rep" << endl;
	p.buildLatteVRepFile();

	cout << "Making dual polymake file" << endl;
	p.buildPolymakeDualFile();
	cout << "making latte dual V-rep" << endl;
	p.buildLatteVRepDualFile();
	cout << "making latte dual H-rep...not implemented." << endl;

	cout << "Polytope " << (p.isSimplicial() ? "is" : "is NOT") << " simplicial" << endl;
	cout << "Polytope " << (p.isSimple() ? "is" : "is NOT") << " simple" << endl;
	cout << "Dual Polytope " << (p.isDualSimplicial() ? "is" : "is NOT") << " simplicial" << endl;
	cout << "Dual Polytope " << (p.isDualSimple() ? "is" : "is NOT") << " simple" << endl;

 	cout << "\nplease remove the pie.* files yourself" << endl;
 	return 0;
 }
 
 
