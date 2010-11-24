 
 
 #include "BuildPolytope.h";
 #include <iostream>
 
 using namespace std;
 
 int main(void)
 {
 
 	cout << "goig to test stuff!!!" << endl;
 	BuildPolytope p;
 	
 	p.forDebugging();
 	p.setBaseFileName("pie");

 	cout << "making polymake" << endl;	
 	p.buildPolymakeFile();
 	cout << "making latte H" << endl;
	p.buildLatteHRepFile();
	cout << "makeing latte V" << endl;
	p.buildLatteVRepFile();
 	
 	return 0;
 }
 
 
