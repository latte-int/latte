#include <gmpxx.h>
#include "interpolation/PolynomialInterpolation.h"
#include <iostream>
#include <stdlib.h>
#include <vector>
#include <fstream>



class EhrhartPolyRootFinder
{
private:
	string matlabCommands, mapleCommands, matlabRoots, mapleRoots;
	vector<mpq_class> poly;
public:
	
	EhrhartPolyRootFinder(vector<mpq_class> p);
	
	void saveMatlab();
	void saveMaple();

};

EhrhartPolyRootFinder::EhrhartPolyRootFinder(vector<mpq_class> p): poly(p)
{
	matlabCommands = "ehrhart3.matlabCommands.txt";
	mapleCommands  = "ehrhart3.mapleCommands.txt";
}

void EhrhartPolyRootFinder::saveMatlab()
{
	ofstream file;
	file.open(matlabCommands.c_str());
	
	file << "p = sym('" << poly[0] << "');\n\n" ;
	for(int i = 1; i < poly.size(); ++i)
		file << "p = p + sym('" << poly[i] << " * t^" << i << "');\n";
		
	file << "syms t;\n";
	file << "fptr = fopen('ehrhart.matlabroots.txt', 'w'); %will delete old data; maybe use 'a'?\n";
	file << "theRoots = solve(p);\n\n";
	file << "for i = 1:length(theRoots)\n";
	file << "	if imag(theRoots(i)) ~= 0\n";
	file << "     fprintf(fptr, '%f  %fi => %f\\n', double(real(theRoots(i))), double(imag(theRoots(i))), double(subs(p, t, theRoots(i))));\n";
	file << "  else\n";
	file << "     fprintf(fptr, '%f => %f\\n', double(theRoots(i)), double(subs(p, t, theRoots(i))));\n";
	file << "  end% if complex number\n";
	file << "end%for i\n";
	file << "fclose(fptr);\n";
	file << "exit\n";
	
	file.close();
}//saveMatlab()


void EhrhartPolyRootFinder::saveMaple()
{
	ofstream file;
	file.open(mapleCommands.c_str());
	
	file << "Digits:=1000;\n";
	file << "p := ";
	for(int i = 0; i < poly.size(); ++i)
		file << " + " << poly[i] << "* t^" << i << " ";
	file << ":\n";
	file << "theRoots := fsolve(p = 0, t, complex):\n";
	file << "writeto(\"ehrhart.mapleroots.txt\"):\n";
	file << "for x in theRoots do printf(\"%Zf = > %Zf\\n\", x, subs(t = x, p)): od:\n";
	
	file.close();
}


using namespace std;

/* this function is called automatically from ehrhart2
 * passing points to polnomial interpolation, it gets 
 * the answer back and displays the result */

int main(int argc, char *argv[]){
	unsigned int deg, count = 0;
	vector<mpq_class> vec;
	mpq_class i = 0, temp = 0;
	cin >> deg;
	//cout << deg << endl;
	PolynomialInterpolation p(deg);
	while(cin >> temp){
		//cout << "adding point " << i << ", " << temp << endl;
		p.addPoint(i++, temp);
	}
	vec = p.solve();
	for(int i = deg; i > 0; i--){
		cout << " + " << vec[i] << "*t^" << i;
		if(vec[i] < 0)
			count++;
	}
	if ( vec[0] != mpq_class(1))
	{
		cout << "ehrhart3:: vec[0] = " << vec[0] << " != 1" << endl;
		exit(1);
	}//check to make sure there was no errors
	cout << " + 1" << endl;
	cout << "Number of negative coefficients: " << count << endl;
	
	EhrhartPolyRootFinder finder(vec);
	finder.saveMatlab();
	finder.saveMaple();
	return 0;
}
