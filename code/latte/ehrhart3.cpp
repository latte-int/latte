#include <gmpxx.h>
#include "interpolation/PolynomialInterpolation.h"
#include <iostream>
#include <stdlib.h>
#include <vector>
#include <fstream>

using namespace std;

/* this function is called automatically from ehrhart2
 * passing points to polnomial interpolation, it gets 
 * the answer back and displays the result */

int main(int argc, char *argv[]){
	unsigned int deg;
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
	cout << vec[deg] << "t^" << deg;
	for(int i = deg - 1; i > -1; i--){
		cout << " + " << vec[i] << "t^" << i;
	}
	cout << endl;
	return 0;
}
