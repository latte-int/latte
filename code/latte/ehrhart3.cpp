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
	cout << " + 1" << endl;
	cout << "Number of negative coefficients: " << count << endl;
	return 0;
}
