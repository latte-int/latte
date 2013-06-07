#include <iostream>
#include <fstream>
#include <vector>
#include "TopKnapsack.h"
#include "latte_ntl.h"
#include "integration/burstTrie.h"
#include "integration/multiply.h"
#include <cstdlib>
#include <ctime>
#include <climits>

using namespace std;

int main(int argc, char *argv[])
{
if(1){
	ifstream file;
	file.open(argv[1]);
	int n;
	file >> n;
	vec_ZZ alpha;
	alpha.SetLength(n);

	for(int i =0; i < n; ++i)
		file >> alpha[i];

	cout << "you said" << endl;
	for(int i = 0; i < n ; ++i)
		cout << alpha[i] << ", ";
	cout << endl;

	TopKnapsack tk;
	tk.set(alpha);
	tk.coeff_NminusK(2);
}
else{
	monomialSum t1, t2, ans;
	t1.varCount = 2;
	t2.varCount = 2;
	ans.varCount=2;

	vector<int> xpow, ypow, coef;

	int n = 5;
	for(int i = 0; i < n; ++i)
	{
		xpow.push_back(i);
		ypow.push_back(-i);
		coef.push_back(i*2);
	}

	srand(time(0));
	for(int last = n; last > 0; --last)
	{
		int index = rand()%last;

		swap(xpow[index], xpow[last-1]);
		swap(ypow[index], ypow[last-1]);
		swap(coef[index], coef[last-1]);
	}

	for(int i = 0; i < n; ++i)
	{
		int e[2];
		e[0] = xpow[i];
		e[1] = ypow[i];
		insertMonomial(RationalNTL(coef[i], 1), e, t1);
		insertMonomial(RationalNTL(coef[i], 1), e, t2);
	}

	BTrieIterator<RationalNTL, int>* it = new BTrieIterator<RationalNTL, int>();
	BTrieIterator<RationalNTL, int>* it2 = new BTrieIterator<RationalNTL, int>();

	it->setTrie(t1.myMonomials, t1.varCount);
	it2->setTrie(t2.myMonomials, t2.varCount);
	int low[2], high[2];
	low[0] = low[1] = INT_MIN;
	high[0] = high[1] = INT_MAX;
	multiply<RationalNTL>(it, it2, ans, low, high);

	cout  << printMonomials(ans) << endl;

	destroyMonomials(t1);
	destroyMonomials(t2);
	destroyMonomials(ans);
	delete it;
	delete it2;
}
}
//main()
