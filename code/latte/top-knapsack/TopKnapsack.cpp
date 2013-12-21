/*
 * TopKnapsack.cpp
 *
 *  Created on: Apr 23, 2013
 *      Author: bedutra
 */

#include "TopKnapsack.h"
#include "preprocess.h"
#include "barvinok/barvinok.h"
#include "barvinok/dec.h"
#include "integration/multiply.h"
#include "print.h"
#include "dual.h"
#include "timing.h"
#include <typeinfo>
#include <climits>
#include <algorithm>




MobiusPair::MobiusPair():gcd(to_ZZ(0)), mu(to_ZZ(0)), mobiusValid(false){}
MobiusPair::MobiusPair(const ZZ& g, const ZZ& m): gcd(g), mu(m), mobiusValid(false) {}


// **************************************************************************************************
// **************************************************************************************************


MobiusList::MobiusList() {}

void MobiusList::insertGCD(const ZZ& v )
{
	bool found = false;
	for(int i = 0; i < list.size(); ++i)
		if ( list[i].gcd == v)
		{
			found = true;
			break;
		}

	if( !found )
		list.push_back(MobiusPair(v, to_ZZ(0)));
}
void MobiusList::computeMobius()
{
	int indexOne =-1;
	for(int i = 0; i < list.size(); ++i)
	{
		list[i].mu = 0;
		list[i].mobiusValid = false;
		if(list[i].gcd == 1)
			indexOne = i;
	}

	if (indexOne == -1)
		THROW_LATTE_MSG(LattException::bug_Unknown, "gcd of entire list is not one");
	computeMobius(indexOne);
}

void MobiusList::computeMobius(int i)
{
	if( list[i].mobiusValid == true)
		return;

	ZZ sum;
	sum = 0;
	for(int j = 0; j < (int) list.size(); ++j)
		if( i != j  &&  divide(list[j].gcd, list[i].gcd) )
		{
			if( list[j].mobiusValid == false)
				computeMobius(j);
			sum += list[j].mu;
		}
	list[i].mu = to_ZZ(1) - sum;
	list[i].mobiusValid = true;
}

void MobiusList::print() const
{
	for(int i = 0; i < (int) list.size(); ++i)
		cout << list[i].mu << ", gcd=" << list[i].gcd << ", isv=" << list[i].mobiusValid << endl;
}


void MobiusSeriesList::computeMobius()
{
	MobiusList::computeMobius();
	unweightedSeries.resize(list.size());
	for(int i = 0; i < (int)unweightedSeries.size(); ++i)
		unweightedSeries[i] = NULL;
}

MobiusSeriesList::~MobiusSeriesList()
{
	for(int i = 0; i < (int)unweightedSeries.size(); ++i)
		delete unweightedSeries[i];
}

// **************************************************************************************************
// **************************************************************************************************
/*
 * Bernoulli of the 2nd kind:
 * Computes B[0], B[1], .... B[k] where
 * t/(1-exp(-t)) = \sum_{m=0}^infinity B[m] (-t)^m/m!
 * See https://en.wikipedia.org/wiki/Bernoulli_number#Generating_function
 *
 * Formula:
 * B[0]= 1;
 * B[m] = 1 - \sum_{j=0}^{m-1} choose(m,j)B[j]/(m-j+1)
 * See https://en.wikipedia.org/wiki/Bernoulli_number#Explicit_definition
 *
 * Then note that Bernoulli of the fist kind is equal to the second kind only b[1] = -1/2, instead of 1/2
 */
void BernoulliFirstKind::setBernoulli(int k)
{
	//compute the 2nd kind
	B.resize(k+1);
	B[0] = 1;

	int m, j;
	for(m = 1; m <=k; m++)
	{

		//TODO: replace the choose function with pascal's triangle type computation to save time.
		//cout << "starting B" << m << endl;
		B[m] = 0;
		if (m > 1 && (m % 2) == 1)
		{
			continue; //B[m] = 0 for odd m larger than 1.
		}
		for( j = 0; j < m; ++j)
		{
			B[m] += B[j]*RationalNTL(TopKnapsack::binomial(m,j), to_ZZ(m-j+1));;
			//cout << " + " << B[j] << "* " << TopKnapsack::binomial(m,j) << "*1/" << m-j+1 << " ";
		}
		//cout << "\nbefore change" << B[m] << endl;
		B[m].changeSign();
		B[m] += RationalNTL(1,1);
	}

	//convert to the first kind.
	if (1 <= k)
		B[1] *= -1;

}


const RationalNTL& BernoulliFirstKind::operator[](int i) const
{
	if(i >= B.size())
		THROW_LATTE(LattException::bug_Unknown);

	return B[i];
}
// **************************************************************************************************
// **************************************************************************************************
// **************************************************************************************************
// **************************************************************************************************
TopKnapsack::TopKnapsack() {
	// TODO Auto-generated constructor stub

}

TopKnapsack::~TopKnapsack() {
	// TODO Auto-generated destructor stub
}

void TopKnapsack::set(const vec_ZZ& list)
{
	alpha = list;
	N = alpha.length()-1;
	bernoulli.setBernoulli(alpha.length());

	//for(int i = 0; i <= N; ++i)
	//	cout << "bernoulli[" << i << "]=" << bernoulli[i] << endl;
}

// compute t^N+t^{N-1}+....+t^{N-k}
void TopKnapsack::coeff_topK(int k)
{
	topKTerms = true;
	coeff(k);
}

//to do: clean at start of run in case user calls this twice or do this in the set() method...?
void TopKnapsack::coeff_NminusK(int k)
{
	topKTerms = false;
	coeff(k);
}

void TopKnapsack::coeff(int k)
{

	assert(0 <= k && k<= N);
	order = k;
	cout << "order=" << order << endl;
	coeffsNminusk.resize(k+1);


	Timer tgcd("Time for gcds");
	tgcd.start();
	for(int i = 0; i <= k; ++i)
		everyGCD(N+1-i);
	gcds.computeMobius();
	tgcd.stop();
	cout << tgcd << endl;

	//cout << "mu found" << endl;
	//gcds.print();


	for(int i = 0; i < (int)gcds.list.size(); ++i)
		if (gcds.list[i].mu != 0)
		{
			E(i);
		}

	packageAnswer();

}

void TopKnapsack::packageAnswer()
{
	for(int i = 0; i < (int)gcds.list.size(); ++i)
		if (gcds.list[i].mu != 0)
		{
			GeneralMonomialSum<PeriodicFunction, int> *series = gcds.unweightedSeries[i];
			if ( series->termCount == 0)
				continue;

			BTrieIterator<PeriodicFunction, int> * itr = new BTrieIterator<PeriodicFunction, int>();
			itr->setTrie(series->myMonomials, series->varCount);
			itr->begin();
			term<PeriodicFunction, int> * oneTerm;
			while( (oneTerm = itr->nextTerm()) )
			{

				PeriodicFunction p(oneTerm->coef);
				int deg = oneTerm->exps[1]; // deg = -N-1 +i for some 0<=i<=k (=order)
				//deg += N+1; //= i
				int h = N-deg; //h= N-i

				ZZ hFactorial;
				hFactorial = 1;
				for(int i = 2; i <= h; ++i)
					hFactorial *= i;

				//cout << "hFactorial" << hFactorial << endl;
				RationalNTL factor;
				if ( h % 2 == 0)
					factor = -1;
				else
					factor = 1;
				factor *= gcds.list[i].mu;
				factor *= gcds.list[i].gcd;
				factor.div(hFactorial);

				//cout << "h=" << h << "FACTOR=" << factor << endl;
				p.times(factor);
				coeffsNminusk[deg].add(p);
			}
			delete itr;
		}

}

//mostly for debugging
void TopKnapsack::printAnswer(ostream & out)
{

	if ( topKTerms == false)
		out << "coeff" << N << "minus" << order << ":= " << coeffsNminusk[order] << ";\n";
	else
	{

		for(int i = 0; i < (int) coeffsNminusk.size(); ++i)
		{
			out << "coeff" << N << "minus" << i << ":= " << coeffsNminusk[i] << ";\n";
		}

		out << "\ntopKPolynomial:=" ;
		for(int i = 0; i < (int) coeffsNminusk.size(); ++i)
		{
			if ( i > 0)
				out << " + ";
			out << "(coeff" << N << "minus" << i << ")*T^(" <<N-i << ")";
		}
		out << ";" << endl;
	}
}

void TopKnapsack::E(int fIndex)
{
	ZZ f = gcds.list[fIndex].gcd;

	assert(gcds.unweightedSeries[fIndex] == NULL);

	gcds.unweightedSeries[fIndex] = new GeneralMonomialSum<PeriodicFunction, int>;
	gcds.unweightedSeries[fIndex]->varCount = 2;


	if (f == 1)
	{
		expandF1Case(*(gcds.unweightedSeries[fIndex]));
		return;
	}
	//cout << "f=" << f << endl;

	vector<ZZ> fDivAlpha, fnDivAlpha; //f (not) divides alpha
	for(int i = 0; i < alpha.length(); ++i)
	{
		if ( divide(alpha[i], f))
			fDivAlpha.push_back(alpha[i]);
		else
			fnDivAlpha.push_back(alpha[i]);
	}

	//cout << "fDivAlpha" << endl;
	//for(int i =0;i < fDivAlpha.size(); ++i)
	//	cout << fDivAlpha[i] << ", ";
	//cout << "\nfnDivAlph" << endl;
	//for(int i = 0; i < fnDivAlpha.size(); ++i)
	//	cout << fnDivAlpha[i] << ", " ;
	//cout << endl;

	mat_ZZ latticeBasis;
	latticeBasis.SetDims(fnDivAlpha.size(), fnDivAlpha.size());
	findLatticeBasis(latticeBasis, fnDivAlpha, f);
	//cols of latticeBasis generate a cone of all x s.t. \sum \alpha_i x_i \in f\Z

	//cout << "lattice basis (in the cols)" << endl;
	//TopKnapsack::printMatrix(latticeBasis);

	mat_ZZ invLatticeBasis, invLatticeBasisScaled;
	ZZ invLatticeBasisD;
	invLatticeBasis.SetDims(fnDivAlpha.size(), fnDivAlpha.size());
	invLatticeBasisScaled.SetDims(fnDivAlpha.size(), fnDivAlpha.size());
	findLatticeBasisInv(invLatticeBasisScaled, invLatticeBasisD, invLatticeBasis, latticeBasis);

	//cout << "cone in basis lambda" << endl;
	//TopKnapsack::printMatrix(invLatticeBasis);


	//find the vertex.
	vec_ZZ tVertex;
	tVertex.SetLength(fnDivAlpha.size());
	findVertex(tVertex, f, fnDivAlpha);

	//cout << "tvertex=" << endl;
	//for(int i = 0; i < tVertex.length(); ++i)
	//	cout << tVertex[i] << ", " ;
	//cout << endl;

	//static Timer ttri("Finding uni cones so far");
	//ttri.start();
	listCone* uniCones = findUnimodularCones(invLatticeBasisScaled); //todp: move the latice scaling to this function
	//ttri.stop();
	//cout << ttri << endl;

	bool finishedResidue = false;

	//static Timer tresidu("residue so far time: ");
	//tresidu.start();
	while ( !finishedResidue)
	{
		try
		{
			findResidue(*(gcds.unweightedSeries[fIndex]), tVertex, uniCones, latticeBasis, invLatticeBasis, invLatticeBasisD, fnDivAlpha, fDivAlpha);
			finishedResidue = true;
		} catch (LattException & e)
		{
			delete gcds.unweightedSeries[fIndex];
			gcds.unweightedSeries[fIndex] = new GeneralMonomialSum<PeriodicFunction, int>;
			gcds.unweightedSeries[fIndex]->varCount = 2;
		}
	}//while.
	//tresidu.stop();
	//cout  << tresidu << endl;
	
	freeListCone(uniCones);

}

/**
 * Find all solutions to Ax:=[fnDivAlpha]x = f\Z where A is the row vector [fnDivAlpha].
 *
 * Then x^TA^T = f.
 * Let y^TUA^T = y^T H = f. H is a col vector with H=(h1, 0,0,....0)^T.
 *    where x^T = y^TU => u^T y = x.
 * If y1 = f/h1, and y2 = ... =yn =0  then Ax=f, but x may not be integral.
 * So we try Ax = 2f, giving y1=2f/h1. Etc.
 *
 * Another way to view this, is to look at U*A^T = H = (h1, 0,0,....0)^T. Then multiplying
 * the first row on both sides by f/gcd(f,h1), we get
 * (U with scaled 1st row)*A^T = (lcm(h1,f), 0,0,....0)^T which gives the smallest positive f\Z we can get.
 *
 * @param latticeBasis: output parameter. The cols contain the the basis such that \sum latticeBasis[:,i]*\Z \in f\Z
 * @param fnDivAlpha: list of numbers
 * @param f: number

 */
void TopKnapsack::findLatticeBasis(mat_ZZ &latticeBasis, const vector<ZZ> & fnDivAlpha, const ZZ & f) const
{
	/*
	  testing

	mat_ZZ A;
	int n = 4;  //row
	int m = 5; //col
	A.SetDims(n,m);
	for(int i = 0; i < n; ++i)
		for(int j = 0; j < m; ++j)
			A[i][j] = 2*i + j;

	vec_ZZ S;
	S.SetLength(n*m);
	for(int i = 0; i < n; ++i)
		for(int j = 0; j < m; ++j)
			S[i*m + j] = A[i][j]; //row order
	vec_ZZ U;
	U.SetLength(n*n);
	vec_ZZ rhs;
	rhs.SetLength(n);
	for(int i =0; i < n; ++i)
		rhs[i] = i+1;
	cout << " rhs " << rhs << endl;



	int r = ihermite(&S, &U, &rhs, m,n);
	cout << "got here" << endl;
	mat_ZZ umat, smat;
	umat.SetDims(n,n);
	smat.SetDims(n,m);

	cout << S << endl;
	cout << U << endl;

		//convert U to matZZ in col-major order.
		for(int i = 0; i < n; ++i)
			for(int j = 0; j < m; ++j)
				smat[i][j] = S[i*m + j];
		for(int i = 0; i < n; ++i)
			for(int j = 0; j < n; ++j)
				umat[i][j] = U[j*n +i];
		cout << "Amat" << endl;
		TopKnapsack::printMatrix(A);

		cout << "umat" << endl;
		TopKnapsack::printMatrix(umat);


		cout << "smat" << endl;
		TopKnapsack::printMatrix(smat);



		cout << "rhs " << rhs << endl;
exit(1);
  */

	vec_ZZ s, u, rhs;

	int n = fnDivAlpha.size();
	int m = 1;

	s.SetLength(n);//column vector
	for(int i = 0; i < n; ++i)
		s[i] = fnDivAlpha[i];


	u.SetLength(n*n);
	rhs.SetLength(n);

	//For us,
	//	s= fnDivAlpha
	//	u = the col vectors give the basis for the lattice (before scaling)
	//	rhs = junk.

	//cout << "before " << s << endl;
	int r = ihermite(&s, &u, &rhs, m,n);

	//U is in row-major order. take the transpose and save taht in latticeBasis
	for(int i = 0; i < n; ++i)
		for(int j = 0; j < n; ++j)
			latticeBasis[i][j] = u[j*n+i];
	/*
	cout << "after ihermite" << endl;
	TopKnapsack::printMatrix(latticeBasis);
    cout << "after " << s << endl;


	cout << "s:" << endl;
	for(int i = 0; i < m; ++i)
	{
		for(int j = 0; j < n; ++j)
			cout << s[m*j+i] << ", ";
		cout << endl;
	}

	cout << "u:" << endl;
	for(int i = 0; i < n; ++i)
	{
		for(int j = 0; j < n; ++j)
			cout << u[n*j+i] << ", ";
		cout << endl;
	}

	cout << "rhs:" << endl;
	for(int i =0; i < m; i++)
		cout << rhs[i] << ", " ;
	cout << endl;

	cout << "r:" << r << endl;
*/


	ZZ newf;
	divide(newf, f, GCD(f,s[0])); //newf = f/gcd(f, u[0][0])
	for(int i = 0; i < n; ++i)
		latticeBasis[i][0] *= newf;


}

void TopKnapsack::findLatticeBasisInv(mat_ZZ & invLatticeBasisScaled, ZZ & invLatticeBasisD, mat_ZZ & invLatticeBasis,const mat_ZZ & latticeBasis) const
{
	ZZ d;
	inv(invLatticeBasisD, invLatticeBasisScaled, latticeBasis);
	invLatticeBasis = invLatticeBasisScaled;
	// (1/d)*invLatticeBasis = latticeBasis^{-1}
	//for latticeBasis^{-1}, scale each col to integer, and divide by gcd. so we can ignore d.
	// However, if d is negative, we do need to mult by -1

	bool changeSign = false;
	if(invLatticeBasisD < 0)
	{
		changeSign = true;
	}

	//scale each col.
	for(int j = 0; j < invLatticeBasisScaled.NumCols(); ++j)
	{
		d=invLatticeBasisScaled[0][j];
		for(int i = 1; i < invLatticeBasisScaled.NumRows(); ++i)
			d = GCD(d, invLatticeBasisScaled[i][j]); //see if we can scale col j
		if ( changeSign)
			d *= -1;

		if ( d != 1 )
		{
			for(int i = 0; i < invLatticeBasisScaled.NumRows(); ++i)
				invLatticeBasisScaled[i][j] /= d;
		}
	}
}

/**
 * Find s1, s2, ..., sI, s_{I+1}, s.s 1 = s1 a_1 + ... s_I a_I + s_{I+1}f = 1
 * where I = | alpha's s.t. f does not divide alpha_i | (call these elements of alpha a_i
 * Solve (x1, ..., x_{I+1})*(a1, a2, .., aI, f)^T = (1, 0, ...., 0)^T
 * y^T UA = y^T H = (1, 0, ...., 0)^T where y^T U = x^T.
 * Because we know the gcd(a1, .., aI, f) = 1, H = (1, 0, .., 0)^T.
 * So y1 = 1, y2=...y_{I+1} = 0.
 * Then x = 1st row of U.
 */
void TopKnapsack::findVertex(vec_ZZ & tVertex, const ZZ &f, const vector<ZZ> &fnDivAlpha) const
{
	vec_ZZ s, u, rhs;
	int n = fnDivAlpha.size()+1; //num rows.
	int m = 1; //num cols.

	s.SetLength(n);
	for(int i = 0; i < n-1; ++i)
		s[i] = fnDivAlpha[i];
	s[n-1]=f;


	u.SetLength(n*n);
	rhs.SetLength(n); //junk.

	int r = ihermite(&s, &u, &rhs, m,n);

	//only need 1st row in U, and in this row, the last col is deleted
	//again, u is in row-major order.
	for(int i = 0; i < n-1; ++i)
		tVertex[i]= u[i];
}


listCone* TopKnapsack::findUnimodularCones(const mat_ZZ & invLatticeBasis)const
{
	listCone * newCone = createListCone();

	mat_ZZ rays = transpose(invLatticeBasis);
	//mat_ZZ rays = invLatticeBasis;


	newCone->vertex = new Vertex;
	newCone->vertex->vertex = new rationalVector(invLatticeBasis.NumRows()); //make a zero vertex.

	newCone->rays = new listVector(rays[0]);

	for(int i = 1; i < rays.NumRows(); ++i)
		newCone->rays = appendVectorToListVector(rays[i], newCone->rays);

	//cout << "newCone before unimodular" << endl;
	//printListCone(newCone, rays.NumRows());
	//cout << "end newCone before unimodular" << endl;

	BarvinokParameters params;
	params.Number_of_Variables = rays.NumRows();
	params.max_determinant = 1;
	listCone *uniTriangulation = decomposeCones(newCone, true, params);


	//cout << "uniTriangulation after unimodular" << endl;
	//printListCone(uniTriangulation, rays.NumRows());

	freeListCone(newCone);

	return uniTriangulation;
}

/*
 * Let xv:=(alpha_I + ae)x where a is a random vector, e is epsilon
 *  exp(<B^T xv, {-TB^{-1}s}_u>)
 *  ---------------------------- ---------------
 *  \prod (1-exp(<B^Txv, g_i>)) \prod (1-e^{a_ixv})
 *   g_i is the ith ray of the unimodular cone u
 *
 *   u should have its rays and facets computed.
 *
 */
void TopKnapsack::findResidue(GeneralMonomialSum<PeriodicFunction, int> & fSeries, const vec_ZZ & s, const listCone *u, const mat_ZZ & B,
		const mat_ZZ & invB, const ZZ & invBd,
		const vector<ZZ> &fnDivAlpha, const vector<ZZ> & fDivAlpha)
{
	int seed = 654354; // for debugging. TODO: change to time(0);
	srand(seed);

	//get random beta vector.
	int I = fnDivAlpha.size();
	int i;
	vec_ZZ beta;
	beta.SetLength(I);
	for(i = 0; i < I; ++i)
		beta[i] = ( rand()%2 ? -1*rand()%500 : rand()%500); //set beta to random vector in [-500,500]^I
	//cout << "beta = " << beta << endl;

	vec_ZZ temp;
	temp.SetLength(I);

	//copy fnDivide to a vec_ZZ
	vec_ZZ alpha;
	alpha.SetLength(I);
	for(i = 0; i < I; ++i)
		alpha[i] = fnDivAlpha[i];

	//store exp(< (alpha  + beta e)x, Bg_i>)...alpha,beta are vectors, x is scalar.
	vector<ZZ> expa, expe;
	expa.resize(I);// expa[i] = alpha * B*g_i *    x
	expe.resize(I);// expe[i] = beta * B*g_i * e * x


	GeneralMonomialSum<PeriodicFunction, int> fDivAlphaExpansion;
	expandNonperiodicPart(fDivAlphaExpansion, fDivAlpha);
	ZZ bottomCoeffNonperiodicPart;
	bottomCoeffNonperiodicPart = 1;
	for(i = 0; i < (int)fDivAlpha.size();++i)
		bottomCoeffNonperiodicPart *= fDivAlpha[i];

	//cout << "non periodic part is (missing sign/X^"<< fDivAlpha.size() << "*" << bottomCoeffNonperiodicPart << ")* ==" << fDivAlphaExpansion.printMonomials().c_str() << endl;


	int minE[2] = {INT_MIN, INT_MIN};
	int maxE[2] = {0, 0};


	//(u->rays->first)[0] = 0;
	//(u->rays->first)[1] = 1;

	//cout << "Ba" << transpose(B)*alpha << endl;

	for(const listCone * oneCone = u; oneCone; oneCone = oneCone->rest)
	{
		static long int abcd = 1;
		abcd++;
		//if ( 414 != abcd)
		//	continue;

		//cout << "#####################################################" << endl;
		//printCone((listCone *) oneCone, I);

		//compute all the inner produces in 1/(1-exp(<(alpha  + beta e)*x, Bg>))
		i=0;
		int numPoles = 0;

		for(const listVector * g = oneCone->rays; g; g = g->rest)
		{
			vec_ZZ Bg;
			Bg.SetLength(I);
			mul(Bg, B, (g->first));
			InnerProduct(expa[i], alpha, Bg);
			InnerProduct(expe[i], beta, Bg);

			if ( IsZero(expa[i]) && IsZero(expe[i]) )
			{
				cout << "beta not random enough";
				THROW_LATTE_MSG(LattException::de_divisionByZero, "trying new random perturbation");
			}
			if( IsZero(expa[i]))
				numPoles++;
			i++;
		}//for g

		//cout << "expa: ";
		//for(i = 0; i < expa.size(); ++i)
		//	cout << expa[i] << ", ";
		//cout << "\nexpe: " ;
		//for(i = 0; i < expe.size(); ++i)
		//	cout << expe[i] << ", ";
		//cout << endl;
		//now,  1/(1-exp(<(alpha  + beta e)*x, Bg>)) = 1/(1-exp(x*expa + x*e*expe))


		GeneralMonomialSum<PeriodicFunction, int> fnDivAlphaExpansion;
		ZZ bottomCoeffPeriodicPart;
		expandPeriodicPart(bottomCoeffPeriodicPart, fnDivAlphaExpansion, numPoles, expa, expe);




		//Goal:compute the numerator exp(<B^T (alpah + beta e)x, {-B^-1 T s}> where {-x} = fractionalPart(x)

		//first, put the rays as the column of the matrix coneRays.
		mat_ZZ coneRays;
		coneRays.SetDims(I, I);
		int row = 0;
		for(const listVector * g = oneCone->rays; g; g = g->rest)
		{
			for(int col = 0; col < I; ++col)
				coneRays[col][row] = g->first[col];
			++row;
		}

		mat_ZZ coneRaysInv;
		ZZ coneRaysInvD;
		inv(coneRaysInvD, coneRaysInv, coneRays);

		temp = coneRaysInv * (invB * s);
		//now, B^-1s = \sum temp[i]/(invBd*coneInvD) * g_i where g_i is the rays of the cone

		//cout << temp << " divided by " << invBd << "*" << coneRaysInvD << endl;

		vector<RationalNTL> fractionalPart;
		fractionalPart.resize(I);
		for(int i = 0; i < I; ++i)
		{
			fractionalPart[i] = RationalNTL(temp[i], invBd*coneRaysInvD);
			if (fractionalPart[i].getDenominator() == 1)
				fractionalPart[i] = 0; //fractionalPart(integer) = 0
			//cout << periodicTerm[i] << ", ";
		}

		//now, fractionalPart[i] = {-B^-1 T s}_i
		vec_ZZ fractionalPartCoeffa, fractionalPartCoeffe;
		fractionalPartCoeffa.SetLength(I);
		fractionalPartCoeffe.SetLength(I);

		//compute < B^T(alpha + beta e)x, \sum fractionalPart[i] * g_i>
		i = 0;
		for(const listVector * g = oneCone->rays; g; g = g->rest)
		{
			temp  = B* (g->first); //B*g_i
			InnerProduct(fractionalPartCoeffa[i], alpha, temp);
			InnerProduct(fractionalPartCoeffe[i], beta, temp);
			++i;
		}

		//now, < B^T(alpha + beta e)x, \sum fractionalPart[i] * g_i> = \sum  ( fractionalPartCoeffa[i] + fractionalPartCoeffe[i]*e)*fractionalPart[i]*x
		//Next: write this as two periodic functions
		PeriodicFunction pa, pe;

		//for(i = 0; i < I; ++i)
		//{
		//	cout << "sum+=" << "(" <<fractionalPartCoeffa[i]  << "+" << fractionalPartCoeffe[i]<<"*e)*{" << fractionalPart[i]<< "}*x" << endl;
		//}




		GeneralMonomialSum<PeriodicFunction, int> exExpansion;
		expandExponentialPart(exExpansion, numPoles, fractionalPartCoeffa, fractionalPartCoeffe, fractionalPart);

		minE[0] = 0;
		minE[1] = 0;

		maxE[0] = numPoles;
		maxE[1] = order;

		//cout << "fnDivAlphaExpansion " << fnDivAlphaExpansion.printMonomials().c_str() << endl;
		//cout << "fDivAlphaExpansion  " << fDivAlphaExpansion.printMonomials().c_str() << endl;
		//cout << "exExpansion         " << exExpansion.printMonomials().c_str() << endl;

		fnDivAlphaExpansion.multiply(fDivAlphaExpansion, minE, maxE);
		fnDivAlphaExpansion.multiply(exExpansion, minE, maxE);
		//fDivAlphaExpansion.destroyMonomials(); //cannot free it, used again on next loop.
		exExpansion.destroyMonomials();



		//cout << "fnDivAlphaExpansion after...." << fnDivAlphaExpansion.printMonomials().c_str() << endl;
		//cout << "fDivAlphaExpansion  " << fDivAlphaExpansion.printMonomials().c_str() << endl;
		//cout << "exExpansion         " << exExpansion.printMonomials().c_str() << endl;


		RationalNTL scale;
		//scale = (-1)^{N+1}
		if ( (N+1) % 2 == 0)
			scale = 1;
		else
			scale = -1;
		scale.div(bottomCoeffNonperiodicPart * bottomCoeffPeriodicPart);

		scale *= oneCone->coefficient;


		PeriodicFunction scaleTerm;
		scaleTerm.setToConstant(scale);
		GeneralMonomialSum<PeriodicFunction, int> finalProductScale;
		finalProductScale.varCount = 2;
		//maxE[0] = -1* numPoles;
		//maxE[1] = -1*(N+1);
		maxE[0] = 0;
		maxE[1] = 0;
		finalProductScale.insertMonomial(scaleTerm, maxE);



		//maxE[0] = 0;
		//minE[0] = 0;

		//maxE[1] = -1*N -1 + order;
		//if ( topKTerms)
		//	minE[1] = -1*N -1;
		//else
		//	minE[1] = -1*N -1 + order;

		maxE[0] = numPoles;
		minE[0] = numPoles;

		maxE[1] = order;
		if ( topKTerms)
			minE[1] = 0;
		else
			minE[1] = order;

		//cout << "before.....full cone expansion" << fnDivAlphaExpansion.printMonomials().c_str() << endl;
		fnDivAlphaExpansion.multiply(finalProductScale, minE, maxE);
		//cout << "after.....full cone expansion" << fnDivAlphaExpansion.printMonomials().c_str() << endl;


		fSeries.add(fnDivAlphaExpansion);

	}//for oneCone

	//cout << "fseries is for all the cones" << fSeries.printMonomials().c_str() << endl;

}


/**
 * Expand 1/(\prod_{i: f does not divide \alpha_i} (1-exp(alaph_i*x)))
 *
 * Note: 1/(1-exp(ax)) = -1/(ax) \sum_{m=0}^inf B_m (ax)^m/m!
 * We compute -1/(ax) \sum_{m=0}^{k} B_m (ax)^m/m!
 *          = \sum_{m=0}^{k} -1 B_m (ax)^{m-1}/m!
 *          = -1B_0 (1/ax) - B_1/1! - B_2 (ax)^1/2! - ... - B_k (ax)^{k-1}/k!
 */
void TopKnapsack::expandNonperiodicPart(GeneralMonomialSum<PeriodicFunction, int> &a, const vector<ZZ>  & fDivAlpha)
{

	a.varCount = 2;
	a.setToConstant(PeriodicFunction(RationalNTL(1,1), true));

	int min[2] = { INT_MIN, INT_MIN}; //does not really matter.
	int max[2] = {0, INT_MAX};
	max[1] = order;

	//max[2] = (int)fDivAlpha.size();

	ZZ mFract;

	int exponents[2];
	exponents[0] = 0;

	for(int j = 0; j <  (int)fDivAlpha.size(); ++j)
	{

		GeneralMonomialSum<PeriodicFunction, int> oneExpansion;
		oneExpansion.varCount = 2;
		mFract = 1;
		RationalNTL am(1, 1);
		for(int m = 0; m <= order; ++m)
		{
			RationalNTL coeff(am);
			coeff.div(mFract);
			coeff *= bernoulli[m];
			exponents[1] = m;

			PeriodicFunction pCoeff;
			pCoeff.setToConstant(coeff);
			oneExpansion.insertMonomial(pCoeff, exponents);

			//cout << "one Expansion " << oneExpansion.printMonomials().c_str() << endl;
			//cout << "coeff=" << coeff << "am=" << am << " bern[" << m << "]=" << bernoulli[m] << endl;

			mFract *= (m+1);
			am *= fDivAlpha[j];
		}//insert one series expression

		//cout << "before multiply" << endl;
		//cout << "one Expansion" << oneExpansion.printMonomials().c_str() << endl;
		//cout << "a " << a.printMonomials().c_str() << endl;
		a.multiply(oneExpansion, min, max);
		//cout << "after multiply" << endl;
		//cout << "a " << a.printMonomials().c_str() << endl;
	}//for j
	//cout << "\n\n\n" << endl;
}

/*
* Note: 1/(1-exp(ax)) = -1/(ax) \sum_{m=0}^inf B_m (ax)^m/m!
* We compute -1/(ax) \sum_{m=0}^{k} B_m (ax)^m/m!
*          = \sum_{m=0}^{k} -1 B_m (ax)^{m-1}/m!
*          = -1B_0 (1/ax) - B_1/1! - B_2 (ax)^1/2! - ... - B_k (ax)^{k-1}/k!
*/
void TopKnapsack::expandPeriodicPart(ZZ & bottomCoeffPeriodicPart, GeneralMonomialSum<PeriodicFunction, int> & a, const int numPoles, const vector<ZZ> & expa, const vector<ZZ> & expe)
{
	a.varCount = 2;
	a.setToConstant(PeriodicFunction(RationalNTL(1,1), true));
	bottomCoeffPeriodicPart = 1;

	RationalNTL coeff;
	int exponents[2];
	int minE[2] = {INT_MIN, INT_MIN};
	int maxE[2];
	maxE[0] = numPoles;
	maxE[1] = order;



	GeneralMonomialSum<PeriodicFunction, int> oneExpansion;
	PeriodicFunction p;

	//expand 1/(1-exp((expa + expe*e)*x))
	for(int i = 0; i < (int) expa.size(); ++i)
	{
		oneExpansion.destroyMonomials();
		oneExpansion.varCount = 2;
		ZZ mFract;
		mFract= 1;

		for(int m = 0; m <= order; ++m)
		{

			coeff = 1;
			coeff.div(mFract);
			coeff *= bernoulli[m];
			//coeff = -B[m]/(m)! for m = 0, 1

			exponents[1] = m;

			//expand coeff*(expa + expe*e)^mx^m
			//insert coeff * x^m ( \sum_{j=0}^{j=m} (m choose j) *(expe^j*e^j * expa^{m-j})  )
			for(int j = 0; j <= min(numPoles,m); ++j)
			{
				RationalNTL newCoeff(coeff);
				newCoeff.mult(TopKnapsack::binomial(m,j));
				newCoeff.mult( power(expe[i],j) * power(expa[i], m-j));
				exponents[0] = j;

				p.setToConstant(newCoeff);
				oneExpansion.insertMonomial(p, exponents);
			}
			mFract *= (m+1);
		}//for m.

		oneExpansion.check();
		a.check();
		a.multiply(oneExpansion, minE, maxE);
	}//for i


	int powerX = 0;
	int powerE = 0;
	coeff = 1;
	for(int i = 0; i < (int)expa.size(); ++i)
	{
		if ( expa[i] == 0 )
		{
			bottomCoeffPeriodicPart *= expe[i];
			powerX++;
			powerE++;
		}//expa = 0, but expe != 0. Then (1/ax) = 1/ (expe*e*x)
		else if ( expe[i] == 0 || numPoles == 0)
		{
			bottomCoeffPeriodicPart *= expa[i];
			powerX++;
		}// expa != 0, and expe = 0. Then (1/ax) = 1/ (expa *x)
		else
		{
			//cout << "good, this happened" << endl;
			//exit(1);
			powerX++;
			oneExpansion.destroyMonomials();
			oneExpansion.varCount = 2;
			PeriodicFunction p;
			exponents[1] = 0;
			for(int m = 0; m <= numPoles; ++m)
			{
				RationalNTL bcoeff;
				if( m % 2 == 0)
					bcoeff = 1;
				else
					bcoeff = -1;
				bcoeff *= power(expe[i], m) ;
				bcoeff *= (RationalNTL(expa[i], 1)).power(-1 - m);
				exponents[0] = m;
				p.setToConstant(bcoeff);
				oneExpansion.insertMonomial(p, exponents);
			}
			oneExpansion.check();
			a.multiply(oneExpansion, minE, maxE);
		}//expa !=0 and expe != 0. Then 1/ax = 1/x * (1/ (expa + expe*e) = sum_{m=0}^{inf} (-1)^m (expe*e)^m * expa^{-1-m}
	}

	a.check();

	assert( powerX == (int)expa.size() && powerE == numPoles);


}

/**
 * Expand exp( \sum (ai + ei*epsilon)*x*{fi} )
 */
void TopKnapsack::expandExponentialPart(GeneralMonomialSum<PeriodicFunction, int> & exExpansion, const int numPoles, const vec_ZZ &a, const vec_ZZ & e, const vector<RationalNTL> & f)
{
	exExpansion.varCount = 2;

	int exponents[2];

	//cout << "going to expand exp(";
	//for(int i = 0; i < a.length(); ++i)
	//	cout << a[i] << "*{" << f[i] << "}+";
	//cout << ")^x \n and exp(";
	//for(int i = 0; i < a.length(); ++i)
	//	cout << e[i] << "*{" << f[i] << "}+";
	//cout << ")^xe" << endl;


	//goal: set pa = \sum ai*{fi}  and pe = \sum ei*{fi} and expand exp(pa*x)*exp(pe*epsilon*x)
	PeriodicFunction pa, pe;

	for(int i = 0; i < a.length(); ++i)
	{
		if ( f[i] == 0)
			continue;
		if ( a[i] != 0)
		{
			PeriodicFunction temp1(f[i], false);
			temp1.times(RationalNTL(a[i],1));
			pa.add(temp1);
		}
		if ( e[i] != 0)
		{
			PeriodicFunction temp2(f[i], false);
			temp2.times(RationalNTL(e[i],1));
			pe.add(temp2);
		}
	}


	//first expand exp(pa*x) = \sum (pa)^m*x^m/m!
	ZZ mFract;
	mFract = 1;
	exponents[0] = 0;
	for(int m = 0; m <= order; ++m)
	{
		exponents[1] = m;
		PeriodicFunction newpa(pa);
		newpa.pow(m);
		//cout << "mfract=" << mFract << endl;
		newpa.div(mFract);

		exExpansion.insertMonomial(newpa, exponents);
		//cout << "inserted " << newpa << "* e^" << exponents[0] << "* x^" << exponents[1] << " + " << endl;
		mFract *= (m+1);
	}

	//next expand exp(pe*x*e) = \sum (pe)^m*epsilon^mx^m/m!
	if ( numPoles == 0 || order == 0)
		return;


	GeneralMonomialSum<PeriodicFunction, int> exExpansion2;
	exExpansion2.varCount = 2;
	mFract = 1;
	for(int m = 0; m <= min(order, numPoles); ++m)
	{
		exponents[1] = m;
		exponents[0] = m;
		PeriodicFunction newpe(pe);
		newpe.pow(m);
		newpe.div(mFract);

		exExpansion2.insertMonomial(newpe, exponents);
		mFract *= (m+1);
	}
	int maxE[2] = {0,0};
	int minE[2] = {INT_MIN, INT_MIN};
	maxE[0] = numPoles;
	maxE[1] = order;
	exExpansion.multiply(exExpansion2, minE, maxE);
}


void TopKnapsack::expandF1Case(GeneralMonomialSum<PeriodicFunction, int> & expansion)
{
	vector<ZZ> alphaCopy;

	alphaCopy.resize(alpha.length());
	for(int i = 0; i < (int) alpha.length(); ++i)
		alphaCopy[i] = alpha[i];


	expandNonperiodicPart(expansion, alphaCopy);
	
	//cout << "f=1 series::" << expansion.printMonomials().c_str() << endl;
	
	ZZ bottomCoeffNonperiodicPart;
	bottomCoeffNonperiodicPart = 1;
	for(int i = 0; i < (int)alphaCopy.size();++i)
		bottomCoeffNonperiodicPart *= alphaCopy[i];


	if( (N+1) % 2 == 1)
		bottomCoeffNonperiodicPart *= -1;
	PeriodicFunction p;
	p.setToConstant(RationalNTL(1, bottomCoeffNonperiodicPart));

	int exponent[2];
	exponent[0] = 0;
	//exponent[1] = -1*(N+1);
	exponent[1] = 0;
	GeneralMonomialSum<PeriodicFunction, int> scaleTerm;
	scaleTerm.varCount = 2;
	scaleTerm.insertMonomial(p, exponent);

	int maxE[2], minE[2];
	maxE[0] = 0;
	minE[0] = 0;

	//maxE[1] = -1*N -1 + order;
	//if ( topKTerms )
	//	minE[1] = -1*N -1;
	//else
	//	minE[1] = -1*N -1 + order;

	//cout << "f=1 series w/o scale" << expansion.printMonomials().c_str() << endl;

	maxE[1] =  order;
	if ( topKTerms )
		minE[1] = 0;
	else
		minE[1] =  order;


	expansion.multiply(scaleTerm, minE, maxE);


}

/*
 * choose(-n,k) = choose(n+k-1,k)*(-1)^k
 *
 * 1/(1-exp(ax +bxe)) = -1/(ax+bxe) \sum_{m=0} B_m (ax+bxe)^m/m!
 *      = -1/(ax+bxe)(B0 + B1 (ax+bxe) + B2(ax+bxe)^2/2! + B3(ax+bxe)^3/3! + b4(ax+bxe)^4/4! +...)
 *      = -1( B0/(ax+bxe) + B1 + B2(ax+bxe)^1/2! + B3(ax+bxe)^2/3! + B4(ax+bxe)^3/4!+...)
 * If a,b!= 0
 *      = -1( (sum_{p=0}  B0(bxe)^p(ax)^{-1-p}choose(-1,p))  + B1 + B2(ax+bxe)^1/2! + B3(ax+bxe)^2/3! + B4(ax+bxe)^3/4!+...)
 *      = -1( (sum_{p=0}  B0(bxe)^p(ax)^{-1-p}(-1)^p)  + B1 + B2(ax+bxe)^1/2! + B3(ax+bxe)^2/3! + B4(ax+bxe)^3/4!+...)
 *      expand to order
 */
/*
void TopKnapsack::expandOneTerm(monomialSum & oneExpansion, const ZZ & a, const ZZ & b)
{
	//[0] = power of x, [1]=power of e
	int monomialExp[2] = {0,0};


	RationalNTL temp;
	ZZ apow, bpow;
	//first insert (sum_{p=0}^{numPoles}  -B0(bxe)^p(ax)^{-1-p}(-1)^p)
	if( a == 0)
	{
		monomialExp[1] = monomialExp[0] = -1;
		insertMonomial(RationalNTL(to_ZZ(-1), b), monomialExp, oneExpansion);
	}//-1B0/(ax+bxe) = -B0/(0+bxe)
	else if ( b == 0 )
	{
		monomialExp[0]=-1;
		monomialExp[1] = 0;
		insertMonomial(RationalNTL(to_ZZ(-1),a), monomialExp, oneExpansion);
	}//-1B0/(ax+bxe) = -B0/(ax+0)
	else
	{
		monomialExp[0] = -1;
		apow = a;
		bpow = 1;
		for(int p = 0; p <= order; ++p)
		{
			temp = bpow;
			temp.div(apow);
			if( p %2 == 0)
				temp.changeSign();
			monomialExp[1] = p;

			insertMonomial(temp, monomialExp, oneExpansion);
			apow *= a;
			bpow *= b;
		}
	}//-1( B0/(ax+bxe) = sum_{p=0}  -B0(bxe)^p(ax)^{-1-p}(-1)^p) =sum_{p=0} -B0(b^p)(e^p)/(x*a^{p+1})(-1)^p)

	//next: -1( + B1 + B2(ax+bxe)^1/2! + B3(ax+bxe)^2/3! + B4(ax+bxe)^3/4!+...)
	monomialExp[0] = monomialExp[1] = 0;
	temp = bernoulli[1];
	temp.changeSign();
	insertMonomial(temp, monomialExp, oneExpansion);

	//next: -(B2(ax+bxe)^1/2! + B3(ax+bxe)^2/3! + B4(ax+bxe)^3/4!+...)
	ZZ fract;
	fract = 1;
	bpow = 1;
	for(int p = 1; p <= order; ++p)
	{
		fract*= (p+1);
		bpow *= b;
		ZZ bpowTemp = bpow; //bpowTemp = b^p
		apow = 1;
		monomialExp[0] = p;
		if ( a == 0)
		{
			temp = bernoulli[p+1];
			temp.div(fract);
			temp.changeSign();
			temp.mult(bpowTemp);
			monomialExp[1] = p;
			insertMonomial(temp, monomialExp, oneExpansion);
		}
		else if ( b == 0)
		{
			temp = bernoulli[p+1];
			temp.div(fract);
			temp.changeSign();
			temp.mult(power(a,p));
			monomialExp[1] = 0;
			insertMonomial(temp, monomialExp, oneExpansion);
		}
		else
		{
			for(int pow = 0; pow <= p; ++pow)
			{
				//-B[p+1](ax+bxe)^p/(p+1)! = -B[p+1]/(p+1)! ( sum (ax)^pow(bxe)^{p-pow} (p choose pow))
				monomialExp[1] = p-pow;
				temp = bernoulli[p+1];
				temp.div(fract);
				temp.changeSign();
				temp.mult(apow*bpowTemp);
				temp.mult(binomial(p,pow));

				insertMonomial(temp, monomialExp, oneExpansion);

				apow *= a;
				bpowTemp /= b;

			}//for pow
		}//else

	}//for p.

}
*/
/**
 * Expand 1/(1-exp(ax)) up to order 'order'
 * 1/(1-exp(ax)) = -1/(ax) \sum_{m=0} B_m (ax)^m/m!
 *               = \sum_{m=0} -B_m (ax)^{m-1}/m!
 */
/*
void TopKnapsack::expandOneNonPerturbedTerm(monomialSum & oneExpansion, const ZZ & a)
{
	int monomialExp = -1;
	insertMonomial(RationalNTL(to_ZZ(-1), a), &monomialExp, oneExpansion);

	ZZ fractorial;
	fractorial = 1;
	RationalNTL temp;
	RationalNTL c;
	c = 1;
	for(int m = 1; m <= order; ++m)
	{

		monomialExp = m-1; //power
		temp = c*bernoulli[m]; //a^{m-1}*B[m]
		fractorial*=to_ZZ(m); // m!
		temp.div(fractorial);
		temp.changeSign(); //temp = -B[m](a)^{m-1}/m!
		insertMonomial(temp, &monomialExp, oneExpansion);
		c *= a;
	}
}
*/
/**
 * @parm xProduct, output after we take the coefficient of e=0.
 * @parm xeProduct, current product of polynomials in x and e
 * @parm B, columns are the basis for the lattice
 * @parm invB B^{-1}
 * @parm current cone. It is assumed its facets and rays are computed.
 */
/*
void TopKnapsack::expandPeriodicExponential(GeneralMonomialSum<PeriodicFunction, int> & pxeProduct,
		const mat_ZZ &B, const mat_ZZ & invB, const ZZ & invBd,
		const listCone * oneCone, const vec_ZZ & s, const vec_ZZ & alpha, const vec_ZZ & beta)
{
	//convert the rays to a matrix.
	int I = alpha.length();
	mat_ZZ raysT;
	raysT.SetDims(I, I);
	int counter= 0;
	for(listVector * ray = oneCone->rays; ray; ray = ray->rest)
	{
		raysT[counter] = ray->first ;
		counter++;
	}

	mat_ZZ BT; //B^T
	BT.SetDims(B.NumRows(), B.NumCols());
	transpose(BT, B);
	//x<av, B[cones]{s}> = <[cones]^TB^Ta, {s}>, likewise for e.
	vec_ZZ av, bv;
	av.SetLength(I);
	bv.SetLength(I);
	av = raysT * (BT*alpha);
	bv = raysT * (BT*beta);

	//now we compute {s} for  < av x + bv xe, {s}>.
	//{s} =  fractionpart( [cones]^-1 * invB s)
	//cones^-1 has already been computed...
	// it is in the facet information of listCone in the form -[cones]^-T
	// however, thoese facets have been scaled.
	// need to unscale them...requires more book keeping or conversion to a rational matrix.
	// todo: do this, for now just recompute the inverse.
	mat_ZZ raysI;
	ZZ d;
	raysI.SetDims(raysT.NumRows(), raysT.NumCols());
	//BT is no longer needed, so I will reuse it.
	transpose(BT, raysT);
	inv(d, raysI, BT); // (raysI/d)* rays = I
	vec_ZZ fractionalS  = raysI*(invB*s);


	//ok, finally have enough data to make a and b.
	//a = <av, {fractionalS/d}>, likewise for b.
	PeriodicFunction a, b;
	RationalNTL f;
	for(int i = 0; i < I; ++i)
	{
		f = fractionalS[i];
		f.div(d*invBd);
		f.changeSign();
		//a.addProduct(RationalNTL(av[i],1), f); //adds the term av[i]*{f}
		//b.addProduct(RationalNTL(bv[i],1), f);
		//cout << "inserted xc" << av[i] << " f=" << f << endl;
		//cout << "inserted ec" << bv[i] << " f=" << f << endl;
	}

	cout << "av=" << av << endl;
	cout << "bv=" << bv << endl;
	cout << "fs=" << fractionalS << endl;
	cout << "d" << d << endl;
	cout << "invb" << invBd << endl;


	//cout << "coeff of x is" << endl;
	//cout << a << endl;
	//cout << "coeff of e is" << endl;
	//cout << b << endl;
	//cout << "end." << endl;

	//insert the expansion of exp(ax+bxe)
	// = sum_{m=0}^order (ax+bxe)^m/m!
	pxeProduct.varCount = 2;


	PeriodicFunction ftemp, ftemp2; //temp function.
	RationalNTL ctemp; //temp coefficient.
	ZZ fract;
	int monomialExp[2] = {0,0};
	fract = 1;


	ftemp.setToConstant(1);
	pxeProduct.insertMonomial(ftemp, monomialExp);//insert 1.

	cout << "coing to insert exp(ax+bxe)\n where a=" << a << "\n and b=" << b << endl;
	//sum_{p=0}^order (ax+bxe)^p/p!
	//= sum_{p=0}^order (1/p!) sum (ax)^pow(bxe)^{p-pow}(p choose pow)
	for(int p = 1; p <= order; ++p)
	{
		fract*= p;
		monomialExp[0] = p;
		if ( a == 0)
		{
			ftemp = b;
			ftemp.pow(p);
			ftemp.div(fract);
			monomialExp[1] = p;
			pxeProduct.insertMonomial(ftemp, monomialExp);
			cout << "inserted " << ftemp << endl;
		}
		else if ( b == 0)
		{
			ftemp = a;
			ftemp.pow(p);
			ftemp.div(fract);
			monomialExp[1] = 0;
			pxeProduct.insertMonomial(ftemp, monomialExp);
			cout << "inserted " << ftemp << endl;
		}
		else
		{ //(1/p!) sum (ax)^pow(bxe)^{p-pow}(p choose pow)
			for(int pow = 0; pow <= p; ++pow)
			{

				monomialExp[1] = p-pow;
				ftemp = a;
				ftemp.pow(p);
				ftemp2  = b;
				ftemp2.pow(p-pow);
				ftemp.times(ftemp2);
				ftemp.times(RationalNTL(binomial(p,pow), fract));//ftemp = (p choose pow)/p! a^pow b^{p-pow}

				pxeProduct.insertMonomial(ftemp, monomialExp);
				cout << "inserted " << ftemp << endl;
			}//for pow
		}//else

	}//for p.

}
*/
/*

void TopKnapsack::expandPerturbationTerms(GeneralMonomialSum<PeriodicFunction, int> & ans,
		const GeneralMonomialSum<PeriodicFunction, int> & p1, const monomialSum & p2)
{
	BTrieIterator<PeriodicFunction, int>* it1 = new BTrieIterator<PeriodicFunction, int>();
	BTrieIterator<RationalNTL, int>* it2 = new BTrieIterator<RationalNTL, int>();

		it1->setTrie(p1.myMonomials, p1.varCount);
		it2->setTrie(p2.myMonomials, p2.varCount);
		int low[2], high[2];
		low[0] = INT_MIN;
		low[1] = 0;
		high[0] = order;
		high[1] = 0;


		ans.myMonomials = new BurstTrie<PeriodicFunction, int>();
		ans.varCount = 1;
		int power[1];

		term<PeriodicFunction, int> *periodicTerm;
		term<RationalNTL, int> *rationalTerm;

		it1->begin();
		it2->begin();

		int i;
		while (periodicTerm = it1->nextTerm())
		{
			while (rationalTerm = it2->nextTerm())
			{
				//[0] = power of x, [1]=power of e
				int xPower = periodicTerm->exps[0] + rationalTerm->exps[0];
				int ePower = periodicTerm->exps[1] + rationalTerm->exps[1];

				if ( ePower != 0)
					continue;
				if (xPower > order)
					continue;

				power[0] = xPower;
				PeriodicFunction temp;
				temp = periodicTerm->coef;
				temp.times(rationalTerm->coef);

				ans.myMonomials->insertTerm(temp, power, 0, ans.varCount, -1);

			}
			it2->begin();
		}

}

*/

//This is a static function
ZZ TopKnapsack::binomial(int n, int k)
{
	if ( k == n || k == 0)
		return to_ZZ(1);
	if ( n/k < 0.5)
		return binomial(n, n-k);

	//n !/ k! (n-k)! = n*...(n-k+1)/k!
	ZZ num, denum;
	num = denum = 1;
	for(int i = n; i >= n-k+1; --i)
		num *= i;
	for(int i = 1; i <= k; ++i)
		denum *= i;
	return num/denum;
}

/**
 * @parm k: integer 1 <= k <= alphaSize
 * Computes every k-subset of alpha and computes the gcd of these sublists.
 */
void TopKnapsack::everyGCD(int k)
{
	if (k == N+1)
	{
		gcds.insertGCD(to_ZZ(1));
		return;
	}

	//sublist is going to be a list of index values to add to the output list.
	//sublist starts off as [k,k-1, ..., 2, 1].
	//We keep adding values to sublist[1] and if sublist[i] > (n - i) we add 1 to the next element
	//Hence the last element added is [n, n-1, n-2, n-3, ..., n-k+1]
	int* sublist = new int[k];
	int i;
	ZZ newGCD;
	for(i = 0; i < k; ++i)
		sublist[i] = k-i;

	ZZ limit, counter;
	limit = binomial(N+1, k);
	counter = 1;
	while (counter <= limit)
	{

		newGCD = alpha[sublist[0]-1];
		for(i = 1; i < k; ++i)
			newGCD = GCD(newGCD, alpha[sublist[i]-1]);
		gcds.insertGCD(newGCD);
		/*
		cout << "("<< counter << "/" << limit << ") List:";
		for(i = 0; i < k; ++i)
			cout << sublist[i] << ", ";
		cout << " gcd=" << newGCD << endl;
		*/

		if (limit == counter)
			break;

		sublist[0] += 1;

		i=0;
		while (sublist[i] > (N+1  - i) )
		{
			sublist[i+1] += 1;
			++i;
		}
		while (i > 0)
		{
			sublist[i-1] = sublist[i]+1;
			--i;
		}
		counter += 1;

	}//while
	delete [] sublist;
}//everyGCD

void TopKnapsack::printMatrix(const mat_ZZ &A)
{
	long int c = A.NumCols();
	long int r = A.NumRows();

	for(long int i = 0; i < r; ++i)
	{
		for(long int j = 0; j < c; ++j)
			cout << A[i][j] << ", ";
		cout << endl;
	}

}

