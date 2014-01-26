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
MobiusList::~MobiusList() {}

/**
 * @param v is inserted at the end of the list if it is not already there.
 */
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

/**
 * Sets things up to compute the mu values. Throws an exception if the gcd of the entire list is not 1.
 */
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


/**
 * Recursively computes the mu values.
 * Formula: mu(f) = 1 - \sum_{j: f divides j} mu(j).
 * @param i, index of the current mu value that will be computed. On first call, this index value must be the gcd=1 case.
 */
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

/**
 * Used for debugging.
 */
void MobiusList::print() const
{
	for(int i = 0; i < (int) list.size(); ++i)
		cout << list[i].mu << ", gcd=" << list[i].gcd << ", isv=" << list[i].mobiusValid << endl;
}

/**
 * Computes the mu values and initializes the unweightedSeries vector.
 */
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
		if (unweightedSeries[i])
			delete unweightedSeries[i];
}

// **************************************************************************************************
// **************************************************************************************************

/*
 * Bernoulli of the 2nd kind:
 * Computes B[0], B[1], .... B[k] where
 * t/(1-exp(-t)) = \sum_{m=0}^infinity B[m] (-t)^m/m!
 * See https://en.wikipedia.org/wiki/Bernoulli_number
 *
 * Formula:
 * B[0]= 1;
 * B[m] = 1 - \sum_{j=0}^{m-1} choose(m,j)B[j]/(m-j+1)
 * See https://en.wikipedia.org/wiki/Bernoulli_number#Explicit_definition
 *
 * Then note that Bernoulli of the fist kind is equal to the second kind only b[1] = -1/2, instead of 1/2
 *
 * @param k starting with B[0], computes up to B[k], where the B's are the Bernoulli of the 1st kind:
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
		B[m] = 0;
		if (m > 1 && (m % 2) == 1)
		{
			continue; //B[m] = 0 for odd m larger than 1.
		}
		for( j = 0; j < m; ++j)
		{
			B[m] += B[j]*RationalNTL(TopKnapsack::binomial(m,j), to_ZZ(m-j+1));;
		}
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

TopKnapsack::TopKnapsack():computeGCDbySubsets(false) {
	
}

TopKnapsack::~TopKnapsack() {
	// TODO Auto-generated destructor stub
}

/**
 * Sets the knapsack list and computes the bernoulli numbers.
 * @param list the knapsack list
 */
void TopKnapsack::set(const vec_ZZ& list)
{
	alpha = list;
	N = alpha.length()-1;
	bernoulli.setBernoulli(alpha.length());
}

/**
 * Sets how to compute the gcds.
 * @param t: if true, uses the every-subsets method, else uses the complete method.
 */
void TopKnapsack::useSubsetsForGCD(bool t)
{
	computeGCDbySubsets = t;
}

/**
 * Sets the seed for rand().
 */
void TopKnapsack::seed(int s)
{
	if (s >= 0)
		srand(s);
	else
		srand(time(0));
}

/**
 * @param k will compute the coefficients of t^N, t^{N-1}, ..., t^{N-k}
 */
void TopKnapsack::coeff_topK(int k)
{
	topKTerms = true;
	coeff(k);
}

/**
 * @param k will compute the coefficient of t^{N-k}
 */
void TopKnapsack::coeff_NminusK(int k)
{
	topKTerms = false;
	coeff(k);
}

/**
 * Start of main computation.
 * Computes the gcd, the mu values, does the residues
 */
void TopKnapsack::coeff(int k)
{

	assert(0 <= k && k<= N);
	order = k;
	cout << "order=" << order << endl;
	coeffsNminusk.resize(k+1);


	Timer tgcd("Time for gcds");
	tgcd.start();
	findGCDs(k);
	gcds.computeMobius();
	tgcd.stop();
	cout << tgcd << endl;

	cout << "mu found" << endl;
	gcds.print();


	for(int i = 0; i < (int)gcds.list.size(); ++i)
		if (gcds.list[i].mu != 0)
		{
			E(i);
		}

	packageAnswer();

}

/**
 * 1) Looks at every term in the series computed so far,
 * 2) figure how for which t^{N-i}-coefficient this term contributes to,
 * 3) compute the final residue value,
 * 4) Factor in f and mu(f)
 */
void TopKnapsack::packageAnswer()
{
	for(int i = 0; i < (int)gcds.list.size(); ++i)
		if (gcds.list[i].mu != 0)
		{
			//series might contain terms that contribute to t^N, .., t^{N-k} or just t^{N-k}
			GeneralMonomialSum<PeriodicFunction, int> *series = gcds.unweightedSeries[i];
			if ( series->termCount == 0)
				continue;

			//loop over every term in the series
			BTrieIterator<PeriodicFunction, int> * itr = new BTrieIterator<PeriodicFunction, int>();
			itr->setTrie(series->myMonomials, series->varCount);
			itr->begin();
			term<PeriodicFunction, int> * oneTerm;
			while( (oneTerm = itr->nextTerm()) )
			{

				PeriodicFunction p(oneTerm->coef);
				//ignore exps[0] as this is the power of epsilon which does not have meaning.
				int deg = oneTerm->exps[1]; // deg = i for some 0<=i<=k (=order)
				//recall that 1/(x^{N+1}}* gcds.unweightedSeries gives you the real series expansion.
				//and so deg = -N -1 +i in the true series expansion.
				//Then this term would be obtained when we computed the residue of ((-x)^{N-i}/(N-i)! * real series expansion)
				//Hence this term contributes to t^{N-i}
				int h = N-deg; //h= N-i

				ZZ hFactorial;
				hFactorial = 1;
				for(int i = 2; i <= h; ++i)
					hFactorial *= i;

				RationalNTL factor;
				if ( h % 2 == 0)
					factor = -1;
				else
					factor = 1;

				factor *= gcds.list[i].mu;
				factor *= gcds.list[i].gcd;
				factor.div(hFactorial);


				p.times(factor);
				coeffsNminusk[deg].add(p); //update coeff of t^{N-deg}
			}
			delete itr;
		}//if

}

/**
 * Prints the answer in a maple-friendly way.
 */
void TopKnapsack::printAnswer(ostream & out)
{

	if ( topKTerms == false)
		out << "coeff" << N << "minus" << order << ":= " << coeffsNminusk[order] << ";\n"; //save the answer in a maple var
	else
	{
		//print each coeff. computed in its own variable.
		for(int i = 0; i < (int) coeffsNminusk.size(); ++i)
		{
			out << "coeff" << N << "minus" << i << ":= " << coeffsNminusk[i] << ";\n";
		}

		//sum the polynomial
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

/**
 * @param k: gets the coeff of t^{N-k}. Assumes coeff_NminusK or coeff_topK was already called
 */
PeriodicFunction TopKnapsack::getCoeffNminusK(int k)
{
	return coeffsNminusk[k];
}

/**
 * Going to compute the expansion of F(\alpha, f, T)(x).
 */
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


	vector<ZZ> fDivAlpha, fnDivAlpha; //f (not) divides alpha
	for(int i = 0; i < alpha.length(); ++i)
	{
		if ( divide(alpha[i], f))
			fDivAlpha.push_back(alpha[i]);
		else
			fnDivAlpha.push_back(alpha[i]);
	}


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

	listCone* uniCones = findUnimodularCones(invLatticeBasisScaled); //todp: move the latice scaling to this function

	bool finishedResidue = false;
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
	
	freeListCone(uniCones);

}


//cols of latticeBasis generate a cone of all x s.t. \sum \alpha_i x_i \in f\Z
/**
 * Find all solutions to Ax:=[fnDivAlpha]x = f\Z where A is the row vector [fnDivAlpha].
 *
 * Then x^TA^T = f.
 * Let y^TUA^T = y^T H = f. H is a col vector with H=(h1, 0,0,....0)^T.
 *    where x^T = y^TU => u^T y = x.
 * If y1 = f/h1, and the other y's anything (say y2 = ... =yn =0)  then Ax=f, but x may not be integral.
 * So we try Ax = 2f, giving y1=2f/h1. Etc.
 *
 * Another way to view this, is to look at U*A^T = H = (h1, 0,0,....0)^T. Then multiplying
 * the first row on both sides by f/gcd(f,h1), we get
 * (U with scaled 1st row)*A^T = (lcm(h1,f), 0,0,....0)^T which gives the smallest positive f\Z we can get.
 *
 * Hence x is in the from U^T y where y1 = f/gcd(f,h1), y2, .., yn free. for \sum alpha * x_i = smallest element in f\Z
 *  Then x is in the form U^T y where y1 = f/gcd(f,h1)z, with z, y2, ..., y free for \sum alpha * x_i in f\Z
 * Hence x is in the from x = (first for of U scaled)^T y were y is free.
 * So the cols of (first rows of U scaled)^T give the basis we want!!!
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

	//U is in row-major order. take the transpose and save that in latticeBasis
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
	divide(newf, f, GCD(f,s[0])); //newf = f/gcd(f, u[0][0]);
	for(int i = 0; i < n; ++i)
		latticeBasis[i][0] *= newf;


}


/**
 * Let L = latticeBasis. We want to write the standard orthant (which is a cone(e1, .., en)) in the L basis.
 * That is, we need vectors c1, c2, .., cn s.t Lci = ei. Hence if C:=[c1, ..., cn], then L*C = I, so C = L^-1.
 * But the problem is that we want integer columns ci (or thinking about the ci as rays, integer rays), so we compute L^-1, but scale each column to integers.
 * Because the columns of C are rays of a cone in a new basis L, we can divide each ci by its gcd, which preserves the cone they span.
 * In summary, we write the cone generated by the standard basis vectors in the L basis while making sure the rays are primitive.
 * @param latticeBasis: if x = latticeBasis*y for any y, then \sum alapha_i x_i \in f\Z
 * @param invLatticeBasisD = determinant(latticeBasis)
 * @param invLatticeBasis matrix such that invLatticeBasisScaled * latticeBasis = determinant(latticeBasis) * I
 * @param invLatticeBasisScaled = invLatticeBasis but each column scaled to have gcd 1
 */
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
 * Find s1, s2, ..., sI, s_{I+1}, such that  s1 a_1 + ... s_I a_I + s_{I+1}f = 1
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

/**
 * @param invLatticeBasis: the columns contains rays for a cone.
 * @return unimodular triangulation of the cone invLatticeBasis
 */
listCone* TopKnapsack::findUnimodularCones(const mat_ZZ & invLatticeBasis)const
{
	listCone * newCone = createListCone();

	//make the rays be in the rows.
	mat_ZZ rays = transpose(invLatticeBasis);


	//build the cone
	newCone->vertex = new Vertex;
	newCone->vertex->vertex = new rationalVector(invLatticeBasis.NumRows()); //make a zero vertex.

	newCone->rays = new listVector(rays[0]);

	for(int i = 1; i < rays.NumRows(); ++i)
		newCone->rays = appendVectorToListVector(rays[i], newCone->rays);


	BarvinokParameters params;
	params.Number_of_Variables = rays.NumRows();
	params.max_determinant = 1;
	listCone *uniTriangulation = decomposeCones(newCone, true, params);

	freeListCone(newCone);

	return uniTriangulation;
}

/*
 * Let xv:=(alpha + beta*e)x
 * 		where alpha is the vector of alpha's that f does not divide
 * 		      beta is a is a random vector
 * 		      e is epsilon and symbolic
 * 		      x is symbolic.
 *
 *  For every unimodular cone  u with rays g_i, want to compute the order^{th} term in x of
 *
 *  exp(<B^T xv, {-TB^{-1}s}_u>)                       (sign of cone u)
 *  ---------------------------- x --------------------------------------------------------
 *  \prod (1-exp(<B^Txv, g_i>))    \prod_{the alpha's that f does divide} (1-e^{alpah_i x})
 *
 * where {-TB^{-1}s}_u = \sum {-T*c_i} g_i
 * where c_i are the unique numbers such that B^{-1}s = \sum c_i g_i
 * Finally,  {-T*c_i} = fractionalPart(T*c_i).
 *
 * @param fSeries: output parameter. Is the series expansion of the above, but with an extra x^{N+1} scaling term.
 *        Also, the polynomial has two variables (epsilon, x), but the epsilon terms do not have meaning (there is no need to set their powers to zero)
 *        So fSeries could be the polynomial 1 + epsilon^3*x^3 + 2epsilon^3, but it should be interpreted as 1 + x^3 + 2
 *        This scaling issue is taken care of in packageAnswer()
 * @param s: the vertex of the cone
 * @param u: list of unimodular cones
 * @param B: colums are a basis for the lattice of all solutions for \sum \alpha_i x_i in f\Z
 * @param invB: matrix such that B*invB/invBd = I
 * @param invBd: det(B)
 * @param fnDivAlpha list of alpha's that f does Not divide
 * @param fDivAlpha list of alpha's  that f does divide
 */
void TopKnapsack::findResidue(GeneralMonomialSum<PeriodicFunction, int> & fSeries, const vec_ZZ & s, const listCone *u, const mat_ZZ & B,
		const mat_ZZ & invB, const ZZ & invBd,
		const vector<ZZ> &fnDivAlpha, const vector<ZZ> & fDivAlpha)
{
	//int seed = 654354; // for debugging. TODO: change to time(0);
	//srand(seed);

	//get random beta vector.
	int I = fnDivAlpha.size();
	int i;
	vec_ZZ beta;
	beta.SetLength(I);
	for(i = 0; i < I; ++i)
		beta[i] = ( rand()%2 ? -1*rand()%500 : rand()%500); //set beta to random vector in [-500,500]^I


	vec_ZZ temp;
	temp.SetLength(I);

	//copy fnDivide to a vec_ZZ
	vec_ZZ alpha;
	alpha.SetLength(I);
	for(i = 0; i < I; ++i)
		alpha[i] = fnDivAlpha[i];

	//goal: to compute <B^Txv, g_i> = < (alpha  + beta e)x, Bg_i>, save the dot products in expa, expe.
	vector<ZZ> expa, expe;
	expa.resize(I);// expa[i] = alpha * B*g_i *    x
	expe.resize(I);// expe[i] = beta * B*g_i * e * x

	//find the expansion of 1/\prod_{the alpha's that f does divide} (1-e^{alpah_i x})
	GeneralMonomialSum<PeriodicFunction, int> fDivAlphaExpansion;
	expandNonperiodicPart(fDivAlphaExpansion, fDivAlpha);
	ZZ bottomCoeffNonperiodicPart;
	bottomCoeffNonperiodicPart = 1;
	for(i = 0; i < (int)fDivAlpha.size();++i)
		bottomCoeffNonperiodicPart *= fDivAlpha[i];

	//cout << "non periodic part is (missing sign/X^"<< fDivAlpha.size() << "*" << bottomCoeffNonperiodicPart << ")* ==" << fDivAlphaExpansion.printMonomials().c_str() << endl;


	int minE[2] = {INT_MIN, INT_MIN}; //min and max exponents.
	int maxE[2] = {0, 0};


	for(const listCone * oneCone = u; oneCone; oneCone = oneCone->rest)
	{
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

		//now,  1/(1-exp(<(alpha  + beta e)*x, Bg>)) = 1/(1-exp(x*expa + x*e*expe))

		//cout << "expa: ";
		//for(i = 0; i < expa.size(); ++i)
		//	cout << expa[i] << ", ";
		//cout << "\nexpe: " ;
		//for(i = 0; i < expe.size(); ++i)
		//	cout << expe[i] << ", ";
		//cout << endl;



		//next, let fnDivAlphaExpansion be the expansion of 1/(1-exp(x*expa + x*e*expe)), up to a scaling.
		GeneralMonomialSum<PeriodicFunction, int> fnDivAlphaExpansion;
		ZZ bottomCoeffPeriodicPart;
		expandPeriodicPart(bottomCoeffPeriodicPart, fnDivAlphaExpansion, numPoles, expa, expe);




		//Next Goal:compute the numerator exp(<B^T (alpah + beta e)x, {-B^-1 T s}_u> where {-x} = fractionalPart(x)
		//sub-goal: find {-B^-1 T s}_u  = \sum {-T temp[i]} g[i] = \sum fractionalPart(T*temp[i])g[i] where B^{-1}s = \sum temp[i] * g[i]

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


		vector<RationalNTL> fractionalPart;
		fractionalPart.resize(I);
		for(int i = 0; i < I; ++i)
		{
			fractionalPart[i] = RationalNTL(temp[i], invBd*coneRaysInvD);
			if (fractionalPart[i].getDenominator() == 1)
				fractionalPart[i] = 0; //fractionalPart(integer) = 0
		}

		//now, fractionalPart[i] = {-B^-1 T s}_i
		vec_ZZ fractionalPartCoeffa, fractionalPartCoeffe;
		fractionalPartCoeffa.SetLength(I);
		fractionalPartCoeffe.SetLength(I);

		//sub-goal: compute < B^T(alpha + beta e)x, \sum fractionalPart[i] * g_i>
		i = 0;
		for(const listVector * g = oneCone->rays; g; g = g->rest)
		{
			temp  = B* (g->first); //B*g_i
			InnerProduct(fractionalPartCoeffa[i], alpha, temp);
			InnerProduct(fractionalPartCoeffe[i], beta, temp);
			++i;
		}

		//now, < B^T(alpha + beta e)x, \sum fractionalPart[i] * g_i> = \sum  ( fractionalPartCoeffa[i] + fractionalPartCoeffe[i]*e)*fractionalPart[i]*x

		//sub-goal: compute the expansion of exp(<B^T (alpah + beta e)x, {-B^-1 T s}_u>
		GeneralMonomialSum<PeriodicFunction, int> exExpansion; //"e to the x" expansion
		expandExponentialPart(exExpansion, numPoles, fractionalPartCoeffa, fractionalPartCoeffe, fractionalPart);

		//goal: multiply every series together and add in the scaling terms.
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

		scale *= oneCone->coefficient; //the sign of the unimodular cone


		//create a constant polynomial equal to scale.
		PeriodicFunction scaleTerm;
		scaleTerm.setToConstant(scale);
		GeneralMonomialSum<PeriodicFunction, int> finalProductScale;
		finalProductScale.varCount = 2;
		maxE[0] = 0;
		maxE[1] = 0;
		finalProductScale.insertMonomial(scaleTerm, maxE);


		maxE[0] = numPoles; //we never divided by 1/epsilon^numPoles, so we want to extract terms with epsilon^numPoles.
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
 * Goal: Expand 1/(\prod_{i: f does not divide \alpha_i} (1-exp(alaph_i*x)))
 * Note: 1/(1-exp(ax)) = -1/(ax) \sum_{m=0}^inf B_m (ax)^m/m! where B is the bernoulli numbers
 * We compute the product of \sum_{m=0}^inf B_m (ax)^m/m! with out the scaling term.
 * @param fDivAlpha: list of numbers
 * @param a: output polynomial. will be cleared before we start. Is equal to 1/(\prod_{i: f does not divide \alpha_i} (1-exp(alaph_i*x))) * alpha_i*x*-1
 */
void TopKnapsack::expandNonperiodicPart(GeneralMonomialSum<PeriodicFunction, int> &a, const vector<ZZ>  & fDivAlpha)
{

	a.varCount = 2;
	a.setToConstant(PeriodicFunction(RationalNTL(1,1), true));

	int min[2] = { INT_MIN, INT_MIN}; //does not really matter. could be 0, 0
	int max[2] = {0, INT_MAX};
	max[1] = order;

	ZZ mFract; //m!
	int exponents[2];
	exponents[0] = 0; //no epsilon term

	for(int j = 0; j <  (int)fDivAlpha.size(); ++j)
	{
		//oneExpansion = \sum_{m=0}^inf B_m (ax)^m/m!
		GeneralMonomialSum<PeriodicFunction, int> oneExpansion;
		oneExpansion.varCount = 2;
		mFract = 1;
		RationalNTL am(1, 1); //a^m
		for(int m = 0; m <= order; ++m)
		{
			RationalNTL coeff(am);
			coeff.div(mFract);
			coeff *= bernoulli[m];
			exponents[1] = m;

			PeriodicFunction pCoeff;  //convert coeff to a periodic function and save it.
			pCoeff.setToConstant(coeff);
			oneExpansion.insertMonomial(pCoeff, exponents);

			mFract *= (m+1);
			am *= fDivAlpha[j];
		}//insert one series expression


		a.multiply(oneExpansion, min, max);
	}//for j

}

/**
* Goal: find the expansion of \prod 1/(1-exp(x*expa[i] + x*e*expe[i]))
* The expansion will be equal to a/bottomCoeffPeriodicPart * 1/x^expa.size() * 1/e^numPoles
*
* Note: 1/(1-exp(ax)) = -1/(ax) \sum_{m=0}^inf B_m (ax)^m/m!
*
* @param bottomCoeffPeriodicPart. output parameter. scaling term.
* @param a. output parameter. the polynomial expansion with missing scaling term.
* @param numPoles is the number of zeros in expa
**/
void TopKnapsack::expandPeriodicPart(ZZ & bottomCoeffPeriodicPart, GeneralMonomialSum<PeriodicFunction, int> & a, const int numPoles, const vector<ZZ> & expa, const vector<ZZ> & expe)
{
	a.varCount = 2;
	a.setToConstant(PeriodicFunction(RationalNTL(1,1), true));
	bottomCoeffPeriodicPart = 1;

	RationalNTL coeff;
	int exponents[2];
	int minE[2] = {INT_MIN, INT_MIN};
	int maxE[2];
	maxE[0] = numPoles; //no need to expand in epsilon past the number of poles.
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
			//coeff = -B[m]/(m)! for m = 0, 1, ...

			exponents[1] = m; //power of x

			//expand coeff*(expa + expe*e)^m x^m via the binomial theorem
			//insert coeff * x^m ( \sum_{j=0}^{j=m} (m choose j) *(expe^j*e^j * expa^{m-j})  )
			for(int j = 0; j <= min(numPoles,m); ++j)
			{
				RationalNTL newCoeff(coeff);
				newCoeff.mult(TopKnapsack::binomial(m,j));
				newCoeff.mult( power(expe[i],j) * power(expa[i], m-j));
				exponents[0] = j; //power of e

				p.setToConstant(newCoeff);
				oneExpansion.insertMonomial(p, exponents);
			}
			mFract *= (m+1);
		}//for m.

		//oneExpansion.check();
		//a.check();
		a.multiply(oneExpansion, minE, maxE);
	}//for i

	//what is left: scale a by \prod 1/((expa[i] + expe[i]*e)*x)

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
			powerX++;
			oneExpansion.destroyMonomials();
			oneExpansion.varCount = 2;
			PeriodicFunction p;
			exponents[1] = 0; //power of x
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
			//oneExpansion.check();
			a.multiply(oneExpansion, minE, maxE);
		}//expa !=0 and expe != 0. Then 1/ax = 1/x * (1/ (expa + expe*e) = 1/x*sum_{m=0}^{inf} (-1)^m (expe*e)^m * expa^{-1-m}
	}

	//a.check();

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


/**
 * Computes the expansion of
 *
 *               x^{N+1}
 *    -------------------------
 *    prod (1 - exp(alpha_i x))
 *
 * @param expansion, an output parameter, assumed to be initialized to zero, else its value will be lost.
 */
void TopKnapsack::expandF1Case(GeneralMonomialSum<PeriodicFunction, int> & expansion)
{
	//we need a vector type to reuse the function expandNonperiodicPart()
	vector<ZZ> alphaCopy;
	alphaCopy.resize(alpha.length());
	for(int i = 0; i < (int) alpha.length(); ++i)
		alphaCopy[i] = alpha[i];


	expandNonperiodicPart(expansion, alphaCopy);
	
	//compute the denominator of the expansion
	ZZ bottomCoeffNonperiodicPart;
	bottomCoeffNonperiodicPart = 1;
	for(int i = 0; i < (int)alphaCopy.size();++i)
		bottomCoeffNonperiodicPart *= alphaCopy[i];

	if( (N+1) % 2 == 1)
		bottomCoeffNonperiodicPart *= -1;
	PeriodicFunction p;
	p.setToConstant(RationalNTL(1, bottomCoeffNonperiodicPart));

	//insert the fraction part first
	int exponent[2];
	exponent[0] = 0; //no epsilon term
	exponent[1] = 0; //insert as a constant
	GeneralMonomialSum<PeriodicFunction, int> scaleTerm;
	scaleTerm.varCount = 2;
	scaleTerm.insertMonomial(p, exponent);

	int maxE[2], minE[2];
	maxE[0] = 0; //no epsilon term
	minE[0] = 0;


	maxE[1] =  order;
	if ( topKTerms )
		minE[1] = 0;
	else
		minE[1] =  order; //ignore the other terms.

	//scale the answer
	expansion.multiply(scaleTerm, minE, maxE);
}




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
 * Finds all the gcds for computing the coeff of t^{N-k}
 */
void TopKnapsack::findGCDs(int k)
{
	cout << "computing gcd using a " << (computeGCDbySubsets ? "" : "non-") << "polynomial time algorithm" << endl;
	if ( computeGCDbySubsets )
		for(int i = 0; i <= k; ++i)
			everyGCDFromSubsets(N+1-i);
	else
		everyGCDFromEntireList(k);
}

/**
 * Finds all gcds from subsets of size N+1-k or larger
 * Algorithm: find the gcd of every subset by dynamic programming.
 */
void TopKnapsack::everyGCDFromEntireList(int k)
{
	vector<ZZ> output; //list of gcds from every subset.
	for(int i= 0; i < alpha.length(); ++i)
	{
		//at the end of this i-loop, output will contain every gcd from subsets of {alpha[0], ..., alpha[i]}
		//add gcd(alpha[i], any element of output) to output
		for(int j = 0; j < output.size(); ++j)
		{
			ZZ g;
			g = GCD(output[j], alpha[i]);
			if ( ! binary_search(output.begin(), output.end(), g))
			{
				output.push_back(g);
				for(int i = ((int) output.size()) -2; i>= 0 && output[i] > output[i+1]; --i)
				{
					ZZ t = output[i+1];
					output[i+1] = output[i];
					output[i] = t;
				}
					
			}//insert new gcd into sorted order
		
		}//for j

		//insert alaph[i] into the output list in sorded order if it is not already in the list.
		if ( ! binary_search(output.begin(), output.end(), alpha[i]))
		{
			output.push_back(alpha[i]);
			for(int i = ((int) output.size()) -2; i>= 0 && output[i] > output[i+1]; --i)
			{
				ZZ t = output[i+1];
				output[i+1] = output[i];
				output[i] = t;
			}		
		}//if
	}
	
	//output now contains the gcd from every subset of the alpha vector.
	//We only care about terms that divide N+1-k many.
	for (int i = 0; i < (int) output.size(); ++i)
	{
		int numDiv = 0;
		for(int j = 0; j < alpha.length(); ++j)
			if ( alpha[j] % output[i] == 0)
				++numDiv;
		if(numDiv >= N+1-k)
			gcds.insertGCD(output[i]);	
	}
	
	
}

/**
 * @parm k: integer 1 <= k <= alphaSize = N+1
 * Computes every k-subset of alpha and computes the gcd of these sublists.
 * Assumes the gcd of the entire list is one.
 */
void TopKnapsack::everyGCDFromSubsets(int k)
{
	if (k == N+1)
	{
		gcds.insertGCD(to_ZZ(1));
		return;
	}

	//sublist is going to be a list of index values to add to the output list.
	//sublist starts off as [k,k-1, ..., 2, 1].
	//We keep adding values to sublist[1] and if sublist[i] > (n - i) we add 1 to the next element
	//Hence the last element checked is [n, n-1, n-2, n-3, ..., n-k+1]
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
		//find the gcd of the current subset index by sublist.
		newGCD = alpha[sublist[0]-1];
		for(i = 1; i < k; ++i)
			newGCD = GCD(newGCD, alpha[sublist[i]-1]);
		gcds.insertGCD(newGCD);

		if (limit == counter)
			break;

		sublist[0] += 1;
		//ex: if N+1=10 and sublist = [11,9,2,1]
		i=0;
		while (sublist[i] > (N+1  - i) )
		{
			sublist[i+1] += 1;
			++i;
		}
		//ex then sublist = [11,10,3,1], i = 2 (0-based)
		while (i > 0)
		{
			sublist[i-1] = sublist[i]+1;
			--i;
		}
		//ex then sublist = [5,4,3,1]
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

