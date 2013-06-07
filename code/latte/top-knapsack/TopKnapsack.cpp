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
#include "print.cpp"
#include "dual.h"
#include <typeinfo>

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

// **************************************************************************************************
// **************************************************************************************************
/*
 * Computes B[0], B[1], .... B[k] where
 * t/(1-exp(-t)) = \sum_{m=0}^infinity B[m] (-t)^m/m!
 * See https://en.wikipedia.org/wiki/Bernoulli_number#Generating_function
 *
 * Formula:
 * B[0]= 1;
 * B[m] = 1 - \sum_{j=0}^{m-1} choose(m,j)B[j]/(m-j+1)
 * See https://en.wikipedia.org/wiki/Bernoulli_number#Explicit_definition
 */
void BernoulliSecondKind::setBernoulli(int k)
{
	B.resize(k+1);
	B[0] = 1;

	int m, j;
	for(m = 1; m <=k; m++)
	{

		//TODO: replace the choose function with pascal's triangle type computation to save time.
		B[m] = 0;
		for( j = 0; j < m; ++j)
		{
			B[m] += B[j]*RationalNTL(TopKnapsack::binomial(m,j), to_ZZ(m-j+1));;
		}
		B[m].changeSign();
		B[m] += RationalNTL(1,1);
	}

}


const RationalNTL& BernoulliSecondKind::operator[](int i) const
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

}
void TopKnapsack::findAll()
{
	for(int i = 0; i <= N; ++i)
		coeff_NminusK(i);
}



void TopKnapsack::coeff_NminusK(int k)
{
	order = k;


	for(int i = 0; i <= k; ++i)
		everyGCD(N+1-i);

	gcds.computeMobius();
	cout << "mu found" << endl;
	gcds.print();

	for(int i = 0; i < (int)gcds.list.size(); ++i)
		if (gcds.list[i].mu != 0)
			E(gcds.list[i].gcd);
}

void TopKnapsack::E(ZZ f)
{
	if (f == 1)
	{
		cout << "skippoing f=1 for now..." << endl;
		return;
	}
	cout << "f=" << f << endl;

	vector<ZZ> fDivAlpha, fnDivAlpha; //f (not) divides alpha
	for(int i = 0; i < gcds.list.size(); ++i)
	{
		if ( divide(gcds.list[i].gcd, f))
			fDivAlpha.push_back(gcds.list[i].gcd);
		else
			fnDivAlpha.push_back(gcds.list[i].gcd);
	}

	mat_ZZ latticeBasis;
	latticeBasis.SetDims(fnDivAlpha.size(), fnDivAlpha.size());
	findLatticeBasis(latticeBasis, fnDivAlpha, f);
	//cols of latticeBasis generate a cone of all x s.t. \sum \alpha_i x_i \in f\Z

	cout << "lattice basis " << endl;
	TopKnapsack::printMatrix(latticeBasis);

	mat_ZZ invLatticeBasis;
	invLatticeBasis.SetDims(fnDivAlpha.size(), fnDivAlpha.size());
	ZZ d;
	inv(d, invLatticeBasis, latticeBasis); // (1/d)*invLatticeBasis = latticeBasis^{-1}
	//for latticeBasis^{-1}, scale each col to integer, and divide by gcd. so we can ignore d.

	//scale each col.
	for(int j = 0; j < fnDivAlpha.size(); ++j)
	{
		d=invLatticeBasis[0][j];
		for(int i = 1; i < fnDivAlpha.size(); ++i)
			d = GCD(d, invLatticeBasis[i][j]);
		if ( d != 1)
			for(int i = 0; i < fnDivAlpha.size(); ++i)
				invLatticeBasis[i][j] /= d;
	}

	cout << "cone in basis lambda" << endl;
	TopKnapsack::printMatrix(invLatticeBasis);

	//find the vertex.
	vec_ZZ tVertex;
	tVertex.SetLength(fnDivAlpha.size());
	findVertex(tVertex, f, fnDivAlpha);

	cout << "tvertex=" << endl;
	for(int i = 0; i < tVertex.length(); ++i)
		cout << tVertex[i] << ", " ;
	cout << endl;

	listCone* uniCones = findUnimodularCones(invLatticeBasis);
	printListCone(uniCones, fnDivAlpha.size());

	findResidue(tVertex, uniCones, latticeBasis, invLatticeBasis, d, fnDivAlpha);

}

void TopKnapsack::findLatticeBasis(mat_ZZ &latticeBasis, const vector<ZZ> & fnDivAlpha, const ZZ & f) const
{
	vec_ZZ s, u, rhs;

	int n = fnDivAlpha.size();
	int m = 1;

	s.SetLength(n);
	for(int i = 0; i < n; ++i)
		s[i] = fnDivAlpha[i];


	u.SetLength(n*n);
	rhs.SetLength(m);

/*
	cout << "fnd=";
	for(int i = 0; i < fnDivAlpha.size(); ++i)
		cout <<fnDivAlpha[i] << ", ";
	cout << endl;
	cout << "calling ihermite with " << f << endl;
*/
	//For us,
	//	s=1 by n matrix = fnDivAlpha
	//	u = the col vectors give the basis for the lattice (before scaling)
	//	rhs = junk.
	int r = ihermite(&s, &u, &rhs, m,n);

	//convert U to matZZ in col-major order.
	for(int i = 0; i < n; ++i)
		for(int j = 0; j < n; ++j)
			latticeBasis[i][j] = u[j*n + i];
/*
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
	divide(newf, f, GCD(f,latticeBasis[0][0])); //newf = f/gcd(f, u[0][0])
	for(int i = 0; i < n; ++i)
		latticeBasis[i][0] *= newf;

}

void TopKnapsack::findVertex(vec_ZZ & tVertex, const ZZ &f, const vector<ZZ> &fnDivAlpha) const
{
	vec_ZZ s, u, rhs;

	cout << "the vertex, f="<< f << endl;
	for(int i = 0; i < fnDivAlpha.size(); ++i)
		cout << fnDivAlpha[i] << ",";
	cout << endl;

	int n = fnDivAlpha.size()+1;
	int m = 1;

	s.SetLength(n);
	for(int i = 0; i < n-1; ++i)
		s[i] = fnDivAlpha[i];
	s[n-1]=f;


	u.SetLength(n*n);
	rhs.SetLength(m);


	//For us,
	//	s=1 by n matrix = [fnDivAlpha, 1]
	//	u = the first col vector gives how to write \sum alpha_rx_r + x_{r+1}f = 1
	//     we want to set tVertex=[x_1, ..., x_r]
	//	rhs = junk.
	int r = ihermite(&s, &u, &rhs, m,n);

	//only need 1st column in U, and in this column, the last row is deleted
	//again, u is in col-major order.
	for(int i = 0; i < n-1; ++i)
		tVertex[i]= u[i];


}


listCone* TopKnapsack::findUnimodularCones(const mat_ZZ & invLatticeBasis)const
{
	listCone * newCone = createListCone();

	mat_ZZ rays = transpose(invLatticeBasis);
	cout << "type: " << typeid(rays[1]).name() << endl;

	newCone->vertex = new Vertex;
	newCone->vertex->vertex = new rationalVector(invLatticeBasis.NumRows()); //make a zero vertex.

	newCone->rays = new listVector(rays[0]);

	for(int i = 1; i < rays.NumRows(); ++i)
		newCone->rays = appendVectorToListVector(rays[i], newCone->rays);

	//cout << "newCone before unimodular" << endl;
	//printListCone(newCone, rays.NumRows());

	BarvinokParameters param;
	param.Number_of_Variables = rays.NumRows();
	param.max_determinant = 1;
	listCone *uniTriangulation = decomposeCones(newCone, true, param);

	//dualizeCones(uniTriangulation, param.Number_of_Variables, &param);
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
void TopKnapsack::findResidue(const vec_ZZ & s, const listCone *u, const mat_ZZ & B, const mat_ZZ & invB, const ZZ & invBd, const vector<ZZ> &fnDivAlpha)
{
	int seed = 0; // for debugging. TODO: change to time(0);
	srand(seed);

	//get random beta vector.
	int I = fnDivAlpha.size();
	int i;
	vec_ZZ beta;
	beta.SetLength(I);
	for(i = 0; i < I; ++i)
		beta[i] = ( rand()%2 ? -1*rand()%500 : rand()%500); //set beta to random vector in [-500,500]^I

	//copy fnDivide to a vec_ZZ
	vec_ZZ alpha;
	alpha.SetLength(I);
	for(i = 0; i < I; ++i)
		alpha[i] = fnDivAlpha[i];

	//store exp(< alpha x + beta e, Bg_i>)
	vector<ZZ> expx, expe;
	expx.resize(I);
	expe.resize(I);

	monomialSum ans;
	//try to compile and
	//start here


	for(const listCone * oneCone = u; oneCone; oneCone = oneCone->rest)
	{
		//compute all the inner produces in 1/(1-exp(<alpha x + beta e, Bg>))
		i=0;
		int numPoles = 0;
		for(const listVector * g = oneCone->rays; g; g = g->rest)
		{
			vec_ZZ Bg;
			Bg.SetLength(I);
			mul(Bg, B, (g->first));
			InnerProduct(expx[i], alpha, Bg);
			InnerProduct(expe[i], beta, Bg);

			if ( IsZero(expx[i]) && IsZero(expe[i]) )
			{
				cout << "beta not random enough";
				exit(1);
			}
			if( IsZero(expx[i]))
				numPoles++;
			i++;
		}//for g


		monomialSum xeProduct; //product of expansions in the x and e varables.
		xeProduct.varCount = 2;

		int low[2], high[2]; //polynomial bounds.
		low[0] = low[1] = INT_MIN;
		high[0] = high[1] = 0;
		insertMonomial(RationalNTL(1,1), high, xeProduct); //insert the polynomial "1"

		high[0] = high[1] = order; //update the upper bound.

		BTrieIterator<RationalNTL, int>* it = new BTrieIterator<RationalNTL, int>();
		BTrieIterator<RationalNTL, int>* it2 = new BTrieIterator<RationalNTL, int>();
		//now we do the expansion of 1/(1-exp{expx[i]*x + expe[i]*e})
		for(i = 0; i < I; ++i)
		{
			monomialSum oneExpansion;
			oneExpansion.varCount = 2;
			expandOneTerm(oneExpansion, expx[i], expe[i]);

			//mult oneExpansion with xeProduct.
			it->setTrie(xeProduct.myMonomials, xeProduct.varCount);
			it2->setTrie(oneExpansion.myMonomials, oneExpansion.varCount);

			monomialSum temp;
			temp.varCount = 2;

			multiply<RationalNTL>(it, it2, temp, low, high);

			//Delete memory and move the answer into xeProducts.
			destroyMonomials(oneExpansion);
			destroyMonomials(xeProduct);
			xeProduct = temp; //pointer copy of data.
		}//for i. for each term

		GeneralMonomialSum<PeriodicFunction, int> pxeProduct, pxProduct;
		expandPeriodicExponential(pxeProduct, B, invB, invBd, oneCone, s, alpha, beta);
		expandPerturbationTerms(pxProduct, pxeProduct, xeProduct);
		//now we do the expansion of exp{<alpha x + beta e, B*[cones]*{s}>) where {s} is a vector and =


	}//for oneCone
}


/*
 * choose(-n,k) = choose(n+k-1,k)*(-1)^k
 *
 * 1/(1-exp(ax +be)) = -1/(ax+be) \sum_{m=0} B_m (ax+be)^m/m!
 *      = -1/(ax+be)(B0 + B1 (ax+be) + B2(ax+be)^2/2! + B3(ax+be)^3/3! + b4(ax+be)^4/4! +...)
 *      = -1( B0/(ax+be) + B1 + B2(ax+be)^1/2! + B3(ax+be)^2/3! + B4(ax+be)^3/4!+...)
 * If a,b!= 0
 *      = -1( (sum_{p=0}  B0(be)^p(ax)^{-1-p}choose(-1,p))  + B1 + B2(ax+be)^1/2! + B3(ax+be)^2/3! + B4(ax+be)^3/4!+...)
 *      = -1( (sum_{p=0}^{numPoles}  B0(be)^p(ax)^{-1-p}(-1)^p)  + B1 + B2(ax+be)^1/2! + B3(ax+be)^2/3! + B4(ax+be)^3/4!+...)
 *      expand p up to numPoles.
 * If a!=0, b=0
 *      = \sum_{m=0} -B_m (a)^{m-1}(x)^{m-1}/m!
 *      expand x up to power order
 * If a=0, b!=0, like the last one.
 */
void TopKnapsack::expandOneTerm(monomialSum & oneExpansion, const ZZ & expx, const ZZ & expe)
{


	//this is an array of size one because this is the exponent "vector" using the BurstTrie.
	//[0] = power of x, [1]=power of y
	int monomialExp[2] = {0,0};
	//int mindeg[2] = {INT_MIN, INT_MIN};
	//int maxdeg[2] = {INT_MAX, INT_MAX};
	//maxdeg[1] = numPoles;

	RationalNTL c, temp;
	c=1;
	ZZ fractorial = to_ZZ(1);
	if ( IsZero(expx) && !IsZero(expe) )
	{

		monomialExp[1] = -1;
		insertMonomial(RationalNTL(to_ZZ(-1), expe), monomialExp, oneExpansion);

		for(int m = 1; m <= order; ++m)
		{

			monomialExp[1] = m-1; //power
			temp = c*bernoulli[m]; //b^{m-1}*B[m]
			fractorial*=to_ZZ(m); // m!
			temp.div(fractorial);
			temp.changeSign(); //temp = -B[m](b)^{m-1}/m!
			insertMonomial(temp, monomialExp, oneExpansion);
			c *= expe;
		}

	}// expansion of 1/(1-e^{ expe*e})
	else if ( !IsZero(expx) ) //&& IsZero(expe) )
	{
		monomialExp[0] = -1;
		insertMonomial(RationalNTL(to_ZZ(-1), expx), monomialExp, oneExpansion);

		for(int m = 1; m <= order; ++m)
		{

			monomialExp[0] = m-1; //power
			temp = c*bernoulli[m]; //b^{m-1}*B[m]
			fractorial*=to_ZZ(m); // m!
			temp.div(fractorial);
			temp.changeSign(); //temp = -B[m](b)^{m-1}/m!
			insertMonomial(temp, monomialExp, oneExpansion);
			c *= expx;
		}
	}//expansion of 1/(1-e^{expx*x})
	/*
	else
	{
		RationalNTL cy(1,1);
		//first insert (sum_{p=0}^{numPoles}  -B0(by)^p(ax)^{-1-p}(-1)^p)
		//when p=0:
		monomialExp[0] = -1;
		monomialExp[1] = 0;
		insertMonomial(RationalNTL(-1, expx), monomialExp, oneExpansion);

		c = RationalNTL(1,expx*expx);
		cy = expe;
		for(int p = 1; p <= numPoles; ++p)
		{

			monomialExp[0] = -1-p; //power
			monomialExp[1] = p;
			temp = c*cy;
			if( p % 2 == 0)
				temp.changeSign(); //temp = -B0(-1)^p a^p (b)^{-1-p}
			insertMonomial(temp, monomialExp, oneExpansion);

			//update the coefficients for the next round.
			c /= expx;
			cy *= expe;
		}

		//secondly, insert -(B1 + B2(ax+by)^1/2! + B3(ax+by)^2/3! + B4(ax+by)^3/4!+...)
		//insert -B1.
		monomialExp[1] = 0;
		monomialExp[0] = 0;
		temp = bernoulli[1];
		temp.changeSign();
		insertMonomial(temp, monomialExp, oneExpansion);

		fractorial = 1;
		//now the rest: -B2(ax+by)^1/2! + -B3(ax+by)^2/3! + -B4(ax+by)^3/4!+...)
		for(int m = 1; m <= order; ++m)
		{
			//insert -B[m+1](ax+by)^m/(m+1)! where the power of y is in [0, min(m, numPoles)]
			//insert -B[m+1]/(m+1)!*( sum_{k=0}^m (ay)^k(bx)^{m-k} choose(m,k) )

			//set c=(b)^(m+1);
			c = expx;
			for(int k = 1; k<= m; ++k)
				c *= expx;
			cy=RationalNTL(1,expe);
			fractorial *= (m+1);
			for(int k = 0; k <= m; ++k)
			{
				c/=expx; //c=b^{m-k}
				cy*=expe; //cy=a^k
				temp = c;
				temp.mult(cy);
				temp.mult(bernoulli[m+1]);
				temp.changeSign();
				temp.mult(binomial(m,k), fractorial); //temp = -B[m+1]/(m+1)! choose(m,k) a^k * b^{m-k}
				monomialExp[0] = m-k;
				monomialExp[1] = k;
				if ( k >= numPoles)
					break; //the largest y term needed is y^{numPoles-1}
				insertMonomial(temp, monomialExp, oneExpansion);

			}//for k

		}//for m. expand x to power order


	}//expansion of 1/(1-e^{expx*x + expe*e})
	*/
}

/**
 * @parm xProduct, output after we take the coefficient of e=0.
 * @parm xeProduct, current product of polynomials in x and e
 * @parm B, columns are the basis for the lattice
 * @parm invB B^{-1}
 * @parm current cone. It is assumed its facets and rays are computed.
 */
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

	//now we compute {s} for  < av x + bv e, {s}>.
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
		a.addProduct(RationalNTL(av[i],1), f); //adds the term av[i]*{f}
		b.addProduct(RationalNTL(bv[i],1), f);
		cout << "inserted xc" << av[i] << " f=" << f << endl;
		cout << "inserted ec" << bv[i] << " f=" << f << endl;
	}

	cout << "coeff of x is" << endl;
	cout << a << endl;
	cout << "coeff of e is" << endl;
	cout << b << endl;
	cout << "end." << endl;

	//insert the expansion of exp(ax+be)
	// = sum_{m=0}^order (ax+be)^m/m!
	pxeProduct.varCount = 2;
	bool isZeroVector = true;
	for(int i = 0; i < I && isZeroVector; ++i)
		if ( ! IsZero(av[i])  )
			isZeroVector = false;

	//do expansion of exp(be) or exp(ax)
	//insert when m = 0.
	int monomialExp[2] = {0,0};
	PeriodicFunction b2;
	b2.setToConstant(1);
	pxeProduct.insertMonomial(b2, monomialExp);
	ZZ fractorial = to_ZZ(1);

	//insert for m=1..order
	for(int m = 1; m <= order; ++m)
	{

		fractorial *= m;
		if(isZeroVector)
			b2 = b;
		else
			b2 = a;
		cout << "m=" <<m << ", b2=" << b2 << ", a=" << a << endl;
		b2.pow(m);
		b2.div(fractorial);
		if(isZeroVector)
			monomialExp[1] = m; //expand exp(be)
		else
			monomialExp[0]= m; //expand exp(ax)
		pxeProduct.insertMonomial(b2, monomialExp);
		cout << "pxeProduct inserted " << monomialExp[0] << ", " << monomialExp[1] << "::" << b2 << endl;
	}
	cout << "pxe monomials" << endl;
	cout << pxeProduct.printMonomials() << endl;

}



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

