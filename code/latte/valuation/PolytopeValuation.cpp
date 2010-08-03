/*
 * PolytopeValuation.cpp
 *
 *  Created on: Jun 25, 2010
 *      Author: Brandon Dutra and Gregory Pinto
 */

#include "PolytopeValuation.h"

using namespace std;

/**
 *
 */
PolytopeValuation::PolytopeValuation(Polyhedron *p, BarvinokParameters &bp) :
	parameters(bp), poly(p), vertexRayCones(NULL), polytopeAsOneCone(NULL),
			triangulatedPoly(NULL), freeVertexRayCones(0),
			freePolytopeAsOneCone(0), freeTriangulatedPoly(0)
{
	numOfVars = parameters.Number_of_Variables; //keep number of origional variables.

	if (p->unbounded)
	{
		cout << "Ops, cannot compute valuation for unbounded polyhedron."
				<< endl;
		exit(1);
	}//check the polytope is bounded.

	if (p->homogenized == false)
		vertexRayCones = p->cones; //Polyhedron is a list of vertex-ray cones.
	else
	{
		polytopeAsOneCone = p->cones;
		cout
				<< "Ops, I have not implemented passing a polyhedron as one cone yet"
				<< endl;
		//if homogenized is true and a one is placed in the first element of each ray/vertex, then
		//we should safely just assign p->cones to polytopeAsOneCone.
		//But if homogenized means the extra one was placed in a different location (say the last element), I'll have
		//to change how I am making and reading such representations of one-cone polytopes in this class to be consistant the latte.
		exit(1);
	}//else the polyhedron is given by one cone.

	srand(time(0));
}//constructor

/**
 * Saves a listCone* depending if the listCone encodes vertex-rays or a triangulation.
 *

 PolytopeValuation::PolytopeValuation(listCone *cones, ConeType coneType,
 int numofvars, BarvinokParameters &bp) :
 vertexRayCones(NULL), parameters(bp), polytopeAsOneCone(NULL),
 triangulatedPoly(NULL), numOfVars(numofvars),
 freeVertexRayCones(0), freePolytopeAsOneCone(0),
 freeTriangulatedPoly(0)
 {
 if (coneType == VertexRayCones)
 vertexRayCones = cones;
 else if (coneType == TriangulatedCones)
 triangulatedPoly = cones;
 srand(time(0));
 }
 */

PolytopeValuation::~PolytopeValuation()
{
	//don't free vertexRayCones, because we did not make them!
	if (polytopeAsOneCone && freePolytopeAsOneCone)
		freeListCone(polytopeAsOneCone);
	if (triangulatedPoly && freeTriangulatedPoly)
		freeListCone(triangulatedPoly);
}

/**
 *  Takes the vertex-ray representation of the polytope, extracts the vertex information,
 *  then then creates one cone with the vertex at the origin and who's rays are integer
 *  multiple of the (dilated) vertex with a leading 1.
 *
 *  Example: if the polytope has vertex { (3, 3/4), (5, 1/2), (1/2, 1/2)} and dilationFactor = 1
 *  then the new cone will have vertex (0 0 0) and integer rays
 *  (1, 3, 3/4)*4, (1, 5, 1/2)*2, (1, 1/2, 1/2)*2
 *
 *  Example: if the polytope has vertex { (3, 3/4), (5, 1/2), (1/2, 1/2)} and dilationFactor = 2
 *  then the new cone will have vertex (0 0 0) and integer rays
 *  (1, 6, 3/2)*2, (1, 10, 1)*1, (1, 1, 1)*1
 *
 */
void PolytopeValuation::convertToOneCone(const RationalNTL &dilationFactor)
{
	if (polytopeAsOneCone)
		return; //already did computation.
	if (triangulatedPoly)
		return; //don't care about converting to one cone for triangulation because already triangulated!
	if (!vertexRayCones)
	{
		cout << "vertexRayCones* is not defined" << endl;
		exit(1);
	}//error.

	listCone * oneCone = new listCone();
	oneCone->coefficient = 1;
	oneCone->determinant = 0;
	oneCone->subspace_generators = NULL;
	oneCone->dual_determinant = 0;
	oneCone->facets = NULL;
	oneCone->equalities = NULL;
	oneCone->latticePoints = NULL;
	oneCone->rest = NULL;

	//set to zero vector of numofvars + 1 size.
	oneCone->vertex = new Vertex();
	oneCone->vertex->vertex = new rationalVector(numOfVars + 1);

	oneCone->rays = new listVector;
	oneCone->rays->rest = 0;

	//now add the vertex to the rays list with a leading 1: (1, old poly cone vertex).
	//The first entry in masterList should be ignored because masterList->first = masterList->rest->first.
	listVector * masterList = new listVector;

	for (listCone * currentCone = vertexRayCones; currentCone; currentCone
			= currentCone->rest)
	{
		vec_ZZ buildRay; //buildRay = [1, old-vertex]
		ZZ nume, denom;
		buildRay.SetLength(numOfVars + 1);

		ZZ scaleFactor; //scaleRationalVectorToInteger() sets scaleFactor.
		vec_ZZ integerVertex;
		if ( dilationFactor == 1)
			integerVertex = scaleRationalVectorToInteger(
				currentCone->vertex->vertex, numOfVars, scaleFactor);
		else
		{
			rationalVector dilatedVector = *currentCone->vertex->vertex * dilationFactor;
			integerVertex = scaleRationalVectorToInteger(
							&(dilatedVector), numOfVars, scaleFactor); //dilate polytope first, then scale to integer.
		}

		buildRay[0] = scaleFactor; // = 1 * scaleFactor.
		for (int i = 0; i < numOfVars; ++i)
		{
			buildRay[i + 1] = integerVertex[i];
		}//for i

		//cout << buildRay << endl;

		masterList->first = buildRay;
		masterList = appendVectorToListVector(buildRay, masterList);
	}//for currentCone

	//cout << "END  BUILDING THE RAYS" << endl;

	oneCone->rest = 0;
	oneCone->rays = masterList->rest; //ignore masterList->first, so just return the rest and NOT masterList.


	polytopeAsOneCone = oneCone; //finally, save the result.
	freePolytopeAsOneCone = 1; //delete this object in the deconstructor.
}//convertToOneCone


/**
 * Converts maple polynomial strings to a vector.  one monomial is saved as [coef., power of x1, power of x2, ... power of xn]. The powers are all rational.
 *
 * This is motly a copy of void parseMonomials(MonomialConsumer<ZZ>* consumer, const string &line)
 */
void PolytopeValuation::breakupPolynomialToMonomials(
		const string & polynomialString,
		vector<vector<RationalNTL> > & monomialList) const
{
	int varCount = 0;
	for (int i = 0; polynomialString[i] != ']'; i++)
	{
		varCount += (polynomialString[i] == ',');
	}
	if (varCount < 1)
	{
		cout << "polynomialString: `" << polynomialString << "'" << endl;
		cout << "There are " << varCount << " variables, bailing." << endl;
		return;
	}
	//consumer->setDimension(varCount);

	int termIndex, lastPos, expIndex, flag;
	termIndex = lastPos = flag = 0; //0 means we expect coefficient, 1 means we expect exponent vector

	int *exponents = new int[varCount];
	vector<RationalNTL> oneMonomial; // coefficient;

	for (int i = 1; i < polynomialString.length() - 1; i++) //ignore outermost square brackets
	{
		if (polynomialString[i] == '[')
		{

			switch (flag)
			{
				case 0: //coefficient
					lastPos = i + 1;
					for (; polynomialString[i] != ','; i++)
						;

					oneMonomial.push_back(RationalNTL(polynomialString.substr(
							lastPos, i - lastPos).c_str()));
					flag = 1;
					break;
				case 1: //exponent vector
					expIndex = 0;
					for (i++; polynomialString[i] != ']'; i++)
					{
						if (polynomialString[i] != ' ')
						{
							lastPos = i;
							for (; polynomialString[i] != ','
									&& polynomialString[i] != ']'; i++)
								;
							oneMonomial.push_back(RationalNTL(
									polynomialString.substr(lastPos, i
											- lastPos).c_str()));
							//exponents[expIndex++] = atoi(polynomialString.substr(lastPos, i
							//		- lastPos).c_str());
						}
					}
					monomialList.push_back(oneMonomial); //consumer->ConsumeMonomial(coefficient, exponents);
					oneMonomial.clear();
					flag = 0;
					break;
				default: //error
					cout << "Flag is " << flag << ", bailing." << endl;
					return;
			}
		}
	}

	delete[] exponents;

}


/**
 * Computes the volume for one simplex.
 *
 * The volume of an n-simplex in n-dimensional space with vertices (v0, ..., vn) is the abs. value of
 * {1\over n!}\det \begin{pmatrix} v_1-v_0 & v_2-v_0& \dots & v_{n-1}-v_0 & v_n-v_0 \end{pmatrix}
 *
 * However, using facts about how the determinat does not change when adding a multiple of a col. to another col,
 * and det(A) = 1/a * det ( times one col of A by a), we can get away with doing this subtraction,
 * and we do not have to project the vertices back down ( i.e. remove the leading one/element).
 *
 * We compute the abs. value of
 * {1\over n!}  {1 \over {v_o[0]*v_1[0]*...v_n[0]}} \det \begin{pmatrix} v_0 & v_1 & v_2 & \dots & v_{n-1} & v_n \end{pmatrix}
 */
RationalNTL PolytopeValuation::findVolumeUsingDeterminant(
		const listCone * oneSimplex) const
{
	int i, numOfRays;
	mat_ZZ mat;

	vec_ZZ head;
	vec_ZZ tail;
	ZZ numerator, denominator;
	numerator = 1;
	denominator = 1;

	numOfRays = lengthListVector(oneSimplex->rays);

	mat.SetDims(numOfRays, numOfVars + 1);

	i = 0;
	for (listVector * currentRay = oneSimplex->rays; currentRay; currentRay
			= currentRay->rest)
	{
		for (int k = 0; k < numOfVars + 1; ++k)
			mat[i][k] = ((currentRay->first)[k]);
		denominator *= (currentRay->first)[0];
		++i;
	}//for currentRay

	numerator = abs(determinant(mat));
	denominator *= factorial(numOfRays - 1);
	//cout << mat << " = " << determinant(mat) << "\n./." << factorial(numOfRays -1) << endl;
	return RationalNTL(numerator, denominator);
}//findDetermiantForVolume


/* computes the volume of a polytope using the lawrence forumla
 * takes into account the coefficient given to a cone when decomposed into unimodular cones
 * thus it works on all inputs
 * @input: a listCone of the cones and the nnumber of variables (dimension of the space)
 * @return RationalNTL: the volume of the polytope
 *
 * Math: For each simple vertex-ray cone, we sum the fractions
 *
 * 	<v, c>^d * det(matrix formed by the rays) * cone's coefficient
 *  --------------------------------------------------------------
 *                <-r_1, c> * <-r_2, c>*...*<-r_d, c>
 * where v is a vertex,
 *       c is a random vector,
 *       d is the dimension
 *       r_i is the ith ray of the cone.
 *
 * We use the cone's coefficient in case the cone decomposition is signed (ex, unimodular decomposition).
 */
RationalNTL PolytopeValuation::findVolumeUsingLawrence()
{
	RationalNTL answer;

	vec_ZZ c = vec_ZZ();
	ZZ scale = ZZ();
	ZZ num = ZZ();
	ZZ denom = ZZ();
	denom = 1;
	ZZ tempNum = ZZ();
	ZZ tempDenom = ZZ();
	vec_ZZ vert = vec_ZZ();
	vec_ZZ ans = vec_ZZ();
	mat_ZZ mat;

	ZZ det = ZZ();
	mat.SetDims(numOfVars, numOfVars);

	c.SetLength(numOfVars);
	for (int i = 0; i < numOfVars; i++)
		c[i] = rand() % 10000;

	for (listCone * simplicialCone = triangulatedPoly; simplicialCone; simplicialCone
			= simplicialCone->rest)
	{

		//find vertex
		vert = scaleRationalVectorToInteger(simplicialCone->vertex->vertex,
				numOfVars, tempDenom);

		//raise f(vertex) to the power of the dimension
		tempNum = vert * c;
		tempNum = power(tempNum, numOfVars);
		tempDenom = power(tempDenom, numOfVars);

		int col = 0;

		for (listVector * currentRay = simplicialCone->rays; currentRay; currentRay
				= currentRay->rest, col++)
		{
			//

			tempDenom *= -1 * c * currentRay->first;

			//generate matrix
			for (int row = 0; row < numOfVars; row++)
			{
				mat[row][col] = currentRay->first[row];
			}//for every component of the ray

		}//for every ray in the simple cone

		//get the determinant
		determinant(det, mat);

		//multiply by the absolute value of the determinant
		tempNum *= abs(det) * simplicialCone->coefficient;

		//add current term to the running total
		//cout << "adding " << tempNum << " / " << tempDenom << " to " << num
		//		<< " / " << denom << endl;
		answer.add(tempNum, tempDenom);
		//add(num, denom, tempNum, tempDenom);
	}//for every simple cone in the cone

	//	}//for every cone
	ZZ one;
	one = 1;
	answer.mult(one, factorial(numOfVars));

	return answer;
}//findVolumeUsingLarence()


/**
 * Converts the input polynomal vertex-ray input to one polytope for trangulation,
 * then sums the volume of each triangulated simplex.
 *
 * Uses the determinant or lawrence method.
 */
RationalNTL PolytopeValuation::findVolume(VolumeType v)
{
	RationalNTL answer;

	if (v == DeterminantVolume)
	{

		convertToOneCone();
		triangulatePolytopeCone();

		for (listCone * oneSimplex = triangulatedPoly; oneSimplex; oneSimplex
				= oneSimplex->rest)
			answer.add(findVolumeUsingDeterminant(oneSimplex));

		//cout << "findVolumeUsingDeterminant(): VOLUME: " << answer << endl;
	} else if (v == LawrenceVolume)
	{
		triangulatePolytopeVertexRayCone();
		answer = findVolumeUsingLawrence();

		//cout << "findVolumeUsingLawrence(): VOLUME: " << answer << endl;
	}

	return answer;

}//findVolume


/**
 * Computes n!, n >= 0.
 */
ZZ PolytopeValuation::factorial(const int n)
{
	ZZ product;
	product = 1;
	for (int i = n; i > 1; --i)
		product *= i;
	return product;
}//factorial


/**
 * Integrates the polynomial over the polytope. The polytope is written in maple syntax.
 * Example: 2*X_0^3*X_1^4 + 5 ==> [ [2, [3, 4], [5, [0, 0]] ]
 *
 * In the non-integer case, we compute
 *  1/a^n * \sum_{i}{coefficient_i * integral( m_i(x1/a, x2/a, ..., xn/a) dX }
 *  Where m_i is a monomial with a coefficient of 1,
 *  and where 1/(a*a*...*a) is the Jacobian of the change of variables to make the polytope have integer vertices.
 *
 *  We use the fact that if P is a polytope with rational vertices, then
 *
 *  integral over P { f(w1, w2, ..., wn) dw } = integral over P' { f(x1/a, x2/a, ..., xn/a)  * |d(w1, .., wn)/d(x1, .., xn)| dX }
 *  where |d(w1, .., wn)/d(x1, .., xn)| is the absolute value of the Jacobian given by the equations xi = wi * a
 *  and a is such that ever vertex becomes integer when mult. by a.
 *  P' is now a dilation of P such that P' has only integer vertices.
 */
RationalNTL PolytopeValuation::integrate(const string& polynomialString)
{
	RationalNTL answer;
	linFormSum linearForms;

	monomialSum originalPolynomial;// polynomial without the updated coefficients.

	loadMonomials(originalPolynomial, polynomialString);

	//find the dilation factor.
	ZZ dilationFactor;
	dilationFactor = 1;

	for (listCone * currentCone = vertexRayCones; currentCone; currentCone
			= currentCone->rest)
		for (int i = 0; i < numOfVars; ++i)
			dilationFactor = lcm(dilationFactor,
					(currentCone->vertex->vertex->denominators())[i]);

	//dilate the polytope and triangulate.
	cout << "dilation factor = " << dilationFactor << endl;
	convertToOneCone(RationalNTL(dilationFactor, 1)); //every vertex should be integer
	Timer t("Triangilate dilated poly");
	t.start();
	triangulatePolytopeCone(); //every tiangulated vertex is now in the form (1, a1, ..., an) such that ai \in Z.
	t.stop();
	cout << t << endl;

	//Goal: find new polynomial.
	//define the new polynomial.
	monomialSum transformedPolynomial; //updated coefficients.
	transformedPolynomial.termCount = 0;
	transformedPolynomial.varCount = originalPolynomial.varCount;

	//make an iterator for the original polynomial.
	MonomialLoadConsumer<RationalNTL>* transformedPolynomialLoader =
			new MonomialLoadConsumer<RationalNTL> ();
	transformedPolynomialLoader->setMonomialSum(transformedPolynomial);
	BTrieIterator<RationalNTL, int>* originalPolynomialIterator =
			new BTrieIterator<RationalNTL, int> ();
	originalPolynomialIterator->setTrie(originalPolynomial.myMonomials,
			originalPolynomial.varCount);
	originalPolynomialIterator->begin();

	term<RationalNTL, int>* originalMonomial;
	RationalNTL coefficient, totalDegree;
	//loop over the original polynomial, and insert the updated monomial into the transformedPolynomial
	for (originalMonomial = originalPolynomialIterator->nextTerm(); originalMonomial; originalMonomial
			= originalPolynomialIterator->nextTerm())
	{
		long totalDegree = 0;
		coefficient = originalMonomial->coef;

		//find the total degree of the monomial.
		for(int currentPower = 0; currentPower < originalMonomial->length; ++currentPower)
			totalDegree += originalMonomial->exps[currentPower];
		//cout << "factor^degree = " << dilationFactor << "^" << originalMonomial->degree << endl;
		//cout << "length = " << originalMonomial->length << endl;
		coefficient.div(power(dilationFactor, totalDegree));

		transformedPolynomialLoader->ConsumeMonomial(coefficient,
				originalMonomial->exps);
	}//for every term in the originalPolynomial

	assert(originalPolynomial.termCount == transformedPolynomial.termCount
			&& originalPolynomial.varCount == transformedPolynomial.varCount );

	//make an iterator for the transformed polynomial and decompost it into linear forms.
	BTrieIterator<RationalNTL, int>* transformedPolynomialIterator =
			new BTrieIterator<RationalNTL, int> ();
	linearForms.termCount = 0;
	linearForms.varCount = transformedPolynomial.varCount;
	transformedPolynomialIterator->setTrie(transformedPolynomial.myMonomials,
			transformedPolynomial.varCount);
	decompose(transformedPolynomialIterator, linearForms);


	answer.add(integratePolytope(linearForms));

	answer.div(power(dilationFactor, transformedPolynomial.varCount));

	destroyLinForms(linearForms);
	destroyMonomials(originalPolynomial);
	destroyMonomials(transformedPolynomial);

	delete transformedPolynomialLoader;

	delete originalPolynomialIterator;
	delete transformedPolynomialIterator;

	return answer;

	/*
	 linFormSum forms;

	 int numberOfPolynomialVariables; //same as the dimension of the polytope.

	 //find each monomial.
	 vector<vector<RationalNTL> > monomialList; // one monomial is saved as [coef., power of x1, power of x2, ... power of xn]. The powers are all rational.
	 breakupPolynomialToMonomials(polynomialString, monomialList); //saves each monomial into the list.

	 numberOfPolynomialVariables = monomialList[0].size() - 1;

	 //cout << "number of monomials =" << monomialList.size() << endl;

	 //find the dilation factor.
	 ZZ dilationFactor;
	 dilationFactor = 1;

	 for (listCone * currentCone = vertexRayCones; currentCone; currentCone
	 = currentCone->rest)
	 for (int i = 0; i < numOfVars; ++i)
	 dilationFactor = lcm(dilationFactor,
	 (currentCone->vertex->vertex->denominators())[i]);

	 cout << "dilation factor = " << dilationFactor << endl;
	 dilatePolytope(RationalNTL(dilationFactor, to_ZZ(1))); //dilate so that every vertex is integer
	 convertToOneCone(); //every vertex should be integer
	 triangulatePolytopeCone(); //every tiangulated vertex is now in the form (1, a1, ..., an) such that ai \in Z.


	 int *exponents;
	 exponents = new int[numberOfPolynomialVariables];

	 //cout << "exponents size = " << numberOfPolynomialVariables1 << endl;

	 for (int i = 0; i < monomialList.size(); ++i)
	 {
	 RationalNTL coefficient;
	 coefficient = monomialList[i][0]; //get the original coefficient.

	 //cout << "original coefficient = " << coefficient << endl;
	 //get the new coefficient.
	 long totalDegree;
	 totalDegree = 0;
	 for (int k = 1; k < monomialList[i].size(); ++k)
	 totalDegree += to_int(monomialList[i][k].getNumerator());
	 coefficient.div(power(dilationFactor, totalDegree));
	 //cout << "new coefficient = " << coefficient << endl;

	 //make a new monomial to integrate.
	 monomialSum oneMonomial;
	 oneMonomial.termCount = 0;
	 MonomialLoadConsumer<ZZ>* myLoader = new MonomialLoadConsumer<ZZ> ();
	 myLoader->setMonomialSum(oneMonomial);
	 myLoader->setDimension(numberOfPolynomialVariables); //number of variables in the monomial.
	 //cout << "My loader made" << endl;

	 //get the exponents as c++ ints.
	 for (int k = 0; k < numberOfPolynomialVariables; ++k)
	 {
	 exponents[k] = to_int(monomialList[i][k + 1].getNumerator());
	 //cout << exponents[k] << ", " ;
	 }
	 //cout << "exponent vector made" << endl;

	 myLoader->ConsumeMonomial(to_ZZ(1), exponents); //the coeff. here is 1 because we extracted the real coefficient in "coefficient"

	 //cout << "consume monomial worked" << endl;

	 //cout << "oneMonomila = " << printMonomials(oneMonomial) << endl;

	 //oneMonomial is now constructed. Next make the linear forms.
	 BTrieIterator<ZZ, int>* iterator = new BTrieIterator<ZZ, int> ();
	 forms.termCount = 0;
	 forms.varCount = oneMonomial.varCount;
	 iterator->setTrie(oneMonomial.myMonomials, oneMonomial.varCount);
	 decompose(iterator, forms);

	 //cout << "linear forms = " << printLinForms(forms) << endl;


	 //integrate every cone over this monomial.
	 answer.add(coefficient * integratePolytope(forms));

	 destroyMonomials(oneMonomial);
	 destroyLinForms(forms);
	 delete myLoader;
	 delete iterator;
	 }//for each monomial.
	 delete[] exponents;

	 answer.div(power(dilationFactor, numberOfPolynomialVariables)); // monomialList[0].size() - 1 = number of variables.


	 */
}//integrate

//integrates the polynomial over the polytope using the rational-vertex integration functions.
RationalNTL PolytopeValuation::integrateRational(const string& polynomialString)
{
	RationalNTL answer;
	linFormSum linearForms;

	monomialSum polynomial;// polynomial without the updated coefficients.

	loadMonomials(polynomial, polynomialString);


	convertToOneCone();
	Timer t2("Triangulate rational poly");
	t2.start();
	triangulatePolytopeCone();
	t2.stop();
	cout << t2 << endl;

	//make an iterator for the polynomial and decompost it into linear forms.
	BTrieIterator<RationalNTL, int>* polynomialIterator =
			new BTrieIterator<RationalNTL, int> ();
	linearForms.termCount = 0;
	linearForms.varCount = polynomial.varCount;
	polynomialIterator->setTrie(polynomial.myMonomials,
			polynomial.varCount);
	decompose(polynomialIterator, linearForms);


	answer.add(integrateRationalPolytope(linearForms));


	destroyLinForms(linearForms);
	destroyMonomials(polynomial);

	delete polynomialIterator;


	return answer;

}//integrate


/**
 * Computes the integral of the linear forms over the polytope.
 * Assumes the polytope has already been triangulated and has integer vertices.
 *
 */
RationalNTL PolytopeValuation::integratePolytope(linFormSum &forms) const
{
	RationalNTL answer;

	BTrieIterator<RationalNTL, ZZ>* linearFormIterator = new BTrieIterator<
			RationalNTL, ZZ> ();
	linearFormIterator->setTrie(forms.myForms, forms.varCount);

	for (listCone * currentCone = triangulatedPoly; currentCone; currentCone
			= currentCone->rest)
	{
		//cout << "starting to integrate over a new cone" << endl;
		//struct simplexZZ
		//	int d;
		//	vec_vec_ZZ s;
		//	ZZ v;

		//First construct a simplex.
		simplexZZ oneSimplex;
		oneSimplex.d = forms.varCount;

		listVector * rays = currentCone->rays;

		int vertexCount = 0; //the current vertex number being processed.
		oneSimplex.s.SetLength(numOfVars + 1);

		for (rays = rays; rays; rays = rays->rest, ++vertexCount)
		{
			oneSimplex.s[vertexCount].SetLength(numOfVars);
			assert( rays->first[0] == 1); //make sure the triangulation is such that the vertices of the original polytope is integer.
			for (int k = 0; k < numOfVars; ++k)
				oneSimplex.s[vertexCount][k] = rays->first[k + 1];

		}//create the simplex. Don't copy the leading 1.

		//comute the volume of the Parallelepiped
		mat_ZZ matt;
		matt.SetDims(oneSimplex.d, oneSimplex.d);
		for (int j = 1; j <= oneSimplex.d; j++)
			matt[j - 1] = oneSimplex.s[j] - oneSimplex.s[0];
		oneSimplex.v = abs(determinant(matt));

		ZZ numerator, denominator;
		//cout << printLinForms(forms) << endl;

		//cout << "print oneSimplex" << endl;
		//oneSimplex.print(cout);

		//Finally, we are ready to integrate the linear form over the simplex!
		integrateLinFormSum(numerator, denominator, linearFormIterator,
				oneSimplex);

		answer.add(numerator, denominator);
		//answer.add(numerator * coefficient.getNumerator(), denominator * coefficient.getDenominator());
	}//for every triangulated simplex.
	return answer;
}//integratePolytope()


/**
 * Computes the integral of the linear forms over the polytope.
 * Assumes the polytope has already been triangulated, does not need integer vertices.
 *
 */
RationalNTL PolytopeValuation::integrateRationalPolytope(linFormSum &forms) const
{
	RationalNTL answer;
	ZZ determinantCorrectionFactor;

	BTrieIterator<RationalNTL, ZZ>* linearFormIterator = new BTrieIterator<
			RationalNTL, ZZ> ();
	linearFormIterator->setTrie(forms.myForms, forms.varCount);

	for (listCone * currentCone = triangulatedPoly; currentCone; currentCone
			= currentCone->rest)
	{
		//cout << "starting to integrate over a new cone" << endl;
		//struct simplexZZ
		//	int d;
		//	vec_vec_ZZ s;
		//	ZZ v;

		determinantCorrectionFactor = 1;

		//First construct a simplex.
		simplexRationalNTL oneSimplex;
		oneSimplex.d = forms.varCount;

		listVector * rays = currentCone->rays;

		int vertexCount = 0; //the current vertex number being processed.
		oneSimplex.s.resize(numOfVars + 1);

		mat_ZZ matt;
		matt.SetDims(oneSimplex.d + 1, oneSimplex.d + 1);
		for (rays = rays; rays; rays = rays->rest, ++vertexCount)
		{
			oneSimplex.s[vertexCount].SetLength(numOfVars);


			for (int k = 0; k < numOfVars; ++k)
			{
				oneSimplex.s[vertexCount][k] = RationalNTL(rays->first[k + 1], rays->first[0]);
				matt[vertexCount][k] = rays->first[k];
			}
			matt[vertexCount][numOfVars] = rays->first[numOfVars];

			determinantCorrectionFactor *= rays->first[0];
		}//create the simplex. Don't copy the leading "1"

		//compute the volume of the Parallelepiped
		//MATH: let v_i = (a_i, v_i') where a_i is integer, and v_i' is the rest of the vector and is integer.
		// Then, det([v_i, i=1..d+1]) = det( [ v_1/a1, v_1/a1 - v_2/a2, ..., v_1/a_1 - v_{d+1}/a_{d+1}]) * (a1*a2*...*a_{d+1}).
		// Recall: det ( one row of A times a) = det(A)*a, and adding a mult. of a row/col to another does not change the det.
		//So then abs(det(v2'/a2 - v1'/a1, ..., v_d+1'/a_d+1 - v1'/a1)) = det([v_i, i=1..d+1])) * 1/ abs( a1*a2*...*a_{d+1}) (N1)
		// because the LHS of (N1) has all zeros except for the first col., which corresponds to v1.
		// also, the LHS of (N1) is the volume of the Parallelepiped!!!
		// Expanding on this almost zero row by minors gives the proof why we do not need to subtract the vertex rays.
		oneSimplex.v = abs(determinant(matt));
		oneSimplex.v.div(abs(determinantCorrectionFactor)); //this factor should already be positive....no need for the abs.....
		//cout << "determinantCorrectionFactor = " << determinantCorrectionFactor << endl;



		RationalNTL integralOneTriangulation;
		//Finally, we are ready to integrate the linear form over the simplex!
		integrateLinFormSum(integralOneTriangulation, linearFormIterator, oneSimplex);
		//void integrateLinFormSum(RationalNTL & answer, PolyIterator<RationalNTL, ZZ>* it, const simplexRationalNTL &mySimplex)



		answer.add(integralOneTriangulation);
	}//for every triangulated simplex.
	return answer;
}//integratePolytope()

/**
 * Returns the lowest common multiple of a and b.
 */
ZZ PolytopeValuation::lcm(const ZZ &a, const ZZ & b)
{
	return (a * b) / GCD(a, b);
}//lcm

/**
 * Find the Lawrence rational function for volume.
 */
void PolytopeValuation::printLawrenceVolumeFunction()
{
	listCone * triangulatedCones;
	vec_ZZ vert = vec_ZZ();
	ZZ temp = ZZ();
	mat_ZZ mat;
	ZZ det;
	mat.SetDims(numOfVars, numOfVars);

	triangulatePolytopeVertexRayCone();

	cout << "( ";
	for (listCone * simplicialCone = triangulatedPoly; simplicialCone; simplicialCone
			= simplicialCone->rest)
	{
		vert = scaleRationalVectorToInteger(simplicialCone->vertex->vertex,
				parameters.Number_of_Variables, temp);
		cout << "( ";

		//dot of c and v raised to the dim power
		for (int i = 0; i < parameters.Number_of_Variables; i++)
		{
			cout << vert[i];
			if (temp != 1)
				cout << " / " << temp;
			cout << " * c" << i;
			if (i != parameters.Number_of_Variables - 1)
				cout << " + ";
		}
		cout << " ) ^ " << parameters.Number_of_Variables << " / ( ";

		//the correct sign on the denominator
		if (parameters.Number_of_Variables % 2 == 1)
			cout << "-";

		//divide by the multiplication of all the dots
		//of the rays dotted with c
		int col = 0;

		for (listVector * currentRay = simplicialCone->rays; currentRay; currentRay
				= currentRay->rest, col++)
		{
			cout << "( ";
			for (int row = 0; row < numOfVars; row++)
			{
				cout << currentRay->first[row] << " * c" << row;
				if (row != parameters.Number_of_Variables - 1)
				{
					cout << " + ";
				}
				mat[row][col] = currentRay->first[row];
			}
			cout << " )";
			if (currentRay->rest != NULL)
				cout << " * ";
		}//for every ray

		//get the determinant
		determinant(det, mat);

		//close up the denominator
		cout << " ) * ";

		//multiply by the coefficient and determinant
		cout << simplicialCone->coefficient;
		if (det != 1)
			cout << " * (" << abs(det) << ')';

		//if more cones, type the +
		if (simplicialCone->rest != NULL)
			cout << " + ";
	}//for every simple cone.

	// divide the sum by the factorial
	cout << ") / ( " << parameters.Number_of_Variables << "!";
	cout << " )" << endl;
}

/**
 * Triangulates 1 cone (which encodes the polytope. See the comments for convertToOneCone()
 * to learn how the polytope is encoded in one cone.)
 */
void PolytopeValuation::triangulatePolytopeCone()
{
	if (triangulatedPoly)
		return; //all ready did computation.
	if (polytopeAsOneCone == NULL)
	{
		cout
				<< "PolytopeValuation::triangulatePolytopeCone(): there is no cone to triangulate"
				<< endl;
		exit(1);
	}

//	cout << "BEFORE TRIANGULATION: one cone" << endl;
//	printListCone(polytopeAsOneCone, numOfVars + 1);

	parameters.Number_of_Variables = numOfVars + 1;
	triangulatedPoly = triangulateCone(polytopeAsOneCone, numOfVars + 1,
			&parameters);
//	cout << "AFTER TRIANGULATION: many cones" << endl;
//	printListCone(triangulatedPoly, numOfVars + 1);
	//parameters.Number_of_Variables = numOfVars; //convert back. This is not really needed, because we never look at this value again.
	freeTriangulatedPoly = 1; //Delete this in the deconstructor.
}//triangulateCone()


/**
 * Triangulate the cones from they vertex-ray cone.
 */
void PolytopeValuation::triangulatePolytopeVertexRayCone()
{
	if (triangulatedPoly)
		return; //already did computation


	triangulatedPoly = decomposeCones(vertexRayCones,
			parameters.Number_of_Variables, parameters.Flags,
			parameters.File_Name, 0, true, //change to false later?!?!?!
			BarvinokParameters::DualDecomposition);
	freeTriangulatedPoly = 1; //Delete this in the deconstructor.
}
