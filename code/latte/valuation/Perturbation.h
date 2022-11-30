#ifndef PERTURBATION_H_
#define PERTURBATION_H_


#include <vector>
#include <NTL/ZZ.h>
#include <NTL/vec_ZZ.h>
#include "rational.h"
#include "cone.h"

using namespace std;

/**
 * For each cone, we perturb the linear form l by l:=(l+e) where e is a vector in the form (a1*epsilon, a2*epsilon, ...).
 * Each <l+e, v> is saved in the rayDotProducts vector.
 * <l+e, v> = <l,v> + <(a1,...,an),v>*epsilon
 *             \           \_-> this value saved in the epsilon variable
 *              \_-> this value saved in the constant variable
 * If there are repeats, the power is updated. The default power is 0.
 *
 *
 * Math: we need to compute
 *  	<v, l>^d+M * abs(det(matrix formed by the rays)) * M!/(d+M)!
 *  --------------------------------------------------------------
 *                <r_1, -l> * <r_2, -l>*...*<r_d, -l>
 *
 * where v is a vertex of a cone
 *       l is the linear form
 *       r is a ray of the simply cone
 *       M is the power of the linear form
 *       d is the dimension.
 *
 * If we divide by zero, we pick a perturbation e s.t
 *  	<v, l+e>^d+M * abs(det(matrix formed by the rays)) * M!/(d+M)!
 *  --------------------------------------------------------------
 *            <r_1, -l -e> * <r_2, -l -e>*...*<r_d, -l-e>
 *
 */
class LinearLawrenceIntegration;

extern void computeResidueLawrence(const int d, const int M, const LinearLawrenceIntegration & coneTerm, ZZ &numerator, ZZ &denominator);


class LinearPerturbationContainer
{
public:
	LinearPerturbationContainer();
	void setLatticeInformation(const mat_ZZ *_latticeInverse, const ZZ * _latticeInverseDilation);
	void setListCones(int dim, listCone * simpleConeList);
	void findPerturbation(const vec_ZZ &l);
	RationalNTL integratePolytope(int m); //integrates the polytope over l.
private:

	bool tryCurrentPerturbation(const vec_ZZ &l);
	bool tryNoPerturbation(const vec_ZZ &l);

	bool divideByZero; //true = currentPerturbation causes a divide by zero (currentPerturbation could be the zero vector).
	int dimension;
	int numOfRays;
	vec_ZZ currentPerturbation;
	vector<LinearLawrenceIntegration> coneTerms;


	const mat_ZZ *latticeInverse;
	const ZZ * latticeInverseDilation;
}; //class LinearPerturbationContainer

class LinearLawrenceIntegration
{
public:
	LinearLawrenceIntegration();
	LinearLawrenceIntegration(listCone * cone);

	void setSimplicialCone(listCone *cone, int numOfRays);
	bool computeDotProducts(const vec_ZZ &e, const vec_ZZ &l); //true = error, we still divide by zero.
	bool computeDotProducts(const vec_ZZ & l, const mat_ZZ * latticeInverse, const ZZ * latticeInverseDilation);//true=we divided by zero. need to try an perturbation.
	void integrateTerm(RationalNTL &totalSum, int m, int dim);
	void printTerm(bool printRayVertex = false) const; //prints out the current term.
	void updatePowers(); //merges powers and places the location of the (0+e) term in index 0 of rayDotProducts.
private:
	struct linearPerturbation
	{
		ZZ constant;	//a number
		ZZ epsilon;		//coeff. of epsilon from the perturbation
		int power;		//0 if not processed yet. power= k where k is the number of times the term repeats or -1 if it is one of the (k-1) repeated terms.
	}; // (constant + epsilon *e )^power if constant!=0. If constant = 0, then the semantic is 0 + epsilon * (e^power)

	bool divideByZero; //true if one of the <l,ray> terms vanish.
	listCone * simplicialCone; //we treat this as a pointer to a cone, not a pointer to a list of cones...but it is a list of cones.
	vector<linearPerturbation> rayDotProducts;
	linearPerturbation numeratorDotProduct;//power term not used.
	RationalNTL determinant;//abs. det. of the rays of the cone.
	friend void computeResidueLawrence(const int d, const int M, const LinearLawrenceIntegration & coneTerm, ZZ &numerator, ZZ &denominator); //does the residue-like calculation.

}; //class LinearLawrenceIntegration



#endif /* PERTURBATION_H_ */
