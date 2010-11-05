#ifndef PERTURBATION_H_
#define PERTURBATION_H_


#include <vector>
#include <NTL/ZZ.h>
#include <NTL/vec_ZZ.h>
#include "cone.h"

using namespace std;

/**
 * For each cone, we perturb the linear form l by l:=(l+e) where e is a vector in the form (epsilon, epsilon, ...).
 * Each <l+e, v> is saved in the rayDotProducts vector.
 * <l+e, v> = <l,v> + <(1,1,1...1),v>*epsilon
 *             \           \_-> this value saved in the epsilon variable
 *              \_-> this value saved in the constant variable
 *If there are repeats, the power is updated. The default power is 1.
 */
class LinearLawrenceIntegration;


class LinearPerturbationContainer
{
public:
	void setListCones(int dim, listCone * simpleConeList);
	bool tryNoPerturbation();
	void findPerturbation(const vec_ZZ &l);

private:
	vec_zz currentPerturbation;
	vector<LinearLawrenceIntegration> coneTerms;
}; //class LinearPerturbationContainer

class LinearLawrenceIntegration
{
public:
	LinearLawrenceIntegration();
	LinearLawrenceIntegration(listCone * cone);

	void setSimplicialCone(listCone *cone);
	bool computeDotProducts(vec_ZZ e); //true = error, we still divide by zero.
	bool computeDotProducts();//true=we divided by zero. need to try an perturbation.
private:


	struct linearPerturbation
	{
		ZZ constant;
		ZZ epsilon;
		ZZ power;
	};

	bool divideByZero; //true if one of the <l,v> terms vanish.
	listCone * simplicialCone; //we treat this as a pointer to a cone, not a pointer to a list of cones...but it is a list of cones.
	vector<linearPerturbation> rayDotProducts;

}; //class LinearLawrenceIntegration



#endif /* PERTURBATION_H_ */
