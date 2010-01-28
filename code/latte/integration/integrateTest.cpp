#include "PolyRep.h"
#include "newIntegration.h"
#include "../timing.h"
#include <iostream>
#include <fstream>
#include <NTL/ZZ.h>

using namespace std;

template <class T>
class FormIntegrateConsumer : public FormSumConsumer<T> {
public:
  FormIntegrateConsumer() { }
  // Take linear form and consume it.
  void ConsumeLinForm(const ZZ& coefficient, int degree, const vec_ZZ& coefs)
  {
	ZZ de = to_ZZ(1);
	for (int i=1;i<=mySimplex->d+degree;i++)
	{
		de=de*i;
	};
	update(coefs,degree,coefficient,de);
  }
  void setFormSum(const string& myForms) { linForms = myForms; }
  string getFormSum() { return linForms; }
  void setDimension(int dimension) { }//also stored in mySimplex
  int getDimension() {}
  void setSimplex(simplexZZ& simplex) { mySimplex = &simplex; numerator = to_ZZ(0); denominator = to_ZZ(0); }
  void getResults(ZZ& num, ZZ& den) {
	if (denominator < 0)
	{
		num = to_ZZ(-1) * numerator; den = to_ZZ(-1) * denominator;
	}
	else
	{
		num = numerator; den = denominator;
	}
   }
  ~FormIntegrateConsumer() {}
private:
  string linForms;
  simplexZZ* mySimplex;
  ZZ numerator, denominator;
  void update(vec_ZZ l,int m, ZZ coe, ZZ de)
  {
	ZZ sum,lcm,total,g,tem;
	int i,j;
	vec_ZZ inner_Pro,sum_Nu,sum_De;
	inner_Pro.SetLength(mySimplex->d+1);
	sum_Nu.SetLength(mySimplex->d+1);
	sum_De.SetLength(mySimplex->d+1);
	total=0;
	lcm=1;
	for (i=0;i<=mySimplex->d;i++)
	{
		sum=0; for (j=0;j<mySimplex->d;j++) {sum=sum+l[j]*mySimplex->s[i][j];};
		inner_Pro[i]=sum;
	};//stores inner product for use
	for (i=0;i<=mySimplex->d;i++)
	{
		sum_Nu[i]=1;for (j=0;j<m+mySimplex->d;j++) sum_Nu[i]=sum_Nu[i]*inner_Pro[i];
		sum_De[i]=1;for (j=0;j<=mySimplex->d;j++) if (i!=j) sum_De[i]=sum_De[i]*(inner_Pro[i]-inner_Pro[j]);
		if (sum_De[i]==0) {cout<<"Warning!"<<l<<" is not regular! Aborted."<<endl; /*exit(1);*/ denominator = to_ZZ(0); return;}; //irregular
		lcm=lcm*sum_De[i]/(GCD(lcm,sum_De[i]));
	};
	for (i=0;i<=mySimplex->d;i++)
	{
		total+=sum_Nu[i]*(lcm/sum_De[i]);
	};
	lcm=lcm*de;
	total=total*mySimplex->v*coe;
	if (numerator==0) {numerator=total;denominator=lcm;}
	else if ((lcm!=0)&&(denominator!=0)) {tem=denominator*lcm/GCD(denominator,lcm);numerator=numerator*tem/denominator+total*tem/lcm;denominator=tem;};	
	g=GCD(numerator,denominator);
	if (g!=0) 
	{
	numerator=numerator/g;
	denominator=denominator/g;};
  }
};

int main(int argc, char *argv[])
{
	if (argc < 3) { cout << "Usage: " << argv[0] << " fileIn fileOut [decompose]" << endl; return 1; }
	bool decomposing = true; //decomposing by default
	bool polynomial = true; //file is assumed to alternate between polynomial and the simplex
	if (argc == 4) { decomposing = (strcmp(argv[3], "1") == 0); };
	string line;
	monomialSum monomials;
	linFormSum forms;
	ifstream myStream (argv[1]);
	ofstream outStream(argv[2]);
	if (!myStream.is_open()) { cout << "Error opening file " << argv[1] << ", please make sure it is spelled correctly." << endl; return 1; }
	int polyCount = 0;
	int dimension;
	int degree = -1;
	int irregularForms = 0;
	float loadTime, decomposeTime, integrateTime = 0.0f;
	Timer myTimer("Integration timer");
	FormIntegrateConsumer<ZZ> *integrator = new FormIntegrateConsumer<ZZ>();
	string testForms;
	while (!myStream.eof())
	{
		getline(myStream, line, '\n');
		if (!line.empty())
		{
			if (polynomial) //reading polynomial
			{
				if (decomposing) //input is sum of monomials that we decompose into sum of linear forms
				{
					myTimer.start();
					loadMonomials(monomials, line);
					myTimer.stop();
					loadTime += myTimer.get_seconds();

					if (monomials.termCount == 0 || monomials.varCount == 0)
					{
						cout << "Error: loaded invalid monomial sum." << endl;
						return 1;
					}

					forms.termCount = 0;
					dimension = forms.varCount = monomials.varCount;
		
					float thisTime = time(NULL);
					cout << "Decomposing " << printMonomials(monomials);
					myTimer.start();
					for (int i = 0; i < monomials.termCount; i++)
					{
						cout << ".";
						decompose(monomials, forms, i);
					}
					myTimer.stop();
					decomposeTime += myTimer.get_seconds();
					cout << endl;
					
					if (forms.termCount == 0)
					{
						cout << "Error: no terms in decomposition to sum of linear forms.";
						return 1;	
					}
					
					outStream << printLinForms(forms) << endl;
					testForms = printLinForms(forms);
					if (degree == -1) //degree is calculated only once
					{
						degree = 0;
						for (int i = 0; i < monomials.varCount; i++)
						{
							degree += monomials.eHead->data[i];
						}
					}
					destroyMonomials(monomials);
				}
				else //input is just linear forms
				{
					myTimer.start();
					loadLinForms(forms, line);
					myTimer.stop();
					loadTime += myTimer.get_seconds();
					if (forms.termCount == 0 || forms.varCount == 0)
					{
						cout << "Error: loaded invalid sum of linear forms.";
						return 1;	
					}
					integrator->setFormSum(line);
				}
				polynomial = false;
				//cout << "Loaded into " << forms.termCount << " linear forms" << endl;
			}
			else //reading simplex
			{
				simplexZZ mySimplex;
				convertToSimplex(mySimplex, line);
				//integrate here
				ZZ numerator, denominator;
				if (decomposing)
				{
					myTimer.start();
					integrateFlatVector(numerator, denominator, forms, mySimplex);
					myTimer.stop();
					integrateTime += myTimer.get_seconds();
					if (IsZero(denominator)) //irregular
					{	
						irregularForms++;
					}
					else //quick and dirty sanity check
					{
						cout << "Verifying by integrating linear forms from string..." << endl;
						integrator->setSimplex(mySimplex);
						parseLinForms(integrator, testForms);
						ZZ a, b;
						integrator->getResults(a, b);
						if (a != numerator || b != denominator)
						{
							cout << "Expected [" << numerator << " / " << denominator << "], ";
							cout << "got [" << a << " / " << b << "]" << endl;
						}
					}
				}
				else
				{
					integrator->setSimplex(mySimplex);
					myTimer.start();
					parseLinForms(integrator, integrator->getFormSum());
					myTimer.stop();
					integrator->getResults(numerator, denominator);
					if (IsZero(denominator)) //irregular
					{	
						irregularForms++;
					}
					integrateTime += myTimer.get_seconds();
				}
				outStream << "[" << numerator << "," << denominator << "]" << endl;
				destroyLinForms(forms);
				polyCount++;
				polynomial = true;
			}
		}
	}
	if (decomposing) { cout << "Dimension " << dimension << ", degree " << degree << ". " << irregularForms << " forms were irregular." << endl; }
	cout << "Total time to load " << polyCount << " polynomials: " << loadTime << ", avg. is " << loadTime / polyCount << endl;
	if (decomposing) { cout << "Total time to decompose " << polyCount << " polynomials: " << decomposeTime << ", avg. is " << decomposeTime / polyCount << endl; }
	cout << "Total time to integrate " << polyCount << " polynomials: " << integrateTime << ", avg. is " << integrateTime / polyCount << endl;
	cout << "Total time is " << (decomposing ? loadTime + integrateTime + decomposeTime : loadTime + integrateTime) << ", avg. is " << (decomposing ? loadTime + integrateTime + decomposeTime : loadTime + integrateTime) / polyCount << endl;

	myStream.close();
	outStream.close();
	return 0; 
}
