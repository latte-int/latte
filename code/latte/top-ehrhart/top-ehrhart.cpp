#include "valuation/valuation.h"
#include "latte_relocatable.h"


Valuation::ValuationContainer Valuation::computeTopEhrhart(Polyhedron *poly,
		BarvinokParameters &myParameters, const IntegrationInput & intInput)
{
	if ( intInput.integrandType != IntegrationInput::inputVolume
	     && intInput.integrandType != IntegrationInput::nothing)
	{
		THROW_LATTE_MSG(LattException::bug_Unknown, "integrand type not supported.");
	}
	if (poly->dualized) {
	  dualizeCones(poly->cones, poly->numOfVars, &myParameters);
	  poly->dualized = false;
	}
	{
	  ofstream maple("compute-top-ehrhart.mpl");
	  // FIXME: Install and use installed.
	  maple << "read(\""
		<< relocated_pathname(MAPLE_SCRIPT_DIR)
		<< "/" << "Conebyconeapproximations_08_11_2010.mpl"
		<< "\"):\n";
	  
	  maple << "Delta := [";
	  {
	    listVector *ray;
	    assert(poly->cones != NULL && poly->cones->rest == NULL);
	    for (ray = poly->cones->rays; ray!=NULL; ray = ray->rest) {
	      int i;
	      maple << "[";
	      for (i = 0; i<poly->numOfVars-1; i++) {
		if (i > 0)
		  maple << ", ";
		maple << ray->first[i];
	      }
	      maple << "]" << "/" << ray->first[poly->numOfVars-1];
	      if (ray->rest != NULL)
		maple << ",\n";
	    }
	    maple << "]" << ":" << endl;
	  }

	  maple << "Form := [";
	  int i;
	  for (i = 1; i<poly->numOfVars; i++) {
	    if (i > 1)
	      maple << ", ";
	    maple << "0";
	  }
	  maple << "]" << ":" << endl;

	  maple << "Exponent := 0" << ":" << endl;

	  if (intInput.numEhrhartCoefficients < 0) {
	    maple << "printIncrementalEhrhartweightedPoly";
	    if (intInput.realDilations)
	      maple << "_real";
	    maple << "(n,Delta,Form,Exponent):\n";
	  }
	  else {
	    maple << "k := " << intInput.numEhrhartCoefficients << ":" << endl;
	    maple << "printTopEhrhartweightedPoly";
	    if (intInput.realDilations)
	      maple << "_real";
	    maple << "(n,Delta,Form,Exponent,k):\n";
	  }
	  
	}

	system_with_error_check(MAPLE_PATH + string(" compute-top-ehrhart.mpl"));
	
}
