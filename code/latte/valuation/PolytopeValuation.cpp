/*
 * PolytopeValuation.cpp
 *
 *  Created on: Jun 25, 2010
 *      Author: bedutra
 */

#include "PolytopeValuation.h"


using namespace std;


PolytopeValuation::PolytopeValuation()
{
	// TODO Auto-generated constructor stub

}

PolytopeValuation::~PolytopeValuation()
{
	// TODO Auto-generated destructor stub
}




void PolytopeValuation::findVolume(listCone *cones, int numOfVars, unsigned int Flags,
	       char *File_Name, int max_determinant,
	       bool dualize,
	       BarvinokParameters::DecompositionType decomposition)
{
	  Collecting_Single_Cone_Parameters parameters;
	  parameters.Flags = Flags;
	  parameters.Number_of_Variables = numOfVars;
	  parameters.max_determinant = max_determinant;
	  parameters.File_Name = File_Name;
	  parameters.decomposition = decomposition;

	  listCone * allDecompCones = decomposeCones(cones, dualize, parameters);

	  printListConeToFile("cones.decomposed", allDecompCones, numOfVars);

	  cout << "findVolumn():: total of simple cones:" << lengthListCone(allDecompCones) << endl;

	  ZZ sum;

	  for(listCone * current = allDecompCones; current; current = current->rest)
		  sum += current->coefficient * current->dual_determinant;
	  cout << "PolytopeValuation::findVolumn:: total = " << sum << endl;
}//findVolumn

void PolytopeValuation::findVolumeSingle(listCone *cones, int numOfVars,
		unsigned int flags, char *File_Name, int max_determinant, bool dualize,
		BarvinokParameters::DecompositionType decomposition)
{
	Standard_Single_Cone_Parameters *Barvinok_Parameters =
			new Standard_Single_Cone_Parameters;

	Barvinok_Parameters->Flags = flags;
	Barvinok_Parameters->Number_of_Variables = numOfVars;
	Barvinok_Parameters->max_determinant = max_determinant;
	Barvinok_Parameters->File_Name = File_Name;
	Barvinok_Parameters->decomposition = decomposition;

	decomposeAndReturnCones(cones, dualize, *Barvinok_Parameters);

	delete Barvinok_Parameters;

	printListConeToFile("cones.decomposedSingle", cones, numOfVars);

}



listCone * PolytopeValuation::decomposeAndReturnCones(listCone *cones, bool dualize,
			   Standard_Single_Cone_Parameters &param)
{
	int numOfAllCones;
	mat_ZZ mat;

	if (dualize)
	{
		dualizeCones(cones, param.Number_of_Variables, &param);
	}


	numOfAllCones = lengthListCone(cones);
	cerr << numOfAllCones << " cones total to be processed";

	//Set Ten_Power to 100 billion
	// FIXME: What is magic about this number? --mkoeppe, Sat Mar  4 21:21:45 PST 2006
	param.Ten_Power = 1;
	for (int i = 0; i < 12; i++)
		param.Ten_Power *= 10;


	param.Taylor_Expansion_Result = new ZZ[1 + 1];

	param.Degree_of_Taylor_Expansion = 1;
	param.Controller = new Node_Controller(param.Number_of_Variables + 1,
			1);

	cerr << "Number of cones: " << numOfAllCones << endl;

	barvinokDecomposition_List(cones, param);

	cerr << endl << "Total Unimodular Cones: " << param.Total_Uni_Cones << endl;




}//decomposeAndReturnCones






