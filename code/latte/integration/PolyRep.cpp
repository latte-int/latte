#include "PolyRep.h"
#include <stdio.h>
#include <sstream>

void loadPolynomial(polynomial &myPoly, const string line)
{
	//cout << "Loading polynomial `" << line << "'" << endl;
	int varCount = line.find("]]") - (line.find(",[") + 2);
	varCount = (varCount + 1) / 2;
	int termCount = 1; //this should be 1 + occurences of "]," in line
	for (int i = 2; i < line.length() - 1; i++) //first 2 characters are "[["
	{
		if (line.at(i) == ']' && line.at(i+1) == ',') { termCount++; }
	}
	myPoly.exponentBlocks = (eBlock*) calloc(termCount / BLOCK_SIZE, sizeof(eBlock));
	myPoly.coefficientBlocks = (cBlock*) calloc(termCount / BLOCK_SIZE, sizeof(cBlock));
	//still need to ensure that terms are sorted by exponent

	int index = 1; string temp; char *exps, *myTok;
	bool exponents = false;
	int termIndex, varIndex;
	termIndex = 0;
	while (index < line.length() - 1)
	{
		if (exponents)
		{
			temp = line.substr(index + 1, line.find("]]", index) - index - 1);
			index += temp.length() + 4;
			exps = (char*) malloc (sizeof(char) * temp.length());
			strcpy(exps, temp.c_str());
			myTok = strtok(exps, ",");
			myPoly.exponentBlocks[termIndex / BLOCK_SIZE].data[termIndex % BLOCK_SIZE].SetLength(varCount);
			varIndex = 0;
			while (myTok != NULL)
			{
				myPoly.exponentBlocks[termIndex / BLOCK_SIZE].data[termIndex % BLOCK_SIZE][varIndex++] = to_ZZ(myTok);
				myTok = strtok (NULL, ",");
			}
			exponents = false;
			termIndex++;
			free (exps);
		}
		else
		{
			temp = line.substr(index + 1, line.find(",", index) - index - 1);
			index += temp.length() + 2;
			myPoly.coefficientBlocks[termIndex / BLOCK_SIZE].data[termIndex % BLOCK_SIZE] = to_ZZ(temp.c_str());
			exponents = true;
		}
	}
	myPoly.termCount = termCount;
	myPoly.varCount = varCount;
}

string printPolynomial(const polynomial &myPoly)
{
	stringstream output (stringstream::in | stringstream::out);
	output << "[";
	for (int i = 0; i < myPoly.termCount; i++)
	{
		output << "[" << myPoly.coefficientBlocks[i / BLOCK_SIZE].data[i % BLOCK_SIZE] << ",[";
		for (int j = 0; j < myPoly.varCount; j++)
		{
			output << myPoly.exponentBlocks[i / BLOCK_SIZE].data[i % BLOCK_SIZE][j];
			if (j + 1 < myPoly.varCount)
			{ output << ","; }
		}
		output << "]]";
		if (i + 1 < myPoly.termCount)
		{ output << ","; }
	}
	output << "]";
	return output.str();
}

void destroyPolynomial(polynomial &myPoly)
{
	for (int i = 0; i < (myPoly.termCount / BLOCK_SIZE + 1); i++)
	{
		for (int j = 0; j < BLOCK_SIZE && (i * BLOCK_SIZE + j) < myPoly.termCount; j++)
		{
			myPoly.exponentBlocks[i].data[j].kill();
		}
	}
	free (myPoly.exponentBlocks);
	free (myPoly.coefficientBlocks);
	myPoly.exponentBlocks = NULL;
	myPoly.coefficientBlocks = NULL;
	myPoly.termCount = myPoly.varCount = 0;
}
