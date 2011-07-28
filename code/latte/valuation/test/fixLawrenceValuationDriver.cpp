/*
 * fixLawrenceValuationDriver.cpp
 *
 *  Created on: Apr 10, 2011
 *      Author: bedutra
 */



#include <string>
#include <cstdio>
#include <vector>
#include <sstream>
#include <iostream>
#include <ctime>
#include "../valuation.h"
#include "../../sqlite/VolumeDB.h"
#include "../../sqlite/IntegrationDB.h"

bool findCorrectionFactor(const string & latteFile, RationalNTL & correctionFactor)
{
	string num, line;
	int rows, cols;
	bool checkForError = 0;

	if ( string::npos == latteFile.find(".vrep."))
		return false; //not a vrep file.

	ifstream file(latteFile.c_str());
	if ( !file.is_open() )
	{
		cout << "could not open file" << endl;
		return false;
	}
	file >> rows >> cols;

	for ( int i = 0; i < rows; ++i)
	{
		file >> num;
		//cout << "read in " << i << ", " << num.c_str() << endl;
		for (int j = 1; j < cols; ++j)
			file >> line; //read in the other numbers (and throw it away);

		if ( checkForError == 0)
		{
			correctionFactor = RationalNTL(num);
			checkForError = 1;
		}
		else
		{
			if(correctionFactor != RationalNTL(num))
			{
				cout << "ERROR: correctionFactor: " << correctionFactor << " num on row " << i << " is "  << num.c_str() << endl;
				return false;
			}
		}
	}
	file.close();

	return true; //no errors.
}//findCorrectionFactor



bool findMonomialDegree(const string & monomialFileName, int & monomialDegree)
{
	int monomialCount = 0;
	string monomialString;
	ifstream file(monomialFileName.c_str());


	//read the file into a string.
	if ( ! file.is_open())
	{
		cout << "could not open file " << endl;
		return false;
	}
	getline(file, monomialString, '\n');
	file.close();

	monomialSum monomials;
	loadMonomials(monomials, monomialString); //get the polynomial from the string.

	assert( monomials.termCount == 1); //if there are more than 1 monomial we cannot "undo" our dilation error. :(

	//make an iterator for the monomial.
	BTrieIterator<RationalNTL, int>* polynomialIterator = new BTrieIterator<RationalNTL, int> ();
	polynomialIterator->setTrie(monomials.myMonomials, monomials.varCount);
	polynomialIterator->begin();

	//set up the iterator
	term<RationalNTL, int>* originalMonomial;
	RationalNTL coefficient, totalDegree;
	//get the degree
	for (originalMonomial = polynomialIterator->nextTerm(); originalMonomial; originalMonomial
			= polynomialIterator->nextTerm())
	{
		long totalDegree = 0;


		//find the total degree of the monomial.
		for (int currentPower = 0; currentPower < originalMonomial->length; ++currentPower)
			totalDegree += originalMonomial->exps[currentPower];

		//finally, we have the degree of the monomial.
		monomialDegree = totalDegree; //save the degree to the ouput parameter.
		monomialCount++;
	}
	assert( monomialCount == 1);

	//free up memory.
	destroyMonomials(monomials);
	delete polynomialIterator;


	return true; //no errors.
}//findMonomialDegree




/**
 * The vrep files contained vertices in the form P:=(a, b1, ..., bn), but
 * the lawrence method incorrectly assumed that this was in the form P:=(a, b1, .., bn)
 *
 * We wanted to find the volume of P, but instead we found the volume of aP.
 * Therefore, we need to find the dilation factor a and then dividy the vol(aP)
 * answer by a^dim(P) to get vol(P).
 *
 * NOTE::I have fixed the volume. Please comment out the db-write commands when done, so
 * no one accidently messes with the db.
 */
void fixTheVolume(const char * dbFile)
{
	SqliteDB db;
	stringstream sql;
	vector<vector<string> > results;

	sql << "select v.rowid, p.latteFilePath, p.dim, v.theVolume, v.flagType, v.flagValue "
		<< " from volume as v "
		<< " join polytope as p on p.rowid = v.polytopeID"
		<< " where v.flagType = 0"
		<< " and p.latteFilePath like '%.vrep.%' "
		<< " and v.theVolume <> 'NA'"
		<< " and v.timeLawrence <> -1"
		<< " order by v.timeLawrence"
		;//<< " and p.dim > 2 ";
		//<< " and not p.latteFilePath like '%dual%' ";
	db.open(dbFile);
	results = db.query(sql.str().c_str());

	for (int i = 0; i < results.size(); ++ i)
	{
		//for (int k = 0; k < results[i].size(); ++k)
		//	cout << results[i][k].c_str() << ", ";

		RationalNTL correctionFactor;
		RationalNTL correctValuation;
		if (! findCorrectionFactor(results[i][1], correctionFactor) )
		{
			cout << "Latte file " << results[i][1].c_str() << " contains a problem\n";
			exit(1);
		}

		cout << "Latte file " << results[i][1].c_str() << " has factor " << correctionFactor << endl;

		if ( correctionFactor == 1)
		{
			//mark this row as "processed." We do not need to update the valuation answer.
			sql.str("");
			sql << "update volume set flagType = 1 where rowid = '" << results[i][0] << "'";
			db.query(sql.str().c_str());
			continue;
		}

		correctionFactor.power(atoi(results[i][2].c_str()));
		correctionFactor = RationalNTL(to_ZZ(1), to_ZZ(1))/correctionFactor;
		RationalNTL oldAns(results[i][3].c_str());
		correctValuation = correctionFactor * oldAns;

		//check this factor is correct.
#if 0
		const char * argv[12];
		int j = 0;
		argv[j++] = "./fixTheVolue";
		argv[j++] = "--valuation=volume";
		argv[j++] = "--lawrence";
		argv[j++] = "--vrep";
		argv[j++] = results[i][1].c_str();
		Valuation::ValuationContainer vc = Valuation::mainValuationDriver(argv, j);
		cout << "valuation ended ok" << endl;
		bool found = false;
		for (int k = 0; k < vc.answers.size(); ++k)
			if ( vc.answers[k].valuationType == Valuation::ValuationData::	volumeLawrence)
			{
				found = true;

				//assert(vc.answers[i].answer == correctValuation);
				if ( vc.answers[k].answer != correctValuation)
				{
					cout << "old ans " << oldAns << "\n"
						 << "new ans " << vc.answers[k].answer << "\n"
						 << "pre ans " << correctValuation << "\n"
						 << "factor  " << correctionFactor << "\n";
					exit(1);
				}
			}
		assert ( found == true);
#endif
		//now save the correct volume and correction factor.
		sql.str("");
		sql << "update volume set"
			<< " theVolume = '" << correctValuation << "'"
			<< " , flagType = '1' "
			<< " , flagValue = '" << correctionFactor << "'"
			<< " where rowid = '" << results[i][0] << "'";
		//cout << sql.str().c_str() <<  endl;
		db.query(sql.str().c_str());
	}//for i

	db.close();
}//fixTheVolume


void fixTheIntegration(const char * dbFile)
{
	SqliteDB db;
	stringstream sql;
	vector<vector<string> > results;

	sql << "select i.rowid, t.latteFilePath, p.filePath, t.dim, i.integral, i.flagType, i.flagValue"
		<< " from integrate as i"
		<< " join polytope as t on t.rowid = i.polytopeID"
		<< " join polynomial as p on p.rowid = i.polynomialID"
		<< " where i.flagType = 0"
		<< " and t.latteFilePath like '%.vrep.%'"
		<< " and i.integral <> 'NA'"
		<< " and i.timeLawrence <> -1"
		<< " order by i.timeLawrence"
		;//<< " and t.dim >= 2"
		//<< " and t.dim <= 6";

	db.open(dbFile);
	results = db.query(sql.str().c_str()); //this is a huge query.

	for (int i = 0; i < results.size(); ++ i)
	{
		RationalNTL correctionFactor;
		RationalNTL correctValuation;
		int monomialDegree;

		if ( ! findMonomialDegree(results[i][2], monomialDegree))
		{
			cout << "Monomial file " << results[i][2].c_str() << " contains a problem\n";
			exit(1);
		}

		if (! findCorrectionFactor(results[i][1], correctionFactor) )
		{
			cout << "Latte file " << results[i][1].c_str() << " contains a problem\n";
			exit(1);
		}

		cout << "Latte file " << results[i][1].c_str() << " has factor " << correctionFactor << endl;

		if ( correctionFactor == 1)
		{
			//mark this row as "processed." We do not need to update the valuation answer.
			sql.str("");
			sql << "update integrate set flagType = 1 where rowid = '" << results[i][0] << "'";
			db.query(sql.str().c_str());
			continue;
		}

		correctionFactor.power(atoi(results[i][3].c_str()) + monomialDegree);
		correctionFactor = RationalNTL(to_ZZ(1), to_ZZ(1))/correctionFactor;
		RationalNTL oldAns(results[i][4].c_str());
		correctValuation = correctionFactor * oldAns;

		//check this factor is correct.
#if 0
		const char * argv[12];
		int j = 0;
		argv[j++] = "./fixTheVolue";
		argv[j++] = "--valuation=integrate";
		//make a buffer for the monomial arg.
		char monomialArg[200];
		assert(results[i][2].length() < 150);
		strcpy(monomialArg, "--monomials=");
		strcat(monomialArg, results[i][2].c_str());
		argv[j++] = monomialArg;
		argv[j++] = "--lawrence";
		argv[j++] = "--vrep";
		argv[j++] = results[i][1].c_str();
		Valuation::ValuationContainer vc = Valuation::mainValuationDriver(argv, j);
		cout << "valuation ended ok" << endl;
		bool found = false;
		for (int k = 0; k < vc.answers.size(); ++k)
			if ( vc.answers[k].valuationType == Valuation::ValuationData::	integrateLawrence)
			{
				found = true;

				//assert(vc.answers[i].answer == correctValuation);
				if ( vc.answers[k].answer != correctValuation)
				{
					cout << "old ans " << oldAns << "\n"
						 << "new ans " << vc.answers[k].answer << "\n"
						 << "pre ans " << correctValuation << "\n"
						 << "factor  " << correctionFactor << "\n"
						 << "monomial degree " << monomialDegree << " dim " << results[i][3].c_str() << "\n";
					exit(1);
				}
			}
		assert ( found == true);
#endif

		//now save the correct volume and correction factor.
		sql.str("");
		sql << "update integrate set"
			<< " integral = '" << correctValuation << "'"
			<< " , flagType = '1' "
			<< " , flagValue = '" << correctionFactor << "'"
			<< " where rowid = '" << results[i][0] << "'";
		//cout << sql.str().c_str() <<  endl;
		db.query(sql.str().c_str());
	}//for i

	db.close();


}//fixTheIntegration



int main(int argc, char *argv[])
{


	if( argc <= 1 )
	{
		cout << "error: usage: " << argv[0] << " database-file [volume | integration ]";
		exit(1);
	}

	if (strcmp(argv[2], "volume") == 0)
	{
		fixTheVolume(argv[1]);
	}
	else if ( strcmp(argv[2], "integration") == 0)
	{
		fixTheIntegration(argv[1]);
	}
	else
	{
		cout << "ops, command " << argv[2] << " is not known\n";
	}
	return 0;
}//main
