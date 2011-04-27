/*
 * VolumeDB.cpp
 *
 *  Created on: Jan 25, 2011
 *      Author: bedutra
 */

#include "VolumeDB.h"
#include <sstream>
#include <cassert>
#include <iostream>

VolumeDB::VolumeDB()
{

}


bool VolumeDB::canVolumeTestFinish(AlgorithemUsed alg, int dim, int vertex, bool useDual, int timeLimit)
{
	stringstream sql;
	vector<vector<string> > q; //query.
	string algStr;
	algStr = (alg == Lawrence ? "timeLawrence" : "timeTriangulate");

	//is there a -2 time for this case?
	if ( useDual == true)
	{
		sql << "select count(*) from volume as v"
			<< " join polytope as dualP on dualP.rowid = v.polytopeID"
			<< " join polytope as orgP on orgP.rowid = dualP.dual"
			<< " where orgP.dim = " << dim
			<< " and orgP.vertexCount <= " << vertex
			<< " and v." << algStr << " = -2";
	}
	else
	{
		sql << "select count(*) from volume as v"
			<< " join polytope as t on t.rowid = v.polytopeID"
			<< " where t.dim = " << dim
			<< " and t.vertexCount <= " << vertex
			<< " and v." << algStr << " = -2"
			<< " and t.dual is null";
	}//not dual.

	//cout << "check -2" << sql.str().c_str() << endl;
	if (queryAsInteger(sql.str().c_str()))
		return false; //


	sql.str("");

	//what is the average? Make sure we at least have 3 tests before we say the average is to large.
	if ( useDual == true)
	{
		sql << "select avg(v." << algStr << "), count(*) from volume as v"
			<< " join polytope as dualP on dualP.rowid = v.polytopeID"
			<< " join polytope as orgP on orgP.rowid = dualP.dual"
			<< " where orgP.dim = " << dim
			<< " and orgP.vertexCount = " << vertex
			<< " and v." << algStr << " >= 0";
	}
	else
	{
		sql << "select avg(v." << algStr << "), count(*) from volume as v"
			<< " join polytope as t on t.rowid = v.polytopeID"
			<< " where t.dim = " << dim
			<< " and t.vertexCount = " << vertex
			<< " and v." << algStr << " >= 0"
			<< " and t.dual is null";

	}//not dual.

	//cout << "avg: " << sql.str().c_str() << endl;
	q = query(sql.str().c_str());
	if ( q[0][0] != "NULL" && atoi(q[0][0].c_str()) > timeLimit && atoi(q[0][1].c_str()) >= 3)
		return false; //we have more than 3 tests and the avg. is larger than timeLimit.

	return true;
}

bool VolumeDB::getLimit(AlgorithemUsed alg, int dim, int vertexCount, bool useDual)
{
	stringstream sql;
	string strAlg;

	if (alg == Lawrence)
		strAlg = "timeLawrence";
	else
		strAlg = "timeTriangulate";


	if (useDual == true)
	{
		sql << "select count(i." << strAlg << ")"
			<< " from volume as i"
			<< " join polytope as dualP on dualP.rowid = i.polytopeID"
			<< " join polytope as orgP on orgP.rowid = dualP.dual"
			<< " where dualP.dual is not null " //and with orgp
			<< " and orgP.dim = " << dim
			<< " and orgP.vertexCount = " << vertexCount
			<< " and i." << strAlg << " = -2 ";
	}//if dual
	else
	{

		sql << "select count(i." << strAlg << ") "
			<< " from volume as i"
			<< " join polytope as t on t.rowid = i.polytopeID"
			<< " where t.dual is null"
			<< " and t.dim = " << dim
			<< " and t.vertexCount = " << vertexCount
			<< " and i." << strAlg << " = -2 ";
	}//regular

	bool ans = (queryAsInteger(sql.str().c_str()) > 0);

	//if (dim == 10)
	//	cerr << "vertexCount=" << vertexCount << endl;
	if (ans)
	{
		cerr << "HEY, I found a -2 in " << dim << ", " << vertexCount << ", " << strAlg.c_str() << endl;
	//	exit(1);
	}

	return ans;


}

int VolumeDB::getNumberVolumeTest(int dim, int vertexCount, bool useDual)
{
	stringstream sql;
	if (useDual == true)
	{
		sql << "select count(*) from volume as v"
			<< " join polytope as dualP on v.polytopeID = dualP.rowid"
			<< " join polytope as orgP on orgP.rowid = dualP.dual"
			<< " where orgP.vertexCount = " << vertexCount
			<< " and orgP.dim = " << dim;
	}
	else
	{
		sql << "select count(*) from volume as v"
			<< " join polytope as t on v.polytopeID = t.rowid"
			<< " where t.vertexCount = " << vertexCount
			<< " and t.dim = " << dim
			<< " and t.dual is null";
	}
	return queryAsInteger(sql.str().c_str());
}//getNumberVolumeTest

vector<vector<string> > VolumeDB::getRowsToFindVolume(int dim, int vertex, bool useDual, int limit)
{
	stringstream sql;

	if ( useDual == true)
	{
		sql << "select dualP.latteFilePath, v.timeLawrence, v.timeTriangulate, v.theVolume, v.rowid"
		    << " from volume as v"
			<< " join polytope as dualP on dualP.rowid = v.polytopeID"
			<< " join polytope as orgP on orgP.rowid = dualP.dual"
			<< " where orgP.dim = " << dim
			<< " and orgP.vertexCount = " << vertex
			<< " limit " << limit;
	}
	else
	{
		sql << "select t.latteFilePath, v.timeLawrence, v.timeTriangulate, v.theVolume, v.rowid"
		    << " from volume as v"
			<< " join polytope as t on t.rowid = v.polytopeID"
			<< " where t.dim = " << dim
			<< " and t.dual is null"
			<< " and t.vertexCount = " << vertex
			<< " limit " << limit;
	}//else not dual

	//cout << "volumeDB::getRowsToFindVolume::" << sql.str() << endl;
	return query(sql.str().c_str());

}//	getRowsToFindVolume

/**
 * Return rows in the form (latte file, time lawrence, time triangulate, volume, rowid)
 */
vector<vector<string> > VolumeDB::getRowsToFindVolumeGivenSpecficFile(const char *polymakeFile, bool useDual)
{
	stringstream sql;

	if ( useDual == true)
	{
		sql << "select dualP.latteFilePath, v.timeLawrence, v.timeTriangulate, v.theVolume, v.rowid"
		    << " from volume as v"
			<< " join polytope as dualP on dualP.rowid = v.polytopeID"
			<< " join polytope as orgP on orgP.rowid = dualP.dual"
			<< " where orgP.polymakeFilePath = '" << polymakeFile << "'";
	}
	else
	{
		sql << "select t.latteFilePath, v.timeLawrence, v.timeTriangulate, v.theVolume, v.rowid"
		    << " from volume as v"
			<< " join polytope as t on t.rowid = v.polytopeID"
			<< " where t.polymakeFilePath = '" << polymakeFile << "'";
	}

	//cout << "getRowsToFindVolumeGivenSpecficFile:: " << sql.str().c_str() << endl;
	return query(sql.str().c_str());
}//getRowsVolumeGivenSpecficFile


/**
 * Returns a table of volume statistics. Rows are dim, column is additional points.
 *
 * I define additional points as the number of additional points added to make a non-simplex.
 * So, additional points = vertexCount - dim -1. (dim+1 = a simplex).
 * TODO: add a text picture.
 */
vector<vector<ValuationDBStatistics> > VolumeDB::getStatistics(bool useDual)
{
	stringstream sql;
	vector<vector<string> > dimList, pointList;
	vector<vector<ValuationDBStatistics> > ans;

	//first, get a list of dim and "additional point" lists.
	sql << "select distinct dim from polytope order by dim asc";
	dimList = query(sql.str().c_str());

	sql.str("");
	sql << "select distinct (vertexCount - dim -1) as addPoint from polytope where dual is null order by addPoint asc";
	pointList = query(sql.str().c_str());

	cerr << "getStatistics:: filling " << dimList.size() << " by " << pointList.size() << "table\n";

	ans.resize(dimList.size());
	//now loop over each dim-point pair
	for(int i = 0; i < dimList.size(); ++i)
	{
		ans[i].resize(pointList.size());
		for(int j = 0; j < pointList.size(); ++j)
		{
			cerr << "row col = " << i << ", " << j << "\n";
			ans[i][j] = getStatisticsByDimAndPointCount(atoi(dimList[i][0].c_str()), atoi(pointList[j][0].c_str()), useDual);
		}
	}

	return ans;

}//getStatistics


ValuationDBStatistics VolumeDB::getStatisticsByDimAndPointCount(int dim, int additionalPoints, bool useDual)
{
	ValuationDBStatistics vdbs;
	vector<double> avgMinMaxCountLawrence, avgMinMaxCountTriangulate;
	int vertexCount = dim + 1 + additionalPoints;

	//save how this function was called.
	vdbs.dim = dim;
	vdbs.vertexCount = dim + 1 + additionalPoints;
	vdbs.degree = 0;
	vdbs.useDual = useDual;

	//get avg, min, man, and number finished
	avgMinMaxCountLawrence    = getStatisticsAvgMinMaxCount(Lawrence, dim, vertexCount, useDual);
	avgMinMaxCountTriangulate = getStatisticsAvgMinMaxCount(Triangulate, dim, vertexCount, useDual);

	vdbs.avgTriangulationTime = avgMinMaxCountTriangulate[0];
	vdbs.avgLawrenceTime      = avgMinMaxCountLawrence[0];

	vdbs.minTriangulationTime = avgMinMaxCountTriangulate[1];
	vdbs.minLawrenceTime      = avgMinMaxCountLawrence[1];

	vdbs.maxTriangulationTime = avgMinMaxCountTriangulate[2];
	vdbs.maxLawrenceTime      = avgMinMaxCountLawrence[2];

	vdbs.totalFinishedTriangulationTestCases = avgMinMaxCountTriangulate[3];
	vdbs.totalFinishedLawrenceTestCases      = avgMinMaxCountLawrence[3];

	vdbs.totalTestCases = getNumberVolumeTest(dim, vertexCount, useDual);

	vdbs.manuallyLimitedLawrence = getLimit(Lawrence, dim, vertexCount, useDual);
	vdbs.manuallyLimitedTriangulation = getLimit(Triangulate, dim, vertexCount, useDual);

	return vdbs;
}//ValuationDBStatistics getStatisticsByDimAndPointCount(int dim, int additionalPoints, bool useDual);

vector<double> VolumeDB::getStatisticsAvgMinMaxCount(AlgorithemUsed alg, int dim, int vertexCount, bool useDual)
{
	vector<double> ans;
	vector<vector<string> > strAns;
	stringstream sql;
	string strAlg;

	strAlg = (alg == Lawrence ? "timeLawrence" : "timeTriangulate");

	if (useDual == true)
	{
		sql << "select avg(v." << strAlg << "), min(v." << strAlg << "), max(v." << strAlg << "), count(*) "
			<< " from volume as v"
			<< " join polytope as dualP on v.polytopeID = dualP.rowid"
			<< " join polytope as orgP on orgP.rowid = dualP.dual"
			<< " where dualP.dual is not null"
			<< " and orgP.dim = " << dim
			<< " and orgP.vertexCount = " << vertexCount
			<< " and v." << strAlg << " >= 0 ";
	}//if dual
	else
	{
		sql << "select avg(v." << strAlg << "), min(v." << strAlg << "), max(v." << strAlg << "), count(*) "
			<< " from volume as v"
			<< " join polytope as t on v.polytopeID = t.rowid"
			<< " where t.dual is null"
			<< " and t.dim = " << dim
			<< " and t.vertexCount = " << vertexCount
			<< " and v." << strAlg << " >= 0 ";

	}//regular

	//get the data and save it
	strAns = query(sql.str().c_str());

	ans.push_back(atof(strAns[0][0].c_str())); //avg
	ans.push_back(atof(strAns[0][1].c_str())); //min
	ans.push_back(atof(strAns[0][2].c_str())); //max
	ans.push_back(atof(strAns[0][3].c_str())); //count

	return ans;
}//getStatisticsAvgMinMaxCount

void VolumeDB::updateVolumeTimeAndValue(AlgorithemUsed alg, double time, RationalNTL theComputedVolume, string previousValueStr, string rowid)
{
	stringstream sql;

	sql << "update volume set ";
	if ( alg == Lawrence)
		sql << " timeLawrence = " << time;
	else //Triangulate};
		sql << " timeTriangulate = " << time;


	if ( previousValueStr == "NA")
	{
		sql << " , theVolume = ' " << theComputedVolume << "' ";
	}
	else
	{
		RationalNTL previousValue(previousValueStr);
		if ( previousValue != theComputedVolume)
		{
			throw SqliteDBexception(string("updateVolumeTimeAndValue::The volumes differ")
					+"\n\tpreviousValue: " + previousValueStr
					+"\n\tcomputedValue: " + theComputedVolume.str()
					+"\n\tcurrent sql stm:" + sql.str()
					+"\n\trowid: " + rowid);
		}
	}
	sql << " , flagValue = '" << -1 << "', flagType = '" << LAWRENCE_VOLUME_VERSION << "'";
	sql << " where rowid = " << rowid << endl;


	//cout << sql.str().c_str() << endl;
	query(sql.str().c_str());

}//updateVolumeTimeAndValue


int VolumeDB::volumeTestsCompleted(AlgorithemUsed alg, int dim, int vertex, bool useDual)
{
	stringstream sql;

	if (useDual == true)
	{
		sql << "select count(*) from volume as v"
			<< " join polytope as dualP on v.polytopeID = dualP.rowid"
			<< " join polytope as orgP on orgP.rowid = dualP.dual"
			<< " where orgP.dim = " << dim
			<< " and orgP.vertexCount = " << vertex
			<< " and v." << (alg == Lawrence ? "timeLawrence" : "timeTriangulate") << " >= 0";
	}
	else
	{
		sql << "select count(*) from volume as v join polytope as t on v.polytopeID = t.rowid"
			<< " where t.dim = " << dim
			<< " and t.vertexCount = " << vertex
			<< " and v." << (alg == Lawrence ? "timeLawrence" : "timeTriangulate") << " >= 0"
			<< " and t.dual is null";
	}

	return queryAsInteger(sql.str().c_str());
}

