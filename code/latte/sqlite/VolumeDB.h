/*
 * VolumeDB.h
 *
 *  Created on: Jan 25, 2011
 *      Author: bedutra
 */

#ifndef VOLUMEDB_H_
#define VOLUMEDB_H_

#include "SqliteDB.h"
#include "../rational.h"
#include "ValuationDBStatistics.h"

using namespace std;
//A B C D E F G H I J K L M N O P Q R S T U V W X Y Z

#define LAWRENCE_VOLUME_VERSION 1


/*
 *
 * The volume table was made from reading in the polytope table via a insert into .... select ... from statement.
 * See the createVolume.sql file.
 *
 */
class VolumeDB:public SqliteDB
{
public:
	enum AlgorithemUsed {Lawrence, Triangulate};
	VolumeDB();


	bool canVolumeTestFinish(AlgorithemUsed alg, int dim, int vertex, bool useDual, int timeLimit);


	vector<vector<string> > getRowsToFindVolume(int dim, int vertex, bool useDual, int limit);
	vector<vector<string> > getRowsToFindVolumeGivenSpecficFile(const char *polymakeFile, bool useDual);

	vector<vector<ValuationDBStatistics> > getStatistics(bool useDual);//main printing function.
	ValuationDBStatistics getStatisticsByDimAndPointCount(int dim, int additionalPoints, bool useDual);
	void getStatistics(const vector<double> &data, double &avg, double &min, double &max, double &sd, int &totalFinished, int &totalExist, bool &manuallyLimited);

	void updateVolumeTimeAndValue(AlgorithemUsed alg, double time, RationalNTL theComputedIntegral, string previousValue, string rowid);
	int volumeTestsCompleted(AlgorithemUsed alg, int dim, int vertex, bool useDual);




};

#endif /* VOLUMEDB_H_ */
