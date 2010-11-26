/*
 * integratePolynomialForPaperDriver.cpp
 *
 *  Created on: Nov 25, 2010
 *      Author: Brandon Dutra and Gregory Pinto
 *
 *  Given a list of directories that contain a bunch of .latte files,
 *  this dirver will find them and integretrate a random monomial over each
 *  and save stats in a log file.
 *
 *  This is used to generate  data for a table in a research paper.
 */

#include "VolumeAndIntegrationTests.h"
#include <string>

using namespace std;

int main(int argc, char *argv[])
{

	int degree[] =
	{ 1, 2, 5, 10, 20, 30, 40, 50, 100, 200, 300 };
	int numberDegrees = 11;

	if (argc == 1)
	{
		cout << "Usage: " << argv[0]
				<< " dim1 dim2 ... where dimi contains *.latte files of dim i" << endl;
		exit(1);
	}

	for (int deg = 0; deg < numberDegrees; ++deg)
	{
		for (int i = 1; i < argc; ++i)
		{
			vector<string> files;
			stringstream logFileName;

			logFileName << argv[i] << "/log_polynomialdDegree_" << degree[deg] << ".log";

			IntegrationPaper::findAllLatteFilesInDirectory(string(argv[i])+'/',
					files);

			IntegrationPaper::integrateFiles(logFileName.str(), files, atoi(argv[i]+3), degree[deg]);

		}
	}

}
