/* count.cpp -- Master program

 Copyright 2002, 2003 Raymond Hemmecke, Ruriko Yoshida
 Copyright 2006, 2007 Matthias Koeppe

 This file is part of LattE.

 LattE is free software; you can redistribute it and/or modify it
 under the terms of the version 2 of the GNU General Public License
 as published by the Free Software Foundation.

 LattE is distributed in the hope that it will be useful, but
 WITHOUT ANY WARRANTY; without even the implied warranty of
 MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 General Public License for more details.

 You should have received a copy of the GNU General Public License
 along with LattE; if not, write to the Free Software Foundation,
 Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307 USA.
 */
#include "count.h"

int main(int argc, char *argv[])
{

	CountAnswerContainer ans;
	cout << "going to call main" << endl;
	ans = mainCountDriver(argc,  argv);

	/*
	cout <<"*%*%*%*%*%**%**%*" << endl;
	cout << "ehrhart_coefficients\n";
	for(int i = 0; i < ans.ehrhart_coefficients.size(); ++i)
		cout << "  " << ans.ehrhart_coefficients[i] << "t^" << i << endl;
	cout << "gen fun name" << ans.multivariateGenFunctionFileName.c_str() << endl;
	cout << "num lattice pt" << ans.numLaticePoints << endl;
	cout << "seriesExpansion" << endl;
	for(int i = 0; i < ans.seriesExpansion.length(); ++i)
		cout << "  " << ans.seriesExpansion[i] << "t^" << i << endl;
	*/

}//main
