/****************************************************************
 * Author: Ruriko Yoshida
 * This code install NTL and latte
 * Date: Oct 22nd, 2003
 * Update: Oct 22nd, 2003
*****************************************************************/
#include <iostream.h>
#include <fstream.h>
#include <ctype.h>
#include <string>
#include <stdlib.h> // exit()

using namespace std;

int main(int argc, char* argv[] ){
  system("mkdir LattE");
  //  system("mv latte LattE");
  system("mv cdd LattE");
  system("mv lrs1 LattE");
  system("mv redcheck_gmp LattE");
  system("mv ComputeAdjacency LattE");
  system("mv maximize LattE"); 
  system("mv minimize LattE");
  system("mv count LattE");
  system("mv ehrhart LattE");
  return 0;
}



