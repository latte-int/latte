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
  char directory[127], makeCommand[250], ntlInstall[127], NTLsrc[127], NTLdirectory[127], NTLlibDirectory[127], NTLlibinclude[127], NTLliblib[127], LattEdirectory[127], command[127];
  string tmpString;
  ifstream in("../../latte_directory");
  in.getline(directory, 127);

  strcpy(NTLdirectory, directory);
  strcat(NTLdirectory, "/NTL");

  strcpy(NTLlibDirectory, NTLdirectory);
  strcat(NTLlibDirectory, "/NTL");

  strcpy(NTLlibinclude, NTLlibDirectory);
  strcat(NTLlibinclude, "/include");

  strcpy(NTLsrc, directory);
  strcat(NTLsrc, "/ntl-5.3.1/src/");

  strcpy(NTLliblib, NTLlibDirectory);
  strcat(NTLliblib, "/lib");

  strcpy(LattEdirectory, directory);
  strcat(LattEdirectory, "/latte");

  strcpy(command, "./configure PREFIX=");
  strcat(command, NTLlibDirectory);
  system(command);
  system("make");
  system("make install");
  
  return 0;
}
