/****************************************************************
 * Author: Ruriko Yoshida
 * This code install NTL and latte
 * Date: Oct 22nd, 2003
 * Update: Oct 22nd, 2003
*****************************************************************/
#include <iostream.h>
#include <fstream.h>
#include <ctype.h>
#include <string.h>
#include <stdio.h>
#include <stdlib.h> // exit()

using namespace std;

int main(int argc, char* argv[] ){
  string tmpString;
  char directory[127], makeCommand[250], ntlInstall[127], NTLdirectory[127], NTLlibDirectory[127], NTLlibinclude[127], NTLliblib[127], LattEdirectory[127], command[127];
  // string tmpString;
  system("pwd > latte_directory");
  ifstream in("latte_directory");
  in.getline(directory, 127);
  //system("rm latte_directory");
  //cout << directory << endl;

  strcpy(NTLdirectory, directory);
  strcat(NTLdirectory, "/NTL");

  strcpy(NTLlibDirectory, NTLdirectory);
  strcat(NTLlibDirectory, "/NTL");

  strcpy(NTLlibinclude, NTLlibDirectory);
  strcat(NTLlibinclude, "/include");

  strcpy(NTLliblib, NTLlibDirectory);
  strcat(NTLliblib, "/lib");

  strcpy(LattEdirectory, directory);
  strcat(LattEdirectory, "/latte");

  strcpy(command, "mkdir ");
  strcat(command, NTLdirectory);
  system(command);

  strcpy(makeCommand, "CFLAGS = -Wall -Wno-deprecated -g -DLINUX -I/");
  strcat(makeCommand, NTLlibinclude);
  strcat(makeCommand, "\n\nLib_Dirs = -L/");
  strcat(makeCommand, NTLliblib);
  strcat(makeCommand, "\n");

  // system("mv ntl-5.3.1.tar.gz NTL/");

  system("gunzip ntl-5.3.1.tar.gz");
  system("tar xvf ntl-5.3.1.tar");
  system("mv ntl-5.3.1.tar NTL/");

 system("mv NTLinstall ntl-5.3.1/src/");

/* strcpy(command, "./ntl-5.3.1/src/NTLinstall ");
  strcat(command, directory);
  system(command);       */


  system("gunzip latte.tar.gz");
  system("tar xvf latte.tar");

  ifstream in2("Makefile.install");
  if(!in2){
    cerr << "Need Makefile.install file."<< endl;
    exit(1);
  }
  ofstream out1("Makefile");


  while(tmpString[0] != 'L'){
  getline(in2, tmpString);   //cout << tmpString << endl;
  out1 << tmpString << endl; }

  out1 << makeCommand << endl;


  while(!in2.eof()){
      getline(in2, tmpString);
      out1 << tmpString << endl;;
  }


  ifstream in3("subMakefile.install");
  if(!in3){
    cerr << "Need subMakefile.install file." << endl;
    exit(2);
  }
  ofstream out2("subMakefile");
  while(tmpString[0] != 'L'){
  getline(in3, tmpString); //  cout << tmpString << endl;
  out2 << tmpString << endl; }

  out2 << makeCommand << endl;


  while(!in3.eof()){
      getline(in3, tmpString);
      out2 << tmpString << endl;;
  }

  system("cp subMakefile barvinok/Makefile");
  system("cp subMakefile vertices/Makefile");
  system("mv subMakefile genFunction/Makefile");
  //system("make");

  return 0;
}



