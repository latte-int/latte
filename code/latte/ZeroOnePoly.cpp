#include <iostream.h>
#include <fstream.h>
#include <ctype.h>
#include <cstring>
#include <stdlib.h> // exit()
#include <string>
#include <math.h>

void AddLinearity(){
  ifstream in("test.lat");
  ofstream out("test.lat.lat");
  string tmpString;
  int dim, numOfConst;
  in >> numOfConst >> dim;
  out << numOfConst << " " << dim << endl;
  getline(in , tmpString);
  for(int i = 0; i < numOfConst; i++){
    getline(in , tmpString);
    out << tmpString << endl;
  }
  out << "linearity 1 1" << endl;
}
void RandamGenerator(ofstream & out, int dim, int seed, int num){
  srand(seed + num);
  ofstream out2("test.ext");
  out2 << "V-representation" << endl << "begin" << endl;
  out2 << 3*dim << " " << dim+1 << " integer" << endl;
  for(int i = 0; i < 3*dim; i++){
    out2 << 1 << " ";
    for(int j = 0; j < dim; j++){
      int tmp = rand() % 2;
      out2 << tmp << " ";
      out << tmp << " ";
    }
    out << endl;
    out2 << endl;
    //  out2 <<"hull" <<endl;
  }
  out2 <<"end"<< endl;
  out2 <<"hull" <<endl;

} 

void ReadExt(){
  ifstream in2("test.ine");
  ifstream in("test.ine");
  ofstream out("test.tmp");
  string tmpString, save;

  save = " ";
  int dim, numOfConst = 0, tmp;
  while (tmpString!="begin"){
    getline(in2,tmpString);
    if(tmpString[0] == 'l') save = tmpString;
  }
  while (tmpString!="end"){
    getline(in2,tmpString);
    numOfConst++;
  }
  numOfConst = numOfConst - 2;

  while (tmpString!="begin") getline(in,tmpString);
  in >> tmpString >> dim >> tmpString;
  out << numOfConst << " " << dim << endl;
  for(int i = 0; i < numOfConst; i++){
    for(int j = 0; j < dim; j++){
      in >> tmp;
      out << tmp << " ";
    }
    //    out << "linearity 1 1" << endl;
    out << endl;
  }
}


void ReadLatte(int numOfDeg){
  ifstream in("test.out");
  ofstream out("MapleIn");
  string tmpString;
  int tmp;
  out << "data:=[";
  while (tmpString!="1") getline(in,tmpString);
  for(int i = 0; i < numOfDeg; i++){
    in >> tmp;
    out << "[" << i + 1<< ", " << tmp << "]";
    if(i != (numOfDeg - 1)) out << ",";
  }
  out << "];" << endl;
  out << "with(CurveFitting):" << endl;
  out << "cubfit:=PolynomialInterpolation(data,x);" << endl;
}

void ReadMaple(ofstream & out){
  ifstream in("MapleOut");
  string tmpString;
  while(tmpString != "> with(CurveFitting):")     getline(in,tmpString);

  while(tmpString != "> quit"){
    getline(in,tmpString);
    out << tmpString << endl;
  }

}

int main(int argc, char *argv[]) {
  int dim, seed, numOfDeg, numOfLoops;
  char command[127];
  dim = atoi(argv[1]);
  numOfLoops = atoi(argv[2]);
  seed = atoi(argv[3]);
  numOfDeg = atoi(argv[4]);
  srand(seed);

  ofstream out("Result");
  for(int i = 0; i < numOfLoops; i++){
    out << "Loop: " << i+1 << endl;
    cout << "Loop: " << i+1 << endl;
    RandamGenerator(out, dim, seed, i);
    system("./lrs1 test.ext > test.ine");
    ReadExt();
    system("./17gon test.tmp test.lat");
    AddLinearity();
    strcpy(command, "./ehrhart ");
    strcat(command, argv[4]);
    strcat(command, " test.lat.lat> test.out");
    system(command);
    ReadLatte(numOfDeg);
    system("maple <MapleIn > MapleOut");
    ReadMaple(out);
    system("rm test.*");
  }

  return 0;
}









