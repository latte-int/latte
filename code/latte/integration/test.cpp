#include <stdlib.h>
#include <fstream>
#include <iostream>

using namespace std;

#include "multipoly.h"

int main()
{   
    int k,n;
    string poly;
    ifstream infile ("multipoly.txt");
    infile>>poly;
    multipoly m1(poly);
    m1.print();
    infile.close();
    return 0;
}
