/****************************************************************
 * Author: Ruriko Yoshida
 * This code install CDD lib
 * Date: July 2nd, 2004
 * Update: July 2nd, 2004
*****************************************************************/
#include <iostream.h>
#include <fstream.h>
#include <ctype.h>
#include <string>
#include <stdlib.h> // exit()


int main(int argc, char* argv[] ){
  char directory[127], cddDirectory[400];
  system("pwd > latte_directory");
  ifstream in("latte_directory");
  in.getline(directory, 127);

  strcpy(cddDirectory, "./configure --prefix=");
  strcat(cddDirectory, directory);
  strcat(cddDirectory, "/CDD");

  system(cddDirectory);
    return 0;
}
