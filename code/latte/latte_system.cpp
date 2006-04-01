#include <cstdlib>
#include <iostream>
using namespace std;

#include "latte_system.h"

void system_with_error_check(const char *command)
{
  int status = system(command);
  if (status != 0) {
    cerr << "Command '" << command << "' returned with exit status "
	 << status << "." << endl;
    exit(1);
  }
}
