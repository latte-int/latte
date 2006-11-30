/* latte_system.cpp -- Interface to system utilities

   Copyright 2006 Matthias Koeppe

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

#include <cstdlib>
#include <iostream>
#include <cassert>
#include <sys/stat.h>
#include <sys/types.h>
#include <unistd.h>
#include <errno.h>

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

void system_with_error_check(const string &command)
{
  system_with_error_check(command.c_str());
}


static bool created_temp_dir = false;
static string temp_dir;

#ifndef HOST_NAME_MAX
#define HOST_NAME_MAX 1023
#endif

void create_temporary_directory()
{
  char hostname[HOST_NAME_MAX + 1];
  if (gethostname(hostname, HOST_NAME_MAX + 1) != 0) {
    perror("create_temporary_directory: gethostname failed");
    exit(1);
  }
  pid_t pid = getpid();
  int i;
  char t[PATH_MAX];
  for (i = 0; i<INT_MAX; i++) {
    sprintf(t, "/tmp/latte-%d@%s-%d", pid, hostname, i);
    int result = mkdir(t, 0700);
    if (result == 0) {
      created_temp_dir = true;
      temp_dir = t;
      temp_dir += "/";
      return;
    }
    else {
      if (errno = EEXIST) continue;
      perror("create_temporary_directory: mkdir failed");
      exit(1);
    }
  }
  cerr << "create_temporary_directory: Unable to create a fresh directory" << endl;
  exit(1);
}
  
const string &temporary_directory_name()
{
  if (!created_temp_dir)
    create_temporary_directory();
  return temp_dir;
}

void remove_temporary_directory()
{
  if (created_temp_dir) {
    char command[PATH_MAX + 10];
    sprintf(command, "rm -rf %s", temp_dir.c_str());
    system_with_error_check(command);
  }
}

string temporary_file_name(const string &name)
{
  return temporary_directory_name() + name;
}
