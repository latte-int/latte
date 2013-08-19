// This is a -*- C++ -*- header file.

/* latte_system.h -- System interface functions

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

#ifndef LATTE_SYSTEM_H
#define LATTE_SYSTEM_H

#include <string>

// Execute COMMAND with the system shell.
// When COMMAND returns with a nonzero error status,
// report and exit the program.
void system_with_error_check(const char *command);
void system_with_error_check(const std::string &command);

// Quote argument for the shell.
std::string shell_quote(const std::string &argument);

// The rename(2) system call with error checking and C++ strings as arguments.
void
rename_with_error_check(const std::string &old_name, const std::string &new_name);

// Functions for storing intermediate data in a secure temporary
// directory.  (More importantly, different runs of LattE are isolated
// from each other.)  The  function `temporary_directory_name' returns
// the directory name with a trailing slash.
void create_temporary_directory();
const std::string &temporary_directory_name();
void remove_temporary_directory();

// Return the pathname of the file NAME in the temporary directory.
std::string temporary_file_name(const std::string &name);

#endif
