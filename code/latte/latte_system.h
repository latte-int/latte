// This is a -*- C++ -*- header file.
// System interface functions

#ifndef LATTE_SYSTEM_H
#define LATTE_SYSTEM_H

// Execute COMMAND with the system shell.
// When COMMAND returns with a nonzero error status,
// report and exit the program.
void system_with_error_check(const char *command);

#endif
