// This is a -*- C++ -*- header file.

/* timing.h -- Measure computation times

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

#ifndef TIMING_H
#define TIMING_H

#include <string>

class Timer {
  std::string name;
  static clock_t ticks_per_second;
  clock_t ticks_elapsed;
  clock_t start_time;
  bool started;
public:
  Timer(const std::string &a_name, bool start_timer = false);
  void start();
  void stop();
  float get_seconds() const;
  friend std::ostream &operator<< (std::ostream &s, const Timer &timer);
};

#endif
