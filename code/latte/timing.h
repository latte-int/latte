// This is a -*- C++ -*- header file.

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
