#include <sys/times.h>
#include <cassert>
#include <unistd.h>
#include "timing.h"
#include <iostream>

using namespace std;

clock_t Timer::ticks_per_second;

Timer::Timer(const std::string &a_name, bool start_timer) :
  name(a_name), ticks_elapsed(0), started(false)
{
  ticks_per_second = sysconf(_SC_CLK_TCK);
  if (start_timer) start();
}

void Timer::start()
{
  assert(!started);
  struct tms buf;
  clock_t t = times(&buf);
  assert(t != -1);
  start_time = buf.tms_utime + buf.tms_stime
    + buf.tms_cutime + buf.tms_cstime;
  started = true;
}

void Timer::stop()
{
  assert(started);
  struct tms buf;
  clock_t t = times(&buf);
  assert(t != -1);
  ticks_elapsed += (buf.tms_utime + buf.tms_stime
		    + buf.tms_cutime + buf.tms_cstime)
    - start_time;
  started = false;
}

float Timer::get_seconds() const
{
  assert(!started);
  return ticks_elapsed / (1.0 * ticks_per_second);
}

ostream &operator<< (ostream &s, const Timer &timer)
{
  s << timer.name << ": " << timer.get_seconds() << " sec" << endl;
  return s;
}

