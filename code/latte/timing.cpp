#include <stdio.h>
#include <sys/times.h>
#include <time.h>
using namespace std;

extern float getcputime(void);

float getcputime (void) {
  struct tms now;
  int sec, hun;
  float total;

  /* returns cpu-time in seconds */

#ifndef CLK_TCK
  long clocks_per_second = 0;
  
//  clocks_per_second = sysconf( _SC_CLK_TCK );
  times( &now );
  
  hun  = ((now.tms_utime % clocks_per_second)*100)/clocks_per_second;
  sec = (now.tms_utime / clocks_per_second);
#else
  times( &now );
  
  hun  = ( ( now.tms_utime % CLK_TCK ) * 100 ) / CLK_TCK;
  sec = ( now.tms_utime / CLK_TCK );
#endif
  
  total = (double)sec + ( (double)(hun) / 100.0 );
  
  return (total); 
}
