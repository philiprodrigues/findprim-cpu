#ifndef TIMER_H
#define TIMER_H

#include <time.h>

//======================================================================
// record wall time using clock_gettime()
class Timer
{
public:
  Timer()
  {
    reset();
  }

  // reset the clock
  void reset()
  {
    clock_gettime(CLOCK_MONOTONIC, &tp_start);
  }

  // Get the number of microseconds since the start, or last reset()
  long long getTime()
  {
    timespec tp_end;
    clock_gettime(CLOCK_MONOTONIC, &tp_end);
    const long long billion=1000000000;
    long long end=billion*tp_end.tv_sec+tp_end.tv_nsec;
    long long start=billion*tp_start.tv_sec+tp_start.tv_nsec;
    return (end-start)/1000;
  }

private:
  timespec tp_start;
};


#endif 
