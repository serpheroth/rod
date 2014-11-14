//**************************************************************************************
// Copyright 2004 Huamin Wang.
//**************************************************************************************
// TIMER classes
//**************************************************************************************
#ifndef __TIMER_H__
#define __TIMER_H__
#include <iostream>
#include <sys/timeb.h>
#include <time.h>

class Timer {
public:
  struct timeb start_time;

  Timer() {Start();}

  ~Timer() {}

  void Start()
  {ftime( &start_time);}

  float Now() {
    struct timeb current_time;
    ftime(&current_time);
    return (float) current_time.time + 0.001f * current_time.millitm;
  }

  float GetTime() {
    struct timeb current_time;
    ftime( &current_time);
    return (float)(current_time.time - start_time.time) + 0.001f * (current_time.millitm - start_time.millitm);
  }
};

#ifdef __linux
bool operator<(const timespec& time1, const timespec& time2);
bool operator==(const timespec& time1, const timespec& time2);
bool operator>(const timespec& time1, const timespec& time2);
const timespec operator+(const timespec& time, const timespec& increment);
const timespec operator-(const timespec& end, const timespec& start);
float operator/(const timespec& dividend, const timespec& divisor);
unsigned long long GetNanoSeconds(const timespec& time);
const timespec GetCurrentTimeSpec(clockid_t clk_id = CLOCK_REALTIME);
std::ostream& operator<<(std::ostream& out, const timespec& time);
class AccurateTimer {
public:
  // CLOCK_REALTIME, a system-wide realtime clock.
  // CLOCK_PROCESS_CPUTIME_ID, high-resolution timer provided by the CPU for each process.
  // CLOCK_THREAD_CPUTIME_ID, high-resolution timer provided by the CPU for each of the threads.
  AccurateTimer(clockid_t clk_id = CLOCK_REALTIME) : clk_id_ (clk_id) {
    Start(clk_id);
  }
  void Start(clockid_t clk_id = CLOCK_REALTIME) {
    clk_id_ = clk_id;
    clock_gettime(clk_id_, &start_time_);
  }

  float GetTime() {
    timespec end_time;
    clock_gettime(clk_id_, &end_time);
    end_time = end_time - start_time_;
    Start(clk_id_);
    return end_time.tv_sec + end_time.tv_nsec * 1e-9;
  }
private:
  timespec start_time_;
  clockid_t clk_id_;
};
#endif

#endif //__TIMER_H__
