#pragma once

#include <cassert>
#include <chrono>
#include <string>

struct Timer {
  using Clock = std::chrono::system_clock;
  using TimePoint = Clock::time_point;
  using Nano = std::chrono::nanoseconds;

  bool isMeasuring;

  Timer() : isMeasuring(false) {}

  void start() {
    assert(!isMeasuring);
    _start = Clock::now();
    isMeasuring = true;
  }

  void stop() {
    assert(isMeasuring);
    _end = Clock::now();
    sum += std::chrono::duration_cast<Nano>(_end - _start).count();
    isMeasuring = false;
  }

  double sec() {
    if (isMeasuring) stop(), start();
    return sum / 1e9;
  }

 private:
  TimePoint _start, _end;
  long long sum = 0;
};