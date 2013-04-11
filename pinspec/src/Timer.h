/**
 * @file Timer.h
 * @brief The Timer static class.
 * @author William Boyd (wboyd@mit.edu)
 * @date January 2, 2012
 */

#ifndef TIMER_H_
#define TIMER_H_

#ifdef __cplusplus
#include <time.h>
#include <sys/time.h>
#include <iostream>
#include <sstream>
#include <iomanip>
#include <utility>
#include <vector>
#include "log.h"
#endif

#ifdef __MACH__		/* For OSX */
#define timespec timeval
#endif


/**
 * @class Timer Timer.h "pinspec/src/Timer.h" 
 * @brief The Timer represents a simulation stopwatch.
 * @details The timer class is for profiling code. It outputs running time in
 *          seconds but has resolution of microseconds on OSX and nanoseconds
 *          on Linux machines.
 */
class Timer {
protected:
    /** The start time of a timer split */
    timespec _start_time;
    /** The end time of a timer split */
    timespec _end_time;
    /** The elapsed time between start and stop for a split */
    double _elapsed_time;
    /** Boolean representing whether or not the time is presently running */
    bool _running;
    /** A container of times corresponding to messages for each split */
    std::vector< std::pair<double, const char*> > _timer_splits;
public:
    Timer();
    virtual ~Timer();
    void start();
    void stop();
    void restart();
    void reset();
    void recordSplit(const char* msg);
    double getTime();
    double diff(timespec start, timespec end);
    void printSplits();
};

#endif /* TIMER_H_ */
