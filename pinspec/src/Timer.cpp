#include "Timer.h"


/**
 * @brief Timer class constructor.
 * @details Sets the default elapsed time to 0.0 and the stopwatch
 *          running attribute to false.
 */
Timer::Timer() {
    this->_running = false;
    this->_elapsed_time = 0.0;
}


/**
 * @brief Default Timer destructor
 */
Timer::~Timer() { }


/**
 * @brief Starts the Timer in a similar fashion to starting a stopwatch.
 */
void Timer::start() {

    if (!this->_running) {
        #ifdef __MACH__ 	                /* For OSX */
        gettimeofday(&this->_start_time, NULL);
        #else   				/* For linux */
        clock_gettime(CLOCK_MONOTONIC, &this->_start_time);
        #endif
        this->_running = true;
    }

    return;
}


/**
 * @brief Stops the timer in a similar fashion to stopping a stopwatch.
 */
void Timer::stop() {
    if (this->_running) {
        #ifdef __MACH__             /* For OSX */
        gettimeofday(&this->_end_time, NULL);
        #else			/* For linux */
        clock_gettime(CLOCK_MONOTONIC, &this->_end_time);
        #endif
        this->_running = false;
	this->_elapsed_time += this->diff(this->_start_time, this->_end_time);
    }

    return;
}


/**
 * @brief Resets the timer in a similar fashion to resetting a stopwatch.
 */
void Timer::reset() {
    this->_elapsed_time = 0;
    this->_running = false;
}


/**
 * @brief Restarts the timer
 * @details The elapsed time will accumulate along with the previous time(s) 
 *          the timer was running. If the timer was already running, this
 *          function does nothing.
 */
void Timer::restart() {
    if (!this->_running) {
        this->_elapsed_time += this->diff(this->_start_time, this->_end_time);
	this->start();
    }
}


/**
 * @brief Records a message corresponding to a given time recorded by the timer.
 * @details When this method is called it assumes that the timer has been 
            stopped and has the current time for the process corresponding to 
            the message.
 * @param msg a msg corresponding to this timer split
 */
void Timer::recordSplit(const char* msg) {
    _timer_splits.push_back(std::pair<double, const char*>(getTime(), msg));
}



/**
 * @brief Returns the amount of time elapsed from start to stop of the timer. 
 * @details If the timer is currently runnning, returns the time from the timer 
 *          start to the present time.
 * @return the elapsed time
 */
double Timer::getTime() {
    /* If the timer is not running */
    if (!this->_running) {
        #ifdef __MACH__		        /* For OSX */
        return this->_elapsed_time * 1.0E-6;
         #else				/* For Linux */
        return this->_elapsed_time * 1.0E-9;
        #endif
    }

    /* If the timer is currently running */
    else {
        timespec temp;
        #ifdef __MACH__ 	        /* For OSX */
	gettimeofday(&temp, NULL);
	#else				/* For Linux */
	clock_gettime(CLOCK_MONOTONIC, &temp);
	#endif

        this->_elapsed_time += this->diff(this->_start_time, temp);
        #ifdef __MACH__		        /* For OSX */
        return this->_elapsed_time * 1.0E-6;
	#else				/* For Linux */
        return this->_elapsed_time * 1.0E-9;
        #endif
    }
}


/**
 * @brief Helper function which computes the time between the values of two
 *        timespec structs.
 * @param start the start time timespec struct
 * @param end the end time timespec struct
 */
double Timer::diff(timespec start, timespec end) {
    double sec, delta;
    #ifdef __MACH__
    double usec;
    delta = end.tv_usec - start.tv_usec;
    #else
    double nsec;
    delta = end.tv_nsec - start.tv_nsec;
    #endif

    if (delta < 0) {
        sec = end.tv_sec - start.tv_sec;
	#ifdef __MACH__
	usec = 1.0E6 + delta;
	#else
	nsec = 1.0E9 + delta;
	#endif

    } 
    else {
        sec = end.tv_sec - start.tv_sec;
	#ifdef __MACH__
        usec = delta;
        #else
        nsec = delta;
        #endif
    }

    #ifdef __MACH__
    return (sec * 1.0E6 + usec);
    #else
    return(sec*1.0E9 + nsec);
    #endif
}


/**
 * @brief Prints the Timer's splits to the console.
 * @details This method will loop through all of the Timer's splits and print a
 *          formatted message string (80 characters in length) to the console
 *          with the message and the time corresponding to that message
 */
void Timer::printSplits() {

    const char* curr_msg;
    double curr_split;
    int msg_length, num_whitespaces;

    for (int i=0; i < (int)_timer_splits.size(); i++) {
        std::stringstream formatted_msg;

        curr_split = _timer_splits.at(i).first;
        curr_msg = _timer_splits.at(i).second;
        msg_length = strlen(curr_msg);

        /* Num whitespaces for formatting is:
         * 80 max char - 13 char for logger - 13 for time - msg length */
        num_whitespaces = 80 - 13 - 11 - msg_length -3;

        formatted_msg << curr_msg;

        /* Add periods to format message to 80 characters length */
        for (int i=0; i < num_whitespaces; i++)
            formatted_msg << ".";

        log_printf(RESULT, "%s%.7f sec", formatted_msg.str().c_str(), 
		   curr_split);
    }
}
