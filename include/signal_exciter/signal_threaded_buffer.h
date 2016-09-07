#ifndef INCLUDED_SIGNAL_THREADED_BUFFER
#define INCLUDED_SIGNAL_THREADED_BUFFER

//functional stuff
#include <vector>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>
#include <math.h>
#include <complex>
#include <stdio.h>


typedef std::complex<float> complexf;
typedef std::complex<double> complexd;
typedef std::complex<long double> complexl;

//thread stuff
#include <boost/thread/thread.hpp>
#include <boost/thread/mutex.hpp>
#include <boost/thread/locks.hpp>
#include <boost/thread/condition_variable.hpp>


typedef boost::thread                     bthread;
typedef boost::mutex                      bmutex;
typedef boost::unique_lock<bmutex>        block;
typedef boost::condition_variable         bcondition_variable;

template <class T>
class signal_threaded_buffer
{
  private:
    size_t d_sys_page_size;
    size_t d_max_items;
    size_t d_item_count;
    size_t d_item_size;
    size_t d_base_pointer;
    size_t d_read__at;
    size_t d_write_at;
    std::vector<T> d_buffer;

    size_t d_min_notify;

    bmutex              d_mutex;
    bcondition_variable d_cond_read;
    bcondition_variable d_cond_write;

    size_t write_room();
    size_t read_room();

    bool check_notify_read();
    bool check_notify_write();
    void read_wait();
    void write_wait();
    

  public:
    signal_threaded_buffer(size_t max_buffer_size = 8192, size_t notify_size = 512);
    ~signal_threaded_buffer();

    size_t size();
    size_t bmemcpy( T* buff, size_t item_count, bool direction );
    size_t bmemcpy( signal_threaded_buffer<T> buff, size_t item_count, bool direction );

    size_t write_item( T item );
    size_t read_item( T& item );

    size_t readable();
    size_t writeable();

    void lock();
    void unlock();

    //void cleanup();
    

};



#endif //INCLUDED_SIGNAL_THREADED_BUFFER
