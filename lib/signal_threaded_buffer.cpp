
#include <signal_exciter/signal_threaded_buffer.h>


template <class T>
signal_threaded_buffer<T>::signal_threaded_buffer(size_t max_buffer_size, size_t notify_size)
  : d_item_count(0),
    d_item_size(0),
    d_read__at(0),
    d_write_at(0),
    d_min_notify(notify_size)
{
  d_item_size = sizeof(T);
  long sz = sysconf(_SC_PAGESIZE);
  d_sys_page_size = (unsigned long) sz;
  d_max_items = (unsigned long) (sz*ceil(float(max_buffer_size)*float(d_item_size)/float(sz))/float(d_item_size));
  d_buffer = std::vector<T>(d_max_items*2,T(0));
  //d_min_notify = 1;
}

template <>
signal_threaded_buffer<complexf>::signal_threaded_buffer(size_t max_buffer_size, size_t notify_size)
  : d_item_count(0),
    d_item_size(0),
    d_read__at(0),
    d_write_at(0),
    d_min_notify(notify_size)
{
  d_item_size = sizeof(complexf);
  long sz = sysconf(_SC_PAGESIZE);
  d_sys_page_size = (unsigned long) sz;
  d_max_items = (unsigned long) (sz*ceil(float(max_buffer_size)*float(d_item_size)/float(sz))/float(d_item_size));
  d_buffer = std::vector<complexf>(d_max_items*2,complexf(0.,0.));
  //d_min_notify = 1;
}

template <>
signal_threaded_buffer<complexd>::signal_threaded_buffer(size_t max_buffer_size, size_t notify_size)
  : d_item_count(0),
    d_item_size(0),
    d_read__at(0),
    d_write_at(0),
    d_min_notify(notify_size)
{
  d_item_size = sizeof(complexd);
  long sz = sysconf(_SC_PAGESIZE);
  d_sys_page_size = (unsigned long) sz;
  d_max_items = (unsigned long) (sz*ceil(float(max_buffer_size)*float(d_item_size)/float(sz))/float(d_item_size));
  d_buffer = std::vector<complexd>(d_max_items*2,complexd(0.,0.));
  //d_min_notify = 1;
}

template <>
signal_threaded_buffer<complexl>::signal_threaded_buffer(size_t max_buffer_size, size_t notify_size)
  : d_item_count(0),
    d_item_size(0),
    d_read__at(0),
    d_write_at(0),
    d_base_pointer(0),
    d_min_notify(notify_size)
{
  d_item_size = sizeof(complexl);
  long sz = sysconf(_SC_PAGESIZE);
  d_sys_page_size = (unsigned long) sz;
  d_max_items = (unsigned long) (sz*ceil(float(max_buffer_size)*float(d_item_size)/float(sz))/float(d_item_size));
  d_buffer = std::vector<complexl>(d_max_items,complexl(0.,0.));
  //d_min_notify = 1;
}

template <class T>
signal_threaded_buffer<T>::~signal_threaded_buffer()
{
}

template <class T>
size_t
signal_threaded_buffer<T>::size()
{
  return (d_buffer.size());
}


template <class T>
size_t
signal_threaded_buffer<T>::bmemcpy( T* buff, size_t item_count, bool direction )
{
  size_t transfered(0);
  lock();
  if(direction){//Write into (*this)
    size_t size = d_buffer.size();
    size_t write_pointer = d_write_at;
    size_t read_pointer = 0;
    size_t copy_count = 0;
    if(item_count <= write_room()){//There is enough space in the buffer to write everything requested
      size_t use_item_count = 0;
      if(write_pointer + item_count <= size){
        //Everything can be copied in 1 go.
        use_item_count = item_count;
        copy_count = 1;
      }
      else{
        //Everything can be copied in 2 goes.
        use_item_count = size-write_pointer;
        copy_count = 2;
      }
      while(copy_count){
        memcpy( &d_buffer[write_pointer], &buff[read_pointer], use_item_count*d_item_size );
        copy_count--;
        if(copy_count){
          write_pointer = 0;
          read_pointer += use_item_count;
          use_item_count = item_count - use_item_count;
        }
      }
      d_item_count += item_count;
      transfered += item_count;
      d_write_at = (write_pointer + use_item_count)%size;
    }
    else{//There is not enough space in the buffer to write everything requested
      size_t new_item_count = write_room();
      size_t use_item_count = 0;
      if(write_pointer + new_item_count <= size){
        //Everything can be copied in 1 go.
        use_item_count = new_item_count;
        copy_count = 1;
      }
      else{
        //Everything can be copied in 2 goes.
        use_item_count = size - write_pointer;
        copy_count = 2;
      }
      while(copy_count){
        memcpy( &d_buffer[write_pointer], &buff[read_pointer], use_item_count*d_item_size );
        copy_count--;
        if(copy_count){
          write_pointer = 0;
          read_pointer += use_item_count;
          use_item_count = new_item_count - use_item_count;
        }
      }
      d_item_count += new_item_count;
      transfered += new_item_count;
      d_write_at = (write_pointer + use_item_count)%size;
    }
  }
  else{//Read from (*this)
    size_t size = d_buffer.size();
    size_t write_pointer = 0;
    size_t read_pointer = d_read__at;
    size_t copy_count = 0;
    if(item_count <= read_room()){//There are enough items in the buffer to read everything requested
      size_t use_item_count = 0;
      if(read_pointer + item_count <= size){
        //Everything can be copied in 1 go.
        use_item_count = item_count;
        copy_count = 1;
      }
      else{
        //Need two goes to copy everything
        use_item_count = size - read_pointer;
        copy_count = 2;
      }
      while(copy_count){
        memcpy( &buff[write_pointer], &d_buffer[read_pointer], use_item_count*d_item_size );
        copy_count--;
        if(copy_count){
          read_pointer = 0;
          write_pointer += use_item_count;
          use_item_count = item_count - use_item_count;
        }
      }
      d_item_count -= item_count;
      transfered += item_count;
      d_read__at = (read_pointer + use_item_count)%size;
    }
    else{//There are not enough items for the requested read amount
      size_t new_item_count = read_room();
      size_t use_item_count = 0;
      if(read_pointer + new_item_count <= size){
        //Everything can be copied in 1 go.
        use_item_count = new_item_count;
        copy_count = 1;
      }
      else{
        //Need two goes to copy everything
        use_item_count = size - read_pointer;
        copy_count = 2;
      }
      while(copy_count){
        memcpy( &buff[write_pointer], &d_buffer[read_pointer], use_item_count*d_item_size );
        copy_count--;
        if(copy_count){
          read_pointer = 0;
          write_pointer += use_item_count;
          use_item_count = new_item_count - use_item_count;
        }
      }
      d_item_count -= new_item_count;
      transfered += new_item_count;
      d_read__at = (read_pointer + use_item_count)%size;
    }
  }
  unlock();
  return transfered;
}


template <class T>
size_t
signal_threaded_buffer<T>::bmemcpy( signal_threaded_buffer<T> buff, size_t item_count, bool direction )
{
  size_t using_count(0);
//FIXME: Come back to this later
/*  if(direction){//write to (this) buffer, read from (buff) buffer
    buff.lock();
    size_t twr = write_room();

    using_count = twr > item_count ? item_count : twr;
    if(using_count+d_write_at > d_max_items){//want to pull over more items than the buffer size will allow, but can wrap the rest
      size_t ender = d_max_items - d_write_at;//How many (this) can accept before wrap
      size_t begin = using_count - ender;//The remaining (this) will accept
      using_count = buff.bmemcpy( &d_buffer[d_write_at], ender, false );
      if(using_count < ender){//Didn't have enough data to read everything, can return now
        d_write_at += using_count;
      }
      else{//There was enough to read, continue reading
        size_t secondary = buff.bmemcpy( &d_buffer[0], begin, false );
        d_write_at = secondary;
        using_count += secondary;
      }
    }
    else{//want to pull over less items than can fit to the end of buffer, just copy them all
      using_count = buff.bmemcpy( &d_buffer[d_write_at], using_count, false );
      d_write_at += using_count;
    }
    boost::lock_guard<boost::mutex> lk(d_mutex);
    d_item_count += using_count;
    if(check_notify_read()){
      d_cond_read.notify_all();
    }
    buff.unlock();
  }
  else{//read from (this) buffer, write to (buff) buffer
    lock();
    size_t trr = read_room();

    using_count = trr > item_count ? item_count : trr;
    if(using_count+d_read__at > d_max_items){//want to send over more items than the buffer size will allow, but can wrap the rest
      size_t ender = d_max_items - d_read__at;//How many (this) can send before wrap
      size_t begin = using_count - ender;//The remaining (this) will send
      using_count = buff.bmemcpy( &d_buffer[d_read__at], ender, true );
      if(using_count < ender){//Didn't have enough data to write everything, can return now
        d_read__at += using_count;
      }
      else{//There was enough to write, continue writing
        size_t secondary = buff.bmemcpy( &d_buffer[0], begin, true );
        d_read__at = secondary;
        using_count += secondary;
      }
    }
    else{//want to send over less items than can fit to the end of buffer, just copy them all
      using_count = buff.bmemcpy( &d_buffer[d_read__at], using_count, true );
      d_read__at += using_count;
    }
    d_item_count -= using_count;
    unlock();
  }*/
  return using_count;
}



template <class T>
size_t
signal_threaded_buffer<T>::write_room()
{
  return d_max_items - d_item_count;
}

template <class T>
size_t
signal_threaded_buffer<T>::read_room()
{
  return d_item_count;
}

template <class T>
void
signal_threaded_buffer<T>::lock()
{
  d_mutex.lock();
}

template <class T>
void
signal_threaded_buffer<T>::unlock()
{
  d_mutex.unlock();
}

template <class T>
size_t
signal_threaded_buffer<T>::write_item( T item )
{
  T item_buff[] = {item};
  return bmemcpy( &item_buff[0], 1, true );
}

template <class T>
size_t
signal_threaded_buffer<T>::read_item( T& item )
{
  T item_buff[1];
  size_t val = bmemcpy( &item_buff[0], 1, false );
  item = item_buff[0];
  return val;
}



template <class T>
size_t
signal_threaded_buffer<T>::readable()
{
  if(check_notify_read())
    return read_room();
  else
    return 0;
}

template <class T>
size_t
signal_threaded_buffer<T>::writeable()
{
  return write_room();
}

template <class T>
bool
signal_threaded_buffer<T>::check_notify_read()
{
  return (read_room() > d_min_notify);
}

template <class T>
bool
signal_threaded_buffer<T>::check_notify_write()
{
  return (write_room() > d_min_notify);
}

template <class T>
void
signal_threaded_buffer<T>::read_wait()
{
  boost::mutex::scoped_lock lock(d_mutex);
  while(!check_notify_read()){
    if(!d_cond_read.timed_wait(lock,boost::posix_time::milliseconds(100))){
      //printf("TIMEOUT WRITE: ENDING\n");
      break;
    }
  }
}

template <class T>
void
signal_threaded_buffer<T>::write_wait()
{
  boost::mutex::scoped_lock lock(d_mutex);
  while(!check_notify_write()){
    if(!d_cond_write.timed_wait(lock,boost::posix_time::milliseconds(100))){
      //printf("TIMEOUT WRITE: ENDING\n");
      break;
    }
  }
}




template class signal_threaded_buffer<float>;
template class signal_threaded_buffer<int>;
template class signal_threaded_buffer<char>;
template class signal_threaded_buffer<unsigned int>;
template class signal_threaded_buffer<unsigned char>;
template class signal_threaded_buffer<complexf>;
template class signal_threaded_buffer<complexd>;
template class signal_threaded_buffer<complexl>;

















