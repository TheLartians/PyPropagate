#pragma once

#include <functional>
#include <thread>
#include <future>
#include <mutex>
#include <algorithm>
#include <lars/iterators.h>

namespace lars {
  
  inline unsigned hardware_thread_count(){ return std::thread::hardware_concurrency(); }
  
  template<typename C1,typename C2,typename F> void parallel_for(C1 start,C2 end,F f,uintptr_t thread_count = hardware_thread_count()){
    if(end-start < thread_count) thread_count = end-start;
    std::vector<std::future<void>> handles(thread_count);
    C2 block_size = (end - start)/thread_count;
    for(uintptr_t i=0;i<thread_count-1;++i){
      handles[i] = std::async(std::launch::async,[=](){ for(C2 j=start+block_size*i;j!=start+block_size*(i+1);++j){ f(j); } });
    }
    handles[thread_count-1] = std::async([&](){ for(C2 j=start+block_size*(thread_count-1);j!=end;++j)f(j); });
    for(auto & handle:handles) handle.wait();
  }
  
  template<typename D,typename I1,typename I2,typename F> void unique_parallel_for(I1 start,I2 end,F f,const D &reference = D(), uintptr_t thread_count = hardware_thread_count()){
    if(end-start < thread_count) thread_count = end-start;
    std::vector<std::future<void>> handles(thread_count);
    uintptr_t block_size = (end - start)/thread_count;
    for(uintptr_t i=0;i<thread_count;++i){
      handles[i] = std::async(std::launch::async,[&,f,i](){
        D unique(reference);
        for(auto j:range(start+block_size*i, (i!=thread_count-1)?start+block_size*(i+1):end)) f(j,unique);
      });
    }
    for(auto & handle:handles) handle.wait();
  }

  
  template<typename I,typename F> void parallel_for_each(I start,I end,F f,uintptr_t thread_count = hardware_thread_count()){
    if(thread_count == 1){ std::for_each(start, end, f); return; }
    if(end-start < thread_count) thread_count = end-start;
    std::vector<std::future<void>> handles(thread_count);
    uintptr_t block_size = (end - start)/thread_count;
    for(uintptr_t i=0;i<thread_count-1;++i){
      handles[i] = std::async(std::launch::async,[=](){ std::for_each(start+block_size*i, start+block_size*(i+1), f); });
    }
    handles[thread_count-1] = std::async([=](){ for_each(start+block_size*(thread_count-1), end, f); });
    for(auto & handle:handles) handle.wait();
  }
  
  template<typename F> void parallel(F f,uintptr_t thread_count = hardware_thread_count()){
    std::vector<std::future<void>> handles(thread_count);
    for(uintptr_t i=0;i<thread_count;++i){ handles[i] = std::async(std::launch::async,[=](){ f(); }); }
    for(auto & handle:handles) handle.wait();
  }
  
  
}
