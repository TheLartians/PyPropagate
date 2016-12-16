#pragma once

#include <ostream>

namespace std{
  
  inline wostream &operator<<(wostream & stream,const string &str){
    stream << str.c_str();
    return stream;
  }
  
}

