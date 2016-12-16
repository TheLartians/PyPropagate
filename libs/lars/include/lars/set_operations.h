
// The MIT License (MIT)
//
// Copyright (c) 2016 Lars Melchior
//
// Permission is hereby granted, free of charge, to any person obtaining a copy
// of this software and associated documentation files (the "Software"), to deal
// in the Software without restriction, including without limitation the rights
// to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
// copies of the Software, and to permit persons to whom the Software is
// furnished to do so, subject to the following conditions:
//
// The above copyright notice and this permission notice shall be included in
// all copies or substantial portions of the Software.
//
// THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
// IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
// FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
// AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
// LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
// OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
// THE SOFTWARE.


#pragma once

#include <unordered_set>
#include <iterator>
#include <algorithm>

namespace lars {
  
  template <typename InIt1, typename InIt2, typename OutIt> OutIt unordered_set_intersection(InIt1 b1, InIt1 e1, InIt2 b2, InIt2 e2, OutIt out) {
    while (!(b1 == e1)) {
      if (!(std::find(b2, e2, *b1) == e2)) {
        *out = *b1;
        ++out;
      }
      ++b1;
    }
    return out;
  }
  
  template <class T> std::unordered_set<T> set_intersection(const std::unordered_set<T> &x,const std::unordered_set<T> &y){
    std::unordered_set<T> z;
    unordered_set_intersection(x.begin(), x.end(),y.begin(), y.end(),std::inserter(z, z.begin()));
    return z;
  }
  
  template <class T> std::unordered_set<T> set_union(const std::unordered_set<T> &x,const std::unordered_set<T> &y){
    std::unordered_set<T> z = x;
    std::copy( std::begin(y), std::end(y), std::inserter( z, std::end(z) ) );
    return z;
  }
  
  template <class T1,class T2> T1 set_intersection(const T1 &x,const T2 &y){
    T1 z;
    std::set_intersection(x.begin(), x.end(),y.begin(), y.end(),std::inserter(z, z.begin()));
    return std::move(z);
  }
  
  template <class T1,class T2> T1 set_difference(const T1 &x,const T2 &y){
    T1 z;
    std::set_difference(x.begin(), x.end(),y.begin(), y.end(),std::inserter(z, z.begin()));
    return std::move(z);
  }
  
  template <class T1,class T2> T1 set_union(const T1 &x,const T2 &y){
    T1 z;
    std::set_union(x.begin(), x.end(),y.begin(), y.end(),std::inserter(z, z.begin()));
    return std::move(z);
  }
  
  
  template <class T1,class T2,class C> T1 set_intersection(const T1 &x,const T2 &y,C c){
    T1 z;
    std::set_intersection(x.begin(), x.end(),y.begin(), y.end(),std::inserter(z, z.begin()),c);
    return std::move(z);
  }
  
  template <class T1,class T2,class C> T1 set_difference(const T1 &x,const T2 &y,C c){
    T1 z;
    std::set_difference(x.begin(), x.end(),y.begin(), y.end(),std::inserter(z, z.begin()),c);
    return std::move(z);
  }
  
  template <class T1,class T2,class C> T1 set_union(const T1 &x,const T2 &y,C c){
    T1 z;
    std::set_union(x.begin(), x.end(),y.begin(), y.end(),std::inserter(z, z.begin()),c);
    return std::move(z);
  }
  
  
}