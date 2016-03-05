///////////////////////////////////////////////////////////
//
//    Created in 2014 by Lars Melchior
//
//    This file is distributed under the GNU
//    General Public License (GPL).
//    See LICENSE.TXT for details.
//

#pragma once

#include <iterator>
#include <type_traits>
#include <vector>
#include <assert.h>
#include <cmath>

namespace lars {
  
#if defined(__GNUC__)
#  define UNUSED __attribute__ ((unused))
#elif defined(_MSC_VER)
#  define UNUSED __pragma(warning(suppress:4100))
#else
#  define UNUSED
#endif
  
  
  template<class T> class range_wrapper {
  
    const T r_start,r_end,r_increment;

  public:
    
    class const_iterator : public std::iterator<std::input_iterator_tag, T>{
    private:
      T current;
      const T increment;
      
    public:
      const_iterator(const T &current,const T &increment) : current(current),increment(increment){ }
      
      T operator*()const{ return current; }
      
      const_iterator &operator++() { current+=increment; return *this; }
      const_iterator &operator++(int) { const_iterator range_copy(*this); current+=increment; return range_copy; }
      bool operator!=(const const_iterator &rhs) { return (increment>=0 && (current < rhs.current)) || (increment<0 && (current > rhs.current)); }
    };
    
    
    range_wrapper(const T &r_start, const T &r_end,const T &r_increment) : r_start(r_start), r_end(r_end), r_increment(r_increment) {

    }

    using const_reverse_iterator = const_iterator;

    const_iterator begin() const { return const_iterator(r_start,r_increment); }
    const_iterator end()   const { return const_iterator(r_end  ,r_increment); }

    const_reverse_iterator rbegin() const { return const_reverse_iterator(r_end-1    ,-r_increment); }
    const_reverse_iterator rend()   const { return const_reverse_iterator(r_start-1  ,-r_increment); }
    
    range_wrapper<T> operator+(const T &d){ return range_wrapper<T>(r_start + d,r_end + d,r_increment); }
    range_wrapper<T> operator*(const T &m){ return range_wrapper<T>(r_start * m,r_end * m,r_increment * m); }
    
  };

  template<class T> range_wrapper<T> range(const T &start, const T &end, const T & increment) { return range_wrapper<T>(start, end, increment); }
  template<class T> range_wrapper<T> range(const T &start, const T &end) { return range_wrapper<T>(start, end, 1); }
  template<class T> range_wrapper<T> range(const T &end) { return range_wrapper<T>(0, end, 1); }
  
  //
  // ---- Reverse ----
  //
  
  template<class T> struct reverse_wrapper{
    T & obj;
    typedef typename T::reverse_iterator iterator;
    typedef typename T::iterator reverse_iterator;
    reverse_wrapper(T &o):obj(o){}
    iterator begin(){ return obj.rbegin(); }
    iterator end(){ return obj.rend(); }
    reverse_iterator rbegin(){ return obj.begin(); }
    reverse_iterator rend(){ return obj.end(); }
  };
  
  template<class T> struct const_reverse_wrapper{
    const T & obj;
    typedef typename T::const_reverse_iterator const_iterator;
    typedef typename T::const_iterator const_reverse_iterator;
    const_reverse_wrapper(const T &o):obj(o){}
    const_iterator begin()const{ return obj.rbegin(); }
    const_iterator end()const{ return obj.rend(); }
    const_reverse_iterator rbegin()const{ return obj.begin(); }
    const_reverse_iterator rend()const{ return obj.end(); }
  };
  
  template<typename T, typename = void> struct is_iterator{
    static constexpr bool value = false;
  };
  
  template<typename T> struct is_iterator<T, typename std::enable_if<!std::is_same<typename std::iterator_traits<T>::value_type, void>::value>::type> {
    static constexpr bool value = true;
  };
  
  template <typename T> struct is_iterable{
    template<typename U,size_t (U::*)() const> struct SFINAE { constexpr static bool value=true; };
    template<typename U> static char Test(SFINAE<U, &U::begin>*);
    template<typename U> static int Test(...);
    static constexpr bool value = sizeof(Test<T>(0)) == sizeof(char);
  };

  template<class T> typename std::enable_if<is_iterable<T>::value,reverse_wrapper<T>>::type reverse(T & obj){
    return reverse_wrapper<T>(obj);
  }
  
  template<class T> typename std::enable_if<is_iterable<T>::value,const_reverse_wrapper<T>>::type reverse(const T & obj){
    return const_reverse_wrapper<T>(obj);
  }

  template<class T> struct const_reverse_wrapper<range_wrapper<T>>{
    const range_wrapper<T> obj;
    typedef typename range_wrapper<T>::const_reverse_iterator const_iterator;
    typedef typename range_wrapper<T>::const_iterator const_reverse_iterator;
    const_reverse_wrapper(const range_wrapper<T> &o):obj(o){}
    const_iterator begin()const{ return obj.rbegin(); }
    const_iterator end()const{ return obj.rend(); }
    const_reverse_iterator rbegin()const{ return obj.begin(); }
    const_reverse_iterator rend()const{ return obj.end(); }
  };

  //
  // ---- Split ----
  //
  
  template<class string_type> class string_ref{
  private:
  
  using value_type = typename string_type::value_type;
  
  value_type const * begin_;
  int             size_;
  
  public:
  int size() const { return size_; }
  value_type const* begin() const { return begin_; }
  value_type const* end() const { return begin_ + size_ ; }
    const value_type &operator[](unsigned i){ return *(begin_ + i); }
  
  string_ref( value_type const* const begin, int const size ) : begin_( begin ) , size_( size ){}
    
    operator string_type()const{
      string_type str;
      for(auto c:*this)str+=c;
      return str;
    }
    
  };
  
  template<class string_type> std::vector<string_ref<string_type>> split( const string_type & str, typename string_type::value_type delimiter = ' ' ){
    
    using value_type = typename string_type::value_type;

  std::vector<string_ref<string_type>>   result;
  
  enum State { inSpace, inToken };
  
  State state = inSpace;
  value_type const*     pTokenBegin = 0;    // Init to satisfy compiler.
  for( auto it = str.begin(); it != str.end(); ++it )
    {
    State const newState = (*it == delimiter? inSpace : inToken);
    if( newState != state )
      {
      switch( newState )
        {
          case inSpace:
          result.push_back( string_ref<string_type>( pTokenBegin, &*it - pTokenBegin ) );
          break;
          case inToken:
          pTokenBegin = &*it;
        }
      }
    state = newState;
  }
  if( state == inToken ) result.push_back( string_ref<string_type>( pTokenBegin, &str.back() - pTokenBegin + 1 ) );
  return result;
  }
  


  
}

