#ifndef _MVECTOR_H_
#define _MVECTOR_H_
/* next line added for AIX */
#define bool int

static const int MYBAD_IND = -1;

#include <stdlib.h> 
typedef int (* compare_func_t) (const void *, const void *);
//  Mvector:
//      Templated array-based list. It resizes
//  itself as necessary. Constructor lets you
//  specify the expected size.
//
//----------------------------------------------
#define cMvector const Mvector
template <class T>
class Mvector {
 protected:
  T      *_array;      // pointer to the data
  int     _size;        // number of elements in the array
  int     _max;        // max elements for currently allocated array

  // needed by LIST (below -- ARRAY of REFptrs): 
  virtual void clear_ele  (int)        {}
  virtual void clear_range(int, int)   {}
  void  append_ele(const T& el) {
    // append element
    if (_size == _max) realloc();
    _array[_size++] = el;
  }
  
 public:
  // ******** MANAGERS ********
  Mvector(int m=0) : _array(0), _size(0), _max(m > 0 ? m : 0) {
    if (_max) _array = new T[_max];
  }
  Mvector(cMvector<T>& l) : _array(0), _size(0), _max(0) { *this = l; }
  virtual ~Mvector() { clear(); delete [] _array;}
  
  // ******** ACCESSORS/CONVENIENCE ********
  int          size()              const { return _size; }
  bool         empty()            const { return (_size<=0); }
  bool         valid_index(int k) const { return (k>=0 && k<_size); }
  
  T*           array()                  { return _array; }
  T&           operator [](int j)       { return _array[j]; }
  T&           last()                   { if (empty()){
    fprintf(stderr,"Mvector::last-list empty"); 
    exit(1);
  }
  return _array[_size-1]; }
  const T& operator [](int j)     const { return _array[j]; }
  const T&           last()       const { if (empty()){
    fprintf(stderr,"Mvector::last-list empty");
    exit(1);
  }
  return _array[_size-1]; }
  
  // ******** MEMORY MANAGEMENT ********
  
  void clear() { clear_range(0,_size); _size=0; }
  virtual void truncate(int n) {
    if (valid_index(n)) {
      clear_range(n,_size);
      _size = n;
    } else {
      fprintf(stderr,"Mvector::truncate: bad index %d", n);
      exit(1);
    };
  }
  
  // realloc() is called automatically as necessary
  // when the array runs out of room, or explicitly
  // when it's known that the array will shortly need
  // to contain a given number of elements. new_max
  // tells how large the array should be -- if the
  // array is already that large nothing happens.
  // otherwise, the array is reallocated and its elements
  // are copied to the new array. if new_max == 0, this
  // is interpreted as a request that _max should be
  // doubled. the exception is when _max is also 0
  // (meaning the array itself is null), in which case
  // _max is set to 1 and the array is allocated to
  // hold a single element.
  virtual void realloc(int new_max=0) {
    if (new_max && new_max <= _max)
      return;
    _max = (new_max == 0) ? (_max ? _max*2 : 1) : new_max;
    T *tmp = new T [_max];
    if (tmp == 0){
      fprintf(stderr,"Mvector: realloc failed");
      exit(1);
    }
    for (int i=0; i<_size; i++) {
      tmp[i] = _array[i];
      clear_ele(i);
    }
    delete [] _array;
    _array = tmp;
  }
  
  // ******** CONTAINMENT ********
  int get_index(const T &el) const {
    // return index of element
    // or MYBAD_IND if element is not found:
    for (int k = _size-1; k >= 0; k--)
      if (_array[k] == el)
	return k;
    return MYBAD_IND;
  }
  
  // ******** ADDING ********
  void operator += (const T& el) {
    append_ele(el);
  }
  
  // append (same as +=):
  void add (const T& p) { *this += p; }     
  

  //append uniquely
  void add_once(const T& p){
    for (int i = 0; i < size(); i++) {
      // Use !(x == y) because == should be available and != may not be
      if ((*this)[i] == p) return;
    }
    *this += p;
  }
  


  // ******** REMOVING ********
  bool remove(int k) {
    // remove element k
    // return 1 on success, 0 on failure:
    if (valid_index(k)) {
      // replace element k with last element and shorten list:
      _array[k] = _array[--_size];
      clear_ele(_size);
      return 1;
    } else if (k != MYBAD_IND) {
      fprintf(stderr,"Mvector::remove: invalid index %d", k);
      return 0;
    } else return 0; // assume the call to get_index() failed
  }
  
  // search for given element, remove it:
  void operator -= (cMvector<T> &l)        { for (int i=0; i < l.size(); i++)
    *this -= l[i]; }
  
  bool operator -= (const T &el)         { return remove(get_index(el)); }
  bool          rem(const T &p)          { return (*this -= p); }
  
  
  // ******** Mvector OPERATORS ********
  Mvector<T>& operator =(cMvector<T>& l) {
    // assignment operator:
    if (&l == this)  // don't do anything if rhs already is lhs
      return *this;
    clear();
    if(!l.empty()) {
      realloc(l._size);
      for (int i=0; i<l._size; i++)
	*this += l[i];
    }
    return *this;
  }
  Mvector<T>& operator +=(cMvector<T>& l) {
    // concatenation operator:
    if(!l.empty()) {
      realloc(_size + l._size);
      for (int i=0; i<l._size; i++)
	*this += l[i];
    }
    return *this;
  }
  bool operator ==(const Mvector<T> &c) const {
    if (size() != c.size())
      return 0;
    for (int i = 0; i < size(); i++) {
      // Use !(x == y) because == should be available and != may not be
      if (!((*this)[i] == c[i]))
	return 0;
    }
    return 1;
  }
};

#endif
