#ifndef __COUPLED_SORT__
#define __COUPLED_SORT__

#include <algorithm>
#include <functional>
#include <iterator>
#include <iostream>
#include <vector>

template <class T> 
class Comp : std::binary_function<size_t,size_t,bool>{
  T * _v;
 public:
    Comp(T * v) : _v(v) {}
    bool operator()(size_t i, size_t j){
     return std::less<T>()(_v[i], _v[j]); 
     }
};

template <class T> 
class FoncTab : std::unary_function<size_t,T>{
  T * _v;
 public:
    FoncTab(T * v) : _v(v) {}
    T operator()(size_t i){
     return _v[i]; 
     }
};


template <class T>
inline void sort_with_perm(T * begin,T * end, std::vector<size_t> & perm)
{
 size_t taille = end - begin;
 perm.clear();
 perm.reserve(taille);
 for(size_t i = 0; i < taille ; i++)
           perm.push_back(i);
 
 std::sort( perm.begin(), perm.end(), Comp<T>(begin) );
 
 T * tmp = new T[taille];
 std::transform( perm.begin(), perm.end(), tmp, FoncTab<T>(begin));
 std::copy( tmp , tmp + taille , begin);

 delete [] tmp;
}

template <class T>
inline void apply_perm(T * begin, T * end, std::vector<size_t> & perm)
{
 size_t taille = end -begin;
 T * tmp = new T[taille];
 if ( (perm.begin() + taille) <= perm.end() ) { 
    std::transform( perm.begin() , perm.begin() + taille, tmp, FoncTab<T>(begin));
    std::copy( tmp , tmp + taille, begin);
 }
 else
    std::cerr << "permutation trop courte" << std::endl;
 delete [] tmp;

}

template <class T1, class T2>
inline void coupled_sort(T1 * begin1,T1 * end1, T2 *begin2)
{
 std::vector<size_t> perm;
 sort_with_perm( begin1, end1, perm);
 apply_perm( begin2, begin2 + (end1 - begin1), perm);
}


#endif

