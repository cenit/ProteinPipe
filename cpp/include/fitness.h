#ifndef FIT_H
#define FIT_H
#include <utils.h>
#ifdef DEBUG
#include <cassert>
#endif

#define fit_rmsd // to comment

struct
{
  template<class Points> float operator()(const Points &a, const Points &b)
  {
#ifdef DEBUG
    assert(a.n == b.n);
#endif

#ifdef fit_rmsd
    float rmsd = 0.f;
    for(int i = 0; i < a.n; ++i) rmsd += distance(a.x[i], b.x[i], a.y[i], b.y[i], a.z[i], b.z[i]);
    rmsd /= a.n;
    return std::sqrt(rmsd);

#elif defined fit_std
    float d, d2,
          dist,
          mean = 0.f,
          var = 0.f;
    int n = 0;
    for(int i = 0; i < a.n; ++i)
    {
      ++n;
      dist  = distance(a.x[i], b.x[i], a.y[i], b.y[i], a.z[i], b.z[i]);
      delta = dist - mean;
      mean += dist / n;
      d2    = dist - mean;
      var  += d * d2;
    }
    return var / n;
#else
#error "Undefined fitness function! Chose between 'fit_rmsd' and 'fit_rmsd' in define variables"
#endif
  }
} fit_func;


#endif // FIT_H
