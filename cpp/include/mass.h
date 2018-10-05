#ifndef MASSES_H
#define MASSES_H
#include <memory>
#include <algorithm>
#include <random>
#include <iterator>
#ifdef DEBUG
#include <cassert>
#endif

#define chunk  // bitmask     // to comment
#define single // multi       // to comment

constexpr float scale_factor = 1.f;
constexpr int seed = 123;
static std::mt19937 engine(seed);
static std::uniform_real_distribution<float> rng;

struct masses
{
  int n;
  std::unique_ptr<float[]> m;
  masses(const int &n)
  {
    this->n = n;
    this->m = std::make_unique<float[]>(this->n);
  };
  // random generator
  masses(const int &n, const unsigned int &seed)
  {
    this->n = n;
    this->m = std::make_unique<float[]>(this->n);
    std::generate_n(this->m.get(), this->n, [&](){return rng(engine) * scale_factor;});
  }
  // crossover function
  masses operator+(const masses &m2)
  {
#ifdef DEBUG
    assert(this->n == m2.n);
#endif
    masses res(this->n);

#ifdef chunk
    int pos = static_cast<int>(rng(engine) * this->n);
    for(int i = 0; i < this->n; ++i) res.m[i] = (i < pos) ? this->m[i] : m2.m[i];
#elif defined bitmask
    for(int i = 0; i < this->n; ++i) res.m[i] = (rng(engine) < .5f) ? this->m[i] : m2.m[i];
#else
#error "Undefined crossover function! Chose between 'chunk' and 'bitmask' in define variables"
#endif
    return res;
  }
  // mutation function
  void operator!(void)
  {
#ifdef single
    int pos = static_cast<int>(rng(engine) * this->n);
    this->m[pos] = rng(engine) * scale_factor;
#elif defined multi
    for(int i = 0; i < multi; ++i)
    {
      int pos = static_cast<int>(rng(engine) * this->n);
      this->m[pos] = rng(engine) * scale_factor;
    }
#else
#error "Undefined mutation function! Chose between 'single' and 'multi' in define variables"
#endif
  }


  friend std::ostream& operator<<(std::ostream &os, const masses &m)
  {
    std::copy(m.m.get(), m.m.get() + m.n - 1, std::ostream_iterator<float>(os, ", "));
    os << m.m[m.n - 1];
    return os;
  }

};


#endif // MASSES_H
