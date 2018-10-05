#ifndef PROTEIN_H
#define PROTEIN_H

#include <utils.h>
#include <iostream>
#include <memory>
#include <fstream>
#include <algorithm>
#include <numeric>
#include <iterator>
#include <cmath>

struct protein
{
  int n;
  std::unique_ptr<float[]> x, y, z;

  explicit protein(const std::string &filename)
  {
    std::ifstream is(filename);
    is.unsetf(std::ios_base::skipws);
    //count the newlines with an algorithm specialized for counting:
    this->n = std::count(std::istream_iterator<char>(is),
                          std::istream_iterator<char>(),
                          '\n');

    is.clear();
    is.seekg(0, std::ios::beg);
    is.setf(std::ios_base::skipws);
    this->x = std::make_unique<float[]>(this->n);
    this->y = std::make_unique<float[]>(this->n);
    this->z = std::make_unique<float[]>(this->n);

    for (int i = 0; i < this->n; ++i)
    {
      is >> this->x[i];
      is >> this->y[i];
      is >> this->z[i];
    }
    is.close();
  }

  auto mat_dist(float thr = 6.f)
  {
    std::unique_ptr<float[]> d(new float[this->n*this->n]);
    float dist;
    for(int i = 0; i < this->n; ++i)
      for(int j = 0; j < i; ++i)
      {
        dist = distance(this->x[i], this->x[j],
                        this->y[i], this->y[j],
                        this->z[i], this->z[j]
                        );
        d[i * this->n + j] = d[j * this->n + i] = (dist < thr) ? dist : 0.f;
      }
    return d;
  }

  void laplacian(std::unique_ptr<float[]> &adj, const int &N)
  {
    float degree;
    std::transform(adj.get(), adj.get() + N*N, adj.get(),
                   [](const float &i)
                   {
                    return i != 0.f ? -i : 0.f;
                   });
    for (int i = 0; i < N; ++i)
    {
      degree = std::accumulate(adj.get() + i*N,
                               adj.get() + i*N + N,
                               0.f);
      adj[i*N + i] = -degree;
    }
  }

  friend std::ostream& operator<<(std::ostream &os, const protein &p)
  {
    int i;
    for(i = 0; i < p.n - 1; ++i) std::cout << p.x[i] << " " << p.y[i] << " " << p.z[i] << std::endl;
    std::cout << p.x[i] << " " << p.y[i] << " " << p.z[i];
    return os;
  }

};


#endif // PROTEIN_H
