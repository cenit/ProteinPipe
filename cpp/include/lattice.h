#ifndef LATTICE_H
#define LATTICE_H
#include <iostream>
#include <algorithm>
#include <numeric>
#include <array>

template<int N> struct lattice2d
{
  int nodes;
  std::array<float, N*N * N*N> adj;
  lattice2d()
  {
    this->nodes = N*N;
    std::fill_n(adj.begin(), this->nodes*this->nodes, 0.f);
    int idx1, idx2;

    for (int j = 1; j < N; ++j)
    {
      idx1 = j * this->nodes;
      this->adj[idx1 + j - 1]           = 1.f;
      this->adj[idx1*N + (j - 1)*N]     = 1.f;
      idx1 = (j - 1)*N;
      this->adj[idx1*N + j]             = 1.f;
      this->adj[idx1*this->nodes + j*N] = 1.f;
    }

    for (int i = 1; i < N; ++i)
      for (int j = 1; j < N; ++j)
      {
        idx1 = (i*N + j);
        idx2 = (i-1)*N + j;
        this->adj[idx1 * this->nodes + idx2] = 1.f;
        this->adj[idx2 * this->nodes + idx1] = 1.f;

        idx2 = (i*N + j - 1);
        this->adj[idx1 * this->nodes + idx2] = 1.f;
        this->adj[idx2 * this->nodes + idx1] = 1.f;
      }
  }

  void laplacian()
  {
    float degree;
    std::transform(this->adj.begin(), this->adj.end(),
                   this->adj.begin(),
                   [](const float &i)
                   {
                    return i != 0.f ? -i : 0.f;
                   });
    for (int i = 0; i < this->nodes; ++i)
    {
      degree = std::accumulate(this->adj.begin() + i*this->nodes,
                               this->adj.begin() + i*this->nodes + this->nodes,
                               0.f);
      this->adj[i*this->nodes + i] = -degree;
    }
  }

  friend std::ostream& operator<<(std::ostream &os, const lattice<N> &l)
  {
    for (int i = 0; i < l.nodes; ++i)
      for (int j = 0; j < l.nodes; ++j)
        if (l.adj[i*l.nodes + j])
          os << "(" << i << ", " << j << ")" << std::endl;
    return os;
  }

};




#endif // LATTICE_H
