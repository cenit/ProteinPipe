#ifndef FIT_H
#define FIT_H
#include <utils.h>
#ifdef DEBUG
#include <cassert>
#endif
#include <memory>
#include <algorithm>
#include <numeric>
#include <eigen3/Eigen/Dense>
#include <eigen3/Eigen/Eigenvalues>

#define fit_rmsd // to comment

struct
{
  float operator()(const Eigen::Matrix3Xf &a, const Eigen::Matrix3Xf &b)
  {
#ifdef DEBUG
    assert(a.n == b.n);
#endif

#ifdef fit_rmsd
    float rmsd = 0.f;
    for(int i = 0; i < a.n; ++i) rmsd += distance(a(0, i), b(0, i),
                                                  a(1, i), b(1, i),
                                                  a(2, i), b(2, i));
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
      dist  = distance(a(0, i), b(0, i),
                       a(1, i), b(1, i),
                       a(2, i), b(2, i));
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


struct
{
  float operator()(const Eigen::MatrixXf &laplacian,
                   const masses &mass,
                   const Eigen::Matrix3Xf &true_coords)
  {
    auto eigen_mass = Eigen::Map<Eigen::Matrix<float,
                                 Eigen::Dynamic,
                                 1,
                                 Eigen::RowMajor>>
                                 (&mass.m[0], laplacian.cols(), 1);
    Eigen::SelfAdjointEigenSolver<Eigen::MatrixXf> eigensolver(laplacian * eigen_mass);
    auto eigvals = eigensolver.eigenvalues();

    std::unique_ptr<int[]> idx(new int[N]);
    std::iota(idx.get(), idx.get() + N, 0);
    std::sort(idx.get(), idx.get() + N, [&](const int &a, const int &b){return eigvals(a) < eigvals(b);});

    Eigen::Matrix3Xf guess_coords(3, N);
    guess_coords.row(0) = eigensolver.eigenvectors().col(idx[1]);
    guess_coords.row(1) = eigensolver.eigenvectors().col(idx[2]);
    guess_coords.row(2) = eigensolver.eigenvectors().col(idx[3]);

    auto transf = kabsch(true_coords, guess_coords);

    guess_coords = transf * guess_coords;

    return fit_func(guess_coords, true_coords);
  }
} overalp;


#endif // FIT_H
