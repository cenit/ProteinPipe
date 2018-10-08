#ifndef EIGENVEC_H
#define EIGENVEC_H

#include <memory>
#include <algorithm>
#include <numeric>
#include <eigen3/Eigen/Dense>
#include <eigen3/Eigen/Eigenvalues>

template<typename matrix> auto compute_eigenvec(const matrix &lap, const int &N)
{
  auto eigen_lap = Eigen::Map<Eigen::Matrix<float,
                              Eigen::Dynamic,
                              Eigen::Dynamic,
                              Eigen::RowMajor>>
                              (&lap[0], N, N);
  Eigen::SelfAdjointEigenSolver<Eigen::MatrixXf> eigensolver(eigen_lap);
  auto eigvals = eigensolver.eigenvalues();

  std::unique_ptr<int[]> idx(new int[N]);
  std::iota(idx.get(), idx.get() + N, 0);
  std::sort(idx.get(), idx.get() + N, [&](const int &a, const int &b){return eigvals(a) < eigvals(b);});

  Eigen::Matrix3Xf coords(3, N);
  coords.row(0) = eigensolver.eigenvectors().col(idx[1]);
  coords.row(1) = eigensolver.eigenvectors().col(idx[2]);
  coords.row(2) = eigensolver.eigenvectors().col(idx[3]);

  return coords;
}


#endif // EIGENVEC_H
