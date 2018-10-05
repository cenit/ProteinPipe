#ifndef EIGENVEC_H
#define EIGENVEC_H

#include <eigen3/Eigen/Dense>

template<typename matrix> auto compute_eigenvec(const matrix &lap, const int &N)
{
  auto eigen_lap = Eigen::Map<Eigen::Matrix<float,
                              Eigen::Dynamic,
                              Eigen::Dynamic,
                              Eigen::RowMajor>>
                              (&lap[0], N, N);
  Eigen::EigenSolver<Eigen::MatrixXf> eigensolver(eigen_lap);

}


#endif // EIGENVEC_H
