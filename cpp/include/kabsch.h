#ifndef KABSCH_H
#define KABSCH_H

#include <eigen3/Eigen/Geometry>

// check if reference passing is possible
auto kabsch(Eigen::Matrix3Xf in, Eigen::Matrix3Xf out)
{
  // Default output
  Eigen::Affine3d A;
  A.linear() = Eigen::Matrix3f::Identity(3, 3);
  A.translation() = Eigen::Vector3f::Zero();

#ifdef DEBUG
  assert(in.cols() == out.cols());
#endif

  // First find the scale, by finding the ratio of sums of some distances,
  // then bring the datasets to the same scale.
  float dist_in = 0.f, dist_out = 0.f;
  float scale;
  for (int i = 0; i < in.cols() - 1; ++i)
  {
    dist_in  += (in.col(i + 1)  - in.col(i) ).norm();
    dist_out += (out.col(i + 1) - out.col(i)).norm();
  }

  if (dist_in <= 0.f || dist_out <= 0.f) return A;
  scale = dist_out / dist_in;
  out /= scale;

  // Find the centroids then shift to the origin
  Eigen::Vector3f in_ctr  = Eigen::Vector3f::Zero();
  Eigen::Vector3f out_ctr = Eigen::Vector3f::Zero();
  for (int i = 0; i < in.cols(); ++i)
  {
    in_ctr  += in.col(i);
    out_ctr += out.col(i);
  }
  in_ctr  /= in.cols();
  out_ctr /= out.cols();

  for (int i = 0; i < in.cols(); ++i)
  {
    in.col(i)  -= in_ctr;
    out.col(i) -= out_ctr;
  }

  // SVD
  Eigen::MatrixXf Cov = in * out.transpose();
  Eigen::JacobiSVD<Eigen::MatrixXf> svd(Cov, Eigen::ComputeThinU | Eigen::ComputeThinV);

  // Find the rotation
  float d = (svd.matrixV() * svd.matrixU().transpose()).determinant();
  d = (d > 0.f) ? 1.f : -1.f;
  Eigen::Matrix3f I = Eigen::Matrix3f::Identity(3, 3);
  I(2, 2) = d;
  Eigen::Matrix3f R = svd.matrixV() * I * svd.matrixU().transpose();

  // The final transform
  A.linear()      = scale * R;
  A.translation() = scale*(out_ctr - R*in_ctr);
  return A;
}


#endif // KABSCH_H
