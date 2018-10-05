#ifndef UTILS_H
#define UTILS_H
#include <cmath>
struct
{
  float operator()(const float &x1, const float &x2, const float &y1, const float &y2, const float &z1, const float &z2)
  {
    return std::sqrt((x2 - x1)*(x2 - x1) + (y2 - y1)*(y2 - y1) + (z2 - z1)*(z2 - z1));
  }
} distance;

#endif // UTILS_H
