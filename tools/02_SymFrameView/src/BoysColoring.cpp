#include "BoysColoring.h"

// Helper functions
static double cc(double na, double nd) {
  return na * std::cos(nd * M_PI / 180.0);
}

static double ss(double na, double nd) {
  return na * std::sin(nd * M_PI / 180.0);
}

Eigen::Vector3d Boys2RGB(const Eigen::Vector3d& v) {
  Eigen::Vector3d color;

  Eigen::Vector3d norm = v.normalized();

  double x = norm[0], y = norm[1], z = norm[2];

  // Calculations similar to Python code
  double x2 = x * x, y2 = y * y, z2 = z * z;
  double x3 = x * x2, y3 = y * y2, z3 = z * z2;
  double z4 = z2 * z2;
  double xy = x * y, xz = x * z, yz = y * z;

  // Constants
  double hh1 = 0.5 * (3 * z2 - 1) / 1.58;
  double hh2 = 3 * xz / 2.745;
  double hh3 = 3 * yz / 2.745;
  double hh4 = 1.5 * (x2 - y2) / 2.745;
  double hh5 = 6 * xy / 5.5;
  double hh6 = (1 / 1.176) * 0.125 * (35 * z4 - 30 * z2 + 3);
  double hh7 = 2.5 * x * (7 * z3 - 3 * z) / 3.737;
  double hh8 = 2.5 * y * (7 * z3 - 3 * z) / 3.737;
  double hh9 = ((x2 - y2) * 7.5 * (7 * z2 - 1)) / 15.85;
  double hh10 = ((2 * xy) * (7.5 * (7 * z2 - 1))) / 15.85;
  double hh11 = 105 * (4 * x3 * z - 3 * xz * (1 - z2)) / 59.32;
  double hh12 = 105 * (-4 * y3 * z + 3 * yz * (1 - z2)) / 59.32;

  double s0 = -23.0;
  double s1 = 227.9;
  double s2 = 251.0;
  double s3 = 125.0;

  double ss23 = ss(2.71, s0);
  double cc23 = cc(2.71, s0);
  double ss45 = ss(2.12, s1);
  double cc45 = cc(2.12, s1);
  double ss67 = ss(0.972, s2);
  double cc67 = cc(0.972, s2);
  double ss89 = ss(0.868, s3);
  double cc89 = cc(0.868, s3);

  double X = 0.0;

  X = X + hh2 * cc23;
  X = X + hh3 * ss23;

  X = X + hh5 * cc45;
  X = X + hh4 * ss45;

  X = X + hh7 * cc67;
  X = X + hh8 * ss67;

  X = X + hh10 * cc89;
  X = X + hh9 * ss89;

  double Y = 0.0;

  Y = Y + hh2 * -ss23;
  Y = Y + hh3 * cc23;

  Y = Y + hh5 * -ss45;
  Y = Y + hh4 * cc45;

  Y = Y + hh7 * -ss67;
  Y = Y + hh8 * cc67;

  Y = Y + hh10 * -ss89;
  Y = Y + hh9 * cc89;

  double Z = 0.0;

  Z = Z + hh1 * -2.8;
  Z = Z + hh6 * -0.5;
  Z = Z + hh11 * 0.3;
  Z = Z + hh12 * -2.5;

  double w_x = 4.1925, trl_x = -2.0425;
  double w_y = 4.0217, trl_y = -1.8541;
  double w_z = 4.0694, trl_z = -2.1899;

  color[0] = 0.9 * std::abs(((X - trl_x) / w_x)) + 0.05;
  color[1] = 0.9 * std::abs(((Y - trl_y) / w_y)) + 0.05;
  color[2] = 0.9 * std::abs(((Z - trl_z) / w_z)) + 0.05;

  return color;
}