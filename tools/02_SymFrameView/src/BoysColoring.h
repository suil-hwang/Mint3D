#pragma once

#include <Eigen/Dense>

/*
 * Given any Line field vector v, this function returns the RGB color that corresponds to the direction using Boys' color scheme.
 * Reference: https://github.com/fury-gl/fury/blob/master/fury/colormap.py#L62
 */
Eigen::Vector3d Boys2RGB(const Eigen::Vector3d& v);