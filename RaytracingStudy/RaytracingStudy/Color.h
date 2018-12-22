#pragma once

#include "Eigen/Core"
#include "Eigen/Dense"

using namespace Eigen;
using Color = Vector3d;

static std::string toString(const Color& c) {
	std::stringstream ss;
	ss << c.x() << "," << c.y() << "," << c.z();
	return ss.str();
}
