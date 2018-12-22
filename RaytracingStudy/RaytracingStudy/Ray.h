#pragma once
#pragma once

#include "Eigen/Core"
#include "Eigen/Dense"
using namespace Eigen;

	class Ray {
	public:
		Ray(const Vector3d &begin_, const Vector3d &dir_)
			: pos(begin_)
			, dir(dir_)
		{
		}

		Vector3d pos;
		Vector3d dir;
	};
