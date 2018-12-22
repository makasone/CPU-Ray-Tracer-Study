#pragma once

#include "Eigen/Core"
#include "Eigen/Dense"

using namespace Eigen;

class HitInformation {
public:
	HitInformation()
		: distance(0)
		, position()
		, normal()
	{
	}

	double distance;
	Vector3d position;
	Vector3d normal;
};