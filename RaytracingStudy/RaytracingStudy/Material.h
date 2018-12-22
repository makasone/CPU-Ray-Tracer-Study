#pragma once

#include <functional>
#include "Color.h"


	class Material {
	public:
		enum REFLECTION_TYPE {
			REFLECTION_TYPE_LAMBERT,	// 
			REFLECTION_TYPE_SPECULAR,	// 
			REFLECTION_TYPE_REFRACTION,	// 
		};

	public:
		Material(const REFLECTION_TYPE type = REFLECTION_TYPE_LAMBERT,
			const Vector3d emission_ = Vector3d(0, 0, 0),
			const Vector3d color_ = Vector3d(0, 0, 0),
			const double refraction_rate = 0.0)
			: reflection_type(type)
			, emission(emission_)
			, color(color_)
			, refraction_rate(refraction_rate)
		{}


		REFLECTION_TYPE reflection_type;
		Color emission;
		Color color;
		double refraction_rate;
	};

	struct MaterialHash {
		std::hash<std::string> sh;
		MaterialHash() : sh() {}

		size_t operator()(const Material &mat) const {
			size_t type = static_cast<size_t>(mat.reflection_type);
			std::stringstream ss;
			ss << type << "_" << toString(mat.color) + "_" + toString(mat.emission) + "_" << mat.refraction_rate;
			return sh(ss.str());
		}
	};

	struct MaterialEq {
		bool operator()(const Material &mat1, const Material &mat2) const {
			return mat1.reflection_type == mat2.reflection_type && mat1.color == mat2.color && mat1.emission == mat2.emission && mat1.refraction_rate == mat2.refraction_rate;
		}
	};
