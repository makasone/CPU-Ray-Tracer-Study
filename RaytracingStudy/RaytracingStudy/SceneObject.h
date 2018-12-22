#pragma once

#include "Color.h"
#include "Material.h"
#include "BoundingBox.h"

	class Ray;
	class HitInformation;

	class SceneObject {
	public:
		SceneObject(const Material &material_)
			: material(material_)
			, position(0, 0, 0)
		{
		}
		virtual ~SceneObject() {}

		virtual bool CheckIntersection(const Ray &ray, HitInformation &hit) const = 0;
		//virtual optional<Hit> Intersect(const Ray& ray, const Vector2d& existArea) const = 0;

		Material material;
		Vector3d position;
		BoundingBox boundingBox;
	};
