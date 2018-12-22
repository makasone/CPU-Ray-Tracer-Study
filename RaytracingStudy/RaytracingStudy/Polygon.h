#pragma once
#pragma once

#include "SceneObject.h"
#include "Ray.h"
#include "Constant.h"
#include "HitInformation.h"

	class Polygon : public SceneObject {
	public:
		Polygon(const Vector3d &pos1, const Vector3d &pos2, const Vector3d &pos3, const Vector3d &normal, const Material &mat, const Vector3d &pos)
			: SceneObject(mat)
			, m_normal(normal)
		{
			m_pos[0] = m_rotatedPos[0] = pos1;
			m_pos[1] = m_rotatedPos[1] = pos2;
			m_pos[2] = m_rotatedPos[2] = pos3;
			position = pos;
			reconstruct_boundingbox();
		}
		Polygon(const Polygon &polygon)
			: SceneObject(polygon.material)
		{
			for (int i = 0; i < 3; i++) {
				m_pos[i] = polygon.m_pos[i];
				m_rotatedPos[i] = polygon.m_rotatedPos[i];
			}
			m_normal = polygon.m_normal;
			reconstruct_boundingbox();
		}
		virtual ~Polygon() {}

		static Vector3d CalculateNormal(const Vector3d &anticlockwise_v0, const Vector3d &anticlockwise_v1, const Vector3d &anticlockwise_v2) {
			Vector3d v((anticlockwise_v1 - anticlockwise_v0).cross(anticlockwise_v2 - anticlockwise_v0));
			v.normalize();
			return v;
		}

		void SetTransform(const Vector3d &pos, const Vector3d &scale = Vector3d::Ones(), const MatrixXd &rot = MatrixXd::Identity(4, 4)) {
			position = pos;
			for (int i = 0; i < 3; i++) {
				//m_rotatedPos[i] = rot.Apply(m_pos[i]);
				m_rotatedPos[i] = m_pos[i] * rot;
				m_rotatedPos[i].x() *= scale.x();
				m_rotatedPos[i].y() *= scale.y();
				m_rotatedPos[i].z() *= scale.z();
			}
			reconstruct_boundingbox();
		}

		bool CheckIntersection(const Ray &ray, HitInformation &hit) const {
			//
			// http://shikousakugo.wordpress.com/2012/07/01/ray-intersection-3/
			Vector3d edge1(m_rotatedPos[1] - m_rotatedPos[0]);
			Vector3d edge2(m_rotatedPos[2] - m_rotatedPos[0]);

			Vector3d P(ray.dir.cross(edge2));
			double det = P.dot(edge1);

			if (det > EPS) {
				// solve u
				Vector3d T(ray.pos - (m_rotatedPos[0] + position));
				double u = P.dot(T);

				if (u >= 0 && u <= det) {
					// solve v
					Vector3d Q(T.cross(edge1));
					double v = Q.dot(ray.dir);

					if (v >= 0 && u + v <= det) {
						double t = Q.dot(edge2) / det;

						if (t >= EPS) {
							hit.distance = t;
							hit.position = ray.pos + ray.dir*t;
							hit.normal = m_normal;

							return true;
						}
					}
				}

			}

			return false;
		}

		Vector3d m_pos[3];
		Vector3d m_normal;

	private:
		void reconstruct_boundingbox() {
			boundingBox.SetBox(Vector3d(
				std::min(std::min(m_rotatedPos[0].x(), m_rotatedPos[1].x()), m_rotatedPos[2].x()),
				std::min(std::min(m_rotatedPos[0].y(), m_rotatedPos[1].y()), m_rotatedPos[2].y()),
				std::min(std::min(m_rotatedPos[0].z(), m_rotatedPos[1].z()), m_rotatedPos[2].z())
			) + position,
				Vector3d(
					std::max(std::max(m_rotatedPos[0].x(), m_rotatedPos[1].x()), m_rotatedPos[2].x()),
					std::max(std::max(m_rotatedPos[0].y(), m_rotatedPos[1].y()), m_rotatedPos[2].y()),
					std::max(std::max(m_rotatedPos[0].z(), m_rotatedPos[1].z()), m_rotatedPos[2].z())
				) + position
			);
		}

	private:
		Vector3d m_rotatedPos[3];
	};
