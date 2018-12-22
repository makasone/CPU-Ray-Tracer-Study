#include <iostream>
#include <fstream>
#include <vector>
#include <optional>
#include <random>
#include <tuple>
#include <omp.h>
#include <memory>

#define _USE_MATH_DEFINES
#include <math.h>

#include "Eigen/Core"
#include "Eigen/Dense"

#include "Model.h"
#include "Ray.h"

//- 設定メモ
//	- C++17を使うためのコマンドライン引数を追加
//	- Eigenでのコンパイルエラー回避のためにSDLチェックをしない、warnning:4996を無視

//- プログラムメモ
//	- Modelという基底クラスを作ろうとしたが、めんどくさくなってやめた

using namespace std;
using namespace Eigen;

//struct Model;
struct Sphere;
using TestModel = Sphere;

struct Hit
{
	double dist;
	Vector3d pos;
	Vector3d normal;
	const TestModel* model;
};

enum class SurfaceType
{
	Diffuse,
	Mirror,
	Fresnel,

	Num
};

struct Sphere : SceneObject
{
	Vector3d m_pos;
	double m_radius;
	Vector3d m_reflectance;
	Vector3d m_illuminance;
	SurfaceType m_type;
	double ior = 1.5168;

	Sphere(const Vector3d& pos, double radius, const Vector3d& reflectance, const Vector3d& illuminance, SurfaceType surfaceType)
		:	SceneObject(Material()),	//#TODO てきとー
			m_pos(pos),
			m_radius(radius),
			m_reflectance(reflectance),
			m_illuminance(illuminance),
			m_type(surfaceType)
	{

	}

	//#TODO 後で消す
	bool CheckIntersection(const Ray &ray, HitInformation &hit) const {
		return false;
	}

	optional<Hit> Intersect(const Ray& ray, const Vector2d& existArea) const/* override*/
	{
		const Vector3d rayToSphere = m_pos - ray.pos;
		const double b = rayToSphere.dot(ray.dir);
		const double det = b * b - rayToSphere.dot(rayToSphere) + m_radius * m_radius;

		if (det < 0) {
			return {};
		}

		const double t1 = b - sqrt(det);
		if (existArea.x() < t1 && t1 < existArea.y()) {
			return Hit{ t1, Vector3d::Zero(), Vector3d::Zero(), this };
		}

		const double t2 = b + sqrt(det);
		if (existArea.x() < t2 && t2 < existArea.y()) {
			return Hit{ t2, Vector3d::Zero(), Vector3d::Zero(), this };
		}

		return {};
	}
};

struct Scene
{
	vector<TestModel> models;

	optional<Hit> InterSect(const Ray& ray, const Vector2d& existArea) const
	{
		optional<Hit> minHit;
		Vector2d searchErea = existArea;

		for (const auto& model : models){
			const auto& hit = (&model)->Intersect(ray, searchErea);
			if (!hit) {
				continue;
			}

			minHit = hit;
			searchErea.y() = minHit->dist;
		}

		if (minHit) {
			minHit->pos = ray.pos + ray.dir * minHit->dist;
			minHit->normal = (minHit->pos - minHit->model->m_pos) / minHit->model->m_radius;
		}

		return minHit;
	}

};

int ToneMap(double v) {
	return clamp(int(pow(v, 1 / 2.2) * 255), 0, 255);
}

struct Random
{
	mt19937 engine;
	uniform_real_distribution<double> dist;

	Random() {};
	Random(int seed) 
	{
		engine.seed(seed);
		dist.reset();
	}

	double next()
	{
		return dist(engine);
	}
};

struct Camera
{
	Vector3d eye;
	Vector3d center;
	Vector3d up;
	Vector3d wE;
	Vector3d uE;
	Vector3d vE;

	Vector2i screenSize;
	double fov;
	double tanFov;
	double aspect;

	Camera()
	{

	}

	Camera(const Vector3d& eye, const Vector3d& center, const Vector3d& up, double fov, const Vector2i& screenSize)
	{
		this->eye = eye;
		this->center = center;
		this->up = up;

		wE = (eye - center).normalized();
		uE = up.cross(wE).normalized();
		vE = wE.cross(uE);

		this->fov = fov * M_PI / 180.0;
		this->screenSize = screenSize;
		this->aspect = (double)screenSize.x() / (double)screenSize.y();
		this->tanFov = tan(this->fov * 0.5);
	}

	Vector3d CalcToPixelDir(const Vector2d& pos) const
	{
		//ピクセル位置を[0,1]へ正規化、さらに[-1,1]へ
		Vector2d normalizedPos(pos.x() / screenSize.x(), pos.y() / screenSize.y());
		normalizedPos = normalizedPos * 2.0 - Vector2d::Ones();

		Vector3d w(aspect * tanFov * normalizedPos.x(), tanFov * normalizedPos.y(), -1.0);
		w.normalize();

		return uE * w.x() + vE * w.y() + wE * w.z();
	}
};

double Sign(double A) {
	return (A > 0) - (A < 0);
}

std::tuple<Vector3d, Vector3d> tangentSpace(const Vector3d& n) {
	const double s = std::copysign(1, n.z());
	const double a = -1 / (s + n.z());
	const double b = n.x()*n.y()*a;
	return {
		Vector3d(1 + s * n.x()*n.x()*a,s*b,-s * n.x()),
		Vector3d(b,s + n.y()*n.y()*a,-n.y())
	};
}

Vector3d CalcMirrorDir(const Vector3d& rayDir, const Vector3d& normal)
{
	const auto wi = -rayDir;
	return 2.0 * wi.dot(normal) * normal - wi;
}

double CalcFresnel(double cos, double ior)
{
	//Schlick's approximation
	const auto r = (1.0 - ior) / (1 + ior);
	return r * r + (1 - r * r) * pow(1.0 - cos, 5.0);
}

Vector3d CalcDir(SurfaceType type, Random& randomNumberGenerator, const optional<Hit>& hit, const Ray& ray)
{
	Vector3d output = Vector3d::Zero();	//こうしておかないとEigenが競合してしまう

	//Sample direction in local coordinates
	if (hit->model->m_type == SurfaceType::Diffuse) {
		const auto n = Sign(hit->normal.dot(-ray.dir)) * hit->normal;
		const auto&[u, v] = tangentSpace(n);
		const auto dir = [&]() {
			const double r = sqrt(randomNumberGenerator.next());
			const double t = 2.0 * M_PI * randomNumberGenerator.next();
			const double x = r * cos(t);
			const double y = r * sin(t);
			return Vector3d(x, y, sqrt(max(0.0, 1.0 - x * x - y * y)));
		}();

		//Convert to world coordinates
		output = u * dir.x() + v * dir.y() + n * dir.z();
	}
	else if (hit->model->m_type == SurfaceType::Mirror) {
		output = CalcMirrorDir(ray.dir, hit->normal);
	}
	else if (hit->model->m_type == SurfaceType::Fresnel) {
		const auto wi = -ray.dir;
		const auto into = wi.dot(hit->normal) > 0.0;
		const auto n = into ? hit->normal : -hit->normal;
		const auto ior = hit->model->ior;
		const auto eta = into ? 1.0 / ior : ior;
		const auto wt = [&]() -> std::optional<Vector3d> {
			//Snell's low (vector form)
			const auto t = wi.dot(n);
			const auto t2 = 1.0 - eta * eta * (1.0 - t * t);
			if (t2 < 0) {
				return {};
			}

			return eta * (n * t - wi) - n * sqrt(t2);
		}();

		if (!wt) {
			//Total internal reflection
			output = CalcMirrorDir(ray.dir, hit->normal);
		}
		else {
			const auto cos = into ? wi.dot(hit->normal) : wt->dot(hit->normal);
			const auto Fr = CalcFresnel(cos, ior);

			//Select reflection or refraction
			//Acocording to the fresnel term
			output = randomNumberGenerator.next() < Fr ? CalcMirrorDir(ray.dir, hit->normal) : *wt;
		}
	}

	return output;
}

void AddObject(SceneObject *obj, bool doDelete = true, bool containedInBVH = true) {
	//m_objects.push_back(SceneObjectInfo(obj, doDelete, containedInBVH));
	//if (containedInBVH) {
	//	m_inBVHObjects.push_back(obj);
	//}
	//else {
	//	m_notInBVHObjects.push_back(obj);
	//}
}

void AddModel(Model *obj, bool doDelete = true, bool containedInBVH = true) {
	//m_models.push_back(ModelObjectInfo(obj, doDelete));

	for (size_t i = 0; i < obj->GetMaterialCount(); i++) {
		const Material &mat = obj->GetMaterial(i);
		const Model::PolygonList &pl = obj->GetPolygonList(mat);
		for (size_t j = 0; j < pl.size(); j++) {
			AddObject(pl[j], false, containedInBVH);
		}
	}
}

void Initialize(Scene& scene, Camera& camera)
{
	camera = Camera(Vector3d(50.0, 52.0, 295.6), Vector3d(50.0, 52.0, 295.6) + Vector3d(0.0, -0.042612, -1.0), Vector3d(0.0, 1.0, 0.0), 30.0, Vector2i(1200, 800));
	//Camera camera(Vector3d(5.0, 5.0, 5.0), Vector3d::Zero(), Vector3d(0.0, 1.0, 0.0), 30.0, Vector2i(1200, 800));

	//scene.models.emplace_back(Sphere{ Vector3d(-0.5, 0.0, 0.0), 1.0, Vector3d(1.0, 0.0, 0.0) });
	//scene.models.emplace_back(Sphere{ Vector3d( 0.5, 0.0, 0.0), 1.0, Vector3d(0.0, 1.0, 0.0) });

	scene.models.emplace_back(Sphere{ Vector3d(1e5 + 1,		40.8, 81.6), 1e5, Vector3d(0.75, 0.25, 0.25),		Vector3d::Zero(),			SurfaceType::Diffuse });
	scene.models.emplace_back(Sphere{ Vector3d(-1e5 + 99,	40.8, 81.6), 1e5, Vector3d(0.25, 0.25, 0.75),		Vector3d::Zero(),			SurfaceType::Diffuse });
	scene.models.emplace_back(Sphere{ Vector3d(50,			40.8,  1e5), 1e5, Vector3d(0.75, 0.75, 0.75),		Vector3d::Zero(),			SurfaceType::Diffuse });

	scene.models.emplace_back(Sphere{ Vector3d(50,			1e5, 81.6),  1e5, Vector3d(0.75, 0.75, 0.75),		Vector3d::Zero(),			SurfaceType::Diffuse });
	scene.models.emplace_back(Sphere{ Vector3d(50,	-1e5 + 81.6, 81.6),  1e5, Vector3d(0.75, 0.75, 0.75),		Vector3d::Zero(),			SurfaceType::Diffuse });
	scene.models.emplace_back(Sphere{ Vector3d(27,		   16.5, 47),	16.5, Vector3d(0.999, 0.999, 0.999),	Vector3d::Zero(),			SurfaceType::Mirror });
	scene.models.emplace_back(Sphere{ Vector3d(73,		   16.5, 78),	16.5, Vector3d(0.999, 0.999, 0.999),	Vector3d::Zero(),			SurfaceType::Fresnel });
	scene.models.emplace_back(Sphere{ Vector3d(50, 681.6 - 0.27, 81.6),	 600, Vector3d::Zero(),					Vector3d(12.0, 12.0, 12.0),	SurfaceType::Diffuse });

	//ポリゴンを読み込んでみる


	//モデルを読み込んでみる
	std::unique_ptr<Model> torii(new Model);
	if (!torii->ReadFromObj("C:\\Users\\makas\\Source\\Repos\\RaytracingStudy\\RaytracingStudy\\torii.obj")) {
		std::cerr << "failed to load cube.obj!!!" << std::endl;
		return;
	}

	torii->SetTransform(Vector3d(50.0, 20, 80), Vector3d(2, 2, 2));

	AddModel(torii.get(), false, false);
}

bool Simulation(const Scene& scene, const Camera& camera, vector<Vector3d>& colors)
{
	//100, 5 で数分でできる　見た目は悪いので動作確認用
	//10000, 10で3~4時間くらいかかる　見た目は非常に良い
	const int samplePixelNum = 100;
	const int depthNum = 5;				//1だと真っ黒になる
	const int pixelNum = (int)colors.size();

	#pragma omp parallel for schedule(dynamic, 1)
	for (int pixelIndx = 0; pixelIndx < pixelNum; pixelIndx++) {
		thread_local Random randomNumberGenerator(42 + omp_get_thread_num());

		for (int sampleIndx = 0; sampleIndx < samplePixelNum; sampleIndx++) {			
			//ピクセル内でのランダムな位置を決定
			Vector2d pos(pixelIndx % camera.screenSize.x(), camera.screenSize.y() - pixelIndx / camera.screenSize.x());	//上下をひっくり返している
			pos += Vector2d(randomNumberGenerator.next(), randomNumberGenerator.next());

			Ray ray = { camera.eye, camera.CalcToPixelDir(pos) };

			//レイの反射を計算
			Vector3d color(Vector3d::Zero()), throughput(Vector3d::Ones());
			Vector2d minMax(1e-4, 1e+10);
			for (int depth = 0; depth < depthNum; depth++) {
				//Intersection
				const auto& hit = scene.InterSect(ray, minMax);	//自分にぶつかることを防ぐため、minに小さな値を入れている
				if (!hit) {
					break;
				}

				//Add contribution
				color += Vector3d(	throughput.x() * hit->model->m_illuminance.x(),
									throughput.y() * hit->model->m_illuminance.y(),
									throughput.z() * hit->model->m_illuminance.z());

				//Update next directino
				ray.pos = hit->pos;
				ray.dir = CalcDir(hit->model->m_type, randomNumberGenerator, hit, ray);

				//Update throughput
				//f/p * cos = R
				throughput = Vector3d(	throughput.x() * hit->model->m_reflectance.x(),
										throughput.y() * hit->model->m_reflectance.y(),
										throughput.z() * hit->model->m_reflectance.z());
				if (throughput.maxCoeff() == 0) {
					break;
				}
			}

			colors[pixelIndx] += color / samplePixelNum;
		}
	}

	return true;
}

bool OutputImage(const Camera& camera, const vector<Vector3d>& colors)
{
	//ファイルへ出力
	ofstream ofs("result.ppm");
	ofs << "P3" << endl;
	ofs << camera.screenSize.x() << " " << camera.screenSize.y() << endl;
	ofs << "255" << endl;

	for (const auto& color : colors) {
		ofs << ToneMap(color.x()) << " "
			<< ToneMap(color.y()) << " "
			<< ToneMap(color.z()) << "\n";
	}

	return true;
}

int main()
{
	Scene scene;
	Camera camera;

	Initialize(scene, camera);

	vector<Vector3d> colors(camera.screenSize.x() * camera.screenSize.y());
	Simulation(scene, camera, colors);

	OutputImage(camera, colors);

	system("result.ppm");

	return 0;
}