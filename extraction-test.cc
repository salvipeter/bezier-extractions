#include "bezier-extractions.hh"

#include <fstream>

#include <Eigen/Dense>

#include <geometry.hh>

using namespace Geometry;

void writeCurve(const BSCurve &curve, double from, double to,
                std::string filename, size_t resolution) {
  std::ofstream f(filename);
  f.exceptions(std::ios::failbit | std::ios::badbit);
  for (size_t i = 0; i <= resolution; ++i) {
    double u = from + (to - from) * ((double)i / resolution);
    f << "v " << curve.eval(u) << std::endl;
  }
  f << "l";
  for (size_t i = 0; i <= resolution; ++i)
    f << " " << i + 1;
  f << std::endl;
}

int main(int argc, char **argv) {
  size_t d = 3, L = 3; // L is the number of segments
  DoubleVector knots = {0,0,0,0,1,3,4,4,4,4};
  PointVector cpts = { {0,0,0},{1,1,0},{2,1,0},{3,2,0},{4,1,0},{5,0,0} };
  size_t spans[] = {3, 4, 5};

  writeCurve(BSCurve(d, knots, cpts), 0, 4, "/tmp/curve.obj", 100);

  std::vector<PointVector> Q(L);

  auto C = bezierExtractionMatrices(d, knots);
  for (size_t i = 0; i < L; ++i) {
    std::cout << "Segment " << i + 1 << ":" << std::endl;
    for (size_t j = 0; j <= d; ++j) {
      Point3D p(0, 0, 0);
      for (size_t k = 0; k <= d; ++k)
        p += cpts[spans[i]+k-d] * C[i](k, j);
      std::cout << p << std::endl;
      Q[i].push_back(p);
    }
    writeCurve(BSCurve(Q[i]), 0, 1, std::string("/tmp/curve") + std::to_string(i) + ".obj", 100);
    std::cout << std::endl;
  }

  std::cout << "Reconstructed curve:" << std::endl;
  for (size_t k = 0; k < cpts.size(); ++k) {
    Point3D p(0, 0, 0);
    size_t i = 0;
    while (spans[i] < k)
      ++i;
    auto Cinv = C[i].inverse();
    for (size_t j = 0; j <= d; ++j)
      p += Q[i][j] * Cinv(j, k + d - spans[i]);
    std::cout << p << " (" << (p - cpts[k]).norm() << ')' << std::endl;
  }

  return 0;
}
