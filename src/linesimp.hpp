#include <boost/heap/fibonacci_heap.hpp>
#include <list>

#include <geoflow/geoflow.hpp>

#include <CGAL/Cartesian.h>

namespace linesimp {

struct point_error;
typedef std::array<float,3> Point;

typedef boost::heap::fibonacci_heap<point_error> Heap;
typedef Heap::handle_type heap_handle;

typedef std::list<std::pair<Point, heap_handle>> PointList;

typedef CGAL::Cartesian<double> K;
typedef K::Point_2 Point2;
typedef K::Line_2 Line;

struct point_error {
  point_error(PointList::iterator i, double e) : it(i), error(e){}
  PointList::iterator it;
  double error;
  
  bool operator<(point_error const & rhs) const
  {
    return error > rhs.error; //smallest error on top
  }
};

std::vector<Point> visvalingam(const std::vector<Point> &line_string, double threshold);

class AproximateLine {
public:
  AproximateLine(Point2 p, Point2 q);
  ~AproximateLine();
  bool canMerge(AproximateLine &l2, double threshold);
  bool merge(AproximateLine &l2);
  std::pair<double, double> calulateC(double A, std::vector<Point2> points);
  static bool compare(const AproximateLine& l1, const AproximateLine& l2);
  std::vector<Point2> points;
  std::vector<Line> lines;
  std::vector<Line> buffer;
  double buffersize;
  double length;
  double A, B, Cmin, Cmax;
};
}