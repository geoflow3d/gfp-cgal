#include "linesimp.hpp"
#include <glm/glm.hpp>
#include <glm/gtc/type_ptr.hpp>

#include <CGAL/min_quadrilateral_2.h>
#include <CGAL/squared_distance_2.h>
#include <CGAL/constructions_d.h>

namespace linesimp {
//--- Visvalingam-Whyatt line simplification for 3D lines, for 2D version see https://bost.ocks.org/mike/simplify/
// inline double compute_error(point& p, point& q, point& r);
inline double compute_error(PointList::iterator qit) {
  auto pit = qit; --pit;
  auto rit = qit; ++qit;
  auto p = glm::make_vec3(pit->first.data());
  auto q = glm::make_vec3(qit->first.data());
  auto r = glm::make_vec3(rit->first.data());
  // the area of the triangle pqr is half of the magnitude of the cross product of q-p and q-r
  return glm::length(glm::cross(q-p,q-r))/2;
}

std::vector<Point> visvalingam(const std::vector<Point>& line_string, double threshold) {
  Heap heap;

  // compute errors for all points except first and last
  PointList point_list;
  for(auto& p : line_string){
    auto it = point_list.insert(point_list.end(), std::make_pair(p,heap_handle()));
  }
  auto start = std::next(point_list.begin()); // skip first point
  auto end = std::prev(point_list.end()); // the last point in the line_string (we'll stop before)
  for (PointList::iterator it=start; it!=end; ++it) {
    auto e = compute_error(it);
    auto handle = heap.push(point_error(it,e));
    it->second = handle;
  }
  
  // insert points, update errors of affected triangles until threshold error is reached
  while (!heap.empty() && heap.top().error < threshold){
    // get top element (with largest error) from heap
    auto maxelement = heap.top();
    auto max_p = *maxelement.it;

    point_list.erase(maxelement.it);
        
    auto it_before = std::prev(maxelement.it);
    auto it_after = std::next(maxelement.it);

    if(it_before != point_list.begin()) {
      (*it_before->second).error = compute_error(it_before);
    }
    if(std::next(it_after) != point_list.end()) {
      (*it_after->second).error = compute_error(it_after);
    }
    // remove the point we just inserted in the triangulation from the heap
    heap.pop();
  }

  std::vector<Point> simplified_lines;
  for (auto& p : point_list) {
    simplified_lines.push_back(p.first);
  }
  return simplified_lines;
}

AproximateLine::AproximateLine(Point2 p, Point2 q) {
  CGAL::line_from_pointsC2(p.x(), p.y(), q.x(), q.y(), A, B, Cmin);
  Cmax = Cmin;
  // Create CGAL line
  Line line = Line(p, q);
  lines.push_back(line);

  // Push points to array for buffer construction
  points.push_back(p);
  points.push_back(q);

  // Set line length
  length = sqrt(CGAL::squared_distance(p, q));
}

AproximateLine::~AproximateLine() {
  lines.clear();
  buffer.clear();
  points.clear();
}

bool AproximateLine::canMerge(AproximateLine &l2, double threshold) {
  // Calculate buffer from lines
  double tempA = tan((atan(A) + atan(l2.A))/2);

  std::vector<Point2> tmpPoints(points);
  tmpPoints.insert(tmpPoints.end(), l2.points.begin(), l2.points.end());
  
  std::pair<double, double> minmaxC = calulateC(tempA, tmpPoints);

  if (abs((minmaxC.second - minmaxC.first) * cos(atan(tempA))) < threshold) {
    return true;
  }
  
  return false;

  //std::vector<Line> tmpBuffer;
  //CGAL::min_strip_2(tmpPoints.begin(), tmpPoints.end(), std::back_inserter(tmpBuffer));
  //double tmpBuffersize = sqrt(CGAL::squared_distance(tmpBuffer[0], tmpBuffer[1]));

  //return tmpBuffersize < threshold;
}

std::pair<double, double> AproximateLine::calulateC(double A, std::vector<Point2> points) {
  double tempC = 0, minC = 999999999999, maxC = -999999999999;

  for (auto &p: points) {
    tempC = p.y() - A * p.x();
    if (tempC < minC) minC = tempC;
    if (tempC > maxC) maxC = tempC;
  }
  return std::make_pair(minC, maxC);
}

bool AproximateLine::merge(AproximateLine &l2) {
  points.insert(points.end(), l2.points.begin(), l2.points.end());
  lines.insert(lines.end(), l2.lines.begin(), l2.lines.end());

  // Calculate buffers for lines
  //buffer.clear();
  CGAL::min_strip_2(points.begin(), points.end(), std::back_inserter(buffer));
  buffersize = sqrt(CGAL::squared_distance(buffer[0], buffer[1]));

  length += l2.length;
  return true;
}

bool AproximateLine::compare(const AproximateLine& l1, const AproximateLine& l2) {
  return l1.length > l2.length;
}
}