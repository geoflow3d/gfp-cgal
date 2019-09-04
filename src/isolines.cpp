#include "isolines.hpp"

namespace isolines {
  // Evaluates whether a vertex is position under, above or on (within range +/- e) the (contour) depth.
  int cntrEvalVertex(Vertex_handle v, double depth) {
    double e = 1e-7;
    double z = v->point().z();

    if (z < depth - e)
      return -1;
    else if (depth + e < z)
      return 1;
    else
      return 0;
  }

  // Calculates and returns the point of intersection between the (contour) depth and the edge formed by two vertices v1 and v2.
  Point cntrIntersectEdge(Vertex_handle v1, Vertex_handle v2, double depth) {
    Point p1 = v1->point();
    Point p2 = v2->point();

    double lambda = (depth - p1.z()) / (p2.z() - p1.z());
    double x = (1 - lambda) * p1.x() + lambda * p2.x();
    double y = (1 - lambda) * p1.y() + lambda * p2.y();

    return Point(x, y, depth);
  }
}