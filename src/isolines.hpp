#include "tinsimp.hpp"

#include <CGAL/Surface_mesh.h>
#include <CGAL/boost/graph/Euler_operations.h>


#include <geoflow/geoflow.hpp>

namespace isolines {
  typedef tinsimp::CDT          CDT;
  typedef CDT::Vertex_handle                                  Vertex_handle;
  typedef CDT::Face_handle                                    Face_handle;
  typedef CDT::Vertex_iterator                                Vertex_iterator;
  typedef CDT::Face_iterator                                  Face_iterator;
  typedef CDT::Point													                Point;

  int cntrEvalVertex(Vertex_handle v, double depth);
  Point cntrIntersectEdge(Vertex_handle v1, Vertex_handle v2, double depth);

  template<class Constrained_Delaunay_triangulation_2, class FG>
  typename boost::graph_traits<FG>::vertex_descriptor
    link_cdt_to_face_graph(const Constrained_Delaunay_triangulation_2& cdt, FG& fg) {
    typedef typename Constrained_Delaunay_triangulation_2::Face_handle Face_handle;
    typedef typename Constrained_Delaunay_triangulation_2::Vertex_handle Vertex_handle;
    typedef typename boost::graph_traits<FG>::vertex_descriptor vertex_descriptor;

    clear(fg);
    vertex_descriptor inf;
    vertex_descriptor nullvertex = boost::graph_traits<FG>::null_vertex();
    fg.clear();
    typedef boost::unordered_map<Vertex_handle, vertex_descriptor> Vertex_map;
    Vertex_map vertex_map;
    std::vector<Face_handle> faces;// (cdt.finite_faces_begin(), cdt.finite_faces_end());
    for (CDT::Finite_faces_iterator fit = cdt.finite_faces_begin();
      fit != cdt.finite_faces_end();
      ++fit) {
      faces.push_back(fit);
    }
    CGAL::cpp11::array<vertex_descriptor, 3> face;

    typename boost::property_map<FG, CGAL::vertex_point_t>::type vpm
      = get(CGAL::vertex_point, fg);

    BOOST_FOREACH(Face_handle f, faces) {
      for (int i = 0; i < 3; i++) {
        Vertex_handle vh = f->vertex(i);
        std::pair<typename Vertex_map::iterator, bool> res
          = vertex_map.insert(std::make_pair(vh, nullvertex));
        if (res.second) {
          res.first->second = add_vertex(fg);
          put(vpm, res.first->second, vh->point());
          if (cdt.is_infinite(vh)) {
            inf = res.first->second;
          }
        }
        face[i] = res.first->second;
      }
      CGAL::Euler::add_face(face, fg);
    }
    return inf;
  }
}