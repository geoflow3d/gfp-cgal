#include <CGAL/algorithm.h>
#include <CGAL/Constrained_triangulation_plus_2.h>
#include <boost/optional.hpp>
#include <glm/glm.hpp>

namespace CGAL {

  template < class Tr >
  class Constrained_triangulation_plus_2;

  namespace Polyline_simplification_2 {
    class Visvalingam_cost {
    public:
      /// Initializes the cost function
      Visvalingam_cost() {
      }

      /*!
      Given a vertex in constraint iterator `viq` computes `vip=std::prev(viq)` and `vir=std::next(vir)`, and the cost of removing vertex `*viq`, replacing edges `(*vip,*viq)` and `(*viq,*vir)` with edge `(*vip,*vir)`.
      \param ct The underlying constrained Delaunay triangulation which embeds the polyline constraints
      \param viq The vertex in constraint iterator of the vertex to remove
      \returns The cost for removing `*viq`. The value `boost::none` can be returned to indicate an infinite or uncomputable cost.
      \tparam CDT must be `CGAL::Constrained_triangulation_plus_2` with a vertex type that
      is model of `PolylineSimplificationVertexBase_2`. `CDT::Geom_traits` must be model of
      the concept `ConstrainedDelaunayTriangulationTraits_2`.
      */
      template <class CDT>
      boost::optional<typename CDT::Geom_traits::FT>
        operator()(const Constrained_triangulation_plus_2<CDT>& pct,
          typename Constrained_triangulation_plus_2<CDT>::Vertices_in_constraint_iterator vicq) const {
        typedef Constrained_triangulation_plus_2<CDT>     CT;
        typedef typename CT::Vertices_in_constraint_iterator Vertices_in_constraint_iterator;
        typedef typename CT::Geom_traits Geom_traits;
        typedef typename Geom_traits::FT                                    FT;
        typedef typename Geom_traits::Point                    Point;

        Vertices_in_constraint_iterator vicp = boost::prior(vicq);
        Vertices_in_constraint_iterator vicr = boost::next(vicq);

        FT d1(0);
        Point const& lP = (*vicp)->point();
        glm::vec3 p = glm::vec3(lP.x(), lP.y(), lP.z());
        Point const& lQ = (*vicq)->point();
        glm::vec3 q = glm::vec3(lQ.x(), lQ.y(), lQ.z());
        Point const& lR = (*vicr)->point();
        glm::vec3 r = glm::vec3(lR.x(), lR.y(), lR.z());

        // the area of the triangle pqr is half of the magnitude of the cross product of q-p and q-r
        d1 = glm::length(glm::cross(q - p, q - r)) / 2;

        return d1;
      }
    };
  }
}