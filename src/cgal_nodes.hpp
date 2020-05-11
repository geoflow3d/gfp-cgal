#pragma once

#include <geoflow/geoflow.hpp>
#include "tinsimp.hpp"


// Alpha Shape

#include <CGAL/Alpha_shape_2.h>
#include <CGAL/Alpha_shape_vertex_base_2.h>
#include <CGAL/Alpha_shape_face_base_2.h>
#include <CGAL/Delaunay_triangulation_2.h>

// Simplification
#include <CGAL/Polyline_simplification_2/simplify.h>

namespace geoflow::nodes::cgal {

    typedef tinsimp::CDT  CDT;

    class CDTNode :public Node {
        bool create_triangles = false;
    public:
        using Node::Node;
        void init() {
            add_input("geometries", { typeid(PointCollection), typeid(LineStringCollection) });
            add_output("cgal_cdt", typeid(CDT));
        }
        void process();
    };

    class DTNode :public Node {
    public:
        using Node::Node;
        void init() {
            add_input("points", typeid(PointCollection));
            add_output("cgal_dt", typeid(CDT));
        }
        void process();
    };

    class ComparePointDistanceNode :public Node {
        std::string las_filepath = "";
        std::string log_filepath = "";
        int thin_nth = 20;
    public:
        using Node::Node;
        void init() {
            add_input("triangles1_vec3f", typeid(vec3f));
            add_input("triangles2_vec3f", typeid(vec3f));
            add_output("points", typeid(vec3f));
            add_output("distances1", typeid(vec1f));
            add_output("distances2", typeid(vec1f));
            add_output("diff", typeid(vec1f));

            add_param(ParamPath(las_filepath, "las_filpath",  "LAS path"));
            add_param(ParamPath(log_filepath, "log_filpath",  "LOG path"));
            add_param(ParamInt(thin_nth, "thin_nth",  "Thin factor"));
        }
        void process();
    };

    class PointDistanceNode :public Node {
        std::string filepath = "";
        int thin_nth = 5;
        bool overwritez = false;
    public:
        using Node::Node;
        void init() {
            add_input("triangles", typeid(TriangleCollection));
            add_output("points", typeid(PointCollection));
            add_output("distances", typeid(vec1f));
            add_output("distance_min", typeid(float));
            add_output("distance_max", typeid(float));

            add_param(ParamPath(filepath, "filepath",  "Filepath"));
            add_param(ParamInt(thin_nth, "thin_nth",  "Thin factor"));
            add_param(ParamBool(overwritez, "overwritez",  "Overwrite Z"));
        }
        void process();
    };

    class CDTDistanceNode :public Node {
    public:
        using Node::Node;
        void init() {
            add_input("cgal_cdt_base", typeid(CDT));
            add_input("cgal_cdt_target", typeid(CDT));
            add_output("points", typeid(PointCollection));
            add_output("distance_min", typeid(float));
            add_output("distance_max", typeid(float));
        }
        void process();
    };

    class DensifyNode :public Node {
        int interval = 2;
    public:
        using Node::Node;
        void init() {
            add_input("geometries", { typeid(LineStringCollection) });
            add_output("dense_linestrings", typeid(LineStringCollection));

            add_param(ParamInt(interval, "interval",   "Interval"));
        }
        void process();
    };

    class TinSimpNode :public Node {
        float thres_error = 2;
        float densify_interval = 2;
        bool create_triangles = true;
    public:
        using Node::Node;
        void init() {
            add_input("geometries", { typeid(PointCollection), typeid(LineStringCollection) });
            add_output("triangles", typeid(TriangleCollection));
            add_output("normals", typeid(vec3f));
            add_output("selected_lines", typeid(LineStringCollection));
            add_output("cgal_cdt", typeid(CDT));
            // add_output("count", typeid(vec1ui));
            // add_output("error", typeid(vec1f));

            add_param(ParamFloat(thres_error, "thres_error", "Error threshold"));
            add_param(ParamFloat(densify_interval, "densify_interval",  "Densify interval"));
            add_param(ParamBool(create_triangles, "create_triangles",  "Create triangles"));
        }
        /*void before_gui() {
            bool is_linestring = input("geometries").is_connected_type(typeid(LineStringCollection));
            auto param = std::get<ParamFloat>(parameters.at("densify_interval"));
            param.set_visible(is_linestring);
        }*/
        //void on_change_parameter(std::string name, ParameterVariant& param) {
        //    if (name == "thres_error")
        //        manager.run(*this);
        //}
        void process();
    };

    class TinSimpLASReaderNode :public Node {
        std::string filepath = "";
        int thin_nth = 5;
        float thres_error = 2;
        bool create_triangles = true;
    public:
        using Node::Node;
        void init() {
            add_output("triangles", typeid(TriangleCollection));
            add_output("cgal_cdt", typeid(CDT));

            add_param(ParamPath(filepath, "filepath", "Filepath"));
            add_param(ParamInt(thin_nth, "thin_nth",  "Thin factor"));
            add_param(ParamFloat(thres_error, "thres_error", "Error threshold"));
            add_param(ParamBool(create_triangles, "create_triangles",  "Create triangles"));
        }
        void process();
    };

    class SimplifyLine3DNode :public Node {
        float area_threshold = 0.1;
    public:
        using Node::Node;
        void init() {
            add_input("lines", typeid(LineStringCollection));
            add_output("lines", typeid(LineStringCollection));

            add_param(ParamFloat(area_threshold, "area_threshold", "Stop cost"));
        }
       /* void on_parameter_change(std::string name, ParameterVariant& param) {
            if (name == "area_threshold")
                manager.run(*this);
        }*/
        void process();
    };

    class SimplifyLineNode :public Node {
        float threshold_stop_cost = 0.1;
    public:
        using Node::Node;
        void init() {
            add_input("lines", typeid(LineStringCollection));
            add_output("lines", typeid(LineStringCollection));

            add_param(ParamFloat(threshold_stop_cost, "threshold_stop_cost",  "Stop cost"));
        }
       /* void on_parameter_change(std::string name, ParameterVariant& param) {
            if (name == "threshold_stop_cost")
                manager.run(*this);
        }*/
        void process();
    };

    class SimplifyLinesNode :public Node {
        float threshold_stop_cost = 0.1;
    public:
        using Node::Node;
        void init() {
            add_input("lines", typeid(LineStringCollection));
            add_output("lines", typeid(LineStringCollection));

            add_param(ParamFloat(threshold_stop_cost, "threshold_stop_cost",  "Stop cost"));
        }
        /*void on_parameter_change(std::string name, ParameterVariant& param) {
            if (name == "threshold_stop_cost")
                manager.run(*this);
        }*/
        void process();
    };

    class SimplifyFootprintsCDTNode :public Node {
        float threshold_stop_cost = 0.1;
    public:
        using Node::Node;
        void init() {
            add_input("polygons", typeid(LinearRingCollection));
            add_output("polygons_simp", typeid(LinearRingCollection));

            add_param(ParamFloat(threshold_stop_cost, "threshold_stop_cost",  "Stop cost"));
        }
        /*void on_parameter_change(std::string name, ParameterVariant& param) {
            if (name == "threshold_stop_cost")
                manager.run(*this);
        }*/
        void process();
    };

    class PLYWriterNode :public Node {
        std::string filepath = "";
        bool write_binary = false;
    public:
        bool multiple_files = true;

        using Node::Node;
        void init() {
            add_input("points", typeid(PointCollection)); //TT_point_collection_list
            add_input("labels", typeid(vec1i));

            add_param(ParamPath(filepath, "filepath",  "Filepath"));
            add_param(ParamBool(write_binary, "write_binary",  "Binary output"));
        }
        void process();
    };

    class PLYReaderNode :public Node {
        std::string filepath = "out.ply";
    public:
        using Node::Node;
        void init() {
            add_output("points", typeid(PointCollection)); //TT_point_collection_list
            add_output("normals", typeid(vec3f));

            add_param(ParamPath(filepath, "filepath",  "Filepath"));
        }
        void process();
    };

    class IsoLineNode :public Node {
        float interval = 1.0;
        std::pair<float, float> exclude_interval = { -0.5, 0.5 };
    public:
        using Node::Node;
        void init() {
            add_input("cgal_cdt", typeid(CDT));
            add_input("min", typeid(float));
            add_input("max", typeid(float));
            add_output("lines", typeid(LineStringCollection));
            add_output("attributes", typeid(vec1i));

            add_param(ParamFloat(interval, "interval",  "Interval"));
            add_param(ParamFloatRange(exclude_interval, "exclude",  "Exclude values"));
        }
        void process();
    };

    class IsoLineSlicerNode :public Node {
    public:
        using Node::Node;
        void init() {
            add_input("cgal_cdt", typeid(CDT));
            add_output("lines", typeid(LineStringCollection));
            add_output("attributes", typeid(AttributeMap));
        }
        void process();
    };

    class LineHeightNode :public Node {
        std::string filepath = "";
        int thin_nth = 5;
    public:
        using Node::Node;
        void init() {
            add_input("lines", typeid(LineStringCollection));
            add_output("lines", typeid(LineStringCollection));

            add_param(ParamPath(filepath, "filepath", "Filepath"));
            add_param(ParamInt(thin_nth, "thin_nth",   "Thin factor"));
        }
        void process();
    };

    class LineHeightCDTNode :public Node {
        bool add_bbox = false;
        float densify_interval = 2;
    public:
        using Node::Node;
        void init() {
            add_input("cgal_cdt", typeid(CDT));
            add_input("lines", typeid(LineStringCollection));
            add_output("lines", typeid(LineStringCollection));

            add_param(ParamBool(add_bbox, "add_bbox",  "Add bounding box to lines"));
            add_param(ParamFloat(densify_interval, "densify_interval",   "Line densify"));
        }
        void process();
    };

    class SimplifyLinesBufferNode :public Node {
        float threshold = 1.0;
    public:
        using Node::Node;
        void init() {
            add_input("polygons", typeid(LinearRingCollection));
            add_output("polygons_simp", typeid(LinearRingCollection));

            add_param(ParamFloat(threshold, "threshold",   "Threshold"));
        }
        void process();
    };

    class CDTAddConstraintNode :public Node {
    public:
        using Node::Node;
        void init() {
            add_input("cgal_cdt", typeid(CDT));
            add_vector_input("lines", typeid(LineString));
            add_output("cgal_cdt", typeid(CDT));
        }
        void process();
    };

    class CDT2TrianglesNode :public Node {
    public:
        using Node::Node;
        void init() {
            add_input("cgal_cdt", typeid(CDT));
            add_output("triangles", typeid(TriangleCollection));
        }
        void process();
    };

    class OBJWriterNode :public Node {
        std::string filepath;
    public:
        using Node::Node;
        void init() {
            add_input("triangles", typeid(TriangleCollection));

            add_param(ParamPath(filepath, "filepath",  "File path"));
        }
        void process();
    };

    class VecOBJWriterNode :public Node {
        std::string filepath;
        bool no_offset = false;
    public:
        using Node::Node;
        void init() {
            add_vector_input("triangles", typeid(TriangleCollection));

            add_param(ParamPath(filepath, "filepath",  "File path"));
            add_param(ParamBool(no_offset, "no_offset",  "Do not apply global offset"));
        }
        void process();
    };

    

    class CGALAlphaShapeR :public Node
    {
    public:
        using Node::Node;
        float alpha_value;
        bool sim_on = true;
        bool write_2_file_on = false;
        void init()
        {
            add_input("points", typeid(PointCollection));
            add_output("boundary_points", typeid(PointCollection));
            add_output("ground_points", typeid(PointCollection));
            add_output("boundary_seg", typeid(SegmentCollection));
            add_output("boundary_rings", typeid(LinearRingCollection));
            add_param(ParamBool(sim_on, "sim_on",  "turn on simplification"));
            add_param(ParamFloat(alpha_value, "alpha_value",  "alpha_value"));
            add_param(ParamBool(write_2_file_on, "write to file",  "Write result to file"));
            add_vector_output("id", typeid(int));
            
        }
        bool if_checkbox = true;
       
        void process();

        LinearRing simplify_footprint(const LinearRing& polygon, float& threshold_stop_cost) {
            namespace PS = CGAL::Polyline_simplification_2;
            typedef CGAL::Exact_predicates_inexact_constructions_kernel K;
            typedef K::Point_2 Point_2;
            typedef CGAL::Polygon_2<K>                   Polygon_2;
            typedef PS::Stop_below_count_ratio_threshold Stop_count_ratio;
            typedef PS::Stop_above_cost_threshold        Stop_cost;
            typedef PS::Squared_distance_cost            Cost;

            if (polygon.size() > 2) {
                Polygon_2 cgal_polygon;
                Cost cost;

                for (auto& p : polygon) {
                    cgal_polygon.push_back(Point_2(p[0], p[1]));
                }
                // cgal_polygon.erase(cgal_polygon.vertices_end()-1); // remove repeated point from the boost polygon

                // polygon = PS::simplify(polygon, cost, Stop_count_ratio(0.5));

                cgal_polygon = PS::simplify(cgal_polygon, cost, Stop_cost(threshold_stop_cost));

                LinearRing footprint_vec3f;
                for (auto v = cgal_polygon.vertices_begin(); v != cgal_polygon.vertices_end(); v++) {
                    footprint_vec3f.push_back({ float(v->x()),float(v->y()),0 });
                }

                // HACK: CGAL does not seem to remove the first point of the input polygon in any case, so we need to check ourselves
                auto p_0 = *(cgal_polygon.vertices_begin());
                auto p_1 = *(cgal_polygon.vertices_begin() + 1);
                auto p_end = *(cgal_polygon.vertices_end() - 1);
                // check the distance between the first vertex and the line between its 2 neighbours
                if (CGAL::squared_distance(Point_2(p_0), K::Segment_2(p_end, p_1)) < threshold_stop_cost) {
                    footprint_vec3f.erase(footprint_vec3f.begin());
                }

                return footprint_vec3f;
            }
            else
                return polygon;
        }

    };

    class OverlapCheckNode: public Node {
    public:
        using Node::Node;
        void init() {
           add_vector_input("linear_rings", typeid(LinearRing));
           add_vector_input("boundary_rings", typeid(LinearRingCollection));
           //add_output("bag_types", typeid(vec1i));
           add_vector_output("bag_type", typeid(int));
           //add_output("bag_types", typeid(std::vector<char>));

        }
        void process();
    };
    class ConvertLinearCollection2LinearRing : public Node {
    public:
      using Node::Node;
      void init() {
        add_vector_input("boundary_rings", typeid(LinearRingCollection));
        add_vector_output("linear_rings", typeid(LinearRing));
      }
      void process();
    };

}