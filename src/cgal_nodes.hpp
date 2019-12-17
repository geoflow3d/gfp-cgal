#pragma once

#include <geoflow/geoflow.hpp>
#include "tinsimp.hpp"

namespace geoflow::nodes::cgal
{

typedef tinsimp::CDT CDT;

class CDTNode : public Node
{
  bool create_triangles = false;

public:
  using Node::Node;
  void init()
  {
    add_input("geometries", {typeid(PointCollection), typeid(LineStringCollection)});
    add_output("cgal_cdt", typeid(CDT));
  }
  void process();
};

class DTNode : public Node
{
public:
  using Node::Node;
  void init()
  {
    add_input("points", typeid(PointCollection));
    add_output("cgal_dt", typeid(CDT));
  }
  void process();
};

class ComparePointDistanceNode : public Node
{
  std::string las_filepath = "";
  std::string log_filepath = "";
  int thin_nth = 20;

public:
  using Node::Node;
  void init()
  {
    add_input("triangles1_vec3f", typeid(vec3f));
    add_input("triangles2_vec3f", typeid(vec3f));
    add_output("points", typeid(vec3f));
    add_output("distances1", typeid(vec1f));
    add_output("distances2", typeid(vec1f));
    add_output("diff", typeid(vec1f));

    add_param("las_filpath", ParamPath(las_filepath, "LAS path"));
    add_param("log_filpath", ParamPath(log_filepath, "LOG path"));
    add_param("thin_nth", ParamBoundedInt(thin_nth, 0, 100, "Thin factor"));
  }
  void process();
};

class PointDistanceNode : public Node
{
  std::string filepath = "";
  int thin_nth = 5;
  bool overwritez = false;

public:
  using Node::Node;
  void init()
  {
    add_input("triangles", typeid(TriangleCollection));
    add_output("points", typeid(PointCollection));
    add_output("distances", typeid(vec1f));
    add_output("distance_min", typeid(float));
    add_output("distance_max", typeid(float));

    add_param("filepath", ParamPath(filepath, "Filepath"));
    add_param("thin_nth", ParamBoundedInt(thin_nth, 0, 100, "Thin factor"));
    add_param("overwritez", ParamBool(overwritez, "Overwrite Z"));
  }
  void process();
};

class CDTDistanceNode : public Node
{
public:
  using Node::Node;
  void init()
  {
    add_input("cgal_cdt_base", typeid(CDT));
    add_input("cgal_cdt_target", typeid(CDT));
    add_output("points", typeid(PointCollection));
    add_output("distance_min", typeid(float));
    add_output("distance_max", typeid(float));
  }
  void process();
};

class DensifyNode : public Node
{
  int interval = 2;

public:
  using Node::Node;
  void init()
  {
    add_input("geometries", {typeid(LineStringCollection)});
    add_output("dense_linestrings", typeid(LineStringCollection));

    add_param("interval", ParamBoundedInt(interval, 0, 100, "Interval"));
  }
  void process();
};

class TinSimpNode : public Node
{
  float thres_error = 2;
  float densify_interval = 2;
  bool create_triangles = true;

public:
  using Node::Node;
  void init()
  {
    add_input("geometries", {typeid(PointCollection), typeid(LineStringCollection)});
    add_output("triangles", typeid(TriangleCollection));
    add_output("normals", typeid(vec3f));
    add_output("selected_lines", typeid(LineStringCollection));
    add_output("cgal_cdt", typeid(CDT));
    // add_output("count", typeid(vec1ui));
    // add_output("error", typeid(vec1f));

    add_param("thres_error", ParamFloat(thres_error, "Error threshold"));
    add_param("densify_interval", ParamFloat(densify_interval, "Densify interval"));
    add_param("create_triangles", ParamBool(create_triangles, "Create triangles"));
  }
  void before_gui()
  {
    bool is_linestring = input("geometries").is_connected_type(typeid(LineStringCollection));
    auto param = std::get<ParamFloat>(parameters.at("densify_interval"));
    param.set_visible(is_linestring);
  }
  void on_change_parameter(std::string name, ParameterVariant &param)
  {
    if (name == "thres_error")
      manager.run(*this);
  }
  void process();
};

class TinSimpLASReaderNode : public Node
{
  std::string filepath = "";
  int thin_nth = 5;
  float thres_error = 2;
  bool create_triangles = true;

public:
  using Node::Node;
  void init()
  {
    add_output("triangles", typeid(TriangleCollection));
    add_output("cgal_cdt", typeid(CDT));

    add_param("filepath", ParamPath(filepath, "Filepath"));
    add_param("thin_nth", ParamBoundedInt(thin_nth, 0, 100, "Thin factor"));
    add_param("thres_error", ParamBoundedFloat(thres_error, 0, 100, "Error threshold"));
    add_param("create_triangles", ParamBool(create_triangles, "Create triangles"));
  }
  void process();
};

class SimplifyLine3DNode : public Node
{
  float area_threshold = 0.1;

public:
  using Node::Node;
  void init()
  {
    add_input("lines", typeid(LineStringCollection));
    add_output("lines", typeid(LineStringCollection));

    add_param("area_threshold", ParamFloat(area_threshold, "Stop cost"));
  }
  void on_parameter_change(std::string name, ParameterVariant &param)
  {
    if (name == "area_threshold")
      manager.run(*this);
  }
  void process();
};

class SimplifyLineNode : public Node
{
  float threshold_stop_cost = 0.1;

public:
  using Node::Node;
  void init()
  {
    add_input("lines", typeid(LineStringCollection));
    add_output("lines", typeid(LineStringCollection));

    add_param("threshold_stop_cost", ParamFloat(threshold_stop_cost, "Stop cost"));
  }
  void on_parameter_change(std::string name, ParameterVariant &param)
  {
    if (name == "threshold_stop_cost")
      manager.run(*this);
  }
  void process();
};

class SimplifyLinesNode : public Node
{
  float threshold_stop_cost = 0.1;

public:
  using Node::Node;
  void init()
  {
    add_input("lines", typeid(LineStringCollection));
    add_output("lines", typeid(LineStringCollection));

    add_param("threshold_stop_cost", ParamFloat(threshold_stop_cost, "Stop cost"));
  }
  void on_parameter_change(std::string name, ParameterVariant &param)
  {
    if (name == "threshold_stop_cost")
      manager.run(*this);
  }
  void process();
};

class SimplifyFootprintsCDTNode : public Node
{
  float threshold_stop_cost = 0.1;

public:
  using Node::Node;
  void init()
  {
    add_input("polygons", typeid(LinearRingCollection));
    add_output("polygons_simp", typeid(LinearRingCollection));

    add_param("threshold_stop_cost", ParamFloat(threshold_stop_cost, "Stop cost"));
  }
  void on_parameter_change(std::string name, ParameterVariant &param)
  {
    if (name == "threshold_stop_cost")
      manager.run(*this);
  }
  void process();
};

class PLYWriterNode : public Node
{
  std::string filepath = "";
  bool write_binary = false;

public:
  bool multiple_files = true;

  using Node::Node;
  void init()
  {
    add_input("points", typeid(PointCollection)); //TT_point_collection_list
    add_input("labels", typeid(vec1i));

    add_param("filepath", ParamPath(filepath, "Filepath"));
    add_param("write_binary", ParamBool(write_binary, "Binary output"));
  }
  void process();
};

class PLYReaderNode : public Node
{
  std::string filepath = "out.ply";

public:
  using Node::Node;
  void init()
  {
    add_output("points", typeid(PointCollection)); //TT_point_collection_list
    add_output("normals", typeid(vec3f));

    add_param("filepath", ParamPath(filepath, "Filepath"));
  }
  void process();
};

class IsoLineNode : public Node
{
  float interval = 1.0;
  std::pair<float, float> exclude_interval = {-0.5, 0.5};

public:
  using Node::Node;
  void init()
  {
    add_input("cgal_cdt", typeid(CDT));
    add_input("min", typeid(float));
    add_input("max", typeid(float));
    add_output("lines", typeid(LineStringCollection));
    add_output("attributes", typeid(vec1i));

    add_param("interval", ParamBoundedFloat(interval, 1, 100, "Interval"));
    add_param("exclude", ParamFloatRange(exclude_interval, "Exclude values"));
  }
  void process();
};

class IsoLineSlicerNode : public Node
{
public:
  using Node::Node;
  void init()
  {
    add_input("cgal_cdt", typeid(CDT));
    add_output("lines", typeid(LineStringCollection));
    add_output("attributes", typeid(AttributeMap));
  }
  void process();
};

class LineHeightNode : public Node
{
  std::string filepath = "";
  int thin_nth = 5;

public:
  using Node::Node;
  void init()
  {
    add_input("lines", typeid(LineStringCollection));
    add_output("lines", typeid(LineStringCollection));

    add_param("filepath", ParamPath(filepath, "Filepath"));
    add_param("thin_nth", ParamBoundedInt(thin_nth, 0, 100, "Thin factor"));
  }
  void process();
};

class LineHeightCDTNode : public Node
{
  bool add_bbox = false;
  float densify_interval = 2;

public:
  using Node::Node;
  void init()
  {
    add_input("cgal_cdt", typeid(CDT));
    add_input("lines", typeid(LineStringCollection));
    add_output("lines", typeid(LineStringCollection));

    add_param("add_bbox", ParamBool(add_bbox, "Add bounding box to lines"));
    add_param("densify_interval", ParamBoundedFloat(densify_interval, 0, 100, "Line densify"));
  }
  void process();
};

class SimplifyLinesBufferNode : public Node
{
  float threshold = 1.0;

public:
  using Node::Node;
  void init()
  {
    add_input("polygons", typeid(LinearRingCollection));
    add_output("polygons_simp", typeid(LinearRingCollection));

    add_param("threshold", ParamBoundedFloat(threshold, 0, 100, "Threshold"));
  }
  void process();
};

class CDTAddConstraintNode : public Node
{
public:
  using Node::Node;
  void init()
  {
    add_input("cgal_cdt", typeid(CDT));
    add_vector_input("lines", typeid(LineString));
    add_output("cgal_cdt", typeid(CDT));
  }
  void process();
};

class CDT2TrianglesNode : public Node
{
public:
  using Node::Node;
  void init()
  {
    add_input("cgal_cdt", typeid(CDT));
    add_output("triangles", typeid(TriangleCollection));
  }
  void process();
};

} // namespace geoflow::nodes::cgal