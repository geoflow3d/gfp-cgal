#include "cgal_nodes.hpp"

using namespace geoflow::nodes::cgal;

void register_nodes(geoflow::NodeRegister& node_register) {
    node_register.register_node<CDTNode>("CDT");
    node_register.register_node<DTNode>("DT");
    node_register.register_node<ComparePointDistanceNode>("ComparePointDistance");
    node_register.register_node<PointDistanceNode>("PointDistance");
    node_register.register_node<CDTDistanceNode>("CDTDistance");
    node_register.register_node<DensifyNode>("Densify");
    node_register.register_node<TinSimpNode>("TinSimp");
    node_register.register_node<TinSimpLASReaderNode>("TinSimpLASReader");
    node_register.register_node<SimplifyLine3DNode>("SimplifyLine3D");
    node_register.register_node<SimplifyLineNode>("SimplifyLine");
    node_register.register_node<SimplifyLinesNode>("SimplifyLines");
    node_register.register_node<SimplifyFootprintsCDTNode>("SimplifyFootprintsCDT");
    node_register.register_node<PLYWriterNode>("PLYWriter");
    node_register.register_node<PLYReaderNode>("PLYReader");
    node_register.register_node<IsoLineNode>("IsoLine");
    node_register.register_node<IsoLineSlicerNode>("IsoLineSlicer");
    node_register.register_node<LineHeightNode>("LineHeight");
    node_register.register_node<LineHeightCDTNode>("LineHeightCDT");
    node_register.register_node<SimplifyLinesBufferNode>("SimplifyLinesBuffer");
    node_register.register_node<CDTAddConstraintNode>("CDTAddConstraint");
    node_register.register_node<CDT2TrianglesNode>("CDT2Triangles");
    node_register.register_node<OBJWriterNode>("OBJWriter");
    node_register.register_node<VecOBJWriterNode>("OBJVecWriter");
    node_register.register_node<CGALAlphaShapeR>("CGALAlphaShape");
    node_register.register_node<OverlapCheckNode>("OverlapCheckNode");
    node_register.register_node<ConvertLinearCollection2LinearRing>("ConvertLinearCollection2LinearRing");
}