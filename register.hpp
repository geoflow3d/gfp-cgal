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
}

namespace geoflow::nodes::las {
  NodeRegisterHandle create_register() {
    auto R = NodeRegister::create(GF_PLUGIN_NAME);
    register_nodes(*R);
    return R;
  }
}