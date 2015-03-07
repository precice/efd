#ifndef FsiSimulation_FluidSimulation_XdmfHdf5Output_XdmfWriter_hpp
#define FsiSimulation_FluidSimulation_XdmfHdf5Output_XdmfWriter_hpp

#include "FluidSimulation/BasicCell.hpp"
#include "FluidSimulation/ParallelDistribution.hpp"

#include <Uni/Logging/macros>

#include <Eigen/Core>

#include <array>
#include <sstream>
#include <string>
#include <vector>

#include <boost/filesystem.hpp>
#include <boost/property_tree/ptree.hpp>
#include <boost/property_tree/xml_parser.hpp>

namespace FsiSimulation {
namespace FluidSimulation {
namespace XdmfHdf5Output {
template <int Dimensions>
class XdmfWriter {
  using VectorDi = Eigen::Matrix<int, Dimensions, 1>;

  using PropertyTree = boost::property_tree::ptree;

  using XmlWriterSettings =
          boost::property_tree::xml_writer_settings<std::string>;

  using Path = boost::filesystem::path;

  using ParallelDistributionType = ParallelDistribution<Dimensions>;

  struct DataItem {
    std::vector<int>          dimensions;
    static  std::string const numberType;
    static  std::string const format;
    std::string               data;

    void
    updateNode(PropertyTree& node) {
      std::stringstream ss;

      for (int d = 0; d < dimensions.size(); ++d) {
        if (ss.str().empty()) {
          ss << dimensions[d];
        } else {
          ss << ' ' << dimensions[d];
        }
      }

      PropertyTree& dat_node = node.add("DataItem", data);
      dat_node.put("<xmlattr>.Dimensions", ss.str());
      dat_node.put("<xmlattr>.NumberType", numberType);
      dat_node.put("<xmlattr>.Format", format);
    }
  };

  struct Topology {
    static  std::string const topologyType;
    std::vector<int>          dimensions;

    void
    updateNode(PropertyTree& node) {
      node.put("Topology.<xmlattr>.TopologyType", topologyType);
      std::stringstream ss;

      for (int d = 0; d < dimensions.size(); ++d) {
        if (ss.str().empty()) {
          ss << dimensions[d];
        } else {
          ss << ' ' << dimensions[d];
        }
      }

      node.put("Topology.<xmlattr>.Dimensions", ss.str());
    }
  };

  struct Geometry {
    static  std::string const        geometryType;
    std::array<DataItem, Dimensions> data_items;

    void
    updateNode(PropertyTree& node) {
      PropertyTree& geometry_node = node.add("Geometry", "");
      geometry_node.put("<xmlattr>.GeometryType", geometryType);

      for (auto& item : data_items) {
        item.updateNode(geometry_node);
      }
    }
  };

  struct Attribute {
    std::string name;
    std::string type;
    DataItem    item;

    static  std::string const center;

    void
    updateNode(PropertyTree& node) {
      PropertyTree& attr_node = node.add("Attribute", "");
      attr_node.put("<xmlattr>.Name", name);
      attr_node.put("<xmlattr>.AttributeType", type);
      attr_node.put("<xmlattr>.Center", center);
      item.updateNode(attr_node);
    }
  };

public:
  void
  initialize(std::string const&                             geometry_hdf_name,
             std::vector<std::string> const&                dimension_names,
             std::vector<FluidSimulation::Attribute> const& attributes,
             VectorDi const&                                dimensions) {
    for (int d = Dimensions - 1; d >= 0; --d) {
      auto& item = geometry.data_items[Dimensions - 1 - d];

      if (d < Dimensions) {
        topology.dimensions.push_back(dimensions(d) + 1);
        item.dimensions.push_back(dimensions(d) + 1);
      } else {
        // topology.dimensions.push_back(1);
        item.dimensions.push_back(1);
      }
      item.data = geometry_hdf_name + ":/" + dimension_names[d];
    }

    for (auto const& attribute : attributes) {
      logInfo("{1}", attribute.name);
      _attributes.emplace_back();
      _attributes.back().name = attribute.name;

      for (int d = Dimensions - 1; d >= 0; --d) {
        _attributes.back().item.dimensions.push_back(dimensions(d));
      }

      int attribute_size = 1;

      if (attribute.type == FluidSimulation::Attribute::Type::Vector) {
        attribute_size = Dimensions;
        _attributes.back().type = "Vector";
        _attributes.back().item.dimensions.push_back(3);
      } else {
        _attributes.back().type = "Scalar";
      }
      std::string name = attribute.name;
    }
  }

  void
  write(Path const&        file_path,
        std::string const& attribute_hdf_name,
        double const&      time) {
    for (auto& attribute : _attributes) {
      attribute.item.data = attribute_hdf_name
                            + ":/"
                            + attribute.name;
    }

    PropertyTree  property_tree;
    PropertyTree& grid_node = property_tree.put("Xdmf.Domain.Grid", "");
    grid_node.put("<xmlattr>.Name", name);
    grid_node.put("<xmlattr>.GridType", gridType);
    grid_node.put("<xmlattr>.xml:id", "gid");
    grid_node.put("Time.<xmlattr>.Value", std::to_string(time));

    topology.updateNode(grid_node);
    geometry.updateNode(grid_node);

    for (auto& attribute : _attributes) {
      attribute.updateNode(grid_node);
    }

    boost::property_tree::xml_parser::write_xml(file_path.string(),
                                                property_tree,
                                                std::locale(),
                                                _settings);
  }

  void
  writeTemporal(Path const&                     file_path,
                std::vector<std::string> const& timesteps) {
    PropertyTree  property_tree;
    PropertyTree& xdmf_node = property_tree.put("Xdmf", "");
    xdmf_node.put("<xmlattr>.xmlns:xi", "http://www.w3.org/2001/XInclude");
    PropertyTree& grid_node = xdmf_node.put("Domain.Grid", "");
    grid_node.put("<xmlattr>.Name", "TimeGrid");
    grid_node.put("<xmlattr>.GridType", "Collection");
    grid_node.put("<xmlattr>.CollectionType", "Temporal");

    for (auto ts : timesteps) {
      PropertyTree& ts_node = grid_node.add("xi:include", "");
      ts_node.put("<xmlattr>.href", ts);
      ts_node.put("<xmlattr>.xpointer", "gid");
    }

    boost::property_tree::xml_parser::write_xml(file_path.string(),
                                                property_tree,
                                                std::locale(),
                                                _settings);
  }

private:
  Topology               topology;
  Geometry               geometry;
  std::vector<Attribute> _attributes;

  static  XmlWriterSettings const _settings;
  static  std::string const       name;
  static  std::string const       gridType;
};

template <int Dimensions>
typename XdmfWriter<Dimensions>::XmlWriterSettings const
XdmfWriter<Dimensions>::_settings(' ', 4, "utf-8");
template <int Dimensions>
std::string const XdmfWriter<Dimensions>::DataItem::numberType = "Float";
template <int Dimensions>
std::string const XdmfWriter<Dimensions>::DataItem::format = "HDF";

template <>
std::string const XdmfWriter<3>::Topology::topologyType = "3DRectMesh";
template <>
std::string const XdmfWriter<2>::Topology::topologyType = "2DRectMesh";
template <>
std::string const XdmfWriter<3>::Geometry::geometryType = "VXVYVZ";
template <>
std::string const XdmfWriter<2>::Geometry::geometryType = "VXVY";
template <int Dimensions>
std::string const XdmfWriter<Dimensions>::Attribute::center = "Cell";
template <int Dimensions>
std::string const XdmfWriter<Dimensions>::name = "Grid";
template <int Dimensions>
std::string const XdmfWriter<Dimensions>::gridType = "Uniform";
}
}
}

#endif
