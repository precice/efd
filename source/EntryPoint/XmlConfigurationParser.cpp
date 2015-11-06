#include "XmlConfigurationParser.hpp"

#include <Uni/ExecutionControl/exception>
#include <Uni/Firewall/Interface>
#include <Uni/Helpers/parsefromstring>
#include <Uni/Logging/macros>

#include <Eigen/Core>

#include <libxml/parser.h>

#include <boost/filesystem/fstream.hpp>
#include <boost/regex.hpp>

using Fluid::EntryPoint::XmlConfigurationParser;
using Fluid::Simulation::Configuration;
using Fluid::Simulation::OutputEnum;
using Fluid::Simulation::ScalarEnum;
using Fluid::Simulation::SolverEnum;

namespace Fluid {
namespace EntryPoint {
//
static inline bool
is_default(xmlChar const* content) {
  std::string content_string
    = reinterpret_cast<char const*>(content);

  static boost::regex default_regex("default", boost::regex::icase);

  if (boost::regex_search(content_string, default_regex)) {
    return true;
  }

  return false;
}

static inline bool
parse_bool(xmlChar const* content) {
  static boost::regex on_regex("on|true|1", boost::regex::icase);
  static boost::regex off_regex("off|false|0", boost::regex::icase);

  std::string content_string
    = reinterpret_cast<char const*>(content);

  if (boost::regex_search(content_string, on_regex)) {
    return true;
  } else if (boost::regex_search(content_string, off_regex)) {
    return false;
  } else {
    throwException("Failed to parse boolean value");

    return false;
  }
}

template <typename T>
static inline void
parse_in(xmlChar const* content, T& number) {
  if (!Uni::Helpers::parse_integer_number(
        std::string(reinterpret_cast<char const*>(content)),
        number)) {
    throwException("Failed to parse integer number '{1}'", content);
  }
}

static inline int
parse_in(xmlChar const* content) {
  int result;
  parse_in(content, result);

  return result;
}

template <typename T>
static inline void
parse_fpn(xmlChar const* content, T& number) {
  if (!Uni::Helpers::parse_floating_point_number(
        std::string(reinterpret_cast<char const*>(content)),
        number)) {
    throwException("Failed to parse floating-point number "
                   "'{1}'", content);
  }
}

static inline long double
parse_fpn(xmlChar const* content) {
  long double result;
  parse_fpn(content, result);

  return result;
}

template <typename T, int D>
static inline int
parse_iv(xmlChar const* content,
         Eigen::Matrix<T, D, 1>& vector) {
  int result;

  if ((result = Uni::Helpers::parse_integer_vector(
         std::string(reinterpret_cast<char const*>(content)),
         vector)) == 0) {
    throwException("Failed to parse integer vector "
                   "'{1}'", content);
  }

  return result;
}

template <typename T, int D>
static inline int
parse_fpv(xmlChar const* content,
          Eigen::Matrix<T, D, 1>& vector) {
  int result;

  if ((result = Uni::Helpers::parse_floating_point_vector(
         std::string(reinterpret_cast<char const*>(content)),
         vector)) == 0) {
    throwException("Failed to parse floating-point vector "
                   "'{1}'", content);
  }

  return result;
}

class XmlConfigurationParserImplementation {
private:
  XmlConfigurationParserImplementation(XmlConfigurationParser* in)
    : _in(in) {}

  void
  parse() const {
    boost::filesystem::ifstream fileStream(filePath, std::ios_base::binary);

    std::filebuf* fileBuf = fileStream.rdbuf();

    auto fileSize = fileBuf->pubseekoff(0, fileStream.end, fileStream.in);
    fileBuf->pubseekpos(0, fileStream.in);

    char* data = new char[fileSize];

    fileBuf->sgetn(data, fileSize);

    fileStream.close();

    xmlDocPtr doc = xmlReadMemory(data,
                                  fileSize,
                                  filePath.c_str(),
                                  "UTF-8",
                                  0);

    if (doc == 0) {
      throwException("Failed to parse configuration");
    }

    xmlNodePtr root = xmlDocGetRootElement(doc);

    if (xmlStrcasecmp(root->name, (xmlChar* const)"scenario") != 0) {
      throwException("Failed to find root element");
    }

    parseScenarioParameters(root);

    parseScenarioChildren(root->children);

    xmlFreeDoc(doc);

    xmlCleanupParser();
  }

  void
  parseScenarioParameters(xmlNodePtr node) const {
    static xmlChar* const re
      = (xmlChar* const)"re";
    static xmlChar* const diffusion_multiplier
      = (xmlChar* const)"diffusion-multiplier";
    static xmlChar* const grad_pressure_multiplier
      = (xmlChar* const)"pressure-gradient-multiplier";
    static xmlChar* const timeLimit             = (xmlChar* const)"timeLimit";
    static xmlChar* const iterationLimit        = (xmlChar* const)"iterationLimit";
    static xmlChar* const plotInterval          = (xmlChar* const)"plotInterval";
    static xmlChar* const tau                   = (xmlChar* const)"tau";
    static xmlChar* const time_step_size_method = (xmlChar* const)"time-step-size-method";
    static xmlChar* const gamma                 = (xmlChar* const)"gamma";
    static xmlChar* const width                 = (xmlChar* const)"width";
    static xmlChar* const size                  = (xmlChar* const)"size";
    static xmlChar* const filename              = (xmlChar* const)"filename";
    static xmlChar* const scalarType            = (xmlChar* const)"scalar";
    static xmlChar* const solverType            = (xmlChar* const)"solver";
    static xmlChar* const outputType            = (xmlChar* const)"output";
    static xmlChar* const parallelizationSize
      = (xmlChar* const)"parallelizationSize";
    static xmlChar* const environment = (xmlChar* const)"environment";

    xmlAttrPtr attr = node->properties;

    int dim = -1;

    while (attr) {
      if (xmlStrcasecmp(attr->name, re) == 0) {
        configuration->set(
          "/Equations/Ins/ReynoldsNumber",
          (long double)(parse_fpn(attr->children->content)));
      } else if (xmlStrcasecmp(attr->name, diffusion_multiplier) == 0) {
        if (is_default(attr->children->content)) {
          configuration->set(
            "/Equations/Ins/DiffusionMultiplier",
            std::string("default"));
        } else {
          configuration->set(
            "/Equations/Ins/DiffusionMultiplier",
            (long double)(parse_fpn(attr->children->content)));
        }
      } else if (xmlStrcasecmp(attr->name, grad_pressure_multiplier) == 0) {
        if (is_default(attr->children->content)) {
          configuration->set(
            "/Equations/Ins/PressureGradientMultiplier",
            std::string("default"));
        } else {
          configuration->set(
            "/Equations/Ins/PressureGradientMultiplier",
            (long double)(parse_fpn(attr->children->content)));
        }
      } else if (xmlStrcasecmp(attr->name, timeLimit) == 0) {
        parse_fpn(attr->children->content,
                  configuration->timeLimit);
      } else if (xmlStrcasecmp(attr->name, iterationLimit) == 0) {
        parse_in(attr->children->content,
                 configuration->iterationLimit);
      } else if (xmlStrcasecmp(attr->name, plotInterval) == 0) {
        parse_fpn(attr->children->content,
                  configuration->plotInterval);
      } else if (xmlStrcasecmp(attr->name, tau) == 0) {
        parse_fpn(attr->children->content,
                  configuration->tau);
      } else if (xmlStrcasecmp(attr->name, time_step_size_method) == 0) {
        parseTimeStepSizeMethod(reinterpret_cast<char const*>(attr->children->content));
      } else if (xmlStrcasecmp(attr->name, gamma) == 0) {
        parse_fpn(attr->children->content,
                  configuration->gamma);
      } else if (xmlStrcasecmp(attr->name, width) == 0) {
        int temp_dim = parse_fpv(attr->children->content,
                                 configuration->width);

        if (dim == -1) {
          dim = temp_dim;
        } else if (temp_dim != dim) {
          throwException("Inconsistent dimension size in width vector,"
                         "'{1}', because there was already found "
                         "dimension size of '{2}'",
                         temp_dim, dim);
        }
      } else if (xmlStrcasecmp(attr->name, size) == 0) {
        int temp_dim = parse_iv(attr->children->content,
                                configuration->size);

        if (dim == -1) {
          dim = temp_dim;
        } else if (temp_dim != dim) {
          throwException("Inconsistent dimension size in size vector,"
                         "'{1}', because there was already found "
                         "dimension size of '{2}'",
                         temp_dim, dim);
        }
      } else if (xmlStrcasecmp(attr->name, parallelizationSize) == 0) {
        int temp_dim = parse_iv(attr->children->content,
                                configuration->parallelizationSize);

        if (dim == -1) {
          dim = temp_dim;
        } else if (temp_dim != dim) {
          throwException("Inconsistent dimension size in parallelization size "
                         "vector, '{1}', because there was already found "
                         "dimension size of '{2}'",
                         temp_dim, dim);
        }
      } else if (xmlStrcasecmp(attr->name, environment) == 0) {
        parse_fpv(attr->children->content,
                  configuration->environment);
      } else if (xmlStrcasecmp(attr->name, filename) == 0) {
        configuration->filename = (const char*)attr->children->content;
      } else if (xmlStrcasecmp(attr->name, scalarType) == 0) {
        parseScalarType(((const char*)attr->children->content));
      } else if (xmlStrcasecmp(attr->name, solverType) == 0) {
        parseSolverType(((const char*)attr->children->content));
      } else if (xmlStrcasecmp(attr->name, outputType) == 0) {
        parseOutputType(((const char*)attr->children->content));
      }
      attr = attr->next;
    }

    if (dim == -1) {
      throwException("Unknown dimension size");
    }

    configuration->dimensions = dim;
  }

  void
  parseTimeStepSizeMethod(std::string type_string) const {
    static boost::regex type1_regex("explicit", boost::regex::icase);

    static boost::regex type2_regex("implicit", boost::regex::icase);

    if (boost::regex_search(type_string, type1_regex)) {
      configuration->set("/Equations/Fsfd/TimeStepSizeMethod",
                         std::string("explicit"));
    } else if (boost::regex_search(type_string, type2_regex)) {
      configuration->set("/Equations/Fsfd/TimeStepSizeMethod",
                         std::string("implicit"));
    }
  }

  void
  parseScalarType(std::string type_string) const {
    static boost::regex type1_regex("Float", boost::regex::icase);
    static boost::regex type2_regex("Double", boost::regex::icase);
    static boost::regex type3_regex("Long[\\-\\s]*Double", boost::regex::icase);

    if (boost::regex_search(type_string, type1_regex)) {
      configuration->scalarType = ScalarEnum::Float;
    } else if (boost::regex_search(type_string, type2_regex)) {
      configuration->scalarType = ScalarEnum::Double;
    } else if (boost::regex_search(type_string, type3_regex)) {
      configuration->scalarType = ScalarEnum::LongDouble;
    }
  }

  void
  parseSolverType(std::string type_string) const {
    static boost::regex type1_regex(
      "Improved\\s*Fractional\\s*Step\\s*Finite\\s*Difference",
      boost::regex::icase);

    static boost::regex type2_regex(
      "Simple\\s*Fractional\\s*Step\\s*Finite\\s*Difference",
      boost::regex::icase);

    if (boost::regex_search(type_string, type1_regex)) {
      configuration->solverType = SolverEnum::Ifsfd;
    } else if (boost::regex_search(type_string, type2_regex)) {
      configuration->solverType = SolverEnum::Sfsfd;
    }
  }

  void
  parseOutputType(std::string type_string) const {
    static boost::regex type1_regex("Xdmf", boost::regex::icase);
    static boost::regex type2_regex("Vtk", boost::regex::icase);

    if (boost::regex_search(type_string, type1_regex)) {
      configuration->outputType = OutputEnum::Xdmf;
    } else if (boost::regex_search(type_string, type2_regex)) {
      configuration->outputType = OutputEnum::Vtk;
    }
  }

  void
  parseScenarioChildren(xmlNodePtr node) const {
    xmlNodePtr currentNode = node;

    while (currentNode) {
      if (xmlStrcasecmp(currentNode->name, (xmlChar* const)"walls") == 0) {
        parseWallsChildren(currentNode);
      } else if (xmlStrcasecmp(currentNode->name,
                               (xmlChar* const)"immersed-boundary") == 0) {
        parseImmersedBoundary(currentNode);
      }
      currentNode = currentNode->next;
    }
  }

  void
  parseWallsChildren(xmlNodePtr node) const {
    xmlNodePtr currentNode = node->children;

    while (currentNode) {
      if (xmlStrcasecmp(currentNode->name, (xmlChar* const)"left") == 0) {
        configuration->walls[0][0] = parseWall(currentNode);
      } else if (xmlStrcasecmp(currentNode->name, (xmlChar* const)"right")
                 == 0) {
        configuration->walls[0][1] = parseWall(currentNode);
      } else if (xmlStrcasecmp(currentNode->name, (xmlChar* const)"bottom")
                 == 0) {
        configuration->walls[1][0] = parseWall(currentNode);
      } else if (xmlStrcasecmp(currentNode->name, (xmlChar* const)"top") == 0) {
        configuration->walls[1][1] = parseWall(currentNode);
      } else if (xmlStrcasecmp(currentNode->name, (xmlChar* const)"back")
                 == 0) {
        configuration->walls[2][0] = parseWall(currentNode);
      } else if (xmlStrcasecmp(currentNode->name, (xmlChar* const)"front")
                 == 0) {
        configuration->walls[2][1] = parseWall(currentNode);
      }
      currentNode = currentNode->next;
    }
  }

  Configuration::UniqueWallType
  parseWall(xmlNodePtr node) const {
    static xmlChar* const _type      = (xmlChar* const)"type";
    static xmlChar* const _velocity  = (xmlChar* const)"velocity";
    static xmlChar* const _typeInput = (xmlChar* const)"Input";
    static xmlChar* const _typeParabolicInput
      = (xmlChar* const)"ParabolicInput";
    static xmlChar* const _typeOutput = (xmlChar* const)"Output";

    xmlAttrPtr                    attr = node->properties;
    xmlChar*                      type = 0;
    Configuration::Wall::VectorDs velocity
      = Configuration::Wall::VectorDs::Zero();

    while (attr) {
      if (xmlStrcasecmp(attr->name, _type) == 0) {
        type = attr->children->content;
      } else if (xmlStrcasecmp(attr->name, _velocity) == 0) {
        parse_fpv(attr->children->content, velocity);
      }
      attr = attr->next;
    }

    Configuration::UniqueWallType wall;

    if (type == 0) {
      throwException("No wall type was provided");
    }

    if (xmlStrcasecmp(type, _typeInput) == 0) {
      wall.reset(new Configuration::Input(velocity));
    } else if (xmlStrcasecmp(type, _typeParabolicInput) == 0) {
      wall.reset(new Configuration::ParabolicInput(velocity));
    } else if (xmlStrcasecmp(type, _typeOutput) == 0) {
      wall.reset(new Configuration::Output());
    } else {
      throwException("Unknown wall type '{1}'", type);
    }

    return wall;
  }

  void
  parseImmersedBoundary(xmlNodePtr node) const {
    static xmlChar const* const start_iteration
      = (xmlChar const*)"start-iteration";
    static xmlChar const* const type
      = (xmlChar const*)"type";
    static xmlChar const* const precice_configuration_path
      = (xmlChar const*)"precice-configuration-path";
    static xmlChar const* const full_prediction
      = (xmlChar const*)"full-prediction";
    static xmlChar const* const developing_structure
      = (xmlChar const*)"developing-structure";
    static xmlChar const* const coupling
      = (xmlChar const*)"coupling";
    static xmlChar const* const outerLayerSize
      = (xmlChar const*)"outerLayerSize";
    static xmlChar const* const innerLayerSize
      = (xmlChar const*)"innerLayerSize";
    static xmlChar const* const support_radius
      = (xmlChar const*)"support-radius";
    static xmlChar const* const imq_shape
      = (xmlChar const*)"imq-shape";
    static xmlChar const* const coupling_forces_name
      = (xmlChar const*)"coupling-forces-name";
    static xmlChar const* const structure_mesh_name
      = (xmlChar const*)"structure-mesh-name";
    static xmlChar const* const ib_structure_mesh_name
      = (xmlChar const*)"ib-structure-mesh-name";
    static xmlChar const* const structure_displacements_name
      = (xmlChar const*)"structure-dispacements-name";
    static xmlChar const* const ib_forces_name
      = (xmlChar const*)"immersed-boundary-forces-name";

    static boost::regex type1_regex("Precice[\\-\\s]*Based",
                                    boost::regex::icase);
    static boost::regex type2_regex("Rbf[\\-\\s]*Based", boost::regex::icase);

    xmlAttrPtr attr = node->properties;

    while (attr) {
      if (xmlStrcasecmp(attr->name, start_iteration) == 0) {
        unsigned value;
        parse_in(attr->children->content, value);
        configuration->set("/Ib/Options/StartIteration", value);
      } else if (xmlStrcasecmp(attr->name, type) == 0) {
        std::string type_string
          = reinterpret_cast<char const*>(attr->children->content);

        if (boost::regex_search(type_string, type1_regex)) {
          configuration->set(
            "/Ib/Schemes/DirectForcing/PreciceBased");
        } else if (boost::regex_search(type_string, type2_regex)) {
          configuration->set("/Ib/Schemes/DirectForcing/RbfBased");
        } else {
          throwException("Unknown type for Immersed Boundary");
        }
      } else if (xmlStrcasecmp(attr->name, precice_configuration_path) == 0) {
        boost::filesystem::path path(
          reinterpret_cast<char const*>(attr->children->content));

        if (path.is_relative()) {
          path = filePath.parent_path() / path;
        }
        configuration->set("/Ib/PreciceConfigurationPath",
                           path);
      } else if (xmlStrcasecmp(attr->name, full_prediction) == 0) {
        configuration->set("/Ib/Features/FullVelocityPrediction",
                           parse_bool(attr->children->content));
      } else if (xmlStrcasecmp(attr->name, developing_structure) == 0) {
        configuration->set("/Ib/Features/DevelopingStructure",
                           parse_bool(attr->children->content));
      } else if (xmlStrcasecmp(attr->name, coupling) == 0) {
        configuration->set("/Ib/Features/Coupling",
                           parse_bool(attr->children->content));
      } else if (xmlStrcasecmp(attr->name, outerLayerSize) == 0) {
        configuration->set(
          "/Ib/Schemes/DirectForcing/PreciceBased/OuterLayerSize",
          (unsigned)(parse_in(attr->children->content)));
      } else if (xmlStrcasecmp(attr->name, innerLayerSize) == 0) {
        configuration->set(
          "/Ib/Schemes/DirectForcing/PreciceBased/InnerLayerSize",
          (unsigned)(parse_in(attr->children->content)));
      } else if (xmlStrcasecmp(attr->name, support_radius) == 0) {
        if (is_default(attr->children->content)) {
          configuration->set(
            "/Ib/Schemes/DirectForcing/RbfBased/SupportRadius",
            std::string("default"));
        } else {
          configuration->set(
            "/Ib/Schemes/DirectForcing/RbfBased/SupportRadius",
            (long double)(parse_fpn(attr->children->content)));
        }
      } else if (xmlStrcasecmp(attr->name, imq_shape) == 0) {
        if (is_default(attr->children->content)) {
          configuration->set(
            "/Ib/Schemes/DirectForcing/RbfBased/ImqShape",
            std::string("default"));
        } else {
          configuration->set(
            "/Ib/Schemes/DirectForcing/RbfBased/ImqShape",
            (long double)(parse_fpn(attr->children->content)));
        }
      } else if (xmlStrcasecmp(attr->name, coupling_forces_name) == 0) {
          configuration->set(
            "/Ib/Options/CouplingForcesName",
            std::string(reinterpret_cast<char const*>(attr->children->content)));
      } else if (xmlStrcasecmp(attr->name, structure_mesh_name) == 0) {
          configuration->set(
            "/Ib/Options/StructureMeshName",
            std::string(reinterpret_cast<char const*>(attr->children->content)));
      } else if (xmlStrcasecmp(attr->name, ib_structure_mesh_name) == 0) {
          configuration->set(
            "/Ib/Options/IbStructureMeshName",
            std::string(reinterpret_cast<char const*>(attr->children->content)));
      } else if (xmlStrcasecmp(attr->name, structure_displacements_name) == 0) {
          configuration->set(
            "/Ib/Options/StructureDisplacementsName",
            std::string(reinterpret_cast<char const*>(attr->children->content)));
      } else if (xmlStrcasecmp(attr->name, ib_forces_name) == 0) {
          configuration->set(
            "/Ib/Options/ForcesName",
            std::string(reinterpret_cast<char const*>(attr->children->content)));
      }
      attr = attr->next;
    }
  }

  void
  parsePoints(xmlNodePtr node) const {
    static boost::regex regexpr("[\\[({<]\\+([^\\])}>]*)[^\\])}>]\\+");
    using Vector = Eigen::Matrix<long double, 3, 1>;
    std::vector<Vector> points;

    std::string content_string = reinterpret_cast<char const*>(node->children->content);

    boost::sregex_iterator it
      = boost::sregex_iterator(content_string.begin(), content_string.end(), regexpr);
    boost::sregex_iterator it_end = boost::sregex_iterator();

    for (; it != it_end; ++it) {
      Vector coordinates;
      Uni::Helpers::parse_floating_point_vector(
        it->
        operator[](1).str(),
        coordinates);
      points.emplace_back(coordinates);
    }
  }

  Simulation::Configuration* configuration;
  boost::filesystem::path         filePath;

  Uni_Firewall_INTERFACE_LINK(XmlConfigurationParser);
};
}
}

XmlConfigurationParser::
XmlConfigurationParser(
  std::unique_ptr<Simulation::Configuration> const& configuration,
  boost::filesystem::path const&                         file_path)
  : _im(new XmlConfigurationParserImplementation(this)) {
  _im->configuration = configuration.get();
  _im->filePath      = file_path;
  _im->parse();
}

XmlConfigurationParser::
~XmlConfigurationParser() {}
