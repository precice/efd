#include "XmlConfigurationParser.hpp"

#include <Uni/ExecutionControl/exception>
#include <Uni/Helpers/parsefromstring>
#include <Uni/Logging/macros>

#include <Eigen/Core>

#include <libxml/parser.h>

#include <boost/filesystem/fstream.hpp>
#include <boost/regex.hpp>

using FsiSimulation::EntryPoint::XmlConfigurationParser;
using FsiSimulation::FluidSimulation::Configuration;
using FsiSimulation::FluidSimulation::OutputEnum;
using FsiSimulation::FluidSimulation::ScalarEnum;
using FsiSimulation::FluidSimulation::SolverEnum;

namespace FsiSimulation {
namespace EntryPoint {
namespace Private {
template <typename T>
inline void
parse_in(xmlChar const* content, T& number) {
  if (!Uni::Helpers::parse_ineger_number(
        std::string(reinterpret_cast<char const*>(content)),
        number)) {
    throwException("Failed to parse integer number "
                   "'{1}'", content);
  }
}

template <typename T>
inline void
parse_fpn(xmlChar const* content, T& number) {
  if (!Uni::Helpers::parse_floating_point_number(
        std::string(reinterpret_cast<char const*>(content)),
        number)) {
    throwException("Failed to parse floating-point number "
                   "'{1}'", content);
  }
}

template <typename T, int D>
inline int
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
inline int
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

static
Configuration::UniqueWallType
parseWall(xmlNodePtr node) {
  static xmlChar* const _type               = (xmlChar* const)"type";
  static xmlChar* const _velocity           = (xmlChar* const)"velocity";
  static xmlChar* const _typeInput          = (xmlChar* const)"Input";
  static xmlChar* const _typeParabolicInput = (xmlChar* const)"ParabolicInput";
  static xmlChar* const _typeOutput         = (xmlChar* const)"Output";

  xmlAttrPtr attr = node->properties;
  xmlChar*   type;
  xmlChar*   velocity_string;

  while (attr) {
    if (xmlStrcasecmp(attr->name, _type) == 0) {
      type = attr->children->content;
    } else if (xmlStrcasecmp(attr->name, _velocity) == 0) {
      velocity_string = attr->children->content;
    }
    attr = attr->next;
  }

  Configuration::UniqueWallType wall;

  if (xmlStrcasecmp(type, _typeInput) == 0) {
    Configuration::Wall::VectorDs velocity;
    Private::parse_fpv(velocity_string, velocity);
    wall.reset(new Configuration::Input(velocity));
  } else if (xmlStrcasecmp(type, _typeParabolicInput) == 0) {
    Configuration::Wall::VectorDs velocity;
    Private::parse_fpv(velocity_string, velocity);
    wall.reset(new Configuration::ParabolicInput(velocity));
  } else if (xmlStrcasecmp(type, _typeOutput) == 0) {
    wall.reset(new Configuration::Output());
  }

  return wall;
}

static void
parseWallsChildren(xmlNodePtr     node,
                   Configuration* configuration) {
  xmlNodePtr currentNode = node->children;

  while (currentNode) {
    if (xmlStrcasecmp(currentNode->name, (xmlChar* const)"left") == 0) {
      configuration->walls[0][0] = parseWall(currentNode);
    } else if (xmlStrcasecmp(currentNode->name, (xmlChar* const)"right") == 0) {
      configuration->walls[0][1] = parseWall(currentNode);
    } else if (xmlStrcasecmp(currentNode->name, (xmlChar* const)"bottom")
               == 0) {
      configuration->walls[1][0] = parseWall(currentNode);
    } else if (xmlStrcasecmp(currentNode->name, (xmlChar* const)"top") == 0) {
      configuration->walls[1][1] = parseWall(currentNode);
    } else if (xmlStrcasecmp(currentNode->name, (xmlChar* const)"back") == 0) {
      configuration->walls[2][0] = parseWall(currentNode);
    } else if (xmlStrcasecmp(currentNode->name, (xmlChar* const)"front") == 0) {
      configuration->walls[2][1] = parseWall(currentNode);
    }
    currentNode = currentNode->next;
  }
}

static void
parseImmersedBoundary(xmlNodePtr     node,
                      Configuration* configuration) {
  static xmlChar const* const outerLayerSize = (xmlChar const*)"outerLayerSize";
  static xmlChar const* const innerLayerSize = (xmlChar const*)"innerLayerSize";

  xmlAttrPtr attr = node->properties;

  while (attr) {
    if (xmlStrcasecmp(attr->name, outerLayerSize) == 0) {
      Private::parse_in(attr->children->content,
                        configuration->outerLayerSize);
    } else if (xmlStrcasecmp(attr->name, innerLayerSize) == 0) {
      Private::parse_in(attr->children->content,
                        configuration->innerLayerSize);
    }
    attr = attr->next;
  }
}

static void
parseScenarioChildren(xmlNodePtr     node,
                      Configuration* configuration) {
  xmlNodePtr currentNode = node;

  while (currentNode) {
    if (xmlStrcasecmp(currentNode->name, (xmlChar* const)"walls") == 0) {
      parseWallsChildren(currentNode, configuration);
    } else if (xmlStrcasecmp(currentNode->name,
                             (xmlChar* const)"immersed-boundary") == 0) {
      parseImmersedBoundary(currentNode, configuration);
    }
    currentNode = currentNode->next;
  }
}

static void
parseScalarType(std::string    type_string,
                Configuration* configuration) {
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

static void
parseSolverType(std::string    type_string,
                Configuration* configuration) {
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

static void
parseOutputType(std::string    type_string,
                Configuration* configuration) {
  static boost::regex type1_regex("Xdmf", boost::regex::icase);
  static boost::regex type2_regex("Vtk", boost::regex::icase);

  if (boost::regex_search(type_string, type1_regex)) {
    configuration->outputType = OutputEnum::Xdmf;
  } else if (boost::regex_search(type_string, type2_regex)) {
    configuration->outputType = OutputEnum::Vtk;
  }
}

static void
parseScenarioParameters(xmlNodePtr     node,
                        Configuration* configuration) {
  static xmlChar* const re             = (xmlChar* const)"re";
  static xmlChar* const timeLimit      = (xmlChar* const)"timeLimit";
  static xmlChar* const iterationLimit = (xmlChar* const)"iterationLimit";
  static xmlChar* const plotInterval   = (xmlChar* const)"plotInterval";
  static xmlChar* const tau            = (xmlChar* const)"tau";
  static xmlChar* const gamma          = (xmlChar* const)"gamma";
  static xmlChar* const width          = (xmlChar* const)"width";
  static xmlChar* const size           = (xmlChar* const)"size";
  static xmlChar* const filename       = (xmlChar* const)"filename";
  static xmlChar* const scalarType     = (xmlChar* const)"scalar";
  static xmlChar* const solverType     = (xmlChar* const)"solver";
  static xmlChar* const outputType     = (xmlChar* const)"output";
  static xmlChar* const parallelizationSize
    = (xmlChar* const)"parallelizationSize";
  static xmlChar* const environment = (xmlChar* const)"environment";

  xmlAttrPtr attr = node->properties;

  int dim = -1;

  while (attr) {
    if (xmlStrcasecmp(attr->name, re) == 0) {
      parse_fpn(attr->children->content,
                configuration->re);
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
      parseScalarType(((const char*)attr->children->content), configuration);
    } else if (xmlStrcasecmp(attr->name, solverType) == 0) {
      parseSolverType(((const char*)attr->children->content), configuration);
    } else if (xmlStrcasecmp(attr->name, outputType) == 0) {
      parseOutputType(((const char*)attr->children->content), configuration);
    }
    attr = attr->next;
  }

  if (dim == -1) {
    throwException("Unknown dimension size");
  }

  configuration->dimensions = dim;
}
}
}
}

void
XmlConfigurationParser::
parse(std::unique_ptr<Configuration> const& configuration,
      boost::filesystem::path const&        filePath) {
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

  Private::parseScenarioParameters(root, configuration.get());

  Private::parseScenarioChildren(root->children, configuration.get());

  xmlFreeDoc(doc);

  xmlCleanupParser();
}
