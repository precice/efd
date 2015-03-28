#include "XmlConfigurationParser.hpp"

#include <Uni/ExecutionControl/exception>
#include <Uni/Logging/macros>

#include <Eigen/Core>

#include <libxml/parser.h>

#include <boost/filesystem/fstream.hpp>
#include <boost/locale.hpp>
#include <boost/regex.hpp>

#include <sstream>

using FsiSimulation::EntryPoint::XmlConfigurationParser;
using FsiSimulation::FluidSimulation::Configuration;
using FsiSimulation::FluidSimulation::OutputEnum;
using FsiSimulation::FluidSimulation::ScalarEnum;
using FsiSimulation::FluidSimulation::SolverEnum;

static
std::locale const locale(boost::locale::generator() ("en_US.UTF-8"));

template <typename TResult>
TResult
parseNumber(xmlChar* const content) {
  std::stringstream ss;
  TResult           result;
  ss.imbue(locale);
  ss << content;
  ss >> result;

  return result;
}

using Scalar = double;

template <typename TType>
using Vector =  Eigen::Matrix<TType, 3, 1>;

template <typename TResult>
Vector<TResult>
parseVector(xmlChar const* const content) {
  xmlChar const* subContent = content;
  std::size_t    size       = xmlStrlen(subContent);

  Vector<TResult> vector;

  int w = 0;

  while (w < vector.size()) {
    xmlChar const* tempSubContent = xmlStrchr(subContent, (xmlChar)(','));

    if (!tempSubContent) {
      tempSubContent = xmlStrchr(subContent, (xmlChar)';');
    }

    if (!tempSubContent) {
      if (size > 0) {
        std::stringstream ss;
        ss.imbue(locale);
        ss << subContent;
        ss >> vector(w);
      }
      break;
    }

    std::stringstream ss;
    ss.imbue(locale);

    std::size_t oldSize = size;
    size = xmlStrlen(tempSubContent);

    for (std::size_t i = 0; i < (oldSize - size); ++i) {
      ss << subContent[i];
    }

    ss >> vector(w);

    ++w;
    subContent = &tempSubContent[1];
    --size;
  }

  return vector;
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
  xmlChar*   velocity;

  while (attr) {
    if (xmlStrcasecmp(attr->name, _type) == 0) {
      type = attr->children->content;
    } else if (xmlStrcasecmp(attr->name, _velocity) == 0) {
      velocity = attr->children->content;
    }
    attr = attr->next;
  }

  Configuration::UniqueWallType wall;

  if (xmlStrcasecmp(type, _typeInput) == 0) {
    wall = Configuration::UniqueWallType(
      new Configuration::Input(parseVector<Scalar>(velocity)));
  } else if (xmlStrcasecmp(type, _typeParabolicInput) == 0) {
    wall = Configuration::UniqueWallType(
      new Configuration::ParabolicInput(parseVector<Scalar>(velocity)));
  } else if (xmlStrcasecmp(type, _typeOutput) == 0) {
    wall = Configuration::UniqueWallType(
      new Configuration::Output());
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
      configuration->outerLayerSize = parseNumber<unsigned>(
        attr->children->content);
    } else if (xmlStrcasecmp(attr->name, innerLayerSize) == 0) {
      configuration->innerLayerSize = parseNumber<unsigned>(
        attr->children->content);
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
  static xmlChar* const dim            = (xmlChar* const)"dim";
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

  while (attr) {
    if (xmlStrcasecmp(attr->name, re) == 0) {
      configuration->re = parseNumber<Scalar>(attr->children->content);
    } else if (xmlStrcasecmp(attr->name, timeLimit) == 0) {
      configuration->timeLimit = parseNumber<Scalar>(attr->children->content);
    } else if (xmlStrcasecmp(attr->name, iterationLimit) == 0) {
      configuration->iterationLimit = parseNumber<int>(attr->children->content);
    } else if (xmlStrcasecmp(attr->name, plotInterval) == 0) {
      configuration->plotInterval = parseNumber<Scalar>(
        attr->children->content);
    } else if (xmlStrcasecmp(attr->name, tau) == 0) {
      configuration->tau = parseNumber<Scalar>(attr->children->content);
    } else if (xmlStrcasecmp(attr->name, gamma) == 0) {
      configuration->gamma = parseNumber<Scalar>(attr->children->content);
    } else if (xmlStrcasecmp(attr->name, dim) == 0) {
      configuration->dimensions = parseNumber<int>(attr->children->content);
    } else if (xmlStrcasecmp(attr->name, width) == 0) {
      configuration->width = parseVector<Scalar>(attr->children->content);
    } else if (xmlStrcasecmp(attr->name, size) == 0) {
      configuration->size = parseVector<int>(attr->children->content);
    } else if (xmlStrcasecmp(attr->name, parallelizationSize) == 0) {
      configuration->parallelizationSize
        = parseVector<int>(attr->children->content);
    } else if (xmlStrcasecmp(attr->name, environment) == 0) {
      configuration->environment = parseVector<Scalar>(attr->children->content);
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

  parseScenarioParameters(root, configuration.get());

  parseScenarioChildren(root->children, configuration.get());

  xmlFreeDoc(doc);

  xmlCleanupParser();
}
