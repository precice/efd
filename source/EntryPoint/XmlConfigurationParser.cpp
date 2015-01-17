#include "XmlConfigurationParser.hpp"

#include <Uni/ExecutionControl/exception>
#include <Uni/Logging/macros>

#include <Eigen/Core>

#include <libxml/parser.h>

#include <boost/filesystem/fstream.hpp>
#include <boost/locale.hpp>

#include <sstream>

using FsiSimulation::EntryPoint::XmlConfigurationParser;
using FsiSimulation::FluidSimulation::Configuration;

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
    if (xmlStrEqual(attr->name, _type)) {
      type = attr->children->content;
    } else if (xmlStrEqual(attr->name, _velocity)) {
      velocity = attr->children->content;
    }
    attr = attr->next;
  }

  Configuration::UniqueWallType wall;

  if (xmlStrEqual(type, _typeInput)) {
    wall = Configuration::UniqueWallType(
      new Configuration::Input(parseVector<Scalar>(velocity)));
  } else if (xmlStrEqual(type, _typeParabolicInput)) {
    wall = Configuration::UniqueWallType(
      new Configuration::ParabolicInput(parseVector<Scalar>(velocity)));
  } else if (xmlStrEqual(type, _typeOutput)) {
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
    if (xmlStrEqual(currentNode->name, (xmlChar* const)"left")) {
      configuration->walls[0][0] = parseWall(currentNode);
    } else if (xmlStrEqual(currentNode->name, (xmlChar* const)"right")) {
      configuration->walls[0][1] = parseWall(currentNode);
    } else if (xmlStrEqual(currentNode->name, (xmlChar* const)"bottom")) {
      configuration->walls[1][0] = parseWall(currentNode);
    } else if (xmlStrEqual(currentNode->name, (xmlChar* const)"top")) {
      configuration->walls[1][1] = parseWall(currentNode);
    } else if (xmlStrEqual(currentNode->name, (xmlChar* const)"back")) {
      configuration->walls[2][0] = parseWall(currentNode);
    } else if (xmlStrEqual(currentNode->name, (xmlChar* const)"front")) {
      configuration->walls[2][1] = parseWall(currentNode);
    }
    currentNode = currentNode->next;
  }
}

static void
parseImmersedBoundary(xmlNodePtr     node,
                      Configuration* configuration) {
  static xmlChar const* const methodAttrName = (xmlChar const*)"method";
  static xmlChar const* const feedbackForcingMethod
    = (xmlChar const*)"FeedbackForcing";
  static xmlChar const* const alphaAttrName = (xmlChar const*)"alpha";

  xmlAttrPtr attr = node->properties;

  while (attr) {
    if (xmlStrEqual(attr->name, methodAttrName)) {
      if (xmlStrEqual(attr->children->content, feedbackForcingMethod)) {
        configuration->immersedBoundaryMethod
          = Configuration::FeedbackForcingMethod;
      }
    } else if (xmlStrEqual(attr->name, alphaAttrName)) {
      configuration->alpha = parseNumber<double>(attr->children->content);
    }
    attr = attr->next;
  }
}

static void
parseScenarioChildren(xmlNodePtr     node,
                      Configuration* configuration) {
  xmlNodePtr currentNode = node;

  while (currentNode) {
    if (xmlStrEqual(currentNode->name, (xmlChar* const)"walls")) {
      parseWallsChildren(currentNode, configuration);
    } else if (xmlStrEqual(currentNode->name, (xmlChar* const)"immersed-boundary")) {
      parseImmersedBoundary(currentNode, configuration);
    }
    currentNode = currentNode->next;
  }
}

static void
parseScenarioParameters(xmlNodePtr     node,
                        Configuration* configuration) {
  static xmlChar* const re                  = (xmlChar* const)"Re";
  static xmlChar* const timeLimit           = (xmlChar* const)"timeLimit";
  static xmlChar* const iterationLimit      = (xmlChar* const)"iterationLimit";
  static xmlChar* const plotInterval        = (xmlChar* const)"plotInterval";
  static xmlChar* const tau                 = (xmlChar* const)"tau";
  static xmlChar* const gamma               = (xmlChar* const)"gamma";
  static xmlChar* const dim                 = (xmlChar* const)"dim";
  static xmlChar* const width               = (xmlChar* const)"width";
  static xmlChar* const size                = (xmlChar* const)"size";
  static xmlChar* const filename            = (xmlChar* const)"filename";
  static xmlChar* const parallelizationSize =
    (xmlChar* const)"parallelizationSize";
  static xmlChar* const environment = (xmlChar* const)"environment";

  xmlAttrPtr attr = node->properties;

  while (attr) {
    if (xmlStrEqual(attr->name, re)) {
      configuration->re = parseNumber<Scalar>(attr->children->content);
    } else if (xmlStrEqual(attr->name, timeLimit)) {
      configuration->timeLimit = parseNumber<Scalar>(attr->children->content);
    } else if (xmlStrEqual(attr->name, iterationLimit)) {
      configuration->iterationLimit = parseNumber<int>(attr->children->content);
    } else if (xmlStrEqual(attr->name, plotInterval)) {
      configuration->plotInterval = parseNumber<Scalar>(
        attr->children->content);
    } else if (xmlStrEqual(attr->name, tau)) {
      configuration->tau = parseNumber<Scalar>(attr->children->content);
    } else if (xmlStrEqual(attr->name, gamma)) {
      configuration->gamma = parseNumber<Scalar>(attr->children->content);
    } else if (xmlStrEqual(attr->name, dim)) {
      configuration->dim = parseNumber<int>(attr->children->content);
    } else if (xmlStrEqual(attr->name, width)) {
      configuration->width = parseVector<Scalar>(attr->children->content);
    } else if (xmlStrEqual(attr->name, size)) {
      configuration->size = parseVector<int>(attr->children->content);
    } else if (xmlStrEqual(attr->name, parallelizationSize)) {
      configuration->parallelizationSize = parseVector<int>(
        attr->children->content);
    } else if (xmlStrEqual(attr->name, environment)) {
      configuration->environment = parseVector<Scalar>(attr->children->content);
    } else if (xmlStrEqual(attr->name, filename)) {
      configuration->filename = (const char*)attr->children->content;
    }
    attr = attr->next;
  }
}

std::unique_ptr<Configuration>
XmlConfigurationParser::
parse(boost::filesystem::path const& filePath) {
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

  if (!xmlStrEqual(root->name, (xmlChar* const)"scenario")) {
    throwException("Failed to find root element");
  }

  std::unique_ptr<Configuration> configuration(new Configuration());

  parseScenarioParameters(root, configuration.get());

  parseScenarioChildren(root->children, configuration.get());

  xmlFreeDoc(doc);

  xmlCleanupParser();

  return configuration;
}
