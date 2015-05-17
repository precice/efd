#include "XmlConfigurationParser.hpp"

#include <Uni/ExecutionControl/exception>
#include <Uni/Firewall/Interface>
#include <Uni/Helpers/parsefromstring>
#include <Uni/Logging/macros>

#include <Eigen/Core>

#include <libxml/parser.h>

#include <boost/filesystem/fstream.hpp>
#include <boost/regex.hpp>

using Structure::XmlConfigurationParser;

namespace Structure {
//
// static inline bool
// is_default(xmlChar const* content) {
// std::string content_string
// = reinterpret_cast<char const*>(content);

// static boost::regex default_regex("default", boost::regex::icase);

// if (boost::regex_search(content_string, default_regex)) {
// return true;
// }

// return false;
// }

// static inline bool
// parse_bool(xmlChar const* content) {
// static boost::regex on_regex("on|true|1", boost::regex::icase);
// static boost::regex off_regex("off|false|0", boost::regex::icase);

// std::string content_string
// = reinterpret_cast<char const*>(content);

// if (boost::regex_search(content_string, on_regex)) {
// return true;
// } else if (boost::regex_search(content_string, off_regex)) {
// return false;
// } else {
// throwException("Failed to parse boolean value");

// return false;
// }
// }

// template <typename T>
// static inline void
// parse_in(xmlChar const* content, T& number) {
// if (!Uni::Helpers::parse_integer_number(
// std::string(reinterpret_cast<char const*>(content)),
// number)) {
// throwException("Failed to parse integer number '{1}'", content);
// }
// }

// static inline int
// parse_in(xmlChar const* content) {
// int result;
// parse_in(content, result);

// return result;
// }

// template <typename T>
// static inline void
// parse_fpn(xmlChar const* content, T& number) {
// if (!Uni::Helpers::parse_floating_point_number(
// std::string(reinterpret_cast<char const*>(content)),
// number)) {
// throwException("Failed to parse floating-point number "
// "'{1}'", content);
// }
// }

// static inline long double
// parse_fpn(xmlChar const* content) {
// long double result;
// parse_fpn(content, result);

// return result;
// }

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

    xmlFreeDoc(doc);

    xmlCleanupParser();
  }

  void
  parseScenarioParameters(xmlNodePtr node) const {
    static xmlChar* const       environment = (xmlChar* const)"environment";
    static xmlChar const* const precice_configuration_path
      = (xmlChar const*)"precice-configuration-path";

    xmlAttrPtr attr = node->properties;

    int dimensions = -1;

    while (attr) {
      if (xmlStrcasecmp(attr->name, environment) == 0) {
        Eigen::Matrix<long double, 3, 1> environment_force;
        dimensions = parse_fpv(attr->children->content, environment_force);
        configuration->set("/EnvironmentForce", environment_force);
      } else if (xmlStrcasecmp(attr->name, precice_configuration_path) == 0) {
        boost::filesystem::path path(
          reinterpret_cast<char const*>(attr->children->content));

        if (path.is_relative()) {
          path = filePath.parent_path() / path;
        }

        configuration->set("/PreciceConfigurationPath", path);
      }
      attr = attr->next;
    }

    if (dimensions < -1) {
      throwException("Unknown dimension size");
    }

    configuration->set("/Dimensions", static_cast<unsigned>(dimensions));
  }

  Configuration*          configuration;
  boost::filesystem::path filePath;

  Uni_Firewall_INTERFACE_LINK(XmlConfigurationParser);
};
}

XmlConfigurationParser::
XmlConfigurationParser(
  std::unique_ptr<Configuration> const& configuration,
  boost::filesystem::path const&        file_path)
  : _im(new XmlConfigurationParserImplementation(this)) {
  _im->configuration = configuration.get();
  _im->filePath      = file_path;
  _im->parse();
}

XmlConfigurationParser::
~XmlConfigurationParser() {}
