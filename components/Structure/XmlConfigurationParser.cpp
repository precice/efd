#include "XmlConfigurationParser.hpp"

#include <Uni/ExecutionControl/exception>
#include <Uni/Helpers/parsefromstring>
#include <Uni/Logging/macros>

#include <Eigen/Core>

#include <libxml/parser.h>

#include <boost/filesystem/fstream.hpp>
#include <boost/regex.hpp>

using Structure::XmlConfigurationParser;
using Structure::Simulation;

namespace Structure {
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

static void
parseScenarioParameters(xmlNodePtr  node,
                        Simulation* configuration) {
  static xmlChar* const type     = (xmlChar* const)"type";
  static xmlChar* const velocity = (xmlChar* const)"velocity";

  xmlAttrPtr attr = node->properties;

  int dim = -1;

  while (attr) {
    if (xmlStrcasecmp(attr->name, type) == 0) {
      parse_in(attr->children->content, configuration->type());
    } else if (xmlStrcasecmp(attr->name, velocity) == 0) {
      int temp_dim = parse_fpv(attr->children->content,
                               configuration->velocity());

      if (dim == -1) {
        dim = temp_dim;
      } else if (temp_dim != dim) {
        throwException("Inconsistent dimension size in velocity vector,"
                       "'{1}', because there was already found "
                       "dimension size of '{2}'",
                       temp_dim, dim);
      }
    }
    attr = attr->next;
  }

  if (dim == -1) {
    throwException("Unknown dimension size");
  }
}
}
}

void
XmlConfigurationParser::
parse(std::unique_ptr<Simulation> const& configuration,
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

  xmlFreeDoc(doc);

  xmlCleanupParser();
}
