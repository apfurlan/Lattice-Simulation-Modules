#ifndef READ_YAML_INPUT_HPP
#define READ_YAML_INPUT_HPP

#include <map>
#include <string>
#include <variant>
#include <vector>
#include <stdexcept>

namespace YAML {
    class Node;
}

class ReadYamlInput {
public:
    using SequenceValue = std::variant<int, double, bool, std::string>;
    using Sequence = std::vector<SequenceValue>;
    using ConfigValue = std::variant<int, double, bool, std::string, Sequence>;
    using ConfigMap = std::map<std::string, ConfigValue>;

    // Existing methods
    ConfigMap readConfig(const std::string& filename);
    
    // New method to specifically get property names
    std::vector<std::string> getPropertyNames(const std::string& filename);
    
private:
    // Helper methods
    ConfigMap readNestedConfig(const YAML::Node& node);
    void handleScalarValue(ConfigMap& config, const std::string& key, const YAML::Node& value);
    void handleSequenceValue(ConfigMap& config, const std::string& key, const YAML::Node& value);
    
    // New helper method for property extraction
    std::vector<std::string> extractPropertyNames(const ConfigMap& config);
};

#endif