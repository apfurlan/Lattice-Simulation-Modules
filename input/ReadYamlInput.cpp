//  bond percolation on a square lattice of size WIDTH x WIDTH = N.  2N bonds
#include "ReadYamlInput.hpp"
#include <yaml-cpp/yaml.h>
#include <map>
#include <string>
#include <variant>
#include <iostream>

ReadYamlInput::ConfigMap ReadYamlInput::readConfig(const std::string& filename) {
    ConfigMap config;
    YAML::Node yaml = YAML::LoadFile(filename);
    
    for (const auto& entry : yaml) {
        std::string key = entry.first.as<std::string>();
        const YAML::Node& value = entry.second;
        
        if (value.IsScalar()) {
            handleScalarValue(config, key, value);
        } else if (value.IsMap()) {
            auto nested = readNestedConfig(value);
            for (auto& [nested_key, nested_value] : nested) {
                config[key + "." + nested_key] = nested_value;
            }
        } else if (value.IsSequence()) {
            handleSequenceValue(config, key, value);
        }
    }
    return config;
}

void ReadYamlInput::handleScalarValue(ConfigMap& config, 
                                    const std::string& key, 
                                    const YAML::Node& value) {
    if (value.Tag() == "!") {
        config[key] = value.as<std::string>();
    } else {
        try { config[key] = value.as<int>(); return; } catch(...) {}
        try { config[key] = value.as<double>(); return; } catch(...) {}
        try { config[key] = value.as<bool>(); return; } catch(...) {}
        config[key] = value.as<std::string>();
    }
}

void ReadYamlInput::handleSequenceValue(ConfigMap& config,
                                      const std::string& key,
                                      const YAML::Node& value) {
    std::vector<std::variant<int, double, bool, std::string>> sequence;
    for (const auto& item : value) {
        if (item.IsScalar()) {
            try { sequence.push_back(item.as<int>()); continue; } catch(...) {}
            try { sequence.push_back(item.as<double>()); continue; } catch(...) {}
            try { sequence.push_back(item.as<bool>()); continue; } catch(...) {}
            sequence.push_back(item.as<std::string>());
        }
        // Handle nested sequences or maps if needed
    }
    config[key] = sequence;
}

ReadYamlInput::ConfigMap ReadYamlInput::readNestedConfig(const YAML::Node& node) {
    ConfigMap nested;
    for (const auto& entry : node) {
        std::string key = entry.first.as<std::string>();
        const YAML::Node& value = entry.second;
        
        if (value.IsScalar()) {
            handleScalarValue(nested, key, value);
        } else if (value.IsSequence()) {
            handleSequenceValue(nested, key, value);
        }
    }
    return nested;
}

#include "ReadYamlInput.hpp"
#include <yaml-cpp/yaml.h>

// [Keep all your existing method implementations here...]

std::vector<std::string> ReadYamlInput::getPropertyNames(const std::string& filename) {
    ConfigMap config = readConfig(filename);
    return extractPropertyNames(config);
}

std::vector<std::string> ReadYamlInput::extractPropertyNames(const ConfigMap& config) {
    std::vector<std::string> properties;
    
    try {
        // Find the properties key
        auto it = config.find("properties");
        if (it == config.end()) {
            throw std::runtime_error("'properties' key not found in YAML");
        }
        
        // Verify it's a sequence
        const Sequence* seq = std::get_if<Sequence>(&it->second);
        if (!seq) {
            throw std::runtime_error("'properties' is not a sequence");
        }
        
        // Reserve space for efficiency
        properties.reserve(seq->size());
        
        // Extract all string values
        for (const auto& item : *seq) {
            if (const auto* str = std::get_if<std::string>(&item)) {
                properties.push_back(*str);
            } else {
                throw std::runtime_error("Non-string value found in properties list");
            }
        }
    } catch (const std::bad_variant_access&) {
        throw std::runtime_error("Type error while reading properties");
    }
    
    return properties;
}