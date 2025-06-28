//  bond percolation on a square lattice of size WIDTH x WIDTH = N.  2N bonds
#include <iostream>
#include <iomanip>
#include <string>
#include <fstream>
#include <yaml-cpp/yaml.h>


int Lx               ;
int Ly               ;
int N                ; 
std::string geometry ;
int nruns            ;

class ReadInput {
public:

    // Constructor that loads and parses the YAML file
    explicit ReadInput(const std::string& filepath) {
        try {
            YAML::Node config = YAML::LoadFile(filepath);
            
            // Parse the YAML content
            SimulationValues.name     = config["name"].as<std::string>();
            SimulationValues.Lx       = config["lattice"]["Lx"].as<int>();
            SimulationValues.Ly       = config["lattice"]["Ly"].as<int>();
            SimulationValues.geometry = config["lattice"]["geometry"].as<std::string>();
            SimulationValues.nruns    = config["simulation"]["number_of_steps"].as<long>();
            
        } catch (const YAML::Exception& e) {
            
          throw std::runtime_error("YAML Error: " + std::string(e.what()));
        
        } catch (const std::exception& e) {
            
          throw std::runtime_error("Configuration Error: " + std::string(e.what()));
        
        }
    }

    // Structure for lattice dimensions
    struct SimulationParameters {
        std::string name 	 ;
		int Lx			 	 ;
        int Ly				 ;
        std::string geometry ;
        int nruns		     ; 
        
        void print() const {
            std::cout 
				<< std::setw(20) << "Simulation Name : " << name << "\n"
				<< std::setw(23) << "Lattice Parameters :\n"
				<< std::setw(18) << "Lx         : " << Lx << "\n"
				<< std::setw(18) << "Ly         : " << Ly << "\n"
                << std::setw(18) << "Geometry   : " << geometry << "\n" 
				<< std::setw(10) << "Runs :\n"  
				<< std::setw(22) << "Number of runs : " << nruns << "\n\n" ; 
        }
    };

    // Getter methods
	const SimulationParameters& getSimulationParameters() const { return SimulationValues; }

    // Print configuration
    void printInputParameters() const {
        std::cout << "Simulation Parameters:\n" ; 
		SimulationValues.print() ;
    }

private:
    //std::string simulation_name;
    SimulationParameters SimulationValues;
};

void ReadYamlInputFile(const YAML::Node& node, const std::string& prefix = "") {
    try {
        if (node.IsMap()) {
            for (const auto& pair : node) {
                std::string key = pair.first.as<std::string>();
                std::cout << prefix << "Key: " << key;
                
                const YAML::Node& value = pair.second;
                if (value.IsScalar()) {
                    std::cout << ", Value: " << value.as<std::string>() << std::endl;
                } else if (value.IsMap() || value.IsSequence()) {
                    std::cout << " (Nested)" << std::endl;
                    ReadYamlInputFile(value, prefix + "  "); // Recursively print nested structures
                }
            }
        } else if (node.IsSequence()) {
            for (size_t i = 0; i < node.size(); ++i) {
                std::cout << prefix << "Item " << i << ":" << std::endl;
                ReadYamlInputFile(node[i], prefix + "  "); // Recursively print sequence items
            }
        } else if (node.IsScalar()) {
            std::cout << prefix << "Value:   " << node.as<std::string>() << std::endl;
        }
    } catch (const YAML::Exception& e) {
            
        throw std::runtime_error("YAML Error: " + std::string(e.what()));
    
    } catch (const std::exception& e) {
        
        throw std::runtime_error("Configuration Error: " + std::string(e.what()));
    
    }
}



int main()
{

    YAML::Node config = YAML::LoadFile("input.yaml");
    ReadYamlInputFile(config);
    //std::cout << config << std::endl;
	//ReadInput input("input.yaml") ;
	//input.printInputParameters()  ;

    // Lx       = input.getSimulationParameters().Lx       ;
	// Ly       = input.getSimulationParameters().Ly       ;
	// N        = Lx * Ly                        			;
	// geometry = input.getSimulationParameters().geometry ;
	// nruns    = input.getSimulationParameters().nruns   	;

	//std::cout << input.getSimulationParameters() << std::endl ; 				

 	return 0;
}
  
