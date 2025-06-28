#include "LatticeCoordinationNumber.hpp"

// Definition of the coordination number map
const std::map<std::string, int> coordinationNumber = {
    {"hexagonal" , 3} ,
    {"square"    , 4} ,
    {"triangular", 6} ,
    {"cubic"     , 6} , 
    {"FCC"       , 8} 
};