#include "LatticeBonds.hpp"
#include "LatticeCoordinationNumber.hpp"

#include <iostream>
#include <map>
#include <string>


/**
 * @brief Main program to demonstrate LatticeBonds functionality
 * Creates a 12x12 lattice with 6 bonds per site and prints first bond vertices
 */
int main() {

    const std::string geometry = "triangular"; 
    const int  printfreq        = 50000000    ;
    const int  printfreq2       = 50000000    ;
    const int  runsmax          = 50000001    ;
    const int  Lx               = 12          ;  // Lattice width
    const int  Ly               = 12          ;  // Lattice height
    const int  N                = Lx*Ly       ;  
    int nCoord                                ;

    try {
        nCoord = coordinationNumber[geometry] ;
    } catch (const std::out_of_range& e) {
        std::cerr << "Invalid Lattice Geometry: " << e.what() << std::endl;
        return 1;
    }

    const int totalBonds       = nCoord/2*N ; 
    long order[nCoord/2 * N]                ; 
    for (int i = 0; i < N; ++i)  order[i] = i;



    // Create and initialize lattice
    LatticeBonds lattice(nCoord, Lx, Ly);
    lattice.setNNBondsList();

    
    
    
    // Query and display bond information
    try {
        auto [v1, v2] = lattice.getBondVertices(2);
        std::cout << "Bond 0 connects vertices: " << v1 << " and " << v2 << std::endl;
    } catch (const std::out_of_range& e) {
        std::cerr << "Error: " << e.what() << std::endl;
        return 1;
    }

    return 0;
}