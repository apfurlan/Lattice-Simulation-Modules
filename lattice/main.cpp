#include "LatticeBondsClass.hpp"
#include "LatticeCoordinationNumber.hpp"

#include <iostream>
#include <map>
#include <string>
#include <random>
#include <chrono>
#include <iterator>
#include <algorithm>

/**
 * @brief Main program to demonstrate LatticeBonds functionality
 * Creates a 12x12 lattice with 6 bonds per site and prints first bond vertices
 */
int main() {

    const std::string geometry = "triangular" ; 
    const int  printfreq        = 50000000    ;
    const int  printfreq2       = 50000000    ;
    const int  runsmax          = 1    ;
    const int  Lx               = 12          ;  // Lattice width
    const int  Ly               = 12          ;  // Lattice height
    const int  N                = Lx*Ly       ;  
    int nCoord                                ;

    static long prt[N] ;
    static long valence[N] ; 
    static long imax[N] ;

    double M1, M2 ; 
    long  nsites, nbonds, big ; // , i, index = 3 * N ; // Index for bonds
    // int index = 3 * N ; // Index for bonds
    int i ;

    try {
        nCoord = coordinationNumber.at(geometry) ;
    } catch (const std::out_of_range& e) {
        std::cerr << "Invalid Lattice Geometry: " << e.what() << std::endl;
        return 1;
    }

    LatticeBonds lattice(nCoord, Lx, Ly);
    lattice.setNNBondsList();

    const int totalBonds = nCoord/2*N         ;  
    long* order = new long[totalBonds]        ; 
    for (int i = 0; i < N; ++i)  order[i] = i ;

    unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();
    std::cout << seed << std::endl;
    std::mt19937_64 engine(seed); 
    
    std::uniform_int_distribution<int> dist6(minval,maxval);
    for(int i=0; i<N; i++)imax[i]= dist6(engine); 

    for (long run=0; run <= runsmax; run++) {

        std::shuffle(order, order + totalBonds,  engine);
        for (int i = 0; i < N; i++){
            prt[i]     = -1 ; 
            valence[i] = 0  ; 
        }

        
        double M1 = N ;  
        double M2 = N ; 
        for (nsites=nbonds=big=i=0; i<totalBonds; i++){

            auto [v1, v2] = lattice.getBondVertices(i);
            ranDir = 2*round(genrand_real2())-1 ;
            
        }
    }
    
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