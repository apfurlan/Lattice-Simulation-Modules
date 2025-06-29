#include "lattice/LatticeBondsClass.hpp"
#include "lattice/LatticeCoordinationNumber.hpp"
#include "random_number/RandomGeneratorClass.hpp"

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
    const int  runsmax          = 1           ;
    const int  Lx               = 4           ;  // Lattice width
    const int  Ly               = 4           ;  // Lattice height
    const int  N                = Lx*Ly       ;  
    int nCoord                                ;

    const int maxValIn =  2 ;   
    const int maxValOut = 2 ;   
    
    static long prt[N]      ;
    static long ptrPlus[N]  ; 
    static long ptrMinus[N] ; 
    
    static long valence[N] ; 
    static long imax[N]    ;
    static int valIn[N]    ;
    static int valOut[N]   ;
    
    int sumbondsIn, sumbondsOut ;

    double M1, M2 ; 
    long  nsites, nbonds, big ; 

    int i ;
    int iBondMinus, iBondPlus ; 

    

    try {
        nCoord = coordinationNumber.at(geometry) ;
    } catch (const std::out_of_range& e) {
        std::cerr << "Invalid Lattice Geometry: " << e.what() << std::endl;
        return 1;
    }

    LatticeBonds lattice(nCoord, Lx, Ly);
    lattice.setNNBondsList();

    const int totalBonds = nCoord/2*N                  ;  
    long * order = new long[totalBonds]                ; 
    for (int i = 0; i < totalBonds; ++i) order[i] = i  ;


    auto& rng    = RandomGenerator::Get();
    auto& engine = rng.GetEngine();

    for (long run=0; run <= runsmax; run++) {
        
        std::shuffle(order, order + totalBonds, engine);
        
        for (int i = 0; i < N; i++){
            prt[i]     = -1 ; 
            valence[i] =  0 ; 
        }

        
        double M1 = N ;  
        double M2 = N ; 
        for (nsites=nbonds=big=i=0; i<totalBonds; i++){

            auto [v1, v2] = lattice.getBondVertices(i);
            
            int ranDir    = rng.Sign() ;
            bool addBond  = false ; 
            int dirBond   = 0 ;
            for (int nAttempts = 0; nAttempts < 2; nAttempts++){
                
                int dirBond = ranDir - nAttempts*(2*ranDir) ;  
                if (dirBond > 0){
                    int sumbondsIn  = valIn[v2]  ; 
                    int sumbondsOut = valOut[v1] ;
                } else {
                    int sumbondsIn  = valIn[v1]  ; 
                    int sumbondsOut = valOut[v2] ;
                }

                if (sumbondsIn < maxValIn && sumbondsOut < maxValOut) {
                    addBond = true ;
                    break ; 
                }

            }

            if (addBond){
                if(dirBond > 0){
                    valIn[v2]++    ; 
                    valOut[v1]++   ;
                    iBondPlus++    ; 
                } else {
                    valIn[v1]++    ; 
                    valOut[v2]++   ;
                    iBondMinus++   ; 
                }
            }

            r1 = findroot(v1);
            r2 = findroot(v2);

        }
    }
    
    // Query and display bond information
    // try {
    //     auto [v1, v2] = lattice.getBondVertices(1);
    //     std::cout << "Bond 2 connects vertices: " << v1 << " and " << v2 << std::endl;
    // } catch (const std::out_of_range& e) {
    //     std::cerr << "Error: " << e.what() << std::endl;
    //     return 1;
    // }

    return 0;
}