#include <vector>
#include <array>
#include <iostream>
#include <stdexcept>

class nearest_neighbors_bonds {
private:
    
    int Lx                    ;
    int Ly                    ; 
    int N                     ;
    int coordination_number   ; 
    
    std::string geometry      ;
    
    std::vector<int> vertex1  ;
    std::vector<int> vertex2  ;
    std::vector<int> nn_bonds ;

    void build_triangular_lattice() {
        
        coordination_number = 6                  ; 
        vertex1.reserve(3 * N)                   ;
        vertex2.reserve(3 * N)                   ;
        nn_bonds.resize(N * coordination_number) ;
        
        int x     ;
        int y     ;
        int index ;

        for (int index = x = 0; x < Lx; ++x) {
            for (int y = 0; y < Ly; ++y) {
                
                int v1 = x + Lx * y;
                int v2 =  ((((x+1)+Lx)%Lx)+(((y)+Ly)%Ly)*Lx);
                vertex1[index]   = v1                       ;
                vertex2[index++] = v2                       ;
                nn_bonds[v1 * coordination_number + 0] = v2 ; 
                
                
                v2 = ((((x)+Lx)%Lx)+(((y+1)+Ly)%Ly)*Lx)     ;
                vertex1[index]   = v1                       ;
                vertex2[index++] = v2                       ;
                nn_bonds[v1 * coordination_number + 1] = v2 ;
                
                v2 = ((((x+1)+Lx)%Lx)+(((y+1)+Ly)%Ly)*Lx)   ;
                vertex1[index]   = v1                       ;
                vertex2[index++] = v2                       ;
                nn_bonds[v1*coordination_number + 2] = v2   ; 

                nn_bonds[v1*coordination_number + 3] = ((((x)+Lx)%Lx)+(((y-1)+Ly)%Ly)*Lx) ;  
                nn_bonds[v1*coordination_number + 4] = ((((x-1)+Lx)%Lx)+(((y)+Ly)%Ly)*Lx)    ;
                nn_bonds[v1*coordination_number + 5] = ((((x-1)+Lx)%Lx)+(((y-1)+Ly)%Ly)*Lx)  ;
                
            }
        }
    }     


 
public:

    nearest_neighbors_bonds(int Lx = 0, int Ly = 0, std::string geometry = "")
        : Lx(Lx), Ly(Ly), N(Lx * Ly), coordination_number(0), geometry(geometry) {
        
        if (!geometry.empty()) {
            build(Lx, Ly, geometry);
        }
    }

    void build(int Lx, int Ly, std::string geometry ) {
        
        this->Lx       = Lx         ;
        this->Ly       = Ly         ;
        this->N        = Lx * Ly    ;
        this->geometry = geometry   ;
        

        if (geometry == "square") {
            std::cout << "to implement" ;
            // build_square_lattice() ;
        } else if (geometry == "triangular") {
            build_triangular_lattice() ; 

        } else {
            throw std::invalid_argument("Unsupported geometry: " + geometry);
        }
    }
        // const std::vector<int>& get_vertex1() const { return vertex1; } 
        // const std::vector<int>& get_vertex2() const { return vertex2; } 
        // const std::vector<int>& get_nn_bonds() const { return nn_bonds; }
        // int get_coordination_number() const { return coordination_number; }
        // std::string get_geometry() const { return geometry; }
};