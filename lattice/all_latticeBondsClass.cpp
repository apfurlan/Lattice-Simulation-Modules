#include <utility>
#include <iostream>
#include <random>

class RandomNumber

/**
 * @class LatticeBonds
 * @brief Simulates bond percolation on a 2D lattice with 
 * periodic boundary conditions.
 *
 * This class constructs a lattice of size `Lx × Ly`, where each 
 * site has `ncoord` neighbors. It manages bond states (open/closed) 
 * and provides methods to access nearest-neighbor (NN) bonds. Used 
 * for studying percolation theory, connectivity, and related statistical mechanics problems.
 */
class LatticeBonds{
private : 
    int nCoord    ;
    int Lx        ;  ///< Lattice size along the x-axis.
    int Ly        ;  ///< Lattice size along the y-axis.
    int N         ;  ///< Total sites in the lattice (`N = Lx * Ly`).
    long* vertex1 ;  ///< Array storing vertex 1 of each bond.
    long* vertex2 ;  ///< Array storing vertex 2 of each bond.
    long* nn      ;  ///< Nearest-neighbor table (`nn[site_index * nnBonds]` gives neighbors).
    int* lattice  ;  ///< Bond states (0 = closed, 1 = open).

public : 
    /**
     * @brief Constructor: Initializes the lattice and allocates memory for bonds.
     * @param nCoord Number of bonds per site (must be even).
     * @param Lx     Lattice width (sites along x-axis).
     * @param Ly     Lattice height (sites along y-axis).
     */
     LatticeBonds(int nCoord, int Lx, int Ly) 
    : nCoord(nCoord), Lx(Lx), Ly(Ly), N(Lx * Ly) {
       
        // Allocate memory for arrays
        vertex1 = new long[nCoord/2 * N];
        vertex2 = new long[nCoord/2 * N];
        nn      = new long[nCoord * N];
        lattice = new int[2 * N];  // 2 bonds per site (adjust if needed)
        
    }
    
    /// @brief Destructor: Frees allocated memory.
    ~LatticeBonds(){
        delete[] vertex1;
        delete[] vertex2;
        delete[] nn;
        delete[] lattice;
    }
    /**
     * @brief Generates the nearest-neighbor bond list.
     * Populates `vertex1`, `vertex2`, and `nn` arrays with 
     * periodic boundary conditions.
     */
    void setNNBondsList(){
        
      int i,j ;
      for (int index = i = 0; i < Lx; ++i){
        for (int j = 0; j < Ly; ++j) {

            int v1 = i + Lx*j;
        
            int v2 =  ((((i+1)+Lx)%Lx)+(((j)+Ly)%Ly)*Lx) ; 
                
            vertex1[index]   = v1 ;
            vertex2[index++] = v2 ;
        
            nn[v1*nCoord + 0]=v2 ;
        
            v2 = ((((i)+Lx)%Lx)+(((j+1)+Ly)%Ly)*Lx) ;

            vertex1[index]   = v1 ;
            vertex2[index++] = v2 ;
        
            nn[v1*nCoord + 2] = v2 ;  
            nn[v1*nCoord + 3] = ((((i)+Lx)%Lx)+(((j-1)+Ly)%Ly)*Lx) ;
            //address(x,y-1);
            nn[v1*nCoord + 4] = ((((i-1)+Lx)%Lx)+(((j)+Ly)%Ly)*Lx) ;
            //address(x-1,y);
            nn[v1*nCoord + 5] = ((((i-1)+Lx)%Lx)+(((j-1)+Ly)%Ly)*Lx) ;
            //address(x-1,y-1);
        }    
    }

    };
    //void setBondState(int bondIndex,int state){lattice[bondIndex] = state ;}
    //setBondIndex

    //int getBondState(int bondIndex) { return lattice[bondIndex]};
    //set

    std::pair<long, long> getBondVertices(int bondIndex) const {
        if (bondIndex < 0 || bondIndex >= (nnBonds / 2) * N) {
            throw std::out_of_range("Invalid bondIndex");
        }
        return {vertex1[bondIndex], vertex2[bondIndex]};
    };
};


int main(int argc, char * argv[])
{
  
    const int nnBonds = 6;  // 6 neighbors (e.g., hexagonal)
    const int Lx = 12;      // Lattice width
    const int Ly = 12;      // Lattice height
    const int N  = Lx*Ly ; 

    static long* order = new long[nnBonds/2 * N];


    LatticeBonds Lattice(nnBonds,Lx,Ly) ; 
    Lattice.setNNBondsList();

    // Get vertices for bondIndex = 0
    auto [v1, v2] = Lattice.getBondVertices(12);     
    std::cout << "Bond 0 connects vertices: " << v1 << " and " << v2 << std::endl;

    const unsigned int seed = time(0) ;
    // mt19937_64 rng(seed);
    // uniform_int_distribution<int> dist6(minval,maxval);
    std::cout << seed << "\n";
    //for (i = 0; i < index; ++i)  order[i] = i;


    return 0;
    

};

