#ifndef LATTICE_BONDS_CLASS_HPP
#define LATTICE_BONDS_CLASS_HPP

#include <utility>
#include <stdexcept>

/**
 * @class LatticeBonds
 * @brief Simulates bond percolation on a 2D lattice with periodic boundary conditions
 * 
 * This class manages a lattice where each site has nnBonds connections (bonds).
 * It provides functionality to initialize the lattice, set bond states,
 * and query bond information.
 */
class LatticeBonds {
private:
    int nCoord    ; ///< Number of bonds per lattice site
    int Lx        ; ///< Lattice size in x-direction
    int Ly        ; ///< Lattice size in y-direction
    int N         ; ///< Total number of sites (Lx * Ly)
    long* vertex1 ; ///< Array storing first vertex of each bond
    long* vertex2 ; ///< Array storing second vertex of each bond
    long* nn      ; ///< Nearest-neighbor table (nn[site_index * nnBonds])
    int* lattice  ; ///< Bond states (0 = closed, 1 = open)

public:
    /**
     * @brief Construct a new LatticeBonds object
     * @param nCoord Number of nearest neighbor bond of a given bond (must be even)
     * @param Lx Lattice width in sites
     * @param Ly Lattice height in sites
     */
    LatticeBonds(int nCoord, int Lx, int Ly);
    
    /**
     * @brief Destroy the LatticeBonds object
     * Frees all allocated memory
     */
    ~LatticeBonds();

    /**
     * @brief Initialize the nearest-neighbor bond list
     * Populates vertex1, vertex2 and nn arrays with periodic boundary conditions
     */
    void setNNBondsList();

    /**
     * @brief Get the vertices connected by a bond
     * @param bondIndex Index of the bond (0 ≤ bondIndex < total bonds)
     * @return std::pair<long, long> Vertices connected by the bond
     * @throws std::out_of_range if bondIndex is invalid
     */
    std::pair<long, long> getBondVertices(int bondIndex) const;
};

#endif // LATTICE_BONDS_HPP