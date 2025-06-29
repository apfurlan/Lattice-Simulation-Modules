#include "LatticeBondsClass.hpp"
#include <iostream>

LatticeBonds::LatticeBonds(int nCoord, int Lx, int Ly) 
    : nCoord(nCoord), Lx(Lx), Ly(Ly), N(Lx * Ly) {
    vertex1 = new long[nCoord/2 * N];
    vertex2 = new long[nCoord/2 * N];
    nn = new long[nCoord * N];
    lattice = new int[2 * N];
}

LatticeBonds::~LatticeBonds() {
    delete[] vertex1;
    delete[] vertex2;
    delete[] nn;
    delete[] lattice;
}

void LatticeBonds::setNNBondsList() {
    int index = 0;
    for (int i = 0; i < Lx; ++i) {
        for (int j = 0; j < Ly; ++j) {
            int v1 = i + Lx * j;
            
            // Right neighbor (periodic boundary)
            int v2 = (((i + 1) % Lx) + (j % Ly) * Lx);
            vertex1[index] = v1;
            vertex2[index++] = v2;
            nn[v1 * nCoord + 0] = v2;

            // Top neighbor (periodic boundary)
            v2 = (i % Lx) + ((j + 1) % Ly) * Lx;
            vertex1[index] = v1;
            vertex2[index++] = v2;
            nn[v1 * nCoord + 2] = v2;

            // Other neighbors
            nn[v1 * nCoord + 3] = (i % Lx) + ((j - 1 + Ly) % Ly) * Lx;  // Bottom
            nn[v1 * nCoord + 4] = ((i - 1 + Lx) % Lx) + (j % Ly) * Lx;  // Left
            nn[v1 * nCoord + 5] = ((i - 1 + Lx) % Lx) + ((j - 1 + Ly) % Ly) * Lx; // Bottom-left
        }
    }
}

std::pair<long, long> LatticeBonds::getBondVertices(int bondIndex) const {
    if (bondIndex < 0 || bondIndex >= (nCoord / 2) * N) {
        throw std::out_of_range("Invalid bondIndex");
    }
    return {vertex1[bondIndex], vertex2[bondIndex]};
}