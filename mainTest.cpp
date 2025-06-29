
#include "RandomGeneratorClass.hpp"
#include <iostream>


// // Thread-local random engine for thread safety
// thread_local std::mt19937_64 engine(std::random_device{}());

// // 1. Fast random integer in [min, max]
// int RandIntOnInterval(int min, int max) noexcept {
//     // Fast modulo method (better than rand() % n)
//     return min + static_cast<int>(engine() % static_cast<unsigned>(max - min + 1));
    
//     /* Alternative with perfect uniformity (slightly slower):
//     static std::uniform_int_distribution<int> dist;
//     return dist(engine, std::uniform_int_distribution<int>::param_type(min, max));
//     */
// }

// // 2. Fast random float in [0, 1)
// float RandFloat() noexcept {
//     // Bit manipulation technique (extremely fast)
//     constexpr float factor = 1.0f / (1ULL << 24);
//     return (engine() >> 40) * factor;
    
//     /* Alternative with perfect uniformity:
//     static std::uniform_real_distribution<float> dist(0.0f, 1.0f);
//     return dist(engine);
//     */
// }


int main() {

    auto & rng = RandomGenerator::Get();
    //std::cout << rng << " " ; 

    // Generate random integers
    std::cout << "Random ints [10, 20]: ";
    for (int i = 0; i < 5; ++i) {
        std::cout << rng.Int(10, 20) << " ";
    }
    
    // Generate random floats
    std::cout << "\nRandom floats [0,1): ";
    for (int i = 0; i < 5; ++i) {
        std::cout << rng.Float() << " ";
    }
    
    // Generate random signs
    std::cout << "\nRandom signs: ";
    for (int i = 0; i < 5; ++i) {
        std::cout << rng.Sign() << " ";
    }

    // Generate random signs
    // std::cout << "\nRandom seed: ";
    // for (int i = 0; i < 5; ++i) {
    //     std::cout << rng.Seed(2) << " ";
    // }
    
    return 0;

}