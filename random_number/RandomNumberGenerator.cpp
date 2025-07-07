#include "RandomNumberGenerator.hpp"
#include <chrono>
#include <iostream>


// Initialize with hardware entropy and time-based fallback
thread_local std::mt19937_64 RandomNumberGenerator::engine = []() {
    std::random_device rd;
    unsigned long seed = rd.entropy() ? rd() : 
        std::chrono::system_clock::now().time_since_epoch().count();
    return std::mt19937_64(seed);
}();

// void RandomGenerator::Seed(unsigned long seed) noexcept {
//     engine.seed(seed);
// }

// Branchless integer generation with optimized modulo
int RandomNumberGenerator::Int(int min, int max) noexcept {
    const unsigned range = static_cast<unsigned>(max - min + 1);
    // Fast modulo reduction using Lemire's method
    return min + static_cast<int>((engine() & 0xFFFFFFFF) * static_cast<uint64_t>(range) >> 32);
    
    /* Alternative: For better uniformity at slight performance cost:
    static thread_local std::uniform_int_distribution<int> dist;
    return dist(engine, std::uniform_int_distribution<int>::param_type(min, max));
    */
}

// Fastest float generation using bit manipulation
float RandomNumberGenerator::Float() noexcept {
    // Uses 23 bits of randomness (IEEE 754 mantissa precision)
    static constexpr float norm = 1.0f / (1UL << 23);
    return (engine() >> (64 - 23)) * norm;
    
    /* Alternative: For perfectly uniform floats:
    static thread_local std::uniform_real_distribution<float> dist(0.0f, 1.0f);
    return dist(engine);
    */
}


// Branchless sign generation
int RandomNumberGenerator::Sign() noexcept {
    // 1 subtraction, 1 bit shift - no branches
    return 1 -  (static_cast<int>(engine() & 1) << 1);
}
