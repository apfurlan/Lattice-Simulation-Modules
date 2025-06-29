#include "RandomGeneratorClass.hpp"
#include <cmath>

thread_local std::mt19937_64 RandomGenerator::engine(std::random_device{}());

int RandomGenerator::Int(int min, int max) noexcept {
    // Fast modulo method with range adjustment
    return min + static_cast<int>(engine() % static_cast<unsigned>(max - min + 1));
    
    /* Uncomment for perfect uniformity (slightly slower):
    static std::uniform_int_distribution<int> dist;
    return dist(engine, std::uniform_int_distribution<int>::param_type(min, max));
    */
}

float RandomGenerator::Float() noexcept {
    // Bit manipulation technique (extremely fast)
    constexpr float factor = 1.0f / (1ULL << 24);
    return (engine() >> 40) * factor;
    
    /* Uncomment for perfect uniformity:
    static std::uniform_real_distribution<float> dist(0.0f, 1.0f);
    return dist(engine);
    */
}

int RandomGenerator::Sign() noexcept {
    // Branchless implementation: 1 AND, 1 shift, 1 subtract
    return 1 - ((engine() & 1) << 1);
}