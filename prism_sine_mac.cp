#include <cmath>
#include <cstdint>
#include <vector>
#include <iostream>
#include <chrono>
#include <algorithm>

static constexpr double PI  = 3.14159265358979323846;
static constexpr double TAU = 2.0 * PI;

// ===============================
// Cubic tile coefficients (Hermite)
// ===============================
struct TileCubic {
    float a0, a1, a2, a3;   // y(k) = a0 + a1*k + a2*k² + a3*k³, k in [0, T-1]
};

// ===============================
// Build truth sine cycle (single period)
// ===============================
std::vector<float> build_truth(uint32_t N) {
    std::vector<float> y(N);
    double dtheta = TAU / double(N);
    for (uint32_t i = 0; i < N; ++i)
        y[i] = std::sin(i * dtheta);
    return y;
}

// ===============================
// Hermite cubic -> tile coefficients
// T = number of samples per tile
// ===============================
TileCubic hermite_to_tile(
    float y0, float y1,
    float dy0, float dy1,
    int T
) {
    // convert derivative w.r.t. continuous index k to derivative w.r.t. t in [0,1]
    float m0 = dy0 * (T - 1);
    float m1 = dy1 * (T - 1);

    // Hermite coefficients in t
    float c0 = y0;
    float c1 = m0;
    float c2 = 3.0f*(y1 - y0) - 2.0f*m0 - m1;
    float c3 = 2.0f*(y0 - y1) + m0 + m1;

    // Convert t = k/(T-1) to polynomial in k
    float s  = 1.0f / float(T - 1);
    float s2 = s * s;
    float s3 = s2 * s;

    return {
        c0,
        c1 * s,
        c2 * s2,
        c3 * s3
    };
}

// ===============================
// Encode a full sine cycle into tiles
// ===============================
std::vector<TileCubic> encode_tiles(
    const std::vector<float>& y,
    int TILE
) {
    const size_t N = y.size();
    const double dtheta = TAU / double(N);

    std::vector<TileCubic> out;
    out.reserve(N / TILE);

    for (size_t i = 0; i < N; i += TILE) {
        size_t j = i + TILE - 1;   // last index in tile

        double t0 = i * dtheta;
        double t1 = j * dtheta;

        // derivatives: dy/dk = cos(theta) * dtheta
        float dy0 = std::cos(t0) * dtheta;
        float dy1 = std::cos(t1) * dtheta;

        out.push_back(
            hermite_to_tile(y[i], y[j], dy0, dy1, TILE)
        );
    }
    return out;
}

// ===============================
// Evaluate a tile at sample index k (0 <= k < TILE)
// ===============================
inline float eval_tile(const TileCubic& t, int k) {
    float x = float(k);
    // Horner's method
    return ((t.a3 * x + t.a2) * x + t.a1) * x + t.a0;
}

// ===============================
// Decode full waveform from tiles
// ===============================
std::vector<float> decode_tiles(
    const std::vector<TileCubic>& tiles,
    int TILE,
    size_t total_samples
) {
    std::vector<float> out(total_samples);
    for (size_t i = 0; i < tiles.size(); ++i) {
        size_t base = i * TILE;
        for (int k = 0; k < TILE; ++k)
            out[base + k] = eval_tile(tiles[i], k);
    }
    return out;
}

// ===============================
// Error statistics
// ===============================
void error_stats(
    const std::vector<float>& a,
    const std::vector<float>& b
) {
    double rms = 0.0, maxe = 0.0;
    for (size_t i = 0; i < a.size(); ++i) {
        double e = a[i] - b[i];
        rms += e*e;
        maxe = std::max(maxe, std::abs(e));
    }
    rms = std::sqrt(rms / a.size());
    std::cout << "RMS error = " << rms
              << "   Max error = " << maxe << "\n";
}

// ===============================
// MAIN
// ===============================
int main() {
    constexpr int SAMPLES = 32768;   // must be divisible by TILE
    constexpr int TILE    = 256;

    std::cout << "Generating " << SAMPLES
              << " samples of sin, tile size = " << TILE << "\n";

    auto truth = build_truth(SAMPLES);
    auto tiles = encode_tiles(truth, TILE);
    auto decoded = decode_tiles(tiles, TILE, SAMPLES);

    std::cout << "--- Accuracy ---\n";
    error_stats(decoded, truth);

    // ===============================
    // Benchmark
    // ===============================
    constexpr int LOOPS = 200;

    auto t0 = std::chrono::high_resolution_clock::now();
    for (int r = 0; r < LOOPS; ++r)
        for (int i = 0; i < SAMPLES; ++i)
            volatile float x = std::sinf(i * TAU / SAMPLES);
    auto t1 = std::chrono::high_resolution_clock::now();

    auto t2 = std::chrono::high_resolution_clock::now();
    for (int r = 0; r < LOOPS; ++r)
        for (size_t i = 0; i < tiles.size(); ++i)
            for (int k = 0; k < TILE; ++k)
                volatile float y = eval_tile(tiles[i], k);
    auto t3 = std::chrono::high_resolution_clock::now();

    double sin_ms =
        std::chrono::duration<double, std::milli>(t1-t0).count();
    double tile_ms =
        std::chrono::duration<double, std::milli>(t3-t2).count();

    std::cout << "\n--- Timing ---\n";
    std::cout << "std::sinf:   " << sin_ms   << " ms\n";
    std::cout << "tile decode: " << tile_ms  << " ms\n";
    std::cout << "Speedup:     " << (sin_ms / tile_ms) << "x\n";

    return 0;
}
