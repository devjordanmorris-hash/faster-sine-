#include <cmath>
#include <cstdint>
#include <vector>
#include <iostream>
#include <chrono>
#include <algorithm>

static constexpr double PI  = 3.14159265358979323846;
static constexpr double TAU = 2.0 * PI;

// ===============================
// Prism tile (cubic in k)
// ===============================
struct PrismTile {
    float a0, a1, a2, a3;
};

// ===============================
// Build truth sine cycle
// ===============================
std::vector<float> build_truth(uint32_t N) {
    std::vector<float> y(N);
    double d = TAU / double(N);
    for (uint32_t i = 0; i < N; ++i)
        y[i] = std::sin(i * d);
    return y;
}

// ===============================
// PRE-LAYER: phase rotation
// Snap phase to nearest prism anchor
// ===============================
void pre_rotate(std::vector<float>& y, int snap) {
    if (snap <= 1) return;

    const size_t N = y.size();
    std::vector<float> tmp(N);

    for (size_t i = 0; i < N; ++i) {
        size_t r = (i / snap) * snap;
        tmp[i] = y[r];
    }
    y.swap(tmp);
}

// ===============================
// Hermite cubic → prism coefficients
// ===============================
PrismTile make_prism(
    float y0, float y1,
    float dy0, float dy1,
    int T
) {
    float m0 = dy0 * (T - 1);
    float m1 = dy1 * (T - 1);

    float c0 = y0;
    float c1 = m0;
    float c2 = 3*(y1-y0) - 2*m0 - m1;
    float c3 = 2*(y0-y1) + m0 + m1;

    float s  = 1.0f / float(T - 1);
    float s2 = s*s;
    float s3 = s2*s;

    return {
        c0,
        c1*s,
        c2*s2,
        c3*s3
    };
}

// ===============================
// Encode cycle → prisms
// ===============================
std::vector<PrismTile> encode_prisms(
    const std::vector<float>& y,
    int TILE
) {
    const size_t N = y.size();
    const double dtheta = TAU / double(N);

    std::vector<PrismTile> out;
    out.reserve(N / TILE);

    for (size_t i = 0; i < N; i += TILE) {
        size_t j = i + TILE - 1;

        double t0 = i * dtheta;
        double t1 = j * dtheta;

        float dy0 = std::cos(t0) * dtheta;
        float dy1 = std::cos(t1) * dtheta;

        out.push_back(
            make_prism(y[i], y[j], dy0, dy1, TILE)
        );
    }
    return out;
}

// ===============================
// Decode prism
// ===============================
inline float eval(const PrismTile& p, int k) {
    float x = float(k);
    return ((p.a3*x + p.a2)*x + p.a1)*x + p.a0;
}

// ===============================
// Error stats
// ===============================
void error_stats(
    const std::vector<float>& a,
    const std::vector<float>& b
) {
    double rms = 0.0;
    double maxe = 0.0;

    for (size_t i = 0; i < a.size(); ++i) {
        double e = a[i] - b[i];
        rms += e*e;
        maxe = std::max(maxe, std::abs(e));
    }
    rms = std::sqrt(rms / a.size());

    std::cout << "RMS = " << rms
              << "   MAX = " << maxe << "\n";
}

// ===============================
// MAIN
// ===============================
int main() {
    constexpr int SAMPLES = 32768;
    constexpr int TILE    = 256;
    constexpr int SNAP    = 16;   // pre-layer rotation

    auto truth = build_truth(SAMPLES);

    auto rotated = truth;
    pre_rotate(rotated, SNAP);

    auto prisms = encode_prisms(rotated, TILE);

    std::vector<float> decoded(SAMPLES);
    for (size_t i = 0; i < prisms.size(); ++i) {
        for (int k = 0; k < TILE; ++k)
            decoded[i*TILE + k] = eval(prisms[i], k);
    }

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
        for (size_t i = 0; i < prisms.size(); ++i)
            for (int k = 0; k < TILE; ++k)
                volatile float y = eval(prisms[i], k);
    auto t3 = std::chrono::high_resolution_clock::now();

    double sin_ms =
        std::chrono::duration<double, std::milli>(t1-t0).count();
    double prism_ms =
        std::chrono::duration<double, std::milli>(t3-t2).count();

    std::cout << "\n--- Timing ---\n";
    std::cout << "sinf:   " << sin_ms   << " ms\n";
    std::cout << "prism:  " << prism_ms << " ms\n";
    std::cout << "speedup: " << (sin_ms / prism_ms) << "x\n";
}