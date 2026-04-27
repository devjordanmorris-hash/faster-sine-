Fast Sine Approximation via Prism Encoding

Overview

This project implements a high-performance approximation of the sine function using a piecewise cubic (Hermite-based) encoding, referred to here as prism tiles.

The goal is to provide a low-cost, high-throughput alternative to sinf, particularly useful in:

* numerical simulations
* signal processing
* large-scale modelling workloads
* performance-critical compute pipelines

The method trades a small, controlled approximation error for a significant speedup.

⸻

Key Idea

Instead of evaluating sin(x) directly for every input, the approach:

1. Samples a full sine cycle at fixed resolution
2. Applies an optional phase snapping (pre-rotation) step
3. Encodes segments into cubic Hermite polynomials (prisms)
4. Reconstructs values via fast polynomial evaluation

Each prism tile represents a local approximation of the sine curve.

⸻

Features

* ~5–6× speedup vs std::sinf (hardware dependent)
* Low memory footprint (compact cubic coefficients)
* Deterministic, branch-free evaluation
* Adjustable accuracy via tile size and sampling resolution

⸻

Implementation Details

Prism Representation

Each tile is encoded as: f(k) = a0 + a1*k + a2*k² + a3*k³

Where:

* k is the local index within a tile
* coefficients are derived from endpoint values and derivatives

* Pre-layer: Phase Rotation

The optional pre_rotate step snaps the signal to local anchors:

* Reduces high-frequency variation within tiles
* Improves approximation stability
* Controlled via SNAP parameter

* Encoding

Each tile uses:

* endpoint values (y0, y1)
* derivatives (cos(x) scaled by step size)

Converted into cubic Hermite form and normalized to tile width.

Decoding

Evaluation is a simple cubic polynomial:
((a3*x + a2)*x + a1)*x + a0

This is:

* fast
* vectorizable
* cache-friendly

* Accuracy

Example configuration:

* Samples: 32,768
* Tile size: 256
* Snap: 16

Typical results:

* RMS error: low (small relative error)
* Max error: bounded and stable

Accuracy can be tuned via:

* tile size (TILE)
* sample resolution (SAMPLES)
* snap granularity (SNAP)

⸻

Performance

Benchmark compares:

* std::sinf
* prism evaluation

Typical outcome: Speedup: ~5–6×

Results depend on:

* CPU architecture
* compiler optimizations
* vectorization

* When to Use

This method is suitable for:

* large-scale repeated sine evaluations
* simulation loops
* DSP-style workloads
* approximate computing scenarios

* Integration Notes

* Designed for fixed-domain evaluation (e.g. one sine cycle)
* Can be extended to arbitrary input via:
    * phase wrapping
    * lookup mapping
* Easily adaptable for SIMD or GPU execution

* xample

Run the provided main() to:

* build truth signal
* encode prisms
* decode approximation
* compute error
* benchmark performance

* License

Free to use, modify, and integrate.
No restrictions.

Notes

This is a practical performance-oriented approximation.
It is not intended to replace high-precision math libraries, but to complement them where speed is critical.

Contact

Open to feedback, improvements, and integration ideas. - dev.jordanmorris@gmail.com
