/*
 * candidacy_manager.hpp — C++ CandidacyManager for the EPD/TF sim pipeline.
 *
 * Mirrors the Python CandidacyManager interface exactly so Phase 3.4 corpus
 * replay validation can confirm bit-level identical CapCandidates outputs.
 *
 * Three-level filtering:
 *   Level 1: center-center distance (brute-force O(P²); cell list later)
 *   Level 2: contact-normal registration via angle arithmetic
 *   Level 3: index-sliding fill ±dj around the registered contact edge
 *
 * CapCandidates: K×E int32 matrix, row-major, 0 = ghost (far/unused).
 * Edges are 1-indexed: edge k of particle p = p*N + e + 1 (1-based global).
 */

#pragma once
#include <cstdint>
#include <cmath>
#include <vector>
#include <algorithm>
#include <stdexcept>
#include <cstring>

class CandidacyManager {
public:
    int P, N, E;
    float R0, skin;
    int   dj;

    /* CapCandidates: (K, E) int32, row-major, 0 = ghost. */
    std::vector<int32_t> CapCandidates;

    /* ── constructor ─────────────────────────────────────────────────────── */
    CandidacyManager(int P, int N, float R0, int E = 32,
                     float skin = 0.3f, int dj = -1)
        : P(P), N(N), E(E), R0(R0), skin(skin)
    {
        if (dj < 0)
            this->dj = std::max(3, N / 8);
        else
            this->dj = dj;

        K = P * N;
        CapCandidates.assign(K * E, 0);
    }

    /* ── update: recompute CapCandidates ─────────────────────────────────── */
    /*
     * x_cm  : flat array of P*(x,y) pairs  — length 2*P
     * theta : flat array of P orientations — length P
     */
    void update(const double* x_cm, const double* theta) {
        std::fill(CapCandidates.begin(), CapCandidates.end(), 0);
        float threshold = 2.0f * R0 + skin;

        for (int pA = 0; pA < P; ++pA) {
            for (int pB = pA + 1; pB < P; ++pB) {
                double dx = x_cm[2*pB]   - x_cm[2*pA];
                double dy = x_cm[2*pB+1] - x_cm[2*pA+1];
                double d  = std::sqrt(dx*dx + dy*dy);
                if (d < threshold)
                    fill_pair(pA, pB, x_cm, theta);
            }
        }
    }

    /* float overload for convenience */
    void update_f(const float* x_cm_f, const float* theta_f) {
        std::vector<double> xd(2*P), td(P);
        for (int i = 0; i < 2*P; ++i) xd[i] = x_cm_f[i];
        for (int i = 0; i <   P; ++i) td[i] = theta_f[i];
        update(xd.data(), td.data());
    }

private:
    int K;

    /* Edge k of a regular N-gon has outward normal at angle:
     *   phi_k_world = 2π*(k+0.5)/N + theta  (midpoint of arc k) */
    int contact_facing_edge(int p_idx, double nx, double ny,
                            const double* theta) const {
        double angle_contact = std::atan2(ny, nx);
        int    best          = 0;
        double best_diff     = 1e30;
        for (int k = 0; k < N; ++k) {
            double phi = 2.0 * M_PI * (k + 0.5) / N + theta[p_idx];
            double diff = angle_contact - phi;
            /* wrap to [-π, π] */
            diff = diff - 2.0 * M_PI * std::round(diff / (2.0 * M_PI));
            double adiff = std::fabs(diff);
            if (adiff < best_diff) { best_diff = adiff; best = k; }
        }
        return best;
    }

    void fill_pair(int pA, int pB,
                   const double* x_cm, const double* theta) {
        double dx  = x_cm[2*pB]   - x_cm[2*pA];
        double dy  = x_cm[2*pB+1] - x_cm[2*pA+1];
        double dist = std::sqrt(dx*dx + dy*dy);
        if (dist < 1e-12) return;
        double nx = dx / dist, ny = dy / dist;

        int jA = contact_facing_edge(pA,  nx,  ny, theta);
        int jB = contact_facing_edge(pB, -nx, -ny, theta);

        for (int de = -dj; de <= dj; ++de) {
            int eA  = ((jA + de) % N + N) % N;
            int row = pA * N + eA;           /* 0-indexed row in CapCandidates */

            /* count slots already used */
            int slots_used = 0;
            for (int s = 0; s < E; ++s)
                if (CapCandidates[row * E + s] != 0) slots_used++;

            for (int dc = -dj; dc <= dj; ++dc) {
                int eB       = ((jB + dc) % N + N) % N;
                int32_t cand = (int32_t)(pB * N + eB + 1);   /* 1-indexed */
                if (slots_used >= E) break;
                /* avoid duplicates */
                bool dup = false;
                for (int s = 0; s < slots_used; ++s)
                    if (CapCandidates[row * E + s] == cand) { dup = true; break; }
                if (!dup) {
                    CapCandidates[row * E + slots_used] = cand;
                    slots_used++;
                }
            }
        }
    }
};
