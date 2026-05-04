/*
 * pair_registration_cm.hpp — C++ port of tmp/pair_registration_cm.py
 *
 * Phase 7E.4. Tmp-only prototype. Mirrors the Python class one-to-one
 * for the surgical-update pipeline. Designed to plug into the same
 * tf_sim_v2 candidacy_step callback.
 *
 * Architecture:
 *   L1 (P × M1)   far neighbor list
 *   L2 (P × M2)   near neighbor list (used for shell registration)
 *   Q  (P × M2)   per-pair registration {n0, m0, gap}
 *   C  (P*N × E)  final per-node candidate matrix (GPU-facing)
 *
 * Update flow: tier-detect → affected set → L_refresh → Q_refresh → C_rebuild.
 */

#pragma once

#include <cstdint>
#include <cmath>
#include <vector>
#include <algorithm>
#include <limits>
#include <numeric>

constexpr int  EMPTY_NEIGHBOR = -1;
constexpr int  GHOST_CAND     =  0;
constexpr long NEVER_CYCLE    = -1;

class PairRegistrationCM {
public:
    int P, N, M1, M2, delta, E, K;
    double Lx, Ly;
    bool periodic_x = false, periodic_y = false;
    double L1_radius_scale, L2_radius_scale;
    double skin;
    double Q_skin, L2_skin, L1_skin;
    double L1_trigger_frac, L2_trigger_frac, Q_trigger_frac;
    double R0_max;

    std::vector<double> R0_arr;     // (P,)
    std::vector<double> r_c_per_p;  // (P,)

    // Neighbor lists
    std::vector<int> L1;            // (P, M1)
    std::vector<int> L2;            // (P, M2)
    std::vector<std::vector<int>> inv_L2;  // P-vec of vecs

    // Pair registration (per i-slot)
    std::vector<int>    Q_n0;       // (P, M2)
    std::vector<int>    Q_m0;       // (P, M2)
    std::vector<long>   Q_cycle;    // (P, M2) — cycle when registered

    // Final candidate matrix (read by force kernel)
    std::vector<int> C;             // (K, E) row-major

    // Per-tier last-state (for trigger checks)
    std::vector<double> L1_last_x_cm, L1_last_theta;
    std::vector<double> L2_last_x_cm, L2_last_theta;
    std::vector<double> Q_last_x_cm,  Q_last_theta;

    // Cycle flags
    long cycle_id = 0;
    std::vector<long> L1_refreshed;       // (P,)
    std::vector<long> L2_refreshed;       // (P,)
    std::vector<long> Q_refreshed;        // (P, M2)
    std::vector<long> C_patched;          // (P,)

    bool initialized = false;

    // L2_full[i]: variable-length list of all peers j within skin radius.
    // Used ONLY for cascade — when particle i triggers, all j in L2_full[i]
    // get marked dirty for rebuild. L2 (the M2-closest) is derived from L2_full.
    std::vector<std::vector<int>> L2_full;
    double D_min;          // min over i of (R0[i] + r_c[i]) — uniform threshold base
    double L2_full_skin;   // peers with cm_dist < (R0_i+R0_j+r_c_i+r_c_j)*L2_full_skin in L2_full

    // Rebuild cadence (units = PHYSICS STEPS, not call count).
    // step_counter is bumped by phys_steps_per_call on every update() call;
    // phys_steps_per_call is set by the Python wrapper from cand_check_interval.
    //
    // Defaults are aligned: cascade_interval and L1L2_rebuild_interval are
    // integer multiples of Q_refresh_interval, so cascades and full rebuilds
    // always land on Q boundaries (no wasted Q-then-rebuild on adjacent steps).
    //
    //   - every Q_refresh_interval phys steps: full Q+C refresh for ALL particles
    //   - every cascade_interval phys steps: also do trigger detect + cascade
    //                                        L1/L2 surgical update
    //   - every L1L2_rebuild_interval phys steps: full L1/L2 rebuild (subsumes Q)
    int  Q_refresh_interval     = 30;
    int  cascade_interval       = 150;  // = 5 × Q_refresh_interval (aligned)
    int  L1L2_rebuild_interval  = 300;  // = 10 × Q (= 2 × cascade) (aligned)
    int  phys_steps_per_call    = 10;   // matches default cand_check_interval
    int  forced_rebuild_interval = 300; // legacy alias for L1L2_rebuild_interval
    long step_counter = 0;              // accumulated physics steps

    // Cascade scratch (P,) — reused across calls to avoid allocation
    std::vector<uint8_t> peer_dirty;

    // Last-cycle diagnostics
    int diag_tier_L1 = 0, diag_tier_L2 = 0, diag_tier_Q = 0;
    int diag_affected = 0, diag_C_rows_patched = 0;
    int diag_cascade_dirty = 0, diag_forced_rebuild = 0;

    PairRegistrationCM(int P_, int N_,
                       const double* R0_arr_, const double* r_c_per_p_,
                       double Lx_, double Ly_,
                       int M1_=20, int M2_=12, int delta_=4,
                       double L1_radius_scale_=1.5, double L2_radius_scale_=1.0,
                       double skin_=0.5,
                       bool periodic_x_=false, bool periodic_y_=false,
                       double L1_trigger_frac_=0.25,
                       double L2_trigger_frac_=0.10,
                       double Q_trigger_frac_=0.25,
                       int L1L2_rebuild_interval_=300,
                       int Q_refresh_interval_=30,
                       int cascade_interval_=150,
                       int phys_steps_per_call_=10,
                       double L2_full_skin_=1.5)
        : P(P_), N(N_), M1(M1_), M2(M2_), delta(delta_),
          Lx(Lx_), Ly(Ly_), periodic_x(periodic_x_), periodic_y(periodic_y_),
          L1_radius_scale(L1_radius_scale_), L2_radius_scale(L2_radius_scale_),
          skin(skin_),
          L1_trigger_frac(L1_trigger_frac_),
          L2_trigger_frac(L2_trigger_frac_),
          Q_trigger_frac(Q_trigger_frac_),
          Q_refresh_interval(Q_refresh_interval_),
          cascade_interval(cascade_interval_),
          L1L2_rebuild_interval(L1L2_rebuild_interval_),
          phys_steps_per_call(phys_steps_per_call_),
          forced_rebuild_interval(L1L2_rebuild_interval_),
          L2_full_skin(L2_full_skin_)
    {
        E = 2*delta + 1;
        K = P * N;
        R0_arr.assign(R0_arr_, R0_arr_ + P);
        r_c_per_p.assign(r_c_per_p_, r_c_per_p_ + P);
        R0_max = *std::max_element(R0_arr.begin(), R0_arr.end());
        // D_min = uniform trigger base (rcpgenerator-style: trigger fires when
        // particle has consumed `frac * D_min` of its closing budget).
        D_min = std::numeric_limits<double>::infinity();
        for (int i = 0; i < P; ++i) {
            double b = R0_arr[i] + r_c_per_p[i];
            if (b < D_min) D_min = b;
        }
        double r_c_min = *std::min_element(r_c_per_p.begin(), r_c_per_p.end());
        double R0_mean = 0.0;
        for (auto v : R0_arr) R0_mean += v;
        R0_mean /= P;
        Q_skin  = Q_trigger_frac  * r_c_min;
        L2_skin = L2_trigger_frac * D_min;
        L1_skin = L1_trigger_frac * D_min;
        L2_full.assign(P, {});
        peer_dirty.assign(P, 0);

        L1.assign(P*M1, EMPTY_NEIGHBOR);
        L2.assign(P*M2, EMPTY_NEIGHBOR);
        inv_L2.assign(P, {});
        Q_n0.assign(P*M2, -1);
        Q_m0.assign(P*M2, -1);
        Q_cycle.assign(P*M2, NEVER_CYCLE);
        C.assign(K*E, GHOST_CAND);

        L1_refreshed.assign(P, NEVER_CYCLE);
        L2_refreshed.assign(P, NEVER_CYCLE);
        Q_refreshed.assign(P*M2, NEVER_CYCLE);
        C_patched.assign(P, NEVER_CYCLE);

        L1_last_x_cm.assign(P*2, 0.0);
        L1_last_theta.assign(P, 0.0);
        L2_last_x_cm.assign(P*2, 0.0);
        L2_last_theta.assign(P, 0.0);
        Q_last_x_cm.assign(P*2, 0.0);
        Q_last_theta.assign(P, 0.0);
    }

    // Force-refresh ALL Q + C every step, skip tier detection. L1/L2 untouched.
    // Used to test whether stale-Q is the failure mode.
    void update_force_Q_refresh(const double* x_all, const double* x_cm,
                                  const double* theta) {
        cycle_id += 1;
        diag_tier_L1 = diag_tier_L2 = 0;
        diag_tier_Q = P;
        diag_affected = 0;
        diag_C_rows_patched = 0;
        if (!initialized) {
            full_build(x_all, x_cm, theta);
            initialized = true;
            return;
        }
        // Refresh every particle's Q + C
        for (int i = 0; i < P; ++i) {
            refresh_Q_for_particle(i, x_all, x_cm, theta);
            Q_last_x_cm[2*i+0] = x_cm[2*i+0];
            Q_last_x_cm[2*i+1] = x_cm[2*i+1];
            Q_last_theta[i]    = theta[i];
        }
        for (int i = 0; i < P; ++i) {
            rebuild_C_for_particle(i, x_all, x_cm);
            C_patched[i] = cycle_id;
            ++diag_affected;
            diag_C_rows_patched += N;
        }
    }

    // Combined: tier-detect L1/L2 (per-particle normalized motion), surgically
    // refresh L1/L2 for triggered particles, then force-refresh Q+C on EVERY
    // particle. Intended cadence: every 3 simulation steps. Q tier check is
    // skipped (force-refresh covers it).
    void update_L1_L2_then_full_Q(const double* x_all, const double* x_cm,
                                     const double* theta) {
        if (!initialized) {
            cycle_id += 1;
            full_build(x_all, x_cm, theta);
            initialized = true;
            step_counter += phys_steps_per_call;
            return;
        }
        cycle_id += 1;
        step_counter += phys_steps_per_call;
        diag_tier_L1 = diag_tier_L2 = diag_tier_Q = 0;
        diag_affected = diag_C_rows_patched = 0;
        diag_cascade_dirty = 0;
        diag_forced_rebuild = 0;

        // 4-layer cadence (units = PHYSICS STEPS):
        //   - every call: no work (early return below if no Q this call)
        //   - every Q_refresh_interval (~30): full Q+C refresh
        //   - every cascade_interval (~120): also do trigger detect + cascade
        //                                    L1/L2 surgical update
        //   - every L1L2_rebuild_interval (~500): full L1/L2 rebuild
        bool do_full_L1L2  = (L1L2_rebuild_interval > 0 &&
                              (step_counter % L1L2_rebuild_interval) < phys_steps_per_call);
        bool do_cascade    = do_full_L1L2 ||
                             (cascade_interval > 0 &&
                              (step_counter % cascade_interval) < phys_steps_per_call);
        bool do_full_Q     = do_full_L1L2 || do_cascade ||
                             (Q_refresh_interval > 0 &&
                              (step_counter % Q_refresh_interval) < phys_steps_per_call);

        // Layer 3: full L1/L2 rebuild every L1L2_rebuild_interval phys steps.
        // Resets all baselines and refreshes everything. Safety net.
        if (do_full_L1L2) {
            full_build(x_all, x_cm, theta);
            diag_forced_rebuild = 1;
            diag_tier_L1 = P; diag_tier_Q = P;
            diag_affected = P; diag_C_rows_patched = P*N;
            return;
        }

        // Per-call (between Q refreshes): no work. Most calls fall through
        // here. Q rebuild fires every Q_refresh_interval phys steps; cascade
        // fires every cascade_interval (typically every 4 Q rebuilds).
        if (!do_full_Q) {
            return;
        }

        // Cascade L1/L2 surgical update — only when do_cascade fires (every
        // cascade_interval phys steps, e.g. every 4 Q rebuilds).
        if (do_cascade) {
            double th_L1 = L1_trigger_frac * D_min;
            std::fill(peer_dirty.begin(), peer_dirty.end(), 0);
            for (int i = 0; i < P; ++i) {
                double dxa = x_cm[2*i+0] - L1_last_x_cm[2*i+0];
                double dya = x_cm[2*i+1] - L1_last_x_cm[2*i+1];
                if (periodic_x) dxa -= Lx * std::round(dxa / Lx);
                if (periodic_y) dya -= Ly * std::round(dya / Ly);
                double slip = std::sqrt(dxa*dxa + dya*dya)
                              + R0_arr[i] * std::abs(theta[i] - L1_last_theta[i]);
                if (slip > th_L1) {
                    ++diag_tier_L1;
                    peer_dirty[i] = 1;
                    for (int j : L2_full[i]) peer_dirty[j] = 1;
                }
            }
            for (int i = 0; i < P; ++i) {
                if (!peer_dirty[i]) continue;
                build_L2_for_particle(i, x_cm);
                L1_refreshed[i] = cycle_id;
                L2_refreshed[i] = cycle_id;
                L1_last_x_cm[2*i+0] = x_cm[2*i+0];
                L1_last_x_cm[2*i+1] = x_cm[2*i+1];
                L1_last_theta[i]    = theta[i];
                L2_last_x_cm[2*i+0] = x_cm[2*i+0];
                L2_last_x_cm[2*i+1] = x_cm[2*i+1];
                L2_last_theta[i]    = theta[i];
                ++diag_cascade_dirty;
            }
            if (diag_cascade_dirty > 0) {
                rebuild_inv_L2();
            }
        }

        // Full Q+C refresh for ALL particles
        cycle_id += 1;
        for (int i = 0; i < P; ++i) {
            refresh_Q_for_particle(i, x_all, x_cm, theta);
            Q_last_x_cm[2*i+0] = x_cm[2*i+0];
            Q_last_x_cm[2*i+1] = x_cm[2*i+1];
            Q_last_theta[i]    = theta[i];
        }
        for (int i = 0; i < P; ++i) {
            rebuild_C_for_particle(i, x_all, x_cm);
            C_patched[i] = cycle_id;
            diag_C_rows_patched += N;
        }
        diag_tier_Q = P;
        diag_affected = P;
    }

    // Main entry. x_all (P*N*2), x_cm (P*2), theta (P,).
    void update(const double* x_all, const double* x_cm, const double* theta) {
        cycle_id += 1;
        diag_tier_L1 = diag_tier_L2 = diag_tier_Q = 0;
        diag_affected = diag_C_rows_patched = 0;

        if (!initialized) {
            full_build(x_all, x_cm, theta);
            initialized = true;
            return;
        }

        // Tier detection
        std::vector<int> tier(P, 0);
        for (int i = 0; i < P; ++i) {
            double dxa = x_cm[2*i+0] - L1_last_x_cm[2*i+0];
            double dya = x_cm[2*i+1] - L1_last_x_cm[2*i+1];
            if (periodic_x) dxa -= Lx * std::round(dxa / Lx);
            if (periodic_y) dya -= Ly * std::round(dya / Ly);
            double d_cm_L1 = std::sqrt(dxa*dxa + dya*dya);
            double d_th_L1 = std::abs(theta[i] - L1_last_theta[i]);
            double slip_L1 = d_cm_L1 + R0_arr[i] * d_th_L1;

            double dxb = x_cm[2*i+0] - L2_last_x_cm[2*i+0];
            double dyb = x_cm[2*i+1] - L2_last_x_cm[2*i+1];
            if (periodic_x) dxb -= Lx * std::round(dxb / Lx);
            if (periodic_y) dyb -= Ly * std::round(dyb / Ly);
            double slip_L2 = std::sqrt(dxb*dxb + dyb*dyb)
                             + R0_arr[i] * std::abs(theta[i] - L2_last_theta[i]);

            double dxc = x_cm[2*i+0] - Q_last_x_cm[2*i+0];
            double dyc = x_cm[2*i+1] - Q_last_x_cm[2*i+1];
            if (periodic_x) dxc -= Lx * std::round(dxc / Lx);
            if (periodic_y) dyc -= Ly * std::round(dyc / Ly);
            double slip_Q = std::sqrt(dxc*dxc + dyc*dyc)
                            + R0_arr[i] * std::abs(theta[i] - Q_last_theta[i]);

            // Per-particle normalized thresholds
            double th_L1 = L1_trigger_frac * R0_arr[i];
            double th_L2 = L2_trigger_frac * R0_arr[i];
            double th_Q  = Q_trigger_frac  * r_c_per_p[i];
            int t = 0;
            if (slip_Q  > th_Q)  t = 1;
            if (slip_L2 > th_L2) t = 2;
            if (slip_L1 > th_L1) t = 3;
            tier[i] = t;
            if (t == 3) ++diag_tier_L1;
            else if (t == 2) ++diag_tier_L2;
            else if (t == 1) ++diag_tier_Q;
        }

        if (diag_tier_L1 + diag_tier_L2 + diag_tier_Q == 0) return;

        // Process tiers in descending order
        for (int i = 0; i < P; ++i) {
            if (tier[i] == 3) L1_refresh_particle(i, x_all, x_cm, theta);
        }
        for (int i = 0; i < P; ++i) {
            if (tier[i] == 2) L2_refresh_particle(i, x_all, x_cm, theta);
        }
        for (int i = 0; i < P; ++i) {
            if (tier[i] == 1) Q_refresh_particle(i, x_all, x_cm, theta);
        }

        // Rebuild C for any particle whose Q got touched this cycle
        for (int i = 0; i < P; ++i) {
            if (C_patched[i] == cycle_id) continue;
            // Was any Q slot of i refreshed this cycle?
            bool touched = false;
            for (int s = 0; s < M2; ++s) {
                if (Q_refreshed[i*M2 + s] == cycle_id) { touched = true; break; }
            }
            if (touched) {
                rebuild_C_for_particle(i, x_all, x_cm);
                C_patched[i] = cycle_id;
                diag_C_rows_patched += N;
                ++diag_affected;
            }
        }
    }

    // Public entry point for forcing a from-scratch rebuild (delegates to the
    // private full_build). Used by adaptive_swell after each compress step
    // since affine CM rescaling invalidates L1/L2 ranking AND Q (n0,m0) basis.
    void force_full_rebuild_public(const double* x_all, const double* x_cm,
                                     const double* theta) {
        cycle_id += 1;
        full_build(x_all, x_cm, theta);
        initialized = true;
        step_counter = 0;   // reset cadence so subsequent calls start clean
    }

private:
    inline int& L_at(std::vector<int>& v, int i, int s) { return v[i*M2 + s]; }
    inline int& L1_at(int i, int s) { return L1[i*M1 + s]; }
    inline int& L2_at(int i, int s) { return L2[i*M2 + s]; }

    inline double pair_threshold(int i, int j, double scale) const {
        return scale * (R0_arr[i] + R0_arr[j] + skin);
    }

    void min_image(double& dx, double& dy) const {
        if (periodic_x) dx -= Lx * std::round(dx / Lx);
        if (periodic_y) dy -= Ly * std::round(dy / Ly);
    }

    int contact_facing_edge(int p, double nx, double ny, const double* theta) const {
        double angle = std::atan2(ny, nx);
        double tp = theta[p];
        double best_d = std::numeric_limits<double>::infinity();
        int best_k = 0;
        for (int k = 0; k < N; ++k) {
            double phi = 2.0 * M_PI * (k + 0.5) / N + tp;
            double d = angle - phi;
            d = std::fmod(d + M_PI, 2.0 * M_PI);
            if (d < 0) d += 2.0 * M_PI;
            d -= M_PI;
            d = std::abs(d);
            if (d < best_d) { best_d = d; best_k = k; }
        }
        return best_k;
    }

    // Closest-node-pair local search (port of sliding_candidacy.closest_node_pair_local)
    // (sx, sy) = periodic shift to apply to xB nodes so they sit in the same
    // image as xA. Caller computes once from CM-CM min-image. Open / hard
    // boundaries pass (0, 0); periodic boundaries pass the wrap offset.
    void closest_node_pair_local(const double* x_all, int pA, int pB,
                                   int i_init, int j_init,
                                   int& i_out, int& j_out,
                                   double sx = 0.0, double sy = 0.0,
                                   int max_iter = 15, int escape_steps = 3) {
        const double* xA = x_all + pA*N*2;
        const double* xB = x_all + pB*N*2;
        auto d2 = [&](int a, int b) {
            int aa = ((a % N) + N) % N;
            int bb = ((b % N) + N) % N;
            double dx = xA[aa*2+0] - (xB[bb*2+0] + sx);
            double dy = xA[aa*2+1] - (xB[bb*2+1] + sy);
            return dx*dx + dy*dy;
        };

        int i = ((i_init % N) + N) % N;
        int j = ((j_init % N) + N) % N;
        double d_best = d2(i, j);

        // Diagonal probe
        double d_fwd = d2(i+1, j-1);
        double d_bwd = d2(i-1, j+1);
        int di = 0, dj = 0;
        if (d_fwd < d_best && d_fwd <= d_bwd) {
            di = +1; dj = -1;
            i = ((i+1) % N + N) % N; j = ((j-1) % N + N) % N; d_best = d_fwd;
        } else if (d_bwd < d_best) {
            di = -1; dj = +1;
            i = ((i-1) % N + N) % N; j = ((j+1) % N + N) % N; d_best = d_bwd;
        }

        if (di != 0) {
            int walks = 0;
            while (walks < max_iter) {
                int ni = ((i+di) % N + N) % N;
                int nj = ((j+dj) % N + N) % N;
                double d_new = d2(ni, nj);
                if (d_new < d_best) {
                    i = ni; j = nj; d_best = d_new;
                    ++walks;
                } else {
                    bool escaped = false;
                    for (int s = 2; s <= escape_steps + 1; ++s) {
                        int fi = ((i+s*di) % N + N) % N;
                        int fj = ((j+s*dj) % N + N) % N;
                        double d_far = d2(fi, fj);
                        if (d_far < d_best) {
                            i = fi; j = fj; d_best = d_far;
                            escaped = true;
                            break;
                        }
                    }
                    if (!escaped) break;
                    walks += escape_steps;
                }
            }
        }

        // 8-direction refinement
        const int moves[8][2] = {{+1,0},{-1,0},{0,+1},{0,-1},
                                  {+1,+1},{+1,-1},{-1,+1},{-1,-1}};
        for (int it = 0; it < max_iter; ++it) {
            bool improved = false;
            int best_ni = i, best_nj = j;
            double best_d_new = d_best;
            for (auto& m : moves) {
                int ni = ((i+m[0]) % N + N) % N;
                int nj = ((j+m[1]) % N + N) % N;
                double d_new = d2(ni, nj);
                if (d_new < best_d_new) {
                    best_d_new = d_new; best_ni = ni; best_nj = nj;
                    improved = true;
                }
            }
            if (improved) {
                i = best_ni; j = best_nj; d_best = best_d_new;
            } else break;
        }
        i_out = i; j_out = j;
    }

    void register_pair(int i, int slot, int j,
                        const double* x_all, const double* x_cm, const double* theta) {
        // Raw CM-CM and min-image CM-CM. Their difference is the periodic
        // shift to apply to xB nodes so they sit in the same image as xA.
        // For open / hard boundaries (periodic_x or _y false) the shift on
        // that axis is exactly zero — min_image is a no-op there.
        double rx_raw = x_cm[2*j+0] - x_cm[2*i+0];
        double ry_raw = x_cm[2*j+1] - x_cm[2*i+1];
        double rx = rx_raw, ry = ry_raw;
        min_image(rx, ry);
        double sx = rx - rx_raw;   // periodic shift for xB → image of xA
        double sy = ry - ry_raw;
        double d = std::sqrt(rx*rx + ry*ry);
        if (d < 1e-12) {
            Q_n0[i*M2+slot] = -1;
            Q_m0[i*M2+slot] = -1;
            return;
        }
        double nx = rx / d, ny = ry / d;
        int jA_hint = contact_facing_edge(i, nx, ny, theta);
        int jB_hint = contact_facing_edge(j, -nx, -ny, theta);
        int n0, m0;
        closest_node_pair_local(x_all, i, j, jA_hint, jB_hint, n0, m0, sx, sy);
        Q_n0[i*M2+slot] = n0;
        Q_m0[i*M2+slot] = m0;
        Q_cycle[i*M2+slot] = cycle_id;
    }

    void build_L2_for_particle(int i, const double* x_cm) {
        // Score every j != i by NORMALIZED CM distance:
        //     score = cm_dist / (R0_i + R0_j + r_c_i + r_c_j)
        // L1 = smallest M1 scores; L2 = smallest M2 scores. No radius cutoff.
        // L2_full = ALL j with score < L2_full_skin (peers within reach for cascade).
        std::vector<std::pair<double,int>> cands;
        cands.reserve(P);
        L2_full[i].clear();
        for (int j = 0; j < P; ++j) {
            if (j == i) continue;
            double dx = x_cm[2*j+0] - x_cm[2*i+0];
            double dy = x_cm[2*j+1] - x_cm[2*i+1];
            min_image(dx, dy);
            double d = std::sqrt(dx*dx + dy*dy);
            double norm = R0_arr[i] + R0_arr[j] + r_c_per_p[i] + r_c_per_p[j];
            double score = d / norm;
            cands.emplace_back(score, j);
            if (score < L2_full_skin) {
                L2_full[i].push_back(j);   // within skin → cascade-relevant
            }
        }
        std::sort(cands.begin(), cands.end());
        for (int s = 0; s < M1; ++s) L1[i*M1+s] = EMPTY_NEIGHBOR;
        for (int s = 0; s < M2; ++s) L2[i*M2+s] = EMPTY_NEIGHBOR;
        int K1 = std::min((int)cands.size(), M1);
        int K2 = std::min((int)cands.size(), M2);
        for (int s = 0; s < K1; ++s) L1[i*M1+s] = cands[s].second;
        for (int s = 0; s < K2; ++s) L2[i*M2+s] = cands[s].second;
    }

    void rebuild_inv_L2() {
        for (auto& v : inv_L2) v.clear();
        for (int i = 0; i < P; ++i) {
            for (int s = 0; s < M2; ++s) {
                int j = L2[i*M2+s];
                if (j == EMPTY_NEIGHBOR) break;
                inv_L2[j].push_back(i);
            }
        }
    }

    void rebuild_C_for_particle(int i, const double* x_all,
                                  const double* x_cm = nullptr) {
        // Collect L2 slots and per-slot periodic shifts (xB → image of xA).
        // For non-periodic axes the shift is zero (open / hard wall behavior).
        // x_cm is optional only for backward compat — periodic correctness
        // requires it. Falls back to (0,0) shift if not supplied.
        int n_slots = 0;
        int slots_j[100], slots_n0[100], slots_m0[100];
        double slots_sx[100], slots_sy[100];
        for (int s = 0; s < M2; ++s) {
            int j = L2[i*M2+s];
            if (j == EMPTY_NEIGHBOR) break;
            slots_j[n_slots]  = j;
            slots_n0[n_slots] = Q_n0[i*M2+s];
            slots_m0[n_slots] = Q_m0[i*M2+s];
            if (x_cm != nullptr) {
                double rx_raw = x_cm[2*j+0] - x_cm[2*i+0];
                double ry_raw = x_cm[2*j+1] - x_cm[2*i+1];
                double rx = rx_raw, ry = ry_raw;
                min_image(rx, ry);
                slots_sx[n_slots] = rx - rx_raw;
                slots_sy[n_slots] = ry - ry_raw;
            } else {
                slots_sx[n_slots] = 0.0;
                slots_sy[n_slots] = 0.0;
            }
            ++n_slots;
        }
        const double* x_A = x_all + i*N*2;
        if (n_slots == 0) {
            for (int n = 0; n < N; ++n) {
                int row = i*N + n;
                for (int e = 0; e < E; ++e) C[row*E + e] = GHOST_CAND;
            }
            return;
        }
        // For each node n, find best j among L2 slots via index transport
        for (int n = 0; n < N; ++n) {
            double best_d2 = std::numeric_limits<double>::infinity();
            int best_j = -1, best_m = -1;
            double xn0 = x_A[n*2+0], xn1 = x_A[n*2+1];
            for (int s = 0; s < n_slots; ++s) {
                int j = slots_j[s];
                int n0 = slots_n0[s], m0 = slots_m0[s];
                int p = ((n - n0) % N + N) % N;
                int m = ((m0 - p) % N + N) % N;
                const double* x_B = x_all + j*N*2;
                double dx = xn0 - (x_B[m*2+0] + slots_sx[s]);
                double dy = xn1 - (x_B[m*2+1] + slots_sy[s]);
                double d2v = dx*dx + dy*dy;
                if (d2v < best_d2) { best_d2 = d2v; best_j = j; best_m = m; }
            }
            int row = i*N + n;
            for (int e = 0; e < E; ++e) {
                int off = e - delta;
                int m = ((best_m + off) % N + N) % N;
                C[row*E + e] = best_j * N + m + 1;
            }
        }
    }

    void refresh_Q_for_particle(int i, const double* x_all,
                                  const double* x_cm, const double* theta) {
        for (int s = 0; s < M2; ++s) {
            int j = L2[i*M2+s];
            if (j == EMPTY_NEIGHBOR) break;
            if (Q_refreshed[i*M2+s] == cycle_id) continue;
            register_pair(i, s, j, x_all, x_cm, theta);
            Q_refreshed[i*M2+s] = cycle_id;
        }
    }

    void Q_refresh_inv_L2_of(int i, const double* x_all,
                               const double* x_cm, const double* theta) {
        for (int k : inv_L2[i]) {
            for (int s = 0; s < M2; ++s) {
                if (L2[k*M2+s] == i) {
                    if (Q_refreshed[k*M2+s] != cycle_id) {
                        register_pair(k, s, i, x_all, x_cm, theta);
                        Q_refreshed[k*M2+s] = cycle_id;
                    }
                    break;
                }
            }
        }
    }

    void Q_refresh_particle(int i, const double* x_all,
                              const double* x_cm, const double* theta) {
        refresh_Q_for_particle(i, x_all, x_cm, theta);
        Q_refresh_inv_L2_of(i, x_all, x_cm, theta);
        Q_last_x_cm[2*i+0] = x_cm[2*i+0];
        Q_last_x_cm[2*i+1] = x_cm[2*i+1];
        Q_last_theta[i]    = theta[i];
    }

    void L2_refresh_particle(int i, const double* x_all,
                               const double* x_cm, const double* theta) {
        build_L2_for_particle(i, x_cm);
        L2_refreshed[i] = cycle_id;
        rebuild_inv_L2();
        refresh_Q_for_particle(i, x_all, x_cm, theta);
        Q_refresh_inv_L2_of(i, x_all, x_cm, theta);
        L2_last_x_cm[2*i+0] = x_cm[2*i+0];
        L2_last_x_cm[2*i+1] = x_cm[2*i+1];
        L2_last_theta[i]    = theta[i];
        Q_last_x_cm[2*i+0]  = x_cm[2*i+0];
        Q_last_x_cm[2*i+1]  = x_cm[2*i+1];
        Q_last_theta[i]     = theta[i];
    }

    void L1_refresh_particle(int i, const double* x_all,
                               const double* x_cm, const double* theta) {
        build_L2_for_particle(i, x_cm);  // builds both L1 and L2
        L1_refreshed[i] = cycle_id;
        L2_refreshed[i] = cycle_id;
        rebuild_inv_L2();
        refresh_Q_for_particle(i, x_all, x_cm, theta);
        Q_refresh_inv_L2_of(i, x_all, x_cm, theta);
        for (int axis = 0; axis < 2; ++axis) {
            L1_last_x_cm[2*i+axis] = x_cm[2*i+axis];
            L2_last_x_cm[2*i+axis] = x_cm[2*i+axis];
            Q_last_x_cm[2*i+axis]  = x_cm[2*i+axis];
        }
        L1_last_theta[i] = theta[i];
        L2_last_theta[i] = theta[i];
        Q_last_theta[i]  = theta[i];
    }

    void full_build(const double* x_all, const double* x_cm, const double* theta) {
        for (int i = 0; i < P; ++i) build_L2_for_particle(i, x_cm);
        rebuild_inv_L2();
        for (int i = 0; i < P; ++i) {
            for (int s = 0; s < M2; ++s) {
                int j = L2[i*M2+s];
                if (j == EMPTY_NEIGHBOR) break;
                register_pair(i, s, j, x_all, x_cm, theta);
                Q_refreshed[i*M2+s] = cycle_id;
            }
        }
        for (int i = 0; i < P; ++i) {
            rebuild_C_for_particle(i, x_all, x_cm);
            C_patched[i] = cycle_id;
        }
        for (int i = 0; i < P; ++i) {
            for (int axis = 0; axis < 2; ++axis) {
                L1_last_x_cm[2*i+axis] = x_cm[2*i+axis];
                L2_last_x_cm[2*i+axis] = x_cm[2*i+axis];
                Q_last_x_cm[2*i+axis]  = x_cm[2*i+axis];
            }
            L1_last_theta[i] = theta[i];
            L2_last_theta[i] = theta[i];
            Q_last_theta[i]  = theta[i];
            L1_refreshed[i]  = cycle_id;
            L2_refreshed[i]  = cycle_id;
        }
    }
};
