/*
 * bindings.cpp — pybind11 bindings for tmp PairRegistrationCM
 */

#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>
#include <pybind11/stl.h>
#include <stdexcept>
#include <cstring>
#include "pair_registration_cm.hpp"

namespace py = pybind11;

PYBIND11_MODULE(_pair_registration_cm, m) {
    m.doc() = "C++ PairRegistrationCM (Phase 7E.4 prototype, tmp/-only)";

    py::class_<PairRegistrationCM>(m, "PairRegistrationCM")
        .def(py::init([](int P, int N,
                          py::array R0_arr, py::array r_c_per_p,
                          double Lx, double Ly,
                          int M1, int M2, int delta,
                          double L1_radius_scale, double L2_radius_scale,
                          double skin,
                          bool periodic_x, bool periodic_y,
                          double L1_trigger_frac,
                          double L2_trigger_frac,
                          double Q_trigger_frac,
                          int L1L2_rebuild_interval,
                          int Q_refresh_interval,
                          int cascade_interval,
                          int phys_steps_per_call,
                          double L2_full_skin) {
            auto rb = R0_arr.request();
            auto cb = r_c_per_p.request();
            if (rb.ndim != 1 || rb.shape[0] != P)
                throw std::runtime_error("R0_arr must be (P,) float64");
            if (cb.ndim != 1 || cb.shape[0] != P)
                throw std::runtime_error("r_c_per_p must be (P,) float64");
            return new PairRegistrationCM(
                P, N,
                static_cast<const double*>(rb.ptr),
                static_cast<const double*>(cb.ptr),
                Lx, Ly, M1, M2, delta,
                L1_radius_scale, L2_radius_scale, skin,
                periodic_x, periodic_y,
                L1_trigger_frac, L2_trigger_frac, Q_trigger_frac,
                L1L2_rebuild_interval, Q_refresh_interval,
                cascade_interval, phys_steps_per_call, L2_full_skin);
        }),
            py::arg("P"), py::arg("N"),
            py::arg("R0_arr"), py::arg("r_c_per_p"),
            py::arg("Lx"), py::arg("Ly"),
            py::arg("M1") = 20, py::arg("M2") = 12, py::arg("delta") = 4,
            py::arg("L1_radius_scale") = 1.5, py::arg("L2_radius_scale") = 1.0,
            py::arg("skin") = 0.5,
            py::arg("periodic_x") = false, py::arg("periodic_y") = false,
            py::arg("L1_trigger_frac") = 0.25,
            py::arg("L2_trigger_frac") = 0.10,
            py::arg("Q_trigger_frac")  = 0.25,
            py::arg("L1L2_rebuild_interval") = 300,
            py::arg("Q_refresh_interval") = 30,
            py::arg("cascade_interval") = 150,
            py::arg("phys_steps_per_call") = 10,
            py::arg("L2_full_skin") = 1.5)

        .def("update", [](PairRegistrationCM& self,
                            py::array x_all, py::array x_cm, py::array theta) {
            auto a = x_all.request();
            auto c = x_cm.request();
            auto t = theta.request();
            if (a.ndim != 3 || a.shape[0] != self.P ||
                a.shape[1] != self.N || a.shape[2] != 2)
                throw std::runtime_error("x_all must be (P, N, 2) float64");
            if (c.ndim != 2 || c.shape[0] != self.P || c.shape[1] != 2)
                throw std::runtime_error("x_cm must be (P, 2) float64");
            if (t.ndim != 1 || t.shape[0] != self.P)
                throw std::runtime_error("theta must be (P,) float64");
            self.update(static_cast<const double*>(a.ptr),
                        static_cast<const double*>(c.ptr),
                        static_cast<const double*>(t.ptr));
        },
            py::arg("x_all"), py::arg("x_cm"), py::arg("theta"))

        .def("update_force_Q_refresh",
             [](PairRegistrationCM& self,
                py::array x_all, py::array x_cm, py::array theta) {
            auto a = x_all.request();
            auto c = x_cm.request();
            auto t = theta.request();
            self.update_force_Q_refresh(
                static_cast<const double*>(a.ptr),
                static_cast<const double*>(c.ptr),
                static_cast<const double*>(t.ptr));
        },
            py::arg("x_all"), py::arg("x_cm"), py::arg("theta"))

        .def("update_L1_L2_then_full_Q",
             [](PairRegistrationCM& self,
                py::array x_all, py::array x_cm, py::array theta) {
            auto a = x_all.request();
            auto c = x_cm.request();
            auto t = theta.request();
            const double* aptr = static_cast<const double*>(a.ptr);
            const double* cptr = static_cast<const double*>(c.ptr);
            const double* tptr = static_cast<const double*>(t.ptr);
            // Release GIL during compute so a worker thread can actually
            // run in parallel with the main TF thread.
            py::gil_scoped_release release;
            self.update_L1_L2_then_full_Q(aptr, cptr, tptr);
        },
            py::arg("x_all"), py::arg("x_cm"), py::arg("theta"))

        // Forced full L1/L2/Q/C rebuild from scratch — for use after events
        // that invalidate the CM-distance ranking (affine box compression in
        // adaptive_swell, large external displacements). Resets step_counter
        // so the cadence starts fresh from this anchor point.
        .def("force_full_rebuild",
             [](PairRegistrationCM& self,
                py::array x_all, py::array x_cm, py::array theta) {
            auto a = x_all.request();
            auto c = x_cm.request();
            auto t = theta.request();
            py::gil_scoped_release release;
            self.force_full_rebuild_public(
                static_cast<const double*>(a.ptr),
                static_cast<const double*>(c.ptr),
                static_cast<const double*>(t.ptr));
        },
            py::arg("x_all"), py::arg("x_cm"), py::arg("theta"))

        .def_property_readonly("CapCandidates", [](PairRegistrationCM& self) {
            py::array_t<int32_t> out(std::vector<py::ssize_t>{self.K, self.E});
            int32_t* dst = out.mutable_data();
            // self.C is std::vector<int>; assume same width as int32_t (true on most platforms)
            for (size_t k = 0; k < self.C.size(); ++k) dst[k] = (int32_t)self.C[k];
            return out;
        })

        .def_readonly("P",     &PairRegistrationCM::P)
        .def_readonly("N",     &PairRegistrationCM::N)
        .def_readonly("E",     &PairRegistrationCM::E)
        .def_readonly("M1",    &PairRegistrationCM::M1)
        .def_readonly("M2",    &PairRegistrationCM::M2)
        .def_readonly("delta", &PairRegistrationCM::delta)
        .def_readonly("K",     &PairRegistrationCM::K)
        .def_readwrite("Lx",   &PairRegistrationCM::Lx)
        .def_readwrite("Ly",   &PairRegistrationCM::Ly)
        // Set both axes' box dims at once — call after any box compression so
        // C++ min_image computations use the current Lx/Ly. Without this,
        // PRCM keeps the construction-time Lx/Ly even after Python sets new
        // values (the readwrite individual properties handle that case too,
        // but set_box is more explicit).
        .def("set_box", [](PairRegistrationCM& self, double Lx, double Ly) {
            self.Lx = Lx;
            self.Ly = Ly;
        }, py::arg("Lx"), py::arg("Ly"))
        .def_readonly("Q_skin", &PairRegistrationCM::Q_skin)
        .def_readonly("L2_skin", &PairRegistrationCM::L2_skin)
        .def_readonly("L1_skin", &PairRegistrationCM::L1_skin)
        .def_readonly("diag_tier_L1", &PairRegistrationCM::diag_tier_L1)
        .def_readonly("diag_tier_L2", &PairRegistrationCM::diag_tier_L2)
        .def_readonly("diag_tier_Q",  &PairRegistrationCM::diag_tier_Q)
        .def_readonly("diag_affected", &PairRegistrationCM::diag_affected)
        .def_readonly("diag_C_rows_patched", &PairRegistrationCM::diag_C_rows_patched)
        .def_readonly("diag_cascade_dirty", &PairRegistrationCM::diag_cascade_dirty)
        .def_readonly("diag_forced_rebuild", &PairRegistrationCM::diag_forced_rebuild)
        .def_readonly("step_counter", &PairRegistrationCM::step_counter)
        .def_readonly("D_min", &PairRegistrationCM::D_min)
        .def_readwrite("forced_rebuild_interval",
                        &PairRegistrationCM::forced_rebuild_interval)
        .def_readwrite("L1L2_rebuild_interval",
                        &PairRegistrationCM::L1L2_rebuild_interval)
        .def_readwrite("Q_refresh_interval",
                        &PairRegistrationCM::Q_refresh_interval)
        .def_readwrite("cascade_interval",
                        &PairRegistrationCM::cascade_interval)
        .def_readwrite("phys_steps_per_call",
                        &PairRegistrationCM::phys_steps_per_call);
}
