/*
 * bindings.cpp — pybind11 bindings for CandidacyManager.
 */

#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>
#include <stdexcept>
#include <cstring>
#include "candidacy_manager.hpp"

namespace py = pybind11;

PYBIND11_MODULE(_candidacy_cpp, m) {
    m.doc() = "C++ CandidacyManager for EPD/TF simulation pipeline";

    py::class_<CandidacyManager>(m, "CandidacyManager")
        .def(py::init<int, int, float, int, float, int>(),
             py::arg("P"),
             py::arg("N"),
             py::arg("R0"),
             py::arg("E")    = 32,
             py::arg("skin") = 0.3f,
             py::arg("dj")   = -1)

        .def("update",
             [](CandidacyManager& self,
                py::array x_cm_arr,
                py::array theta_arr) {
                 /* Accept any contiguous double array (P,2) and (P,) */
                 auto x_buf = x_cm_arr.request();
                 auto t_buf = theta_arr.request();
                 if (x_buf.ndim != 2 || x_buf.shape[0] != self.P || x_buf.shape[1] != 2)
                     throw std::runtime_error("x_cm must be shape (P, 2) float64");
                 if (t_buf.ndim != 1 || t_buf.shape[0] != self.P)
                     throw std::runtime_error("theta must be shape (P,) float64");
                 self.update(static_cast<const double*>(x_buf.ptr),
                             static_cast<const double*>(t_buf.ptr));
             },
             py::arg("x_cm"), py::arg("theta"),
             "Recompute CapCandidates from x_cm (P,2) float64 and theta (P,) float64.")

        /* Return CapCandidates as a numpy int32 array (K, E). */
        .def_property_readonly("CapCandidates",
             [](const CandidacyManager& self) {
                 int K = self.P * self.N;
                 py::array_t<int32_t> out(std::vector<py::ssize_t>{K, self.E});
                 std::memcpy(out.mutable_data(),
                             self.CapCandidates.data(),
                             (size_t)K * self.E * sizeof(int32_t));
                 return out;
             })

        .def_readonly("P",    &CandidacyManager::P)
        .def_readonly("N",    &CandidacyManager::N)
        .def_readonly("E",    &CandidacyManager::E)
        .def_readonly("R0",   &CandidacyManager::R0)
        .def_readonly("skin", &CandidacyManager::skin)
        .def_readonly("dj",   &CandidacyManager::dj)
        .def_property_readonly("K",
             [](const CandidacyManager& self){ return self.P * self.N; });
}
