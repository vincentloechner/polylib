#include <pybind11/pybind11.h>

extern "C" {
#include "polylib/polylibgmp.h"
}

namespace py = pybind11;

PYBIND11MODULE(pypolylib, m) { 
    
    py::class_<Matrix>(m, "Matrix")
        .def_readonly("nbrows", &Matrix::NbRows)
        .def_readonly("nbcolumns", &Matrix::NbColumns);
    m.def("MatrixAlloc", &Matrix_Alloc);


    py::class<Polyhedron>(m, "Polyhedron")
        .def_readonly("dimension", &Polyhedron::Dimension)
        .def_readonly("nbconstraints", &Polyhedron::NbConstraints)
        .def_readonly("nbrays", &Polyhedron::NbRays);
    m.def("Constraints2Polyhedron", &Constraints2Polyhedron);


    py::class<LBL>(m, "LBL")
        .def(py::init<>())
        .def_readwrite("Lat", &LBL::Lat)
        .def_readwrite("P", &LBL::P);
        .def_readwrite("next", &LBL::next);
        

    m.def("LBLAlloc", &LBLAlloc);

    
}

