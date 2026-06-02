
// Commande de compilation de ce fichier source en un module python :
// (à partir du répertoire Test/Zpolytest/PyPolyLib/)
// c++ -Wall -shared -std=c++11 -fPIC pypolylib.cpp ../../../.libs/libpolylibgmp.so $(python3 -m pybind11 --includes) -o pypolylib$(python3-config --extension-suffix) -I../../../include
// cette commande génère un fichier 'pypolylib.cpython-312-x86_64-linux-gnu.so' qui peut être utilisé comme un module python.
// 
// Ensuite pour lancer python3 il faut ajouter l'emplacement où trouver libpolylibgmp.so dans cette variable d'environnement :
// export LD_LIBRARY_PATH=$HOME/polylib/.libs/
// 
// Puis on peut lancer :
// python3
// >>> import pypolylib
// >>> a = pypolylib.LBL()


#include <pybind11/pybind11.h>

#include <gmpxx.h>
extern "C" {
#include "polylib/polylibgmp.h"
}

namespace py = pybind11;

PYBIND11_MODULE(pypolylib, m) { 
    
    py::class_<Matrix>(m, "Matrix")
        .def_readonly("nbrows", &Matrix::NbRows)
        .def_readonly("nbcolumns", &Matrix::NbColumns);
    m.def("MatrixAlloc", &Matrix_Alloc);


    py::class_<Polyhedron>(m, "Polyhedron")
        .def_readonly("dimension", &Polyhedron::Dimension)
        .def_readonly("nbconstraints", &Polyhedron::NbConstraints)
        .def_readonly("nbrays", &Polyhedron::NbRays);
    m.def("Constraints2Polyhedron", &Constraints2Polyhedron);


    py::class_<LBL>(m, "LBL")
        .def(py::init<>())
        .def_readwrite("Lat", &LBL::Lat)
        .def_readwrite("P", &LBL::P)
        .def_readwrite("next", &LBL::next);
        

    m.def("LBLAlloc", &LBLAlloc);

    
}

