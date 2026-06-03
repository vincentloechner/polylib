// Commande de compilation de ce fichier source en un module python :
// (à partir du répertoire Test/Zpolytest/PyPolyLib/)
// c++ -Wall -shared -std=c++11 -fPIC pypolylib.cpp ../../../.libs/libpolylibgmp.so /
//     $(python3 -m pybind11 --includes) /
//     -o pypolylib$(python3-config --extension-suffix) /
//     -I../../../include
//
// export LD_LIBRARY_PATH=$HOME/polylib/.libs/
// python3
// >>> import pypolylib
// >>> a = pypolylib.LBL()



//import pypolylib
//# Créer une matrice 2x3
//m = pypolylib.MatrixAlloc(2, 3)
//print(m.nbrows, m.nbcolumns)  # -> 2 3




#include <gmpxx.h>   // DOIT être avant polylibgmp.h
#include <pybind11/pybind11.h>

extern "C" {
#include "polylib/polylibgmp.h"
}

namespace py = pybind11;

// --- Destructeur personnalisé pour LBL ---
struct LBLDeleter {
    void operator()(LBL* p) const {
        if (p) LBLFree(p);
    }
};
using LBLPtr = std::unique_ptr<LBL, LBLDeleter>;

// --- Idem pour Matrix et Polyhedron ---
struct MatrixDeleter {
    void operator()(Matrix* p) const {
        if (p) Matrix_Free(p);
    }
};
using MatrixPtr = std::unique_ptr<Matrix, MatrixDeleter>;

struct PolyhedronDeleter {
    void operator()(Polyhedron* p) const {
        if (p) Polyhedron_Free(p);
    }
};
using PolyhedronPtr = std::unique_ptr<Polyhedron, PolyhedronDeleter>;


PYBIND11_MODULE(pypolylib, m) {

    // --- Matrix ---
    py::class_<Matrix, MatrixPtr>(m, "Matrix")
        .def_readonly("nbrows",    &Matrix::NbRows)
        .def_readonly("nbcolumns", &Matrix::NbColumns);

    m.def("MatrixAlloc", [](unsigned nbrows, unsigned nbcols) {
        return MatrixPtr(Matrix_Alloc(nbrows, nbcols));
    }, py::arg("nbrows"), py::arg("nbcols"));

    m.def("MatrixReadFromString", [](const std::string &s) {
        FILE *old_stdin = stdin;
        stdin = fmemopen((void*)s.c_str(), s.size(), "r");
        Matrix *mat = Matrix_Read();
        fclose(stdin);
        stdin = old_stdin;
        return MatrixPtr(mat);
    });

    // --- Polyhedron ---
    py::class_<Polyhedron, PolyhedronPtr>(m, "Polyhedron")
        .def_readonly("dimension",      &Polyhedron::Dimension)
        .def_readonly("nbconstraints",  &Polyhedron::NbConstraints)
        .def_readonly("nbrays",         &Polyhedron::NbRays);

    m.def("Constraints2Polyhedron", [](Matrix *m, unsigned flags) {
        return PolyhedronPtr(Constraints2Polyhedron(m, flags));
    }, py::arg("matrix"), py::arg("flags") = 0);


    // --- LBL ---
    py::class_<LBL, LBLPtr>(m, "LBL")
        .def_property_readonly("Lat",
            [](LBL &l) { return l.Lat; },
            py::return_value_policy::reference)
        .def_property_readonly("P",
            [](LBL &l) { return l.P; },
            py::return_value_policy::reference)
        .def_property_readonly("next",
            [](LBL &l) { return l.next; },
            py::return_value_policy::reference);

    m.def("LBLAlloc", [](Matrix *lat, Polyhedron *domain) {
        return LBLPtr(LBLAlloc(lat, domain));
    }, py::arg("lat"), py::arg("domain"));




    // Setter pour remplir une matrice case par case
    m.def("MatrixSetValue", [](Matrix *mat, int i, int j, long val) {
        mpz_set_si(mat->p[i][j], val);   // ← mpz_set_si au lieu de value_assign
    }, py::arg("mat"), py::arg("i"), py::arg("j"), py::arg("val"));

    m.def("MatrixGetValue", [](Matrix *mat, int i, int j) -> long {
        return mpz_get_si(mat->p[i][j]);
    }, py::arg("mat"), py::arg("i"), py::arg("j"));

    // Affichage d'un LBL
    m.def("LBLPrint", [](LBL *l) {
        LBLPrint(stdout, "%Zd", l);   // ← format GMP correct
    });
 }