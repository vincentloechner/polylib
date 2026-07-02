// Commande de compilation de ce fichier source en un module python :
// (à partir du répertoire Test/Zpolytest/PyPolyLib/)
// c++ -Wall -shared -std=c++11 -fPIC pypolylib_core.cpp ../../../.libs/libpolylibgmp.so    $(python3 -m pybind11 --includes) -o pypolylib_core$(python3-config --extension-suffix)      -I../../../include
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
//# Create Matrix 2x3
//m = pypolylib.MatrixAlloc(2, 3)
//print(m.nbrows, m.nbcolumns)  # -> 2 3




#include <gmpxx.h>   // MUST come before polylibgmp.h
#include <pybind11/pybind11.h>

extern "C" {
#include "polylib/polylibgmp.h"
}

namespace py = pybind11;

// --- Custom Shredder for LBL ---
struct LBLDeleter {
    void operator()(LBL* p) const {
        if (p) LBLFree(p);
    }
};
using LBLPtr = std::unique_ptr<LBL, LBLDeleter>;

// --- The same goes for Matrix and Polyhedron ---
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


PYBIND11_MODULE(pypolylib_core, m) {

    // ----------------------- Matrix ----------------------
    py::class_<Matrix, MatrixPtr>(m, "Matrix")
        // Matrix properties
        .def_readonly("nbrows",    &Matrix::NbRows)
        .def_readonly("nbcolumns", &Matrix::NbColumns);

    // Setter to fill a matrix cell by cell
    m.def("MatrixSetValue", [](Matrix *mat, int i, int j, long val) {
        mpz_set_si(mat->p[i][j], val);   // ← mpz_set_si instead of value_assign
    }, py::arg("mat"), py::arg("i"), py::arg("j"), py::arg("val"));

    m.def("MatrixGetValue", [](Matrix *mat, int i, int j) -> long {
        return mpz_get_si(mat->p[i][j]);
    }, py::arg("mat"), py::arg("i"), py::arg("j"));

    m.def("MatrixPrint", [](Matrix *mat) {
        Matrix_Print(stdout, " %s", mat);
    }), py::arg("mat");

    
    m.def("MatrixAlloc", [](unsigned nbrows, unsigned nbcols) {
        return MatrixPtr(Matrix_Alloc(nbrows, nbcols));
    }, py::arg("nbrows"), py::arg("nbcols"));

    m.def("MatrixReadFromString", [](const std::string &s) {
        FILE *f = tmpfile();
        fputs(s.c_str(), f);
        rewind(f);
        // read dimensions
        unsigned nb_rows, nb_cols;
        fscanf(f, "%u %u", &nb_rows, &nb_cols);
        // Allocate and Fill
        Matrix *mat = Matrix_Alloc(nb_rows, nb_cols);
        Matrix_Read_InputFile(mat, f);
        fclose(f);
        return MatrixPtr(mat);
    });
    // Product of two matrices
    m.def("MatrixProduct", [](Matrix *a, Matrix *b) {
        Matrix *result = Matrix_Alloc(a->NbRows, b->NbColumns);
        Matrix_Product(a, b, result);
        return MatrixPtr(result);
    }, py::arg("a"), py::arg("b"));

    // Inverse of a matrix
    m.def("MatrixInverse", [](Matrix *mat) {
        Matrix *result = Matrix_Alloc(mat->NbRows, mat->NbColumns);
        int ok = Matrix_Inverse(mat, result);
        if (!ok) {
            Matrix_Free(result);
            throw std::runtime_error("The matrix is not invertible");
        }
        return MatrixPtr(result);
    }, py::arg("mat"));


    // --------------------- Polyhedron ---------------------------
    py::class_<Polyhedron, PolyhedronPtr>(m, "Polyhedron")
        // Polyhedron properties
        .def_readonly("dimension",      &Polyhedron::Dimension)
        .def_readonly("nbconstraints",  &Polyhedron::NbConstraints)
        .def_readonly("nbrays",         &Polyhedron::NbRays)
        .def_readonly("nbbid",          &Polyhedron::NbBid)
        .def_property_readonly("next",
            [](Polyhedron &p) -> Polyhedron* { return p.next; },
            py::return_value_policy::reference)
        .def_property_readonly("constraints", [](Polyhedron &p) {
            return MatrixPtr(Polyhedron2Constraints(&p));
        })
        .def_property_readonly("rays", [](Polyhedron &p) {
            return MatrixPtr(Polyhedron2Rays(&p));
        });

    // -- Polyhedron methods --
    m.def("PolyhedronScan", [](Polyhedron *D, Polyhedron *C, unsigned NbMaxRays) {
        // single polyhedron scan only, unset/restore next:
        Polyhedron *next = D->next;
        D->next = NULL;
        Polyhedron *res = Polyhedron_Scan(D, C, NbMaxRays);
        D->next = next;
        return PolyhedronPtr(res);
    }, py::return_value_policy::reference);

    m.def("Constraints2Polyhedron", [](Matrix *m, unsigned flags) {
        return PolyhedronPtr(Constraints2Polyhedron(m, flags));
    }, py::arg("matrix"), py::arg("flags") = 0);

    m.def("PolyhedronImage", [](Matrix *m, Polyhedron *P, unsigned flags) {
        return PolyhedronPtr(Polyhedron_Image(P, m, flags));
    }, py::arg("matrix"), py::arg("polyhedron"), py::arg("flags") = 0);



    // ---------------------------- LBL ------------------------------
    py::class_<LBL, LBLPtr>(m, "LBL")
        .def_property_readonly("Lat",
            [](LBL &l) { return l.Lat; },
            py::return_value_policy::reference)
        .def_property_readonly("P",
            [](LBL &l) { return l.P; },
            py::return_value_policy::reference)
        .def_property_readonly("next",
            [](LBL &l) { return l.next; },
            py::return_value_policy::reference)
        .def("__repr__", [](LBL *l) {
            char *buf = nullptr;
            size_t size = 0;
            FILE *f = open_memstream(&buf, &size);
            LBLPrint(f, " %s", l);
            fclose(f);
            std::string result(buf, size);
            free(buf);
            return result;
        })
        .def("intersection", [](LBL *a, LBL *b) {
            return LBLPtr(LBLIntersection(a, b));
        })
        .def("difference", [](LBL *a, LBL *b) {
            return LBLPtr(LBLDifference(a, b));
        })
        .def("union", [](LBL *a, LBL *b) {
            return LBLPtr(LBLUnion(a, b));  
        })
        .def("included", [](LBL *a, LBL *b) {
            return (bool)LBLIncluded(a, b);
        })
        .def("zdomain", [](LBL *a) {
            return LBLPtr(LBL2ZDomain(a));
        });


    m.def("LBLAlloc", [](Matrix *lat, Polyhedron *domain) {
        return LBLPtr(LBLAlloc(lat, domain));
    }, py::arg("lat"), py::arg("domain"));

    m.def("LBLImage", [](LBL *a, Matrix *func) {
        return LBLPtr(LBLImage(a, func));
    }, py::arg("lbl"), py::arg("func"));

    m.def("LBLPreimage", [](LBL *a, Matrix *func) {
        return LBLPtr(LBLPreimage(a, func));
    }, py::arg("lbl"), py::arg("func"));

    // Displaying an LBL
    m.def("LBLPrint", [](LBL *l) {
        LBLPrint(stdout, " %s", l);
    });

 }