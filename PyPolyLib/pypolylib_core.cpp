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

// just another name for the polyhedra Constraints and Rays accessors:
struct ConstraintsView {
    Polyhedron *poly;
};
struct RaysView {
    Polyhedron *poly;
};

PYBIND11_MODULE(pypolylib_core, m) {
    // ----------------------- Memory management ----------------------
    m.def("PolylibClose", []() {
        polylib_close();
    });


    // ----------------------- Matrix ----------------------
    py::class_<Matrix, MatrixPtr>(m, "Matrix")
        // Matrix properties
        .def_readonly("nbrows",    &Matrix::NbRows)
        .def_readonly("nbcolumns", &Matrix::NbColumns)
        // .def_readonly("value",     &Matrix::p)
        // --- accessors to gmp values ---
        .def("__getitem__", [](Matrix &m, py::tuple idx) -> py::object {
            int i = idx[0].cast<int>();
            int j = idx[1].cast<int>();
    
            char *s = mpz_get_str(nullptr, 10, m.p[i][j]);
            PyObject *obj = PyLong_FromString(s, nullptr, 10);
            free(s);
            return py::reinterpret_steal<py::int_>(obj);
            // py::object val = py::module_::import("builtins").attr("int")(s);
            // free(s);
            // return val;
        })
        .def("__setitem__", [](Matrix &m, py::tuple idx, py::object val) {
            int i = idx[0].cast<int>();
            int j = idx[1].cast<int>();
    
            std::string s = py::str(val);
            mpz_set_str(m.p[i][j], s.c_str(), 10);
        });

    // --- methods on Matrices ---
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
    // Product
    m.def("MatrixProduct", [](Matrix *a, Matrix *b) {
        Matrix *result = Matrix_Alloc(a->NbRows, b->NbColumns);
        Matrix_Product(a, b, result);
        return MatrixPtr(result);
    }, py::arg("a"), py::arg("b"));

    // Inverse
    m.def("MatrixInverse", [](Matrix *mat) {
        // do a copy, Matrix_Inverse writes to its arg matrix
        mat = Matrix_Copy(mat);
        Matrix *result = Matrix_Alloc(mat->NbRows, mat->NbColumns);
        int ok = Matrix_Inverse(mat, result);
        Matrix_Free(mat);
        if (!ok) {
            Matrix_Free(result);
            throw std::runtime_error("The matrix is not invertible");
        }
        return MatrixPtr(result);
    }, py::arg("mat"));


    // --------------------- Polyhedron ---------------------------
    py::class_<ConstraintsView>(m, "ConstraintsView")
        .def("__getitem__", [](ConstraintsView &p, py::tuple idx) {
            int i = idx[0].cast<int>();
            int j = idx[1].cast<int>();

            char *s = mpz_get_str(nullptr, 10, p.poly->Constraint[i][j]);
            // py::object v = py::module_::import("builtins").attr("int")(s);
            PyObject *obj = PyLong_FromString(s, nullptr, 10);
            free(s);

            // return v;
            return py::reinterpret_steal<py::int_>(obj);
        });
    py::class_<RaysView>(m, "RaysView")
        .def("__getitem__", [](RaysView &p, py::tuple idx) {
            int i = idx[0].cast<int>();
            int j = idx[1].cast<int>();

            char *s = mpz_get_str(nullptr, 10, p.poly->Ray[i][j]);
            // py::object v = py::module_::import("builtins").attr("int")(s);
            PyObject *obj = PyLong_FromString(s, nullptr, 10);
            free(s);

            // return v;
            return py::reinterpret_steal<py::int_>(obj);
        });
    py::class_<Polyhedron, PolyhedronPtr>(m, "Polyhedron")
        // Polyhedron properties
        .def_readonly("dimension",      &Polyhedron::Dimension)
        .def_readonly("nbconstraints",  &Polyhedron::NbConstraints)
        .def_readonly("nbrays",         &Polyhedron::NbRays)
        .def_readonly("nbbid",          &Polyhedron::NbBid)
        .def_property_readonly("next",
            [](Polyhedron &p) -> Polyhedron* { return p.next; },
            py::return_value_policy::reference)
        // constraint and ray accessors
        .def_property_readonly("constraint", [](Polyhedron &p) {
            return ConstraintsView{&p};   // return a constraints view of the polyhedron
        }, py::return_value_policy::move)
        .def_property_readonly("ray", [](Polyhedron &p) {
            return RaysView{&p};   // return a rays view of the polyhedron
        }, py::return_value_policy::move);

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

    m.def("PolyhedronPrint", [](Polyhedron *pol) {
        Polyhedron_Print(stdout, " %s", pol);
    }), py::arg("pol");
    
    m.def("isBoundedPolyhedron", [](Polyhedron *pol) {
        if(pol->NbBid != 0)
            return false;
        for(int r = 0; r < pol->NbRays; r++)
            if(value_zero_p(pol->Ray[r][pol->Dimension + 1]))
                return false;
        return true;
    }), py::arg("pol");



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