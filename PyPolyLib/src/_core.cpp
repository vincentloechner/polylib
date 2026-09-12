/*
 * This file is the pypolylib core python interface to the C polylib.
 *
 * copyright 2026 Vincent Loechner
 */

#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>

#define MAX_RAYS 100

#include <gmpxx.h>   // MUST come before polylibgmp.h
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


// just a view for the polyhedra Constraints and Rays accessors:
struct Poly_Matrix_View {
    Polyhedron *p;
    Value **value_ptr;
};
py::object get_value_from_matrix(Value **mat_ptr, int i, int j) {
    char *s = mpz_get_str(nullptr, 10, mat_ptr[i][j]);
    PyObject *obj = PyLong_FromString(s, nullptr, 10);
    free(s);
    return py::reinterpret_steal<py::object>(obj);
}
// // build buffer matrix and access it
// py::array_t<py::object> make_buffer_matrix(int rows, int cols)
// {
//     py::array_t<py::object> M({rows, cols});
//     auto b = M.mutable_unchecked<2>();
//     for (ssize_t i = 0; i < rows; ++i)
//         for (ssize_t j = 0; j < cols; ++j)
//             b(i, j) = py::none();

//     return M;
// }
// py::object get_value_from_buffer(Value **poly_ptr, py::array_t<py::object>& buf,
//                                  int i, int j) {
//     auto b = buf.mutable_unchecked<2>();
//     py::object value = b(i, j);
//     if (value.is_none()) {
//         char *s = mpz_get_str(nullptr, 10, poly_ptr[i][j]);
//         PyObject *obj = PyLong_FromString(s, nullptr, 10);
//         free(s);
//         value = py::reinterpret_steal<py::object>(obj);
//         b(i, j) = value;
//     }

//     return value;
// }


PYBIND11_MODULE(_core, m) {
    // ----------------------- Memory management ----------------------
    m.def("PolylibClose", []() {
        polylib_close();
    });


    // ----------------------- Matrix ----------------------
    py::class_<Matrix, MatrixPtr>(m, "Matrix")
        // Matrix properties
        .def_readonly("nbrows",    &Matrix::NbRows)
        .def_readonly("nbcolumns", &Matrix::NbColumns)
        // --- accessors to gmp values ---
        .def("__getitem__", [](Matrix &self, py::tuple idx) -> py::object {
            int i = idx[0].cast<int>();
            int j = idx[1].cast<int>();
            return get_value_from_matrix(self.p, i, j);
        }, py::arg("idx:tuple"))
        .def("__setitem__", [](Matrix &self, py::tuple idx, py::object val) {
            int i = idx[0].cast<int>();
            int j = idx[1].cast<int>();

            std::string s = py::str(val);
            mpz_set_str(self.p[i][j], s.c_str(), 10);
        }, py::arg("idx:tuple"), py::arg("value"))
        .def("print", [](Matrix &self) {
            Matrix_Print(stdout, " %s", &self);
        })
        // --- methods on Matrices ---
        // Inverse
        .def("inverse", [](Matrix *self) {
            // do a copy, Matrix_Inverse writes to its arg matrix
            Matrix *mat = Matrix_Copy(self);
            Matrix *result = Matrix_Alloc(mat->NbRows, mat->NbColumns);
            int ok = Matrix_Inverse(mat, result);
            Matrix_Free(mat);
            if (!ok) {
                Matrix_Free(result);
                throw std::domain_error("The matrix is not invertible");
            }
            return MatrixPtr(result);
        })
        // Product
        .def("multiply", [](Matrix *self, Matrix *b) {
            Matrix *result = Matrix_Alloc(self->NbRows, b->NbColumns);
            Matrix_Product(self, b, result);
            return MatrixPtr(result);
        }, py::arg("b"))
        ;

    // --- Matrix creation ---
    m.def("matrix_alloc", [](unsigned nbrows, unsigned nbcols) {
        return MatrixPtr(Matrix_Alloc(nbrows, nbcols));
    }, py::arg("nbrows"), py::arg("nbcols"));

    m.def("matrix_read_from_string", [](const std::string &s) {
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
    }, py::arg("string"));


    // --------------------- Polyhedron ---------------------------
    // Specific class to get read access to Constraint and to Ray matrices:
    py::class_<Poly_Matrix_View>(m, "Poly_Matrix_View")
        .def("__getitem__", [](Poly_Matrix_View &p, py::tuple idx) {
            int i = idx[0].cast<int>();
            int j = idx[1].cast<int>();
            return get_value_from_matrix(p.value_ptr, i, j);
        });
    // Main Polyhedron class: polyhedral domain
    py::class_<Polyhedron, PolyhedronPtr>(m, "Polyhedron")
        // Polyhedron properties
        .def_readonly("dimension",      &Polyhedron::Dimension)
        .def_readonly("nbconstraints",  &Polyhedron::NbConstraints)
        .def_readonly("nbrays",         &Polyhedron::NbRays)
        .def_readonly("nbbid",          &Polyhedron::NbBid)
        .def_property_readonly("next",
            [](Polyhedron &p) -> Polyhedron* { return p.next; },
            py::return_value_policy::reference)

        // constraint and ray accessors with a buffer
        .def_property_readonly("constraint", [](Polyhedron &p) {
            return Poly_Matrix_View{
                &p,
                p.Constraint
            };
        }, py::return_value_policy::move)
        .def_property_readonly("ray", [](Polyhedron &p) {
            return Poly_Matrix_View{
                &p,
                p.Ray
            };
        }, py::return_value_policy::move)

        // polylib operators:
        .def("print", [](Polyhedron &self) {
            Polyhedron_Print(stdout, " %s", &self);
        })
        .def("copy_single_pol", [](Polyhedron &self) {
            return PolyhedronPtr(Polyhedron_Copy(&self));
        })
        .def("copy", [](Polyhedron &self) {
            return PolyhedronPtr(Domain_Copy(&self));
        })
        .def("image", [](Polyhedron &self, Matrix &m) {
            return PolyhedronPtr(Polyhedron_Image(&self, &m, MAX_RAYS));
        }, py::arg("matrix"))
        .def("preimage", [](Polyhedron &self, Matrix &m) {
            return PolyhedronPtr(Polyhedron_Preimage(&self, &m, MAX_RAYS));
        }, py::arg("matrix"))
        .def("add_constraints", [](Polyhedron &self, Matrix &m) {
            return PolyhedronPtr(AddConstraints(m.p[0], m.NbRows, &self, MAX_RAYS));
        }, py::arg("matrix"))
        .def("scan", [](Polyhedron &self, Polyhedron &C) {
            // single polyhedron scan only, unset/restore next:
            Polyhedron *next = self.next;
            self.next = NULL;
            Polyhedron *res = Polyhedron_Scan(&self, &C, MAX_RAYS);
            self.next = next;
            return PolyhedronPtr(res);
        }, py::arg("context_polyhedron"))
        .def("is_bounded", [](Polyhedron &self) {
            if(self.NbBid != 0)
                return false;
            for(unsigned r = 0; r < self.NbRays; r++)
                if(value_zero_p(self.Ray[r][self.Dimension + 1]))
                    return false;
            return true;
        })
        .def("intersect", [](Polyhedron *self, Polyhedron *P2) {
            return PolyhedronPtr(AddConstraints(self->Constraint[0], self->NbConstraints, P2, MAX_RAYS));
        }, py::arg("polyhedron"));
        ;

    // --- Polyhedron creation ---
    m.def("constraints2polyhedron", [](Matrix *m) {
        return PolyhedronPtr(Constraints2Polyhedron(m, MAX_RAYS));
    }, py::arg("constraint matrix"));


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

        // LBL operations
        .def("intersection", [](LBL *self, LBL *b) {
            return LBLPtr(LBLIntersection(self, b));
        }, py::arg("LBL"))
        .def("difference", [](LBL *self, LBL *b) {
            return LBLPtr(LBLDifference(self, b));
        }, py::arg("LBL"))
        .def("union", [](LBL *self, LBL *b) {
            return LBLPtr(LBLUnion(self, b));
        }, py::arg("LBL"))
        .def("contains_point", [](LBL *self, py::sequence point) {
            Matrix *mat = Matrix_Alloc(1, py::len(point));
            for (size_t i = 0; i < py::len(point); i++) {
                std::string s = py::str(point[i]);
                mpz_set_str(mat->p[0][i], s.c_str(), 10);
            }
            bool res = LBLContainsPoint(self, mat->p[0]);
            Matrix_Free(mat);
            return res;
        }, py::arg("sequence"))
        .def("included", [](LBL *self, LBL *b) {
            return (bool)LBLIncluded(self, b);
        }, py::arg("LBL"))
        .def("z_domain", [](LBL *self) {
            return LBLPtr(LBL2ZDomain(self));
        })
        .def("disjoint", [](LBL *self) {
            return LBLPtr(LBLDisjointUnion(self));
        })
        .def("image", [](LBL *self, Matrix *func) {
            return LBLPtr(LBLImage(self, func));
        }, py::arg("matrix"))
        .def("preimage", [](LBL *self, Matrix *func) {
            return LBLPtr(LBLPreimage(self, func));
        }, py::arg("matrix"))
        ;


    // --- LBL creation ---
    m.def("LBLAlloc", [](Matrix *lat, Polyhedron *domain) {
        return LBLPtr(LBLAlloc(lat, domain));
    }, py::arg("lat"), py::arg("domain"));

}

