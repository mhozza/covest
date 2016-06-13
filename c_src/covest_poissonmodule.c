#include <Python.h>
#include <math.h>
#include <stdio.h>

#define MAX_EXP 200

static PyObject *truncated_poisson( PyObject *self, PyObject *args )
{
    double l;
    int j;
    if (!PyArg_ParseTuple(args, "di:truncated_poisson", &l, &j)) {
        return NULL;
    }

    if (l == 0 || l != l) {
        return Py_BuildValue("d", 0);
    }

    long double p1 = 1;
    long double p3 = l;

    for (int i = 1; i<=j; i++) {
        p1 *= l / i;
    }
    while (l > MAX_EXP && p1 > 0) {
        p1 /= expl(MAX_EXP);
        l -= MAX_EXP;
    }
    if (l > 1e-8 && p1 > 0) {
        p3 = expl(l) - 1;
    }
    double res = (double)(p1 / p3);

    return PyFloat_FromDouble(res);
}


static PyObject *poisson( PyObject *self, PyObject *args )
{
    double l;
    int j;
    if (!PyArg_ParseTuple(args, "di:poisson", &l, &j)) {
        return NULL;
    }

    if (l == 0 || l != l) {
        return Py_BuildValue("d", 0);
    }

    long double p1 = 1;

    for (int i = 1; i<=j; i++) {
        p1 *= l / i;
    }
    while (l > MAX_EXP && p1 > 0) {
        p1 /= expl(MAX_EXP);
        l -= MAX_EXP;
    }
    double res = (double)(p1/expl(l));
    return PyFloat_FromDouble(res);
}


static PyObject *poisson_dist( PyObject *self, PyObject *args )
{
    double l;
    int max_j;
    if (!PyArg_ParseTuple(args, "di:poisson_dist", &l, &max_j)) {
        return NULL;
    }

    PyObject *res = PyList_New(max_j);
    if (!res) { return NULL; }

    if (l == 0 || l != l) {
        PyObject *value = PyFloat_FromDouble(0);
        if (!value) {
            Py_DECREF(res);
            return NULL;
        }
        for (int j = 1; j <= max_j; j++) {
            Py_INCREF(value);
            PyList_SetItem(res, j - 1, value);
        }
        Py_DECREF(value);
        return res;
    }

    long double p1 = 1;
    long double d1 = expl(MAX_EXP);
    long double d2 = expl(l);
    for (int j = 1; j <= max_j; j++) {
        p1 *= l / j;
        long double p1c = p1;
        while (l > MAX_EXP && p1c > 0) {
            p1c /= d1;
            l -= MAX_EXP;
        }
        double v = (double)(p1c/d2);
        PyObject *value = PyFloat_FromDouble(v);
        if (!value) {
            Py_DECREF(res);
            return NULL;
        }
        PyList_SetItem(res, j - 1, value);
    }
    return res;
}


static PyMethodDef PoissonMethods[] = {
    {"truncated_poisson",  truncated_poisson, METH_VARARGS,
     "Compute truncated poisson pmf value."},
    {"poisson",  poisson, METH_VARARGS,
     "Compute poisson pmf value."},
    {"poisson_dist",  poisson_dist, METH_VARARGS,
     "Compute poisson pmf value for whole distribution."},
    {NULL, NULL, 0, NULL}        /* Sentinel */
};

static struct PyModuleDef covest_poissonmodule = {
   PyModuleDef_HEAD_INIT,
   "covest_poisson",   /* name of module */
   NULL, /* module documentation, may be NULL */
   -1,       /* size of per-interpreter state of the module,
                or -1 if the module keeps state in global variables. */
   PoissonMethods
};

PyMODINIT_FUNC
PyInit_covest_poisson(void)
{
    return PyModule_Create(&covest_poissonmodule);
}
