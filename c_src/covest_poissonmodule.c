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

    return Py_BuildValue("d", res);
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
    return Py_BuildValue("d", res);
}


static PyMethodDef PoissonMethods[] = {
    {"truncated_poisson",  truncated_poisson, METH_VARARGS,
     "Compute truncated poisson pmf value."},
    {"poisson",  poisson, METH_VARARGS,
     "Compute poisson pmf value."},
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
