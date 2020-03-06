/*
 * pyheliosat2.cxx
 *
 *  Created on: 15 janv. 2020
 *      Author: benoit.gschwind
 */

#include "../solar_geometry.h"

#include <tuple>

#define PY_SSIZE_T_CLEAN
#include <Python.h>

#define NPY_NO_DEPRECATED_API NPY_1_7_API_VERSION
#include <numpy/arrayobject.h>
#include <numpy/npy_common.h>

// short cut for reinterpret_cast
template<typename T>
static inline T * _as(void * x) {
	return reinterpret_cast<T*>(x);
}

template<typename T>
static inline T & _getarg(PyArrayObject ** arr, int const i, int const arg)
{
	return reinterpret_cast<T*>(PyArray_DATA(arr[arg]))[i];
};

static PyObject * py_sunset(PyObject * self, PyObject * args)
{
	int const PHI = 0, DELTA = 1, MAXARGS = 2;

	PyObject * arg[MAXARGS] = {NULL};
	PyArrayObject * arr[MAXARGS] = {NULL};

	if (PyTuple_GET_SIZE(args) != MAXARGS)
		return NULL;

	for (int i = 0; i < MAXARGS; ++i) {
		arg[i] = PyTuple_GET_ITEM(args, i);
	}

	for (int i = 0; i < MAXARGS; ++i) {
		arr[i] = _as<PyArrayObject>(PyArray_FROM_OTF(arg[i], NPY_DOUBLE, NPY_ARRAY_IN_ARRAY));
		if (arr[i] == NULL)
			return NULL;
	}

	auto nd = PyArray_SIZE(arr[0]);

	for (int i = 1; i < MAXARGS; ++i) {
		if (PyArray_SIZE(arr[i]) != nd)
			return NULL;
	}

	npy_intp out_dims[] = {nd};
	PyArrayObject * out_arr = reinterpret_cast<PyArrayObject*>(PyArray_SimpleNew(1, out_dims, NPY_DOUBLE));

	for (int i = 0; i < nd; ++i) {
		*_as<int>(PyArray_GETPTR1(out_arr,i)) = sg1::sunset(
				_getarg<double>(arr, i, PHI),
				_getarg<double>(arr, i, DELTA)
				);
	}

	return _as<PyObject>(out_arr);

}

#define TPL_FUNCTION(name) {#name, py_##name, METH_VARARGS, "Not documented"}

static PyMethodDef methods[] =
{
	TPL_FUNCTION(sunset),
	{NULL, NULL, 0, NULL}
};

struct module_state {
	PyObject * Error;
};

static int solar_geometry_traverse(PyObject *m, visitproc visit, void *arg) {
	Py_VISIT(reinterpret_cast<module_state*>(PyModule_GetState(m))->Error);
	return 0;
}

static int solar_geometry_clear(PyObject *m) {
	Py_CLEAR(reinterpret_cast<module_state*>(PyModule_GetState(m))->Error);
	return 0;
}


static struct PyModuleDef moduledef = {
	PyModuleDef_HEAD_INIT,
	"solar_geometry",
	NULL,
	sizeof(module_state),
	methods,
	NULL,
	solar_geometry_traverse,
	solar_geometry_clear,
	NULL
};

PyMODINIT_FUNC
PyInit_solar_geometry(void)
{

	PyObject * m = PyModule_Create(&moduledef);
	if (m == NULL)
		return NULL;

	module_state * state = reinterpret_cast<module_state*>(PyModule_GetState(m));
	state->Error = PyErr_NewException("solar_geometry.error", NULL, NULL);
	Py_INCREF(state->Error);
	PyModule_AddObject(m, "error", state->Error);

	/** init numpy **/
	import_array();

	return m;

}





