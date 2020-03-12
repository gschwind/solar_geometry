/*
 * pyheliosat2.cxx
 *
 *  Created on: 15 janv. 2020
 *      Author: benoit.gschwind
 */

// Helper to easy the python binding
#include "python-bind-helper.hxx"

#include "solar_geometry.h"

#include <tuple>

#define PY_SSIZE_T_CLEAN
#include <Python.h>

#define NPY_NO_DEPRECATED_API NPY_1_7_API_VERSION
#include <numpy/arrayobject.h>
#include <numpy/npy_common.h>
#include <numpy/ufuncobject.h>

#include <iostream>
#include <string>



#define make_binding(name) \
static python_bind_helper::build_ufunc<decltype(sg1::name), sg1::name> ufunc_##name(#name);

make_binding(sunset)
make_binding(gamma_sun)
make_binding(day_angle)
make_binding(corr_distance)
make_binding(nbdays_month)
make_binding(declination_sun)
make_binding(solar_hour_angle)
make_binding(omega_to_LAT)
make_binding(ymd_to_day_of_year)
make_binding(geogr_to_geoce)
make_binding(azimuth_sun)
make_binding(julian_day_to_ymd)
make_binding(day_of_year_to_ymd)


static PyMethodDef methods[] =
{
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

	/** init numpy **/
	import_array();
	import_umath();

	module_state * state = reinterpret_cast<module_state*>(PyModule_GetState(m));
	state->Error = PyErr_NewException("solar_geometry.error", NULL, NULL);
	Py_INCREF(state->Error);
	PyModule_AddObject(m, "error", state->Error);

#define register_ufunc(name) \
ufunc_##name.register_to(m)

	register_ufunc(sunset);
	register_ufunc(gamma_sun);
	register_ufunc(day_angle);
	register_ufunc(corr_distance);
	register_ufunc(nbdays_month);
	register_ufunc(declination_sun);
	register_ufunc(solar_hour_angle);
	register_ufunc(omega_to_LAT);
	register_ufunc(ymd_to_day_of_year);
	register_ufunc(geogr_to_geoce);
	register_ufunc(azimuth_sun);
	register_ufunc(julian_day_to_ymd);
	register_ufunc(day_of_year_to_ymd);

	return m;

}





