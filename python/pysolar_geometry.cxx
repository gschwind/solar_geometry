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

template<typename>
struct _python_bind_type_info;

template<>
struct _python_bind_type_info<double> {
	enum : int { npy_type = NPY_DOUBLE };
};

template<>
struct _python_bind_type_info<float> {
	enum : int { npy_type = NPY_FLOAT };
};

template<>
struct _python_bind_type_info<int> {
	enum : int { npy_type = NPY_INT };
};

template<>
struct _python_bind_type_info<char> {
	enum : int { npy_type = NPY_INT8 };
};


template<typename T>
struct handle_numpy_1d_array {
	PyArrayObject * _arr;
	handle_numpy_1d_array(PyObject * arr) : _arr{nullptr}
	{
		_arr = _as<PyArrayObject>(PyArray_FROM_OTF(arr, _python_bind_type_info<T>::npy_type, NPY_ARRAY_IN_ARRAY));
	}

	T const & operator[](int i) const
	{
		return reinterpret_cast<T*>(PyArray_DATA(_arr))[i];
	};

	T & operator[](int i)
	{
		return reinterpret_cast<T*>(PyArray_DATA(_arr))[i];
	};

	int size() const
	{
		return PyArray_SIZE(_arr);
	}

};

template<typename>
struct _function_signature;

template<typename R, typename ... ARGS>
struct _function_signature<R(ARGS...)> {
	using return_type = R;
};

// T0, TN expected to be handle_numpy_1d_array<X>
template<typename F, typename T0, typename ... TN>
static PyObject * _final_vectorized_call (F * func, T0 && args0,  TN && ... args)
{
	using return_type = typename _function_signature<F>::return_type;

	int nd = args0.size();
	npy_intp out_dims[] = {nd};
	PyArrayObject * out_arr = reinterpret_cast<PyArrayObject*>(PyArray_SimpleNew(1, out_dims,
			_python_bind_type_info<return_type>::npy_type));

	for (int i = 0; i < nd; ++i) {
		*_as<return_type>(PyArray_GETPTR1(out_arr,i)) = func(args0[i], args[i]...);
	}

	return _as<PyObject>(out_arr);
}

template<typename, int, typename ...>
struct _build_vectorized_function_call_pass1;

template<typename R, typename ... ARGS, int I, typename HEAD, typename ... TAIL>
struct _build_vectorized_function_call_pass1<R(ARGS...), I, HEAD, TAIL...>
{
	using func_type = R(ARGS...);
	template<typename ...XARGS>
	static PyObject * call (func_type * func, PyObject * args, XARGS && ... xargs)
	{
		return _build_vectorized_function_call_pass1<R(ARGS...), I+1, TAIL...>::call(func, args,
				handle_numpy_1d_array<HEAD>(PyTuple_GET_ITEM(args, I)), xargs...);
	}
};

template<typename F, int I>
struct _build_vectorized_function_call_pass1<F, I>
{
	template<typename ...ARGS>
	static PyObject * call (F * func, PyObject * args, ARGS && ... xargs)
	{
		/* check for length of argument */
		if (PyTuple_GET_SIZE(args) != I)
			return NULL;
		return _final_vectorized_call(func, xargs...);
	}
};

template<typename>
struct _build_vectorized_function_call_pass0;

template<typename R, typename ... ARGS>
struct _build_vectorized_function_call_pass0<R(ARGS...)>
{
	using func_type = R(ARGS...);
	static PyObject * call (func_type * func, PyObject * self, PyObject * args)
	{
		return _build_vectorized_function_call_pass1<R(ARGS...), 0, ARGS...>::call(func, args);
	}
};


/**
 * This function use fonction argment type to build the corresponging vectorized function.
 **/
template<typename F>
PyObject * call_vectorized_function_with_python_args(F * func, PyObject * self, PyObject * args)
{
	return _build_vectorized_function_call_pass0<F>::call(func, self, args);
}

static PyObject * py_sunset(PyObject * self, PyObject * args)
{
	return call_vectorized_function_with_python_args(sg1::sunset, self, args);
}

static PyObject * py_gamma_sun(PyObject * self, PyObject * args)
{
	return call_vectorized_function_with_python_args(sg1::gamma_sun, self, args);
}

static PyObject * py_day_angle(PyObject * self, PyObject * args)
{
	return call_vectorized_function_with_python_args(sg1::day_angle, self, args);
}

static PyObject * py_corr_distance(PyObject * self, PyObject * args)
{
	return call_vectorized_function_with_python_args(sg1::corr_distance, self, args);
}

static PyObject * py_nbdays_month(PyObject * self, PyObject * args)
{
	return call_vectorized_function_with_python_args(sg1::nbdays_month, self, args);
}

static PyObject * py_declination_sun(PyObject * self, PyObject * args)
{
	return call_vectorized_function_with_python_args(sg1::declination_sun, self, args);
}

static PyObject * py_solar_hour_angle(PyObject * self, PyObject * args)
{
	return call_vectorized_function_with_python_args(sg1::solar_hour_angle, self, args);
}

static PyObject * py_omega_to_LAT(PyObject * self, PyObject * args)
{
	return call_vectorized_function_with_python_args(sg1::omega_to_LAT, self, args);
}

#define TPL_FUNCTION(name) {#name, py_##name, METH_VARARGS, "Not documented"}

static PyMethodDef methods[] =
{
	TPL_FUNCTION(sunset),
	TPL_FUNCTION(gamma_sun),
	TPL_FUNCTION(day_angle),
	TPL_FUNCTION(corr_distance),
	TPL_FUNCTION(nbdays_month),
	TPL_FUNCTION(declination_sun),
	TPL_FUNCTION(solar_hour_angle),
	TPL_FUNCTION(omega_to_LAT),
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





