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

//#include <iostream>
#include <string>

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

// convenient function for debuging.
static std::string _python_repr(PyObject *obj)
{
    PyObject* repr = PyObject_Repr(obj);
    PyObject* str = PyUnicode_AsEncodedString(repr, "utf-8", "strict");
    auto ret = std::string{PyBytes_AS_STRING(str)};
    Py_XDECREF(repr);
    Py_XDECREF(str);
    return ret;
}

// convenient function for debuging.
static std::string _python_str(PyObject *obj)
{
    PyObject* repr = PyObject_Str(obj);
    PyObject* str = PyUnicode_AsEncodedString(repr, "utf-8", "strict");
    auto ret = std::string{PyBytes_AS_STRING(str)};
    Py_XDECREF(repr);
    Py_XDECREF(str);
    return ret;
}

template<typename>
struct _python_bind_type_info;

template<>
struct _python_bind_type_info<double> {
	enum : int { npy_type = NPY_DOUBLE };

	static PyObject * format() {
		return Py_BuildValue("s", "f8");
	}

};

template<>
struct _python_bind_type_info<float> {
	enum : int { npy_type = NPY_FLOAT };
	static PyObject * format() {
		return Py_BuildValue("s", "f8");
	}
};

template<typename T>
static PyObject * _int_format()
{
	switch(sizeof(T)) {
	case 1:
		return Py_BuildValue("s", "i1");
	case 2:
		return Py_BuildValue("s", "i2");
	case 4:
		return Py_BuildValue("s", "i4");
	case 8:
		return Py_BuildValue("s", "i8");

	}
}

template<>
struct _python_bind_type_info<int> {
	enum : int { npy_type = NPY_INT };
	static PyObject * format() {
		return _int_format<int>();
	}
};

template<>
struct _python_bind_type_info<char> {
	enum : int { npy_type = NPY_INT8 };
	static PyObject * format() {
		return _int_format<char>();
	}
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

template<typename T, int I>
inline ptrdiff_t get_field_offset() {
	return reinterpret_cast<char*>(&std::get<I>(*reinterpret_cast<T*>(0)))-reinterpret_cast<char*>(0);
}

template<typename, int, typename...>
struct _create_desc_for_tuple;

template<typename T, typename HEAD, int I>
void _update_items(PyObject * names, PyObject * formats, PyObject * offsets)
{
	PyList_SetItem(formats, I, _python_bind_type_info<HEAD>::format());
	PyList_SetItem(offsets, I, PyLong_FromLong(get_field_offset<T, I>()));
}

template<typename T, int I, typename HEAD, typename ... TAIL>
struct _create_desc_for_tuple<T, I, HEAD, TAIL...> {
	static void update(PyObject * names, PyObject * formats, PyObject * offsets)
	{
		_update_items<T, HEAD, I>(names, formats, offsets);
		_create_desc_for_tuple<T, I+1, TAIL...>::update(names, formats, offsets);

	}
};

template<typename T, int I, typename HEAD>
struct _create_desc_for_tuple<T, I, HEAD> {
	static void update(PyObject * names, PyObject * formats, PyObject * offsets)
	{
		_update_items<T, HEAD, I>(names, formats, offsets);
	}
};

template<typename>
struct _create_desc_for_tuple_pass0;

template<typename ... ARGS>
struct _create_desc_for_tuple_pass0<std::tuple<ARGS...>> {
	static void update(PyObject * names, PyObject * formats, PyObject * offsets) {
		_create_desc_for_tuple<std::tuple<ARGS...>, 0, ARGS...>::update(names, formats, offsets);
	}
};



template<typename T>
PyArray_Descr * create_desc(std::array<std::string, std::tuple_size<T>::value> const & field_names) {
	static auto const tsize = std::tuple_size<T>::value;
	auto spec = PyDict_New();
	auto names = PyList_New(tsize);
	auto formats = PyList_New(tsize);
	auto offsets = PyList_New(tsize);

	PyDict_SetItemString(spec, "names", names);
	PyDict_SetItemString(spec, "formats", formats);
	PyDict_SetItemString(spec, "offsets", offsets);

	_create_desc_for_tuple_pass0<T>::update(names, formats, offsets);

	/* must be handled by spec */
	Py_DECREF(names);
	Py_DECREF(formats);
	Py_DECREF(offsets);

	for (int i = 0; i < std::tuple_size<T>::value; ++i) {
			PyList_SetItem(names, i, Py_BuildValue("s", field_names[i].c_str()));
	}

	PyArray_Descr * final_desc;
	if (PyArray_DescrConverter(spec, &final_desc) != NPY_SUCCEED)
		throw std::runtime_error("SHOULD NEVER APPEND");

	Py_DECREF(spec);
	return final_desc;

}


template<typename, typename...>
struct _final_vectorized;


// T0, TN expected to be handle_numpy_1d_array<X>
template<typename ... TARGS, typename ... ARGS, typename T0, typename ... TN>
struct _final_vectorized<std::tuple<TARGS...>(ARGS...), T0, TN...> {
	using F = typename std::tuple<TARGS...>(ARGS...);
	using R = typename std::tuple<TARGS...>;
	static PyObject * call (F * func, std::array<std::string, sizeof...(TARGS)> const & field_names, T0 && args0,  TN && ... args)
	{
		int nd = args0.size();

		/* create the output */
		PyObject * buffer = PyByteArray_FromStringAndSize(nullptr, nd*sizeof(R));
		auto output_array_data = reinterpret_cast<R*>(PyByteArray_AsString(buffer));

		for (int i = 0; i < nd; ++i) {
			new (&output_array_data[i]) R; // call tuple inplace constructor
			output_array_data[i] = func(args0[i], args[i]...);
		}

		npy_intp out_dims[] = {nd};
		npy_intp strides[] = {sizeof(R)};

		auto desc = create_desc<R>(field_names);

		PyArrayObject * out_arr = reinterpret_cast<PyArrayObject*>(PyArray_NewFromDescr(&PyArray_Type, desc, 1, out_dims, strides, output_array_data, 0, nullptr));

		return _as<PyObject>(out_arr);
	}
};


// T0, TN expected to be handle_numpy_1d_array<X>
template<typename F, typename T0, typename ... TN>
struct _final_vectorized<F, T0, TN...> {
	static PyObject * call (F * func, T0 && args0,  TN && ... args)
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
};

template<typename F, typename ... TN>
static PyObject * _final_vectorized_call(F * func, TN && ... args)
{
	return _final_vectorized<F, TN...>::call(func, args...);
}

template<typename F, typename ... TN>
static PyObject * _final_vectorized_call(F * func, std::array<std::string, std::tuple_size<typename _function_signature<F>::return_type>::value> const & field_names, TN && ... args)
{
	return _final_vectorized<F, TN...>::call(func, field_names, args...);
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
		return _build_vectorized_function_call_pass1<func_type, I+1, TAIL...>::call(func, args,
				handle_numpy_1d_array<HEAD>(PyTuple_GET_ITEM(args, I)), xargs...);
	}
};

template<typename ...TARGS, typename ... ARGS, int I, typename HEAD, typename ... TAIL>
struct _build_vectorized_function_call_pass1<std::tuple<TARGS...>(ARGS...), I, HEAD, TAIL...>
{
	using func_type = std::tuple<TARGS...>(ARGS...);

	template<typename ...XARGS>
	static PyObject * call (func_type * func, std::array<std::string, sizeof...(TARGS)> const & field_names, PyObject * args, XARGS && ... xargs)
	{
		return _build_vectorized_function_call_pass1<func_type, I+1, TAIL...>::call(func, field_names, args,
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

template<typename ...TARGS, typename ... ARGS, int I>
struct _build_vectorized_function_call_pass1<std::tuple<TARGS...>(ARGS...), I>
{
	using func_type = std::tuple<TARGS...>(ARGS...);

	template<typename ... XARGS>
	static PyObject * call (func_type * func, std::array<std::string, sizeof...(TARGS)> const & field_names, PyObject * args, XARGS && ... xargs)
	{
		/* check for length of argument */
		if (PyTuple_GET_SIZE(args) != I)
			return NULL;
		return _final_vectorized_call(func, field_names, xargs...);
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
		return _build_vectorized_function_call_pass1<func_type, 0, ARGS...>::call(func, args);
	}
};

/* If the function to binf return tuple */
template<typename ... TARGS, typename ... ARGS>
struct _build_vectorized_function_call_pass0<std::tuple<TARGS...>(ARGS...)>
{
	using func_type = std::tuple<TARGS...>(ARGS...);
	static PyObject * call (func_type * func, std::array<std::string, sizeof...(TARGS)> const & field_names, PyObject * self, PyObject * args)
	{
		return _build_vectorized_function_call_pass1<func_type, 0, ARGS...>::call(func, field_names, args);
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

template<typename F>
PyObject * call_vectorized_function_with_python_args(F * func, std::array<std::string, std::tuple_size<typename _function_signature<F>::return_type>::value> const & field_names, PyObject * self, PyObject * args)
{
	return _build_vectorized_function_call_pass0<F>::call(func, field_names, self, args);
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

static PyObject * py_julian_day_to_ymd(PyObject * self, PyObject * args)
{
	return call_vectorized_function_with_python_args(sg1::julian_day_to_ymd, {"year", "month", "day_of_month"}, self, args);
}

static PyObject * py_day_of_year_to_ymd(PyObject * self, PyObject * args)
{
	return call_vectorized_function_with_python_args(sg1::day_of_year_to_ymd, {"month", "day_of_month"}, self, args);
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
	TPL_FUNCTION(julian_day_to_ymd),
	TPL_FUNCTION(day_of_year_to_ymd),
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





