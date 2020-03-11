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
#include <numpy/ufuncobject.h>

#include <iostream>
#include <string>

template<typename ... ARGS>
void fold(ARGS && ... args) { }


template <std::size_t ...>
struct index_sequence { };

template <std::size_t N, std::size_t ... TAIL>
struct _index_sequence : public _index_sequence<N-1, N-1, TAIL...> { };

template <std::size_t ... TAIL>
struct _index_sequence<0U, TAIL ... > {
	using type = index_sequence<TAIL ... >;
};

template <std::size_t N>
using make_index_sequence = typename _index_sequence<N>::type;


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

template<typename F, F &FUNC>
struct _my_build_ufunc;

// If return type is a tuple
template<typename O_ARGS, typename ... I_ARGS, O_ARGS(&FUNC)(I_ARGS...)>
struct _my_build_ufunc<O_ARGS(I_ARGS...), FUNC>
{
	char types[sizeof...(I_ARGS) + 1]; // handle function signature
	PyUFuncGenericFunction func[1]; // handle vectorized function
	void * data[1]; // handle extra data (not used by us currently

	std::string name;
	std::string doc;

	enum : int { ISIZE = sizeof...(I_ARGS) };

	using ISEQ_TYPE = make_index_sequence<ISIZE>;

	template<typename T>
	static inline T & _assign(T & dst, T && src) { return dst = src; }

	template<typename>
	struct _update_types;

	template<std::size_t ... ISEQ>
	struct _update_types<index_sequence<ISEQ...>>
	{
		static void update(char * types)
		{
			fold(_assign<char>(types[ISEQ], _python_bind_type_info<I_ARGS>::npy_type)...);
			fold(_assign<char>(types[ISIZE], _python_bind_type_info<O_ARGS>::npy_type));
		}
	};

	_my_build_ufunc(std::string const & name, std::string const & doc = "") :
		data{nullptr}, func{&ufunc}, name{name}, doc{doc}
	{
		_update_types<ISEQ_TYPE>::update(types);
	}

	template<typename T>
	struct data_handler {
		char * const   _base;
		npy_intp const _step;

		data_handler(char * base, npy_intp step) : _base{base}, _step{step} { }

		T & operator[](int i)
		{
			return *reinterpret_cast<T*>(_base+_step*i);
		};
	};

	template<typename>
	struct _final;

	template<std::size_t ... ISEQ>
	struct _final<index_sequence<ISEQ...>>
	{
		static void call(char **args, npy_intp *dimensions, npy_intp *steps, void *extra)
		{
		    auto inputs  = std::make_tuple(data_handler<I_ARGS>{args[ISEQ], steps[ISEQ]}...);
		    auto outputs = data_handler<O_ARGS>{args[ISIZE], steps[ISIZE]};
		    npy_intp const n = dimensions[0];
		    for (int i = 0; i < n; i++) {
				outputs[i] = FUNC(std::get<ISEQ>(inputs)[i]...);
		     }
		}
	};

	static void ufunc(char **args, npy_intp *dimensions, npy_intp *steps, void *extra)
	{
		_final<ISEQ_TYPE>::call(args, dimensions, steps, extra);
	}

	PyObject * create_ufunc()
	{
		return PyUFunc_FromFuncAndData(func, data, types, 1, ISIZE, 1, PyUFunc_None, name.c_str(), doc.c_str(), 0);
	}

	void register_to(PyObject * module)
	{
	    auto ufunc = create_ufunc();
	    auto d = PyModule_GetDict(module);
	    PyDict_SetItemString(d, name.c_str(), ufunc);
	    Py_DECREF(ufunc);
	}

};

// If return type is a tuple
template<typename ... O_ARGS, typename ... I_ARGS, std::tuple<O_ARGS...>(&FUNC)(I_ARGS...)>
struct _my_build_ufunc<std::tuple<O_ARGS...>(I_ARGS...), FUNC>
{
	char types[sizeof...(I_ARGS) + sizeof...(O_ARGS)]; // handle function signature
	PyUFuncGenericFunction func[1]; // handle vectorized function
	void * data[1]; // handle extra data (not used by us currently

	std::string name;
	std::string doc;

	enum : int { ISIZE = sizeof...(I_ARGS) };
	enum : int { OSIZE = sizeof...(O_ARGS) };

	using ISEQ_TYPE = make_index_sequence<ISIZE>;
	using OSEQ_TYPE = make_index_sequence<OSIZE>;

	template<typename T>
	static inline T & _assign(T & dst, T && src) { return dst = src; }

	template<typename, typename>
	struct _update_types;

	template<std::size_t ... ISEQ, std::size_t ... OSEQ>
	struct _update_types<index_sequence<ISEQ...>, index_sequence<OSEQ...>>
	{
		static void update(char * types)
		{
			fold(_assign<char>(types[ISEQ], _python_bind_type_info<I_ARGS>::npy_type)...);
			fold(_assign<char>(types[ISIZE+OSEQ], _python_bind_type_info<O_ARGS>::npy_type)...);
		}
	};

	_my_build_ufunc(std::string const & name, std::string const & doc = "") :
		data{nullptr}, func{&ufunc}, name{name}, doc{doc}
	{
		_update_types<ISEQ_TYPE, OSEQ_TYPE>::update(types);
	}

	template<typename T>
	struct data_handler {
		char * const   _base;
		npy_intp const _step;

		data_handler(char * base, npy_intp step) : _base{base}, _step{step} { }

		T & operator[](int i)
		{
			return *reinterpret_cast<T*>(_base+_step*i);
		};
	};

	template<typename, typename>
	struct _final;

	template<std::size_t ... ISEQ, std::size_t ... OSEQ>
	struct _final<index_sequence<ISEQ...>, index_sequence<OSEQ...>>
	{
		static void call(char **args, npy_intp *dimensions, npy_intp *steps, void *extra)
		{
		    auto inputs  = std::make_tuple(data_handler<I_ARGS>{args[ISEQ], steps[ISEQ]}...);
		    auto outputs = std::make_tuple(data_handler<O_ARGS>{args[ISIZE+OSEQ], steps[ISIZE+OSEQ]}...);
		    npy_intp const n = dimensions[0];
		    for (int i = 0; i < n; i++) {
				std::tie(std::get<OSEQ>(outputs)[i]...) = FUNC(std::get<ISEQ>(inputs)[i]...);
		     }
		}
	};

	static void ufunc(char **args, npy_intp *dimensions, npy_intp *steps, void *extra)
	{
		_final<ISEQ_TYPE, OSEQ_TYPE>::call(args, dimensions, steps, extra);
	}

	PyObject * create_ufunc()
	{
		return PyUFunc_FromFuncAndData(func, data, types, 1, ISIZE, OSIZE, PyUFunc_None, name.c_str(), doc.c_str(), 0);
	}

	void register_to(PyObject * module)
	{
	    auto ufunc = create_ufunc();
	    auto d = PyModule_GetDict(module);
	    PyDict_SetItemString(d, name.c_str(), ufunc);
	    Py_DECREF(ufunc);
	}

};


#define make_binding(name) \
static _my_build_ufunc<decltype(sg1::name), sg1::name> ufunc_##name(#name);

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





