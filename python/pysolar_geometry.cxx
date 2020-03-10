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

#include <iostream>
#include <string>

struct buffer {
	virtual ~buffer() { }
};

template<typename T>
struct buffer_impl : public buffer {
	T * data;

	buffer_impl(size_t n) {
		data = new T[n];
	}

	virtual ~buffer_impl() override
	{
		delete[] data;
	}
};

struct buffer_pyobject
{
    PyObject_HEAD
	buffer * data;
};

static void buffer_dealloc(buffer_pyobject * ths)
{
	delete ths->data;
	Py_TYPE(ths)->tp_free((PyObject*)ths);
}

static PyObject * buffer_new(PyTypeObject *type, PyObject *args, PyObject *kwds)
{
	buffer_pyobject * ths = (buffer_pyobject*)type->tp_alloc(type, 0);
	return (PyObject*)ths;
}

static int buffer_init(buffer_pyobject * ths, PyObject *args, PyObject *kwds)
{
	ths->data = nullptr;
	return 0;
}

static PyTypeObject buffer_pytypeobject = {
	    PyVarObject_HEAD_INIT(NULL, 0)
	    "solar_geometry._data_buffer", /*tp_name*/
	    sizeof(buffer_pyobject), /*tp_basicsize*/
	    0,                         /*tp_itemsize*/
	    (destructor)buffer_dealloc,  /*tp_dealloc*/
	    0,                         /*tp_print*/
	    0,                         /*tp_getattr*/
	    0,                         /*tp_setattr*/
	    0,                         /*tp_reserved*/
	    0,                         /*tp_repr*/
	    0,                         /*tp_as_number*/
	    0,                         /*tp_as_sequence*/
	    0,                         /*tp_as_mapping*/
	    0,                         /*tp_hash */
	    0,                         /*tp_call*/
	    0,                         /*tp_str*/
	    0,                         /*tp_getattro*/
	    0,                         /*tp_setattro*/
	    0,                         /*tp_as_buffer*/
	    Py_TPFLAGS_DEFAULT,        /*tp_flags*/
	    "Buffer for internal use only", /* tp_doc */
	    0,		               /* tp_traverse */
	    0,		               /* tp_clear */
	    0,		               /* tp_richcompare */
	    0,		               /* tp_weaklistoffset */
	    0,		               /* tp_iter */
	    0,		               /* tp_iternext */
		0,   /* tp_methods */
		0,   /* tp_members */
	    0,                         /* tp_getset */
	    0,                         /* tp_base */
	    0,                         /* tp_dict */
	    0,                         /* tp_descr_get */
	    0,                         /* tp_descr_set */
	    0,                         /* tp_dictoffset */
	    (initproc)buffer_init, /* tp_init */
	    0,                         /* tp_alloc */
		buffer_new,              /* tp_new */
};


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

template<typename T>
struct handle_numpy_1d_array {
	PyArrayObject * _arr;
	handle_numpy_1d_array(PyObject * arr) : _arr{nullptr}
	{
		_arr = reinterpret_cast<PyArrayObject*>(PyArray_FROM_OTF(arr, _python_bind_type_info<T>::npy_type, NPY_ARRAY_IN_ARRAY));
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


template<typename F, F &FUNC>
struct _my_build_vectorized_function;

// If return type is not a tuple
template<typename R, typename ... ARGS, R(&FUNC)(ARGS...)>
struct _my_build_vectorized_function<R(ARGS...), FUNC> {
	using func_type = R(ARGS...);

	// T0, TN expected to be handle_numpy_1d_array<X>
	template<typename T0, typename ... TN>
	static PyObject * _final_call(T0 args0,  TN ... args)
	{
		int nd = args0.size();
		npy_intp out_dims[] = {nd};
		PyArrayObject * out_arr = reinterpret_cast<PyArrayObject*>(PyArray_SimpleNew(1, out_dims,
				_python_bind_type_info<R>::npy_type));

		for (int i = 0; i < nd; ++i) {
			*reinterpret_cast<R*>(PyArray_GETPTR1(out_arr,i)) = FUNC(args0[i], args[i]...);
		}

		return reinterpret_cast<PyObject *>(out_arr);
	}

	template<typename>
	struct _pass1;

	template<std::size_t ... I>
	struct _pass1<index_sequence<I...>>
	{
		static PyObject * call(PyObject * args)
		{
			/* check for length of argument */
			if (PyTuple_GET_SIZE(args) != sizeof...(ARGS))
				return NULL;
			return _final_call(handle_numpy_1d_array<ARGS>(PyTuple_GET_ITEM(args, I))...);
		}
	};

	static PyObject * call(PyObject * self, PyObject * args)
	{
		return _pass1<make_index_sequence<sizeof...(ARGS)>>::call(args);
	}

};



// If return type is a tuple
template<typename ... TARGS, typename ... ARGS, std::tuple<TARGS...>(&FUNC)(ARGS...)>
struct _my_build_vectorized_function<std::tuple<TARGS...>(ARGS...), FUNC>
{
	using R = std::tuple<TARGS...>;
	using func_type = R(ARGS...);
	using field_name_type = std::array<std::string, sizeof...(TARGS)>;

	template<int I>
	static inline ptrdiff_t get_field_offset() {
		return reinterpret_cast<char*>(&std::get<I>(*reinterpret_cast<R*>(0)))-reinterpret_cast<char*>(0);
	}

	template<typename>
	struct _create_desc_for_tuple;

	template<std::size_t ... I>
	struct _create_desc_for_tuple<index_sequence<I...>> {
		static void update(PyObject * names, PyObject * formats, PyObject * offsets)
		{
			std::make_tuple(PyList_SetItem(formats, I, _python_bind_type_info<TARGS>::format())...);
			std::make_tuple(PyList_SetItem(offsets, I, PyLong_FromLong(get_field_offset<I>()))...);
		}
	};

	static inline PyArray_Descr * create_desc(field_name_type const & field_names) {
		static auto const tsize = sizeof...(TARGS);
		auto spec = PyDict_New();
		auto names = PyList_New(tsize);
		auto formats = PyList_New(tsize);
		auto offsets = PyList_New(tsize);

		PyDict_SetItemString(spec, "names", names);
		PyDict_SetItemString(spec, "formats", formats);
		PyDict_SetItemString(spec, "offsets", offsets);

		_create_desc_for_tuple<make_index_sequence<sizeof...(TARGS)>>::update(names, formats, offsets);

		/* must be handled by spec */
		Py_DECREF(names);
		Py_DECREF(formats);
		Py_DECREF(offsets);

		for (int i = 0; i < sizeof...(TARGS); ++i) {
				PyList_SetItem(names, i, Py_BuildValue("s", field_names[i].c_str()));
		}

		PyArray_Descr * final_desc;
		if (PyArray_DescrConverter(spec, &final_desc) != NPY_SUCCEED)
			throw std::runtime_error("SHOULD NEVER APPEND");

		Py_DECREF(spec);
		return final_desc;

	}

	// T0, TN expected to be handle_numpy_1d_array<X>
	template<typename T0, typename ... TN>
	static PyObject * final_call(field_name_type const & field_names, T0 args0,  TN ... args)
	{
		int nd = args0.size();

		auto buffer = new buffer_impl<R>(nd);
		auto output_array_data = buffer->data;

		for (int i = 0; i < nd; ++i) {
			new (&output_array_data[i]) R; // call tuple inplace constructor
			output_array_data[i] = FUNC(args0[i], args[i]...);
		}

		npy_intp out_dims[] = {nd};
		npy_intp strides[] = {sizeof(R)};

		auto desc = create_desc(field_names);

		/* create the output */
		PyObject *argList = Py_BuildValue("()");
		auto pybuffer = PyObject_CallObject(reinterpret_cast<PyObject*>(&buffer_pytypeobject), argList);
		if (pybuffer == nullptr)
			return nullptr;
		reinterpret_cast<buffer_pyobject*>(pybuffer)->data = buffer;

		PyArrayObject * out_arr = reinterpret_cast<PyArrayObject*>(PyArray_NewFromDescr(
				&PyArray_Type, desc, 1, out_dims, strides, output_array_data, 0, nullptr));
		PyArray_SetBaseObject(out_arr, pybuffer);

		return reinterpret_cast<PyObject*>(out_arr);
	}

	template<typename>
	struct _pass1;

	template<std::size_t ... I>
	struct _pass1<index_sequence<I...>> {
		static PyObject * call(field_name_type const & field_names, PyObject * args)
		{
			/* check for length of argument */
			if (PyTuple_GET_SIZE(args) != sizeof...(ARGS))
				return NULL;

			return final_call(field_names, handle_numpy_1d_array<ARGS>(PyTuple_GET_ITEM(args, I))...);
		}
	};

	static PyObject * call(field_name_type const & field_names, PyObject * self, PyObject * args)
	{
		return _pass1<make_index_sequence<sizeof...(ARGS)>>::call(field_names, args);
	}

};



#define make_binding(name) \
static PyObject * py_##name(PyObject * self, PyObject * args) \
{ \
	return _my_build_vectorized_function<decltype(sg1::name), sg1::name>::call(self, args); \
}

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

static PyObject * py_julian_day_to_ymd(PyObject * self, PyObject * args)
{
	return _my_build_vectorized_function<decltype(sg1::julian_day_to_ymd), sg1::julian_day_to_ymd>::call({"year", "month", "day_of_month"}, self, args);
}

static PyObject * py_day_of_year_to_ymd(PyObject * self, PyObject * args)
{
	return _my_build_vectorized_function<decltype(sg1::day_of_year_to_ymd), sg1::day_of_year_to_ymd>::call({"month", "day_of_month"}, self, args);
}


#define TPL_FUNCTION(name) {#name, py_##name, METH_VARARGS, "TODO: documment"}

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
	TPL_FUNCTION(ymd_to_day_of_year),
	TPL_FUNCTION(geogr_to_geoce),
	TPL_FUNCTION(azimuth_sun),
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

	if (PyType_Ready(&buffer_pytypeobject) != 0)
		return NULL;
	Py_INCREF(&buffer_pytypeobject);

	/** init numpy **/
	import_array();

	return m;

}





