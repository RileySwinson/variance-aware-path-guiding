#ifdef NANOGUI_PYTHON

#include "python.h"

typedef floatBox<double> DoubleBox;
typedef IntBox<int64_t> Int64Box;

DECLARE_WIDGET(TextBox);
DECLARE_WIDGET(DoubleBox);
DECLARE_WIDGET(Int64Box);

void register_textbox(py::module &m) {
    py::class_<TextBox, Widget, ref<TextBox>, PyTextBox> tbox(m, "TextBox", D(TextBox));
    tbox
        .def(py::init<Widget *, const std::string &>(), py::arg("parent"),
            py::arg("value") = std::string("Untitled"), D(TextBox, TextBox))
        .def("editable", &TextBox::editable, D(TextBox, editable))
        .def("setEditable", &TextBox::setEditable, D(TextBox, setEditable))
        .def("spinnable", &TextBox::spinnable, D(TextBox, spinnable))
        .def("setSpinnable", &TextBox::setSpinnable, D(TextBox, setSpinnable))
        .def("value", &TextBox::value, D(TextBox, value))
        .def("setValue", &TextBox::setValue, D(TextBox, setValue))
        .def("defaultValue", &TextBox::defaultValue, D(TextBox, defaultValue))
        .def("setDefaultValue", &TextBox::setDefaultValue, D(TextBox, setDefaultValue))
        .def("alignment", &TextBox::alignment, D(TextBox, alignment))
        .def("setAlignment", &TextBox::setAlignment, D(TextBox, setAlignment))
        .def("units", &TextBox::units, D(TextBox, units))
        .def("setUnits", &TextBox::setUnits, D(TextBox, setUnits))
        .def("unitsImage", &TextBox::unitsImage, D(TextBox, unitsImage))
        .def("setUnitsImage", &TextBox::setUnitsImage, D(TextBox, setUnitsImage))
        .def("format", &TextBox::format, D(TextBox, format))
        .def("setFormat", &TextBox::setFormat, D(TextBox, setFormat))
        .def("callback", &TextBox::callback, D(TextBox, callback))
        .def("setCallback", &TextBox::setCallback, D(TextBox, setCallback));

    py::enum_<TextBox::Alignment>(tbox, "Alignment", D(TextBox, Alignment))
        .value("Left", TextBox::Alignment::Left)
        .value("Center", TextBox::Alignment::Center)
        .value("Right", TextBox::Alignment::Right);

    py::class_<Int64Box, TextBox, ref<Int64Box>, PyInt64Box>(m, "IntBox", D(IntBox))
        .def(py::init<Widget *, int64_t>(), py::arg("parent"), py::arg("value") = (int64_t) 0, D(IntBox, IntBox))
        .def("value", &Int64Box::value, D(IntBox, value))
        .def("setValue", (void (Int64Box::*)(int64_t)) &Int64Box::setValue, D(IntBox, setValue))
        .def("setCallback", (void (Int64Box::*)(const std::function<void(int64_t)>&))
                &Int64Box::setCallback, D(IntBox, setCallback))
        .def("setValueIncrement", &Int64Box::setValueIncrement, D(IntBox, setValueIncrement))
        .def("setMinValue", &Int64Box::setMinValue, D(IntBox, setMinValue))
        .def("setMaxValue", &Int64Box::setMaxValue, D(IntBox, setMaxValue))
        .def("setMinValue", &Int64Box::setMinMaxValues, D(IntBox, setMinMaxValues));

    py::class_<DoubleBox, TextBox, ref<DoubleBox>, PyDoubleBox>(m, "floatBox", D(floatBox))
        .def(py::init<Widget *, double>(), py::arg("parent"), py::arg("value") = 0.0)
        .def("value", &DoubleBox::value, D(floatBox, value))
        .def("setValue", (void (DoubleBox::*)(double)) &DoubleBox::setValue, D(floatBox, setValue))
        .def("setCallback", (void (DoubleBox::*)(const std::function<void(double)>&))
                &DoubleBox::setCallback, D(floatBox, setCallback))
        .def("setValueIncrement", &DoubleBox::setValueIncrement, D(floatBox, setValueIncrement))
        .def("setMinValue", &DoubleBox::setMinValue, D(floatBox, setMinValue))
        .def("setMaxValue", &DoubleBox::setMaxValue, D(floatBox, setMaxValue))
        .def("setMinValue", &DoubleBox::setMinMaxValues, D(floatBox, setMinMaxValues));
}

#endif
