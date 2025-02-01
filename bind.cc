#include <sstream>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/stl_bind.h>
#include <pybind11/functional.h>

#include "rigid_body.h"
#include "rigid_body_system.h"
#include "Elements/Beam.h"
#include "Elements/Connector.h"
#include "Elements/Spring.h"
#include "equation_rbs.h"

namespace py = pybind11;

PYBIND11_MODULE(SR, rbd) {

    // adds the BLA bindings as submodule, accessible as rigid_body.bla.Matrix etc.
    // https://github.com/pybind/pybind11/discussions/4027
    auto m = rbd.def_submodule("bla", "basic linear algebra");
    #include "bind_bla_obj.h"



    // the main bindings:
    rbd.doc() = "rigid body Shake and Rattle simulator";

    py::class_<Connector>(rbd,"Connector", py::module_local())
      .def(py::init<>())
      .def(py::init<Vector<double>>())
      .def(py::init<Vector<double>, size_t>())
      .def_property("RefPos",[](Connector& c){return py::make_tuple(c.RefPosition(0),c.RefPosition(1),c.RefPosition(2));},
                             [](Connector&c, std::array<double,3> vals){c.RefPosition() = {vals[0], vals[1], vals[2]};})
      .def_property("body_index", &Connector::BodyIndex, &Connector::BodyIndex)
      .def_property("type", &Connector::Fix, &Connector::Fix);

    py::class_<Beam>(rbd,"Beam", py::module_local())
      .def(py::init<>([](Connector a, Connector b){return Beam{a,b};}))
      .def_property_readonly("length", [](Beam& b){return b.Length();})
      .def_property_readonly("connectorA", [](Beam& b){return b.ConnectorA();})
      .def_property_readonly("connectorB",[](Beam& b){return b.ConnectorB();});

    py::class_<Spring>(rbd,"Spring", py::module_local())
      .def(py::init<>([](Connector a, Connector b){return Spring{a,b};})) // stiffness should be positive! (-k)
      .def_property_readonly("stiffness", [](Spring& b){return b.Stiffness();}) // stiffness should be positive! (-k)
      .def_property_readonly("connectorA", [](Spring& b){return b.ConnectorA();})
      .def_property_readonly("connectorB",[](Spring& b){return b.ConnectorB();});

    py::class_<RigidBody> (rbd, "RigidBody_FEM")
      .def(py::init<>())
      //.def(py::init<Vector<double>, Vector<double>>())
      .def_property_readonly("index", [](RigidBody& rb){ return rb.Index();});
    
    py::class_<RigidBodySystem> (rbd, "RBS_FEM_SR")
      .def(py::init<>())
      .def("add", py::overload_cast<RigidBody&>(&RigidBodySystem::add))
      .def("add", py::overload_cast<Beam&>(&RigidBodySystem::add))
      .def("add", py::overload_cast<Spring&>(&RigidBodySystem::add))
      .def("saveState", &RigidBodySystem::SaveState);

    rbd.def("simulate",[](RigidBodySystem& rbs, double tend, double steps) {simulate_rbs(rbs, tend, steps);});
}
