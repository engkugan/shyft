#include "boostpython_pch.h"
#include <boost/python/docstring_options.hpp>
#include "core/utctime_utilities.h"
#include "core/priestley_taylor.h"
#include "core/actual_evapotranspiration.h"
#include "core/precipitation_correction.h"
#include "core/hbv_snow.h"
#include "core/hbv_soil.h"
#include "core/hbv_tank.h"
#include "core/hbv_stack.h"
#include "api/api.h"
#include "api/hbv_stack.h"
#include "core/hbv_stack_cell_model.h"
#include "core/region_model.h"
#include "core/model_calibration.h"
#include "expose_statistics.h"
#include "expose.h"

static char const* version() {
	return "v1.0";
}

namespace expose {
	namespace hbv_stack {
		using namespace boost::python;
		using namespace shyft::core;
		using namespace shyft::core::hbv_stack;

		static void parameter_state_response() {
			class_<parameter, bases<>, std::shared_ptr<parameter>>("HbvParameter",
				"Contains the parameters to the methods used in the HBV assembly\n"
				"priestley_taylor,hbv_snow,hbv_actual_evapotranspiration,hbv_soil, hbv_tank,precipitation_correction\n"
				)
				.def(init<const priestley_taylor::parameter&, const hbv_snow::parameter&, const hbv_actual_evapotranspiration::parameter&, const hbv_soil::parameter&, const hbv_tank::parameter&, const precipitation_correction::parameter&>(args("pt", "snow", "ae", "soil", "tank", "p_corr"), "create object with specified parameters"))
				.def(init<const parameter&>(args("p"), "clone a parameter"))
				.def_readwrite("pt", &parameter::pt, "priestley_taylor parameter")
				.def_readwrite("ae", &parameter::ae, "actual evapotranspiration parameter")
				.def_readwrite("snow", &parameter::snow, "hbv_snow parameter")
				.def_readwrite("soil", &parameter::soil, "hbv_soil parameter")
				.def_readwrite("tank", &parameter::soil, "hbv_tank parameter")
				.def_readwrite("p_corr", &parameter::p_corr, "precipitation correction parameter")
				.def("size", &parameter::size, "returns total number of calibration parameters")
				.def("set", &parameter::set, args("p"), "set parameters from vector/list of float, ordered as by get_name(i)")
				.def("get", &parameter::get, args("i"), "return the value of the i'th parameter, name given by .get_name(i)")
				.def("get_name", &parameter::get_name, args("i"), "returns the i'th parameter name, see also .get()/.set() and .size()")
				;
			typedef std::map<int, parameter> HbvParameterMap;
			class_<HbvParameterMap>("HbvParameterMap", "dict (int,parameter)  where the int is 0-based catchment_id")
				.def(map_indexing_suite<HbvParameterMap>())
				;

			class_<state>("HbvState")
				.def(init<hbv_snow::state, hbv_soil::state, hbv_tank::state>(args("snow", "soil", "tank"), "initializes state with hbv_snow, hbv_soil and hbv_tank"))
				.def_readwrite("snow", &state::snow, "hbv_snow state")
				.def_readwrite("soil", &state::soil, "soil state")
				.def_readwrite("tank", &state::soil, "tank state")
				;

			typedef std::vector<state> HbvStateVector;
			class_<HbvStateVector, bases<>, std::shared_ptr<HbvStateVector> >("HbvStateVector")
				.def(vector_indexing_suite<HbvStateVector>())
				;
			class_<response>("HbvResponse", "This struct contains the responses of the methods used in the Hbv assembly")
				.def_readwrite("pt", &response::pt, "priestley_taylor response")
				.def_readwrite("snow", &response::snow, "hbb_snow response")
				.def_readwrite("ae", &response::ae, "hbv_actual evapotranspiration response")
				.def_readwrite("soil", &response::soil, "hbv_soil response")
				.def_readwrite("tank", &response::tank, "hbv_tank response")
				.def_readwrite("total_discharge", &response::total_discharge, "total stack response")
				;
		}

		static void
			collectors() {
			typedef shyft::core::hbv_stack::all_response_collector HbvAllCollector;
			class_<HbvAllCollector>("HbvAllCollector", "collect all cell response from a run")
				.def_readonly("destination_area", &HbvAllCollector::destination_area, "a copy of cell area [m2]")
				.def_readonly("pe_output", &HbvAllCollector::pe_output, "pot evap mm/h")
				//.def_readonly("snow_sca",&HbvAllCollector::snow_sca," hbv snow covered area fraction, sca.. 0..1 - at the end of timestep (state)")
				//.def_readonly("snow_swe",&HbvAllCollector::snow_swe,"hbv snow swe, [mm] over the cell sca.. area, - at the end of timestep")
				.def_readonly("snow_outflow", &HbvAllCollector::snow_outflow, " hbv snow output [m3/s] for the timestep")
				.def_readonly("ae_output", &HbvAllCollector::ae_output, "actual evap mm/h")
				.def_readonly("soil_outflow", &HbvAllCollector::soil_outflow, " hbv soil output [m3/s] for the timestep") // Check??
				.def_readonly("avg_discharge", &HbvAllCollector::avg_discharge, "Tank outflow given in[m3/s] for the timestep")
				.def_readonly("end_reponse", &HbvAllCollector::end_reponse, "end_response, at the end of collected")
				;

			typedef shyft::core::hbv_stack::discharge_collector HbvDischargeCollector;
			class_<HbvDischargeCollector>("HbvDischargeCollector", "collect all cell response from a run")
				.def_readonly("destination_area", &HbvDischargeCollector::destination_area, "a copy of cell area [m2]")
				.def_readonly("snow_sca", &HbvDischargeCollector::snow_sca, " hbv snow covered area fraction, sca.. 0..1 - at the end of timestep (state)")
				.def_readonly("snow_swe", &HbvDischargeCollector::snow_swe, "hbv snow swe, [mm] over the cell sca.. area, - at the end of timestep")
				.def_readonly("avg_discharge", &HbvDischargeCollector::avg_discharge, "Tank Outflow given in [m3/s] for the timestep")
				.def_readonly("end_reponse", &HbvDischargeCollector::end_response, "end_response, at the end of collected")
				.def_readwrite("collect_snow", &HbvDischargeCollector::collect_snow, "controls collection of snow routine");

			typedef shyft::core::hbv_stack::null_collector HbvNullCollector;
			class_<HbvNullCollector>("HbvNullCollector", "collector that does not collect anything, useful during calibration to minimize memory&maximize speed")
				;

			typedef shyft::core::hbv_stack::state_collector HbvStateCollector;
			class_<HbvStateCollector>("HbvStateCollector", "collects state, if collect_state flag is set to true")
				.def_readwrite("collect_state", &HbvStateCollector::collect_state, "if true, collect state, otherwise ignore (and the state of time-series are undefined/zero)")
				.def_readonly("snow_swe", &HbvStateCollector::snow_swe, "")
				.def_readonly("snow_sca", &HbvStateCollector::snow_sca, "")
				.def_readonly("soil_moisture", &HbvStateCollector::soil_moisture, "")
				.def_readonly("tank_uz", &HbvStateCollector::tank_uz, "")
				.def_readonly("tank_lz", &HbvStateCollector::tank_lz, "");
		}

		static void
			cells() {
			typedef shyft::core::cell<parameter, environment_t, state, state_collector, all_response_collector> HbvCellAll;
			typedef shyft::core::cell<parameter, environment_t, state, null_collector, discharge_collector> HbvCellOpt;
			expose::cell<HbvCellAll>("HbvCellAll", "tbd: HbvCellAll doc");
			expose::cell<HbvCellOpt>("HbvCellOpt", "tbd: HbvCellOpt doc");
			expose::statistics::priestley_taylor<HbvCellAll>("HbvCell");
			expose::statistics::hbv_snow<HbvCellAll>("HbvCell");//it only gives meaning to expose the *All collect cell-type
			expose::statistics::hbv_actual_evapotranspiration<HbvCellAll>("HbvCell");
			expose::statistics::hbv_soil<HbvCellAll>("HbvCell");
			expose::statistics::hbv_tank<HbvCellAll>("HbvCell");
		}

		static void
			models() {
			typedef shyft::core::region_model<hbv_stack::cell_discharge_response_t> HbvOptModel;
			typedef shyft::core::region_model<hbv_stack::cell_complete_response_t> HbvModel;
			expose::model<HbvModel>("HbvModel", "Hbv_stack");
			expose::model<HbvOptModel>("HbvOptModel", "Hbv_stack");
		}

		static void
			state_io() {
			expose::state_io<shyft::api::hbv_stack_state_io, shyft::core::hbv_stack::state>("HbvStateIo");
		}


		static void
			model_calibrator() {
			expose::model_calibrator<shyft::core::region_model<hbv_stack::cell_discharge_response_t>>("HbvOptimizer");
		}
	}
}

BOOST_PYTHON_MODULE(_hbv_stack)
{

	boost::python::scope().attr("__doc__") = "SHyFT python api for the hbv_stack model";
	boost::python::def("version", version);
	boost::python::docstring_options doc_options(true, true, false);// all except c++ signatures
	expose::hbv_stack::state_io();
	expose::hbv_stack::parameter_state_response();
	expose::hbv_stack::cells();
	expose::hbv_stack::models();
	expose::hbv_stack::collectors();
	expose::hbv_stack::model_calibrator();
}

