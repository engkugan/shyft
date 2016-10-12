#include "test_pch.h"
#include "hbv_stack_test.h"
#include "core/hbv_stack.h"
#include "mocks.h"
#include "core/timeseries.h"
#include "core/utctime_utilities.h"

// Some typedefs for clarity
using namespace shyft::core;
using namespace shyft::timeseries;
using namespace shyft::core::hbv_stack;

using namespace shyfttest::mock;
using namespace shyfttest;

namespace pt = shyft::core::priestley_taylor;
namespace snow = shyft::core::hbv_snow;
namespace soil = shyft::core::hbv_soil;
namespace tank = shyft::core::hbv_tank;
namespace ae = shyft::core::hbv_actual_evapotranspiration;
namespace pc = shyft::core::precipitation_correction;

typedef TSPointTarget<point_timeaxis> catchment_t;

namespace shyfttest {
	namespace mock {
		// need specialization for hbv_stack_response_t above
		template<> template<>
		void ResponseCollector<timeaxis>::collect<response>(size_t idx, const response& response) {
			_snow_output.set(idx, response.snow.outflow);
		}
		template <> template <>
		void DischargeCollector<timeaxis>::collect<response>(size_t idx, const response& response) {
			// hbv_outflow is given in mm, so compute the totals
			avg_discharge.set(idx, destination_area*response.tank.outflow / 1000.0 / 3600.0);
		}
	};
}; // End namespace shyfttest

void hbv_stack_test::test_call_stack() {
	xpts_t temp;
	xpts_t prec;
	xpts_t rel_hum;
	xpts_t wind_speed;
	xpts_t radiation;

	calendar cal;
	utctime t0 = cal.time(YMDhms(2014, 8, 1, 0, 0, 0));
	size_t n_ts_points = 3 * 24;
	utctimespan dt = deltahours(1);
	utctime t1 = t0 + n_ts_points*dt;
	shyfttest::create_time_series(temp, prec, rel_hum, wind_speed, radiation, t0, dt, n_ts_points);

	utctime model_dt = deltahours(24);
	vector<utctime> times;
	for (utctime i = t0; i <= t1; i += model_dt)
		times.emplace_back(i);
	timeaxis time_axis(t0, dt, n_ts_points);
	timeaxis state_time_axis(t0, dt, n_ts_points + 1);
	// Initialize parameters
	std::vector<double> s = { 1.0, 1.0, 1.0, 1.0, 1.0 }; // Zero cv distribution of snow (i.e. even)
	std::vector<double> a = { 0.0, 0.25, 0.5, 0.75, 1.0 };
	pt::parameter pt_param;
	snow::parameter snow_param(s, a);
	ae::parameter ae_param;
	soil::parameter soil_param;
	tank::parameter tank_param;
	pc::parameter p_corr_param;

	// Initialize the state vectors
	soil::state soil_state = {50.0};
	tank::state tank_state = {20.0, 10.0 };  // Check , I follow kirchner
	snow::state snow_state(10.0, 0.5);

	// Initialize response
	response run_response;

	// Initialize collectors
	shyfttest::mock::ResponseCollector<timeaxis> response_collector(1000 * 1000, time_axis);
	shyfttest::mock::StateCollector<timeaxis> state_collector(state_time_axis);

	state state{snow_state,soil_state, tank_state};
	parameter parameter(pt_param, snow_param, ae_param, soil_param, tank_param, p_corr_param);
	geo_cell_data geo_cell_data;
	hbv_stack::run_hbv_stack<direct_accessor, response>(geo_cell_data, parameter, time_axis, temp,  //What is the difference between ptgsk & pthsk??
		prec, wind_speed, rel_hum, radiation, state,
		state_collector, response_collector);

	auto snow_swe = response_collector.snow_swe();
	for (size_t i = 0; i < snow_swe.size(); ++i)
		TS_ASSERT(std::isfinite(snow_swe.get(i).v) && snow_swe.get(i).v >= 0);
}