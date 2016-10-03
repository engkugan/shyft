#pragma once

#include "priestley_taylor.h"
#include "hbv_tank.h"
#include "hbv_snow.h"
#include "actual_evapotranspiration.h"
#include "precipitation_correction.h"
#include "hbv_soil.h"

namespace shyft {
	namespace core {
		namespace hbv_stack {
			using namespace std;

			/** \brief Simple parameter struct for the hbv_stack method stack
			*
			* This struct contains the parameters to the methods used in the HBV assembly.
			*
			* \tparam PTParameter PriestleyTaylor parameter type that implements the interface:
			*    - PTParameter.albedo const --> double, land albedo parameter in PriestleyTaylor.
			*    - PTParameter.alpha const --> double, alpha parameter in PriestleyTaylor.
			* \tparam SnowState HBVSnow parameter type that implements the parameter interface for HBVSnow.
			* \tparam KState hbv_tank parameter type that implements the parameter interface for hbv_tank.
			* \sa HBVSnowParameter \sa hbv_tank \sa hbv_stack  \sa PriestleyTaylor
			*/
			struct parameter {
				typedef priestley_taylor::parameter pt_parameter_t;
				typedef hbv_snow::parameter snow_parameter_t;
				typedef actual_evapotranspiration::parameter ae_parameter_t;
				typedef hbv_tank::parameter tank_parameter_t;
				typedef precipitation_correction::parameter precipitation_correction_parameter_t;
				typedef hbv_soil::parameter soil_parameter_t;

				pt_parameter_t pt;
				snow_parameter_t snow;
				ae_parameter_t ae;
				tank_parameter_t tank;
				precipitation_correction_parameter_t p_corr;
				soil_parameter_t soil;
				
				parameter(pt_parameter_t pt,
					snow_parameter_t snow,
					ae_parameter_t ae,
					tank_parameter_t tank,
					precipitation_correction_parameter_t p_corr,
					soil_parameter_t soil)
					: pt(pt), snow(snow), ae(ae), tank(tank), p_corr(p_corr), soil(soil) { /*Do nothing */
				}

				parameter() {}
				parameter(const parameter& other)
					: pt(other.pt), snow(other.snow), ae(other.ae),
					tank(other.tank), p_corr(other.p_corr), soil(other.soil) { /*Do nothing */
				}

				pt_parameter_t pt;
				snow_parameter_t snow;
				ae_parameter_t ae;
				tank_parameter_t  tank;
				precipitation_correction_parameter_t p_corr;
				soil_parameter_t soil;

				///<calibration support, needs vector interface to params, size is the total count
				size_t size() const { return 16 + 1; }
				///<calibration support, need to set values from ordered vector
				void set(const vector<double>& p) {
					if (p.size() != size())
						throw runtime_error("HBV_Stack Parameter Accessor: .set size missmatch");
					int i = 0;
					tank.uz1 = p[i++];
					tank.kuz2 = p[i++];
					tank.kuz1 = p[i++];
					tank.perc = p[i++];
					tank.klz = p[i++];
					ae.ae_scale_factor = p[i++];
					snow.lw = p[i++];
					snow.tx = p[i++];
					snow.cx = p[i++];
					snow.ts = p[i++];
					snow.cfr = p[i++];
					p_corr.scale_factor = p[i++];
					pt.albedo = p[i++];
					pt.alpha = p[i++];
					soil.fc = p[i++];
					soil.beta = p[i++];
					soil.lp = p[i++];
				}

				///< calibration support, get the value of i'th parameter
				double get(size_t i) const {
					switch (i) {
					case  0:return tank.uz1;
					case  1:return tank.kuz2;
					case  2:return tank.kuz1;
					case  3:return tank.perc;
					case  4:return tank.klz;
					case  5:return ae.ae_scale_factor;
					case  6:return snow.lw;
					case  7:return snow.tx;
					case  8:return snow.cx;
					case  9:return snow.ts;
					case 10:return snow.cfr;
					case 11:return p_corr.scale_factor;
					case 12:return pt.albedo;
					case 13:return pt.alpha;
					case 14: return soil.fc;
					case 15: return soil.beta;
					case 16: return soil.lp;
					default:
						throw runtime_error("HBV_stack Parameter Accessor:.get(i) Out of range.");
					}
					return 0;
				}

				///< calibration and python support, get the i'th parameter name
				string get_name(size_t i) const {
					static const char *names[] = {
						"tank.uz1",
						"tank.kuz2",
						"tank.kuz1",
						"tank.perc",
						"tank.klz",
						"ae.ae_scale_factor",
						"snow.lw",
						"snow.tx",
						"snow.cx",
						"snow.ts",
						"snow.cfr",
						"p_corr.scale_factor",
						"pt.albedo",
						"pt.alpha",
						"soil.fc",
						"soil.beta",
						"soil.lp"
					};
					if (i >= size())
						throw runtime_error("hbv_stack Parameter Accessor:.get_name(i) Out of range.");
					return names[i];
				}
			};

			/** \brief Simple state struct for the hbv_stack method stack
			*
			* This struct contains the states of the methods used in the hbv_stack assembly.
			*
			* \tparam SnowState HBVSnow state type that implements the state interface for GammaSnow.
			* \tparam TankState hbv_tank state type that implements the state interface for hbv_tank.
			* \sa HBVSnowState \sa hbv_tank \sa hbv_stack\sa PriestleyTaylor \sa HBVSnow
			*/
			struct state {
				typedef hbv_snow::state snow_state_t;
				typedef hbv_tank::state tank_state_t;
				typedef hbv_soil::state soil_state_t;
				state() {}
				state(snow_state_t snow, tank_state_t tank, soil_state_t soil) : snow(snow), tank(tank), soil(soil) {}
				snow_state_t snow;
				tank_state_t tank;
				soil_state_t soil;
				bool operator==(const state& x) const { return snow == x.snow && tank == x.tank && soil == x.soil; }
			};


			/** \brief Simple response struct for the hbv_stack method stack
			*
			* This struct contains the responses of the methods used in the HBV_stack assembly.
			*/
			struct response {
				// Model responses
				typedef priestley_taylor::response  pt_response_t;
				typedef hbv_snow::response snow_response_t;
				typedef actual_evapotranspiration::response  ae_response_t;
				typedef hbv_tank::response tank_response_t;
				typedef hbv_soil::response soil_response_t;
				pt_response_t pt;
				snow_response_t snow;
				ae_response_t ae;
				tank_response_t tank;
				soil_response_t soil;


				// Stack response
				double total_discharge;
			};

			/** \brief Calculation Model using assembly of PriestleyTaylor, GammaSnow, Soil and hbv_tank
			*
			* This model first uses PriestleyTaylor for calculating the potential
			* evapotranspiration based on time series data for temperature,
			* radiation and relative humidity. Then it uses the GammaSnow method
			* to calculate the snow/ice adjusted runoff using time series data for
			* precipitation and wind speed in addition to the time series used in
			* the PriestleyTaylor method. The PriestleyTaylor potential evaporation is
			* is used to calculate the actual evapotranspiration that is passed on to the
			* last step, hbv_tank.
			* hbv_tank is run with the gamma snow output
			* and actual evaporation response from the two methods above to
			* calculate the discharge.
			*
			* TODO: This method should construct an internal time stepping hbv_stack struct.
			* This stepper should be used as an iterator in time integration loops. This
			* paves the way for inter-cell communication (e.g. snow transport) without
			* touching this simple interface.
			*
			* \tparam TS Time series type that implements:
			*    - TS::source_accessor_type --> Type of accessor used to retrieve time series data.
			*    - TS.accessor(const T& time_axis) const --> TS::source_accessor_type. Accessor object
			*      for this time series.
			*
			* \tparam TS::source_accessor_type Time series source accessor type that implements:
			*    - TS::source_accessor_type(const TS& time_series, const T& time_axis) --> Construct
			*      accessor for the given time series and time axis.
			*    - TS::source_accessor_type.value(size_t i) --> double, -value of the time series for
			*      period i in the time axis.
			* \tparam T Time axis type that implements:
			*    - T.size() const --> Number of periods in the time axis.
			*    - T(size_t i) const --> shyft::core::utcperiod, - Start and end as shyft::core::utctime
			*      of the i'th period.
			* \tparam S State type that implements:
			*    - S::snow_state_type --> State type for the GammaSnow method.
			*    - S::hbv_tank_state_type --> State type for the hbv_tank method.
			*    - S.snow --> S::snow_state_type, - State variables for the GammaSnow method
			*    - S.hbv_tank --> S::hbv_tank_state_type, - State variables for the hbv_tank method
			* \tparam R Response type that implements:
			*    - R::snow_response_type --> Response type for the GammaSnow routine.
			*    - R.snow --> R::snow_response_type, -Response object passed to the GammaSnow routine.
			* \tparam P Parameter type that implements:
			*    - P::pt_parameter_type --> Parameter type for the PriestleyTaylor method
			*    - P::snow_parameter_type --> Parameter type for the GammaSnow method.
			*    - P::ae_parameter_type --> Parameter type for the ActualEvapotranspiration method.
			*    - P::tank_parameter_type --> Parameter type for the hbv_tank method.
			*    - P.pt --> P::pt_parameter_type --> Parameters for the PriestleyTaylor method.
			*    - P.snow --> P::snow_parameter_type --> Parameters for thge GammaSnow method.
			*    - P.ae --> P::ae_parameter_type --> Parameters for thge ActualEvapotranspiration method.
			*    - P.tank --> P::hbv_tank_parameter_type --> Parameters for the hbv_tank method.
			* \tparam SC State collector type that implements:
			*    - SC.collect(utctime t, const S& state) --> Possibly save some states at time t.
			* \tparam RC Response collector type that implements:
			*    - RC.collect(utctime t, const R& response) --> Possibly save some responses at time t.
			*/

			template<template <typename, typename> class A, class R, class T_TS, class P_TS, class WS_TS, class RH_TS, class RAD_TS, class T,
				class S, class GCD, class P, class SC, class RC>
				void run_hbv_stack(const GCD& geo_cell_data,
					const P& parameter,
					const T& time_axis,
					const T_TS& temp,
					const P_TS& prec,
					const WS_TS& wind_speed,
					const RH_TS& rel_hum,
					const RAD_TS& rad,
					S& state,
					SC& state_collector,
					RC& response_collector
				) {
				// Access time series input data through accessors of template A (typically a direct accessor).
				using temp_accessor_t = A<T_TS, T>;
				using prec_accessor_t = A<P_TS, T>;
				using wind_speed_accessor_t = A<WS_TS, T>;
				using rel_hum_accessor_t = A<RH_TS, T>;
				using rad_accessor_t = A<RAD_TS, T>;

				auto temp_accessor = temp_accessor_t(temp, time_axis);
				auto prec_accessor = prec_accessor_t(prec, time_axis);
				auto wind_speed_accessor = wind_speed_accessor_t(wind_speed, time_axis);
				auto rel_hum_accessor = rel_hum_accessor_t(rel_hum, time_axis);
				auto rad_accessor = rad_accessor_t(rad, time_axis);

				// Initialize the method stack
				precipitation_correction::calculator p_corr(parameter.p_corr.scale_factor);
				priestley_taylor::calculator pt(parameter.pt.albedo, parameter.pt.alpha);
				hbv_snow::calculator<typename P::snow_parameter_t, typename S::snow_state_t> hbv_snow(parameter.snow, state.snow);
				// -- 
				// hbv_soil::calculater<typename P::hbv_soil_paramer_t, typename S::hbv_soil_state, ...> hbv_soil(parameter.hbv_soil);
				hbv_soil::calculator<typename P::soil_parameter_t> soil(parameter.soil);
				hbv_tank::calculator<typename P::tank_parameter_t> tank(parameter.tank);
				// hbv_tank::calculator<hbv_tank::trapezoidal_average, typename P::hbv_tank_parameter_t> hbv_tank(parameter.hbv_tank);
				//
				snow.set_glacier_fraction(geo_cell_data.land_type_fractions_info().glacier());
				// Get the initial states
				auto &snow_state = state.snow;
				//double q = state.tank.q;
				double uz = state.tank.uz; // Check with Prof
				double lz = state.tank.lz; // Cjeck with Prof
				R response;
				const double forest_fraction = geo_cell_data.land_type_fractions_info().forest();
				const double altitude = geo_cell_data.mid_point().z;
				// Step through times in axis
				for (size_t i = 0; i < time_axis.size(); ++i) {
					utcperiod period = time_axis.period(i);
					double temp = temp_accessor.value(i);
					double rad = rad_accessor.value(i);
					double rel_hum = rel_hum_accessor.value(i);
					double prec = p_corr.calc(prec_accessor.value(i));
					state_collector.collect(i, state);///< \note collect the state at the beginning of each period (the end state is saved anyway)
													  //
													  // Land response:
													  //

													  // PriestleyTaylor (scale by timespan since it delivers results in mm/s)
					double pot_evap = pt.potential_evapotranspiration(temp, rad, rel_hum)*calendar::HOUR; // mm/h
					response.pt.pot_evapotranspiration = pot_evap;

					// Gamma Snow
					//snow.step(snow_state, response.snow, period.start, period.timespan(), parameter.snow,
					//	temp, rad, prec, wind_speed_accessor.value(i), rel_hum, forest_fraction, altitude);
					//HBV Snow
					hbv_snow.step(state.snow, response.snow, period.start, period.end, parameter.snow, prec, temp);

					// TODO: Communicate snow
					// At my pos xx mm of snow moves in direction d.

					// Actual Evapotranspiration //Check with Prof
						/*double act_evap = actual_evapotranspiration::calculate_step(q, pot_evap,
						parameter.ae.ae_scale_factor, response.gs.sca, period.timespan());
					response.ae.ae = act_evap;
					*/
					// soil layer goes here:
					// if( geo_cell_data.something), decide action.
					soil.step(state.soil, response.soil, period.start, period.end, response.snow.outflow, pot_evap, act_evap);
					// q_soil=hbv_soil.step(state.hbv_soil, response.snow.outflow);
					// Use responses from PriestleyTaylor and HBVSnow in hbv_tank
					//double q_avg;                                    // q_soil
					//hbv_tank.step(period.start, period.end, q, q_avg, response.soil.outflow, 0.0, geo_cell_data.urban_area);
					//state.hbv_tank.q = q; // Save discharge state variable
					//response.hbv_tank.q_avg = q_avg;

					//
					// Adjust land response for lakes and reservoirs (Treat them the same way for now)
					double total_lake_fraction = geo_cell_data.land_type_fractions_info().lake() +
						geo_cell_data.land_type_fractions_info().reservoir();
					double corrected_discharge = prec*total_lake_fraction +
						response.snow.outflow*geo_cell_data.land_type_fractions_info().glacier() +
						(1.0 - total_lake_fraction - geo_cell_data.land_type_fractions_info().glacier())*q_avg;
					response.total_discharge = corrected_discharge;

					// Possibly save the calculated values using the collector callbacks.
					response_collector.collect(i, response);///< \note collect the response valid for the i'th period (current state is now at the end of period)
					if (i + 1 == time_axis.size())
						state_collector.collect(i + 1, state);///< \note last iteration,collect the  final state as well.
				}
				response_collector.set_end_response(response);
			}
		} // hbv_stack
	} // core
} // shyft
