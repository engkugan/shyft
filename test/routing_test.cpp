#include "test_pch.h"
#include "routing_test.h"
#include "core/timeseries.h"
#include "core/utctime_utilities.h"
#include "core/geo_cell_data.h"
#include "api/timeseries.h"
#include <boost/math/distributions/gamma.hpp>
#include <iomanip>

namespace shyft {
    namespace core {
        namespace routing {

			/** make_uhg_from_gamma a simple function to create a uhg (unit hydro graph) weight vector
			 * containing n_steps, given the gamma shape factor alpha and beta.
			 * ensuring the sum of the weight vector is 1.0
			 * and that it has a min-size of one element (1.0)
			 */

			std::vector<double>  make_uhg_from_gamma(int n_steps, double alpha, double beta) {
				using boost::math::gamma_distribution;
				gamma_distribution<double> gdf(alpha, beta);
				std::vector<double> r;r.reserve(n_steps);
				double s = 0.0;
				double d = 1.0 / double(n_steps);
				for (double q = d;q < 1.0; q += d) {
					double x = quantile(gdf, q);
					double y = pdf(gdf, x);
					s += y;
					r.push_back(y);
				}
				for (auto& y : r) y /= s;
				if (r.size() == 0) r.push_back(1.0);// at a minimum 1.0, no delay
				return std::move(r);
			};

            /// just for emulating a cell-node that do have
            /// the needed properties that we will require
            /// later when promoting the stuff to cell_model
            /// either by explicit requirement, or by concept
            template <class ts_t>
            struct cell_node {
				typedef ts_t ts_t;
				typedef typename  timeseries::convolve_w_ts<ts_t> output_m3s_t;
                geo_cell_data geo;
                ts_t discharge_m3s;
				/**create a unit hydro graph weight vector for the cell-node utilizing the
				 *  cell geo.specific routing information, distance&speed.
				 */
                std::vector<double> uhg(utctimespan dt) const {
					double steps = (geo.routing.distance / geo.routing.velocity)/dt;// time = distance / velocity[s] // dt[s]
					int n_steps = int(steps + 0.5);
					return std::move(make_uhg_from_gamma(n_steps, 3.0, 0.7));//std::vector<double>{0.1,0.5,0.2,0.1,0.05,0.030,0.020};
                }
                timeseries::convolve_w_ts<ts_t> output_m3s() const {
                    // return discharge, notice that this function assumes that time_axis() do have a uniform delta() (requirement)
                    return std::move(timeseries::convolve_w_ts<ts_t>(discharge_m3s,uhg(discharge_m3s.time_axis().delta()),timeseries::convolve_policy::USE_ZERO));
                }
            };

			/**\brief A routing node
			 *
			 * The routing node have flow from 
			 * -# zero or more 'cell_nodes',  typically a cell_model type, lateral flow (.output_m3s())
			 * -# zero or more upstream connected routing nodes, taking their inputs (.output_m3s())
			 * then a routing node can be connected to a down-stream node,
			 * providing a routing function (currently just a convolution of a uhg).
			 *
			 * This definition is recursive, and we use other means to ensure the routing graph
			 * is directed and with no cycles.
			 */
            struct node {
                // binding information
                //  not really needed at core level, as we could to only with ref's in the core
                //  but we plan to expose to python and external persistence models
                //  so providing some handles do make sense for now.
                int id;// self.id, >0 means valid id, 0= null
				shyft::core::routing_info downstream;
				std::vector<double> uhg(utctimespan dt) const {
					double steps = (downstream.distance / downstream.velocity) / dt;// time = distance / velocity[s] // dt[s]
					int n_steps = int(steps + 0.5);
					return std::move(make_uhg_from_gamma(n_steps, 3.0, 0.7));//std::vector<double>{0.1,0.5,0.2,0.1,0.05,0.030,0.020};
				}
            };

            // routing network
            // vector<rn> rn
            //  1) a map<int,node> of identified nodes,
            //     with identified down-streams
            //     cells have routing id to this network.
            //  2) bind to cells
            //
            //
            template<class C>
            struct model {
				typedef typename C::ts_t rts_t; // the result computed in the cell.rc.avg_discharge [m3/s]
				typedef typename C::output_m3s_t ots_t;
				typedef typename shyft::timeseries::uniform_sum_ts<ots_t> sum_ts_t;
				std::map<int, node> node_map; // keeps structure and routing properties
				std::shared_ptr<std::vector<C>> cells;

				/** compute the local lateral inflow from connected shyft-cells into given node-id */
				sum_ts_t local_inflow(int node_id) const {
					std::vector<ots_t> tsv;tsv.reserve(cells->size());// make room for all
					for (const auto& c : *cells) {
						if (c.geo.routing.id == node_id)
							tsv.push_back(c.output_m3s());
					}
					if (tsv.size() == 0) { // would be nice with a default convolve_w_ts(0.0) here,. for now we hack a 0.0 response, later we can avoid this entirely
						//ts_t nts((*cells)[0].discharge_ts.time_axis(), 0.0, shyft::timeseries::POINT_AVERAGE_VALUE);
						//std::vector<double> w;w.push_back(1.0);
						//tsv.push_back(ots_t(nts, w));
						throw std::runtime_error("routing::model: requested node_id has no local inflow");
					}
					return std::move(sum_ts_t(std::move(tsv)));
				}
				bool has_local_inflow(int node_id) const {
					for (const auto& c : *cells) if (c.geo.routing.id == node_id) return true;
					return false;
				}
				bool has_upstream_inflow(int node_id) const {
					for (const auto&i : node_map)
						if ( i.second.downstream.id == node_id) return true;
					return false;
				}
				/** compute the local lateral inflow from connected shyft-cells into given node-id */
				sum_ts_t upstream_inflow(int node_id) const {
					std::vector<ots_t> tsv;tsv.reserve(node_map.size());// make room for all
					for (const auto& i : node_map) {
						if ( i.second.downstream.id == node_id)
							tsv.push_back(output_m3s(i.first));
					}
					if (tsv.size() == 0) { // would be nice with a default convolve_w_ts(0.0) here,. for now we hack a 0.0 response, later we can avoid this entirely
						throw std::runtime_error("routing::model: requested node_id has no upstream inflow");
					}
					return std::move(sum_ts_t(std::move(tsv)));
				}
#if 0
				// will not go through compiler, since the expressions have different types for 
				// cases
				//  local & upstream
				//  local only
				//  upstream only
				// additionally, the recursive dimension creates separate type-trees for each level upstream.
				// so we either need to create 3 methods (one for each case)
				//   plus using result_of ..
				// or flatten out the result by extracting ts_values 
				// or go to dynamic dispatch for the ts.

				ots_t output_m3s(int node_id) const {

					if (has_local_inflow(node_id) && has_upstream_inflow(node_id)) {
						auto sum_inflow = local_inflow(node_id) + upstream_inflow(node_id);
						utctimespan dt=sum_inflow.time_axis().delta();
						return timeseries::convolve_w_ts(sum_inflow, node_map[node_id].uhg(dt),timeseries::POINT_AVERAGE...)
					} else if (has_local_inflow(node_id)) {
						auto sum_inflow = local_inflow(node_id) ;
						utctimespan dt = sum_inflow.time_axis().delta();
						return timeseries::convolve_w_ts(sum_inflow, node_map[node_id].uhg(dt), timeseries::POINT_AVERAGE...)
					} else if (has_upstream_inflow(node_id)) {

					} else {
						throw std::runtime_error("routing::model: requested node has no upstream or local cell inflow");
					}
				}



											  // input cells..
                rts_t input_ts(int node_id /*cell begin..end */) {
                    // find all nodes with this as down_stream
                    // ask those nodes for
                    // .response_ts(..)
                    // return us1.response_ts()+..us2.response_ts();
                    return rts_t();
                }
                rts_t output_ts(int node_id) {
                    // return (local_ts(node_id) + input_ts(node_id)).convolve(w)
                    return rts_t();

                }
#endif
            };

        }
    }
}

namespace shyft {namespace timeseries {


// maybe a more convenient function to get the frozen values out:
template <class Ts>
std::vector<double> ts_values(const Ts& ts) {
	std::vector<double> r; r.reserve(ts.size());
	for (size_t i = 0;i < ts.size();++i) r.push_back(ts.value(i));
	return std::move(r);
}
}
}
void routing_test::test_hydrograph() {
	//
	//  Just for playing around with sum of routing node response
	// 
    using ta_t =shyft::time_axis::fixed_dt;
    using ts_t = shyft::timeseries::point_ts<ta_t>;
	
	//using ats_t = shyft::api::apoint_ts;
    using namespace shyft::core;
    //using namespace shyft::timeseries;
	using node_t = routing::cell_node<ts_t>;
	node_t c1;
	c1.geo.routing.id = 1; // target routing point 1
	c1.geo.routing.distance = 10000;
	c1.geo.routing.velocity = c1.geo.routing.distance / (10 * 3600.0);// takes 10 hours to propagate the distance
	calendar utc;
    ta_t ta(utc.time(2016,1,1),deltahours(1),24);
    c1.discharge_m3s= ts_t(ta,0.0,shyft::timeseries::POINT_AVERAGE_VALUE);
    c1.discharge_m3s.set(0,10.0);

    auto c2 =c1;c2.geo.routing.distance = 5 * 1000;
    c2.discharge_m3s.set(0,20.0);
	auto c3 = c1; c3.geo.routing.distance = 1 * 1000;
	c3.discharge_m3s.set(0, 5.0);
	c3.discharge_m3s.set(1, 35.0);
	std::vector<node_t::output_m3s_t> responses;
	responses.push_back(c1.output_m3s());
	responses.push_back(c2.output_m3s());
	responses.push_back(c3.output_m3s());
	shyft::timeseries::uniform_sum_ts<node_t::output_m3s_t> sum_output_m3s(responses);
	std::cout << "\nresult:\n";
    for(size_t i=0;i<ta.size();++i) {
        std::cout<<i<<"\t"<<std::setprecision(2)<<c1.output_m3s().value(i)<<"\t"<<c2.output_m3s().value(i)<<"\t"<<c3.output_m3s().value(i)<<"\t"<<sum_output_m3s.value(i)<<"\n";
		double exepected_value = c1.output_m3s().value(i) + c2.output_m3s().value(i) + c3.output_m3s().value(i);
		TS_ASSERT_DELTA(exepected_value, sum_output_m3s.value(i), 0.00001);
    }

}
#include <dlib/graph_utils.h>
#include <dlib/graph.h>
#include <dlib/directed_graph.h>

void routing_test::test_routing_model() {
    using namespace shyft::core;
	using ta_t = shyft::time_axis::fixed_dt;
	using ts_t = shyft::timeseries::point_ts<ta_t>;

    routing::model<routing::cell_node<ts_t>> m;
    // simple network
    // a->b-> d
    //    c-/

	routing::node a{ 1,routing_info(2) };
	routing::node b{ 2,routing_info(4) };
	routing::node c{ 3,routing_info(4) };
	routing::node d{ 4,routing_info(0) };
    m.node_map[a.id]=a;
    m.node_map[b.id]=b;
    m.node_map[c.id]=c;
    m.node_map[d.id]=d;
	TS_ASSERT(m.has_upstream_inflow(2));
	TS_ASSERT(!m.has_upstream_inflow(1));
    dlib::directed_graph<routing::node>::kernel_1a_c md;
    md.set_number_of_nodes(4);
    md.add_edge(0,1);
    md.add_edge(1,3);
    md.add_edge(2,3);
    md.node(0).data=a;
    md.node(1).data=b;
    md.node(2).data=c;
    md.node(3).data=d;

    TS_ASSERT_EQUALS( md.node(3).number_of_parents(),2);
    TS_ASSERT_EQUALS( md.node(0).number_of_parents(),0);
    TS_ASSERT_EQUALS( md.node(0).number_of_children(),1);
    TS_ASSERT_EQUALS( md.node(3).number_of_children(),0);
}
