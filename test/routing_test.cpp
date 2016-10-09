#include "test_pch.h"
#include "routing_test.h"
#include "core/timeseries.h"
#include "core/utctime_utilities.h"
#include "core/geo_cell_data.h"
#include "api/timeseries.h"
#include <boost/math/distributions/gamma.hpp>

namespace shyft {
    namespace core {
        namespace routing {

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
				typedef typename  timeseries::convolve_w_ts<ts_t> output_m3s_t;
                geo_cell_data geo;
                ts_t discharge_m3s;
                std::vector<double> uhg(utctimespan dt) const {
					double steps = (geo.routing.distance / geo.routing.velocity)/dt;// time = distance / velocity[s] // dt[s]
					int n_steps = int(steps + 0.5);
					return make_uhg_from_gamma(n_steps, 3.0, 0.7);//std::vector<double>{0.1,0.5,0.2,0.1,0.05,0.030,0.020};
                }
                timeseries::convolve_w_ts<ts_t> output_m3s() const {
                    // return discharge
                    return timeseries::convolve_w_ts<ts_t>(discharge_m3s,uhg(discharge_m3s.time_axis().delta()),timeseries::convolve_policy::USE_ZERO);
                }
            };

            struct node {
                // binding information
                //  not really needed at core level, as we could to only with ref's in the core
                //  but we plan to expose to python and external persistence models
                //  so providing some handles do make sense for now.
                int id;// self.id
                int downstream_id;
                // ct(const vector<rn>& upstreams,vector<c> inflow_cells, transfer_fx)
                // upstreams .. computed
                //  a list of refs to routing_node's
                // each providing output_m3s
                // local-inflows, ..,
                // state ?

                //downstream node
                // downstream filter..
                //
                // output_m3s() = (locals+upstreams).convolve(uhg)
            };
            // routing network
            // vector<rn> rn
            //  1) a map<int,node> of identified nodes,
            //     with identified down-streams
            //     cells have routing id to this network.
            //  2) bind to cells
            //
            //
            template <class ts_t=int>
            struct model {
                std::map<int,node> node_map;
                // input cells..
                ts_t input_ts(int node_id /*cell begin..end */) {
                    // find all nodes with this as down_stream
                    // ask those nodes for
                    // .response_ts(..)
                    // return us1.response_ts()+..us2.response_ts();
                    return ts_t();
                }
                ts_t local_ts(int node_id/*cell begin..end */) {
                    // find all cells. with routing.id == node_id
                    // return c1.routing.output_ts() + c2.routing.output_ts()
                    return ts_t();
                }
                ts_t output_ts(int node_id) {
                    // return (local_ts(node_id) + input_ts(node_id)).convolve(w)
                    return ts_t();

                }
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
        std::cout<<i<<"\t"<<c1.output_m3s().value(i)<<"\t"<<c2.output_m3s().value(i)<<"\t"<<c3.output_m3s().value(i)<<"\t"<<sum_output_m3s.value(i)<<"\n";
		double exepected_value = c1.output_m3s().value(i) + c2.output_m3s().value(i) + c3.output_m3s().value(i);
		TS_ASSERT_DELTA(exepected_value, sum_output_m3s.value(i), 0.00001);
    }

}
#include <dlib/graph_utils.h>
#include <dlib/graph.h>
#include <dlib/directed_graph.h>

void routing_test::test_routing_model() {
    using namespace shyft::core;
    routing::model<int> m;
    // simple network
    // a->b-> d
    //    c-/

    routing::node a{1,2};
    routing::node b{2,4};
    routing::node c{3,4};
    routing::node d{4,0};
    m.node_map[a.id]=a;
    m.node_map[b.id]=b;
    m.node_map[c.id]=c;
    m.node_map[d.id]=d;
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
