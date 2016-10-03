#include "test_pch.h"
#include "routing_test.h"
#include "core/timeseries.h"
#include "core/utctime_utilities.h"

namespace shyft {
    namespace core {
        namespace routing {

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
void routing_test::test_hydrograph() {

// moved

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
