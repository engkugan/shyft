#include "test_pch.h"
#include "routing_test.h"


void routing_test::test_hydrograph() {
 // arrange
 //  hydrograph = vector<double>, weights, .size()
 //   each for one timestep,
 //   .size()*dt is the length of the hydrograph
 //
 //  hydrograph_ts ( hydrograph, point_ts<ta_fixed> ts )
 //   .ts
 //   .hg (aka. vector<double>
 //  where
 //  .value(i) =
 //        double v_i = 0;
 //        size_t w = hg.size();
 //        for (size_t k= 0;k<w;++k)
 //               v_i += hg[k]* ts[i-k]; // this step, plus weights of
 //
 //        return v_it
 //
}
