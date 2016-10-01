#include "test_pch.h"
#include "routing_test.h"
#include "core/timeseries.h"
#include "core/utctime_utilities.h"

namespace shyft {
    namespace timeseries {

        /* READY TO BE PROMOTED to timeseries */

        /** The convolve_policy determines how the convolve_w_ts functions deals with the
         *  initial boundary condition, i.e. when trying to convolve with values before
         *  the first value, value(0)
         * \sa convolve_w_ts
         */
        enum convolve_policy {
            USE_FIRST, ///< ts.value(0) is used for all values before value(0): 'mass preserving'
            USE_ZERO, ///< fill in zero for all values before value(0):shape preserving
            USE_NAN ///< nan filled in for the first length of the filter
        };
        /** \brief convolve_w convolves a time-series with weights w
         *
         * The resulting time-series value(i) is the result of convolution (ts*w)|w.size()
         *  value(i) => ts(i-k)*w(k), k is [0..w.size()>
         *
         * the convolve_policy determines how to resolve start-condition-zone
         * where i < w.size(), \ref convolve_policy
         *
         * \tparam Ts any time-series
         * \tparam W any container with .size() and [i] operator
         */
        template<class Ts,class W>
        struct convolve_w_ts {
            typedef typename Ts::ta_t ta_t;
            Ts ts;
            point_interpretation_policy fx_policy; // inherited from ts
            W w;
            convolve_policy policy=convolve_policy::USE_FIRST;
            //-- default stuff, ct/copy etc goes here
            convolve_w_ts():fx_policy(POINT_AVERAGE_VALUE) {}
            convolve_w_ts(const convolve_w_ts& c):ts(c.ts),fx_policy(c.fx_policy),policy(c.policy) {}
            convolve_w_ts(convolve_w_ts&&c):ts(std::move(c.ts)),fx_policy(c.fx_policy),policy(c.policy) {}

            convolve_w_ts& operator=(const convolve_w_ts& o) {
                if(this != &o) {
                    ts=o.ts;
                    fx_policy=o.fx_policy;
                    w=o.w;
                    policy=o.policy;
                }
                return *this;
            }

            convolve_w_ts& operator=(convolve_w_ts&& o) {
                ts=std::move(o.ts);
                fx_policy=o.fx_policy;
                w=o.w;
                policy=o.policy;
                return *this;
            }

            //-- useful ct goes here
            template<class A_,class W_>
            convolve_w_ts(A_ && ts,W_ && w, convolve_policy policy=convolve_policy::USE_FIRST)
                :ts(std::forward<A_>(ts)),
                 fx_policy(ts.fx_policy),
                 w(std::forward<W_>(w)),
                 policy(policy)
                  {}

            const ta_t& time_axis() const { return ts.ta;}
            point_interpretation_policy point_interpretation() const { return fx_policy; }
            void set_point_interpretation(point_interpretation_policy point_interpretation) { fx_policy=point_interpretation;}

            point get(size_t i) const {return point(ts.time(i),value(i));}

            // BW compatiblity ?
            size_t size() const { return ts.size();}
            size_t index_of(utctime t) const {return ts.index_of(t);}
            utcperiod total_period() const {return ts.total_period();}
            utctime time(size_t i ) const {return ts.time(i);}

            //--
            double value(size_t i) const {
                double v=0.0;
                for(size_t j=0;j<w.size();++j)
                    v+=  j<=i ?
                            w[j]*ts.value(i-j)
                            :
                            ( policy == convolve_policy::USE_FIRST ? w[j]*ts.value(0) :
                             (policy == convolve_policy::USE_ZERO  ? 0.0 : shyft::nan)); // or nan ? policy based ?
                return  v;
            }
            double operator()(utctime t) const {
                return value(ts.index_of(t));
            }
        };
    }
}
void routing_test::test_hydrograph() {
    using namespace shyft::core;
    using namespace shyft;
    calendar utc;
    utctime t0=utc.time(2016,1,1);
    utctimespan dt=deltahours(1);
    time_axis::fixed_dt ta(t0,dt,24);

    timeseries::point_ts<decltype(ta)> ts(ta,10.0);
    for(size_t i=0;i<5;++i) ts.set(10+i,i);
    std::vector<double> w{0.1,0.15,0.5,0.15,0.1};
    timeseries::convolve_w_ts<decltype(ts),decltype(w)> cts_first(ts,w,timeseries::convolve_policy::USE_FIRST);
    timeseries::convolve_w_ts<decltype(ts),decltype(w)> cts_zero(ts,w,timeseries::convolve_policy::USE_ZERO);
    timeseries::convolve_w_ts<decltype(ts),decltype(w)> cts_nan(ts,w,timeseries::convolve_policy::USE_NAN);
    //std::cout<<"\ni\tts\tfirst\tzero\tnan\n";
    //for(size_t i=0;i<cts_first.size();++i) {
    //    std::cout<<"i:"<<i<<"\t"<<ts.value(i)<<"\t"<<cts_first.value(i) <<"\t"<<cts_zero.value(i)<<"\t"<<cts_nan.value(i)<<"\n";
    //}
    // first policy will just repeat the first value through the filter, thus equal first 4 steps.
    TS_ASSERT_DELTA(ts.value(0),cts_first.value(0),0.0001);
    TS_ASSERT_DELTA(ts.value(1),cts_first.value(1),0.0001);
    TS_ASSERT_DELTA(ts.value(2),cts_first.value(2),0.0001);
    TS_ASSERT_DELTA(ts.value(3),cts_first.value(3),0.0001);

    // zero policy will fill in 0 for the values before 0, -loosing some mass.
    TS_ASSERT_DELTA(cts_zero.value(0),1.0,0.0001);
    TS_ASSERT_DELTA(cts_zero.value(1),2.5,0.0001);
    TS_ASSERT_DELTA(cts_zero.value(2),7.5,0.0001);
    TS_ASSERT_DELTA(cts_zero.value(3),9.0,0.0001);
    TS_ASSERT_DELTA(cts_zero.value(4),10.0,0.0001);

    // nan policy will fill in nan for the values before 0, -loosing some mass. inserting nan on the output
    for (size_t i=0;i+1 <w.size();++i)
        TS_ASSERT(!std::isfinite(cts_nan.value(i)));
    std::vector<double> expected{10,10,10,10,10,10,10,10,10,10,9,7.6,2.85,2.1,2.0,3.5,5.15,8.4,9.4,10,10,10,10,10};
    for(size_t i=4;i<w.size();++i) {
        TS_ASSERT_DELTA(expected[i],cts_first.value(i),0.0001);
        TS_ASSERT_DELTA(expected[i],cts_zero.value(i),0.0001);
        TS_ASSERT_DELTA(expected[i],cts_nan.value(i),0.0001);
    }


}
