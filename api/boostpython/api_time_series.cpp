#include "boostpython_pch.h"

#include "core/utctime_utilities.h"
#include "core/time_axis.h"
#include "core/timeseries.h"
#include "api/api.h"
#include "api/timeseries.h"

namespace expose {
    using namespace shyft;
    using namespace shyft::core;
    using namespace boost::python;
    using namespace std;

    #define DEF_STD_TS_STUFF() \
            .def("point_interpretation",&pts_t::point_interpretation,"returns the point interpretation policy")\
            .def("set_point_interpretation",&pts_t::set_point_interpretation,args("policy"),"set new policy")\
            .def("value",&pts_t::value,args("i"),"returns the value at the i'th time point")\
            .def("time",&pts_t::time,args("i"),"returns the time at the i'th point")\
            .def("get",&pts_t::get,args("t"),"returns the i'th Point ")\
            .def("set",&pts_t::set,args("i","v"),"set the i'th value")\
            .def("fill",&pts_t::fill,args("v"),"all values with v")\
            .def("scale_by",&pts_t::scale_by,args("v"),"scale all values by specified factor")\
            .def("size",&pts_t::size,"returns number of points")\
            .def("index_of",&pts_t::index_of,args("t"),"return the index of the intervall that contains t, or npos if not found")\
            .def("total_period",&pts_t::total_period,"returns the total period covered by the time-axis of this time-series")\
            .def("__call__",&pts_t::operator(),args("t"),"return the f(t) value for the time-series")


    template <class TA>
    static void point_ts(const char *ts_type_name,const char *doc) {
        typedef timeseries::point_ts<TA> pts_t;
        class_<pts_t,bases<>,shared_ptr<pts_t>,boost::noncopyable>(ts_type_name, doc)
            .def(init<const TA&,const vector<double>&,optional<timeseries::point_interpretation_policy>>(args("ta","v","policy"),"constructs a new timeseries from timeaxis and points"))
            .def(init<const TA&,double,optional<timeseries::point_interpretation_policy>>(args("ta","fill_value","policy"),"constructs a new timeseries from timeaxis and fill-value"))
            DEF_STD_TS_STUFF()
            .def_readonly("v",&pts_t::v,"the point vector<double>, same as .values, kept around for backward compatibility")
			.def("get_time_axis", &pts_t::time_axis, "returns the time-axis", return_internal_reference<>()) // have to use func plus init.py fixup due to boost py policy
            ;
    }


    BOOST_PYTHON_MEMBER_FUNCTION_OVERLOADS(point_ts_overloads     ,shyft::api::TsFactory::create_point_ts,4,5);
    BOOST_PYTHON_MEMBER_FUNCTION_OVERLOADS(time_point_ts_overloads,shyft::api::TsFactory::create_time_point_ts,3,4);

    static void TsFactory() {
        class_<shyft::api::TsFactory>("TsFactory","TsFactory is used to create point time-series that exposes the ITimeSeriesOfPoint interface, using the internal ts-implementations")
            .def("create_point_ts",&shyft::api::TsFactory::create_point_ts,point_ts_overloads())//args("n","tStart","dt","values","interpretation"),"returns a new fixed interval ts from specified arguments")
            .def("create_time_point_ts",&shyft::api::TsFactory::create_time_point_ts,time_point_ts_overloads())//args("period","times","values","interpretation"),"return a point ts from specified arguments")
            ;
    }

    static void expose_apoint_ts() {
        typedef shyft::api::apoint_ts pts_t;
        typedef pts_t (pts_t::*self_dbl_t)(double) const;
        typedef pts_t (pts_t::*self_ts_t)(const pts_t &)const;
        typedef  pts_t ( *static_ts_ts_t)(const pts_t&,const pts_t&);
        static_ts_ts_t min_stat_ts_ts_f=&pts_t::min;
        static_ts_ts_t max_stat_ts_ts_f=&pts_t::max;

        self_dbl_t min_double_f=&pts_t::min;
        self_ts_t  min_ts_f =&pts_t::min;
        self_dbl_t max_double_f=&pts_t::max;
        self_ts_t  max_ts_f =&pts_t::max;


		class_<shyft::api::apoint_ts>("Timeseries", "A timeseries providing mathematical and statistical operations and functionality")
			.def(init<const time_axis::generic_dt&, double, optional<timeseries::point_interpretation_policy> >(args("ta", "fill_value", "point_fx"), "construct a timeseries with timeaxis ta and specified fill-value, default point_fx=POINT_INSTANT_VALUE"))
			.def(init<const time_axis::generic_dt&, const std::vector<double>&, optional<timeseries::point_interpretation_policy> >(args("ta", "values", "point_fx"), "construct a timeseries timeaxis ta and corresponding values, default point_fx=POINT_INSTANT_VALUE"))

			.def(init<const time_axis::fixed_dt&, double, optional<timeseries::point_interpretation_policy> >(args("ta", "fill_value", "point_fx"), "construct a timeseries with timeaxis ta and specified fill-value, default point_fx=POINT_INSTANT_VALUE"))
			.def(init<const time_axis::fixed_dt&, const std::vector<double>&, optional<timeseries::point_interpretation_policy> >(args("ta", "values", "point_fx"), "construct a timeseries timeaxis ta and corresponding values, default point_fx=POINT_INSTANT_VALUE"))

			.def(init<const time_axis::point_dt&, double, optional<timeseries::point_interpretation_policy> >(args("ta", "fill_value", "point_fx"), "construct a timeseries with timeaxis ta and specified fill-value, default point_fx=POINT_INSTANT_VALUE"))
			.def(init<const time_axis::point_dt&, const std::vector<double>&, optional<timeseries::point_interpretation_policy> >(args("ta", "values", "point_fx"), "construct a timeseries timeaxis ta and corresponding values, default point_fx=POINT_INSTANT_VALUE"))
            .def(init<const shyft::api::rts_t &>(args("core_result_ts"),"construct a timeseries from a shyft core time-series, to allow full ts-functionality in python"))

			.def(init<const shyft::api::apoint_ts&>(args("clone"), "creates a shallow copy of clone"))

			.def(init<const vector<double>&, utctimespan, const time_axis::generic_dt&>(args("pattern", "dt", "ta"), "construct a timeseries given a equally spaced dt pattern and a timeaxis ta"))
			.def(init<const vector<double>&, utctimespan,utctime, const time_axis::generic_dt&>(args("pattern", "dt","t0", "ta"), "construct a timeseries given a equally spaced dt pattern, starting at t0, and a timeaxis ta"))
			DEF_STD_TS_STUFF()
			// expose time_axis sih: would like to use property, but no return value policy, so we use get_ + fixup in init.py
			.def("get_time_axis", &shyft::api::apoint_ts::time_axis, "returns the time-axis", return_internal_reference<>())
			.add_property("values", &shyft::api::apoint_ts::values, "return the values (possibly calculated on the fly)")
			// operators
			.def(self * self)
			.def(double() * self)
			.def(self * double())

			.def(self + self)
			.def(double() + self)
			.def(self + double())

			.def(self / self)
			.def(double() / self)
			.def(self / double())

			.def(self - self)
			.def(double() - self)
			.def(self - double())

			.def(-self)
			.def("average", &shyft::api::apoint_ts::average, args("ta"), "create a new ts that is the true average of self\n over the specified time-axis ta")
			.def("accumulate", &shyft::api::apoint_ts::accumulate, args("ta"), "create a new ts that is\n the integral f(t) *dt, t0..ti,\n the specified time-axis ta")
			.def("time_shift", &shyft::api::apoint_ts::time_shift,args("delta_t"),
				"create a new ts that is a the time-shift'ed  version of self\n"
				"Parameters\n"
				"----------\n"
				"delta_t : number\n"
				"\t number of seconds to time-shift\n"
				"\t e.g. to move a time-series from 2015 to 2016,\n"
				"\t dt should be number of seconds in 2015\n"
				"Returns\n"
				"-------\n"
				"a new time-series, time-shifted version of self\n"
			)
			.def("convolve_w",&shyft::api::apoint_ts::convolve_w,args("weights","policy"),
                "create a new ts that is the convolved ts with the supporteds weights list"
                "Parameters\n"
                "----------\n"
                "weights : DoubleVector\n"
                "\t the weights profile, use DoubleVector.from_numpy(...) to create these.\n"
                "\t it's the callers responsibility to ensure the sum of weights are 1.0\n"
                "policy : convolve_policy(.USE_FIRST|USE_ZERO|USE_NAN)\n"
                "\t Specifies how to handle initial weight.size()-1 values\n"
                "\t  see also ConvolvePolicy\n"

            )
            .def("min",min_double_f,args("number"),"create a new ts that contains the min of self and number for each time-step")
            .def("min",min_ts_f,args("ts_other"),"create a new ts that contains the min of self and ts_other")
            .def("max",max_double_f,args("number"),"create a new ts that contains the max of self and number for each time-step")
            .def("max",max_ts_f,args("ts_other"),"create a new ts that contains the max of self and ts_other")
            .def("max",max_stat_ts_ts_f,args("ts_a","ts_b"),"create a new ts that is the max(ts_a,ts_b)").staticmethod("max")
            .def("min",min_stat_ts_ts_f,args("ts_a","ts_b"),"create a new ts that is the max(ts_a,ts_b)").staticmethod("min")
			.def("partition_by",&shyft::api::apoint_ts::partition_by,args("calendar","t", "partition_interval", "n_partitions","common_t0"),
				"convert ts to a list of n_partitions partition-ts\n"
				"each partition covers partition_interval, starting from utctime t\n"
				"Parameters\n"
				"----------\n"
				"cal : Calendar\n"
				"\t The calendar to use, typically utc\n"
				"t : utctime\n"
				"\tspecifies where to pick the first partition\n"
				"partition_interval : utctimespan\n"
				"\tthe length of each partition, Calendar.YEAR,Calendar.DAY etc.\n"
				"n_partitions : int\n"
				"\tnumber of partitions\n"
				"common_t0 : utctime\n"
				"\tspecifies the time to correlate all the partitions\n"
				"Returns\n"
				"-------\n"
				"TsVector with len n_partitions"
				)

        ;
        typedef shyft::api::apoint_ts (*avg_func_t)(const shyft::api::apoint_ts&,const shyft::time_axis::generic_dt&);
        avg_func_t avg=shyft::api::average;
		avg_func_t acc = shyft::api::accumulate;
        def("average",avg,args("ts","time_axis"),"creates a true average time-series of ts for intervals as specified by time_axis");
		def("accumulate", acc, args("ts", "time_axis"), "create a new ts that is the integral f(t) *dt, t0..ti, the specified time-axis");
        //def("max",shyft::api::max,(boost::python::arg("ts_a"),boost::python::arg("ts_b")),"creates a new time-series that is the max of the supplied ts_a and ts_b");

        typedef shyft::api::apoint_ts (*ts_op_ts_t)(const shyft::api::apoint_ts&a, const shyft::api::apoint_ts&b);
        typedef shyft::api::apoint_ts (*double_op_ts_t)(double, const shyft::api::apoint_ts&b);
        typedef shyft::api::apoint_ts (*ts_op_double_t)(const shyft::api::apoint_ts&a, double);

        ts_op_ts_t max_ts_ts         = shyft::api::max;
        double_op_ts_t max_double_ts = shyft::api::max;
        ts_op_double_t max_ts_double = shyft::api::max;
        def("max",max_ts_ts    ,args("ts_a","ts_b"),"returns a new ts as max(ts_a,ts_b)");
        def("max",max_double_ts,args("a"   ,"ts_b"),"returns a new ts as max(a,ts_b)");
        def("max",max_ts_double,args("ts_a","b"   ),"returns a new ts as max(ts_a,b)");

        ts_op_ts_t min_ts_ts         = shyft::api::min;
        double_op_ts_t min_double_ts = shyft::api::min;
        ts_op_double_t min_ts_double = shyft::api::min;
        def("min",min_ts_ts    ,args("ts_a","ts_b"),"returns a new ts as min(ts_a,ts_b)");
        def("min",min_double_ts,args("a"   ,"ts_b"),"returns a new ts as min(a,ts_b)");
        def("min",min_ts_double,args("ts_a","b"   ),"returns a new ts as min(ts_a,b)");


        typedef vector<pts_t> TsVector;
        class_<TsVector>("TsVector","A vector of Ts")
            .def(vector_indexing_suite<TsVector>())
            ;

        TsVector (*percentile_1)(const TsVector&,const time_axis::generic_dt&, const std::vector<int>&)=shyft::api::percentiles;
        TsVector (*percentile_2)(const TsVector&,const time_axis::fixed_dt&, const std::vector<int>&)=shyft::api::percentiles;
        const char *percentile_doc="return the percentiles (as TsVector type) (NIST R7, excel,R definition) of the timeseries\n"
            " over the specified time_axis.\n"
            " the time-series point_fx interpretation is used when performing \n"
            " the true-average over the time_axis periods\n"
            " percentiles: 0..100, -1 means arithmetic average,e.g.\n"
            "   [ 0, 25,50,-1,75,100] will return 6 time-series\n"
            "";
        def("percentiles",percentile_1,args("timeseries","time_axis","percentiles"),percentile_doc);
        def("percentiles",percentile_2,args("timeseries","time_axis","percentiles"),percentile_doc);

        def("time_shift", shyft::api::time_shift,args("timeseries","delta_t"),
            "returns a delta_t time-shifted time-series\n"
            " the values are the same as the original,\n"
            " but the time_axis equals the original + delta_t\n");


		/* local scope */ {

			typedef shyft::time_axis::fixed_dt ta_t;
			typedef shyft::timeseries::average_accessor<pts_t, ta_t> AverageAccessorTs;
			class_<AverageAccessorTs>("AverageAccessorTs", "Accessor to get out true average for the time-axis intervals for a point time-series", no_init)
				.def(init<const pts_t&, const ta_t&>(args("ts", "ta"), "construct accessor from ts and time-axis ta"))
				.def(init<shared_ptr<pts_t>, const ta_t&>(args("ts", "ta"), "constructor from ref ts and time-axis ta"))
				.def("value", &AverageAccessorTs::value, args("i"), "returns the i'th true average value")
				.def("size", &AverageAccessorTs::size, "returns number of intervals in the time-axis for this accessor")
				;
		}
    }
	static void expose_correlation_functions() {
		const char * kg_doc =
			"Computes the kling-gupta KGEs correlation for the two time-series over the specified time_axis\n"
			"Parameters\n"
			"----------\n"
			"observed_ts : Timeseries\n"
			"\tthe observed time-series\n"
			"model_ts : Timeseries\n"
			"\t the time-series that is the model simulated / calculated ts\n"
			"time_axis : Timeaxis2\n"
			"\tthe time-axis that is used for the computation\n"
			"s_r : float\n"
			"\tthe kling gupta scale r factor(weight the correlation of goal function)\n"
			"s_a : float\n"
			"\tthe kling gupta scale a factor(weight the relative average of the goal function)\n"
			"s_b : float\n"
			"\tthe kling gupta scale b factor(weight the relative standard deviation of the goal function)\n"
			"Return\n"
			"------\n"
			" float: The  KGEs= 1-EDs that have a maximum at 1.0\n";

		def("kling_gupta", shyft::api::kling_gupta, args("observation_ts", "model_ts", "time_axis", "s_r", "s_a", "s_b"),
			kg_doc
		);

		const char *ns_doc = "Computes the Nash-Sutcliffe model effiency coefficient (n.s) \n"
			" for the two time-series over the specified time_axis\n"
			" Ref:  http://en.wikipedia.org/wiki/Nash%E2%80%93Sutcliffe_model_efficiency_coefficient \n"

			"Parameters\n"
			"----------\n"
			"observed_ts : Timeseries\n"
			"\tthe observed time-series\n"
			"model_ts : Timeseries\n"
			"\t the time-series that is the model simulated / calculated ts\n"
			"time_axis : Timeaxis2\n"
			"\tthe time-axis that is used for the computation\n"
			"Return\n"
			"------\n"
			" float: The  n.s performance, that have a maximum at 1.0\n";

		def("nash_sutcliffe", shyft::api::nash_sutcliffe, args("observation_ts", "model_ts", "time_axis"),
			ns_doc
		);


	}
	static void expose_periodic_ts() {
		const char *docstr =
			"Create a Timeseries by repeating the pattern-specification\n"
			"Parameters\n"
			"----------\n"
			"pattern : DoubleVector\n"
			"\tthe value-pattern as a sequence of values\n"
			"dt : int\n"
			"\tnumber of seconds between the pattern values, e.g. deltahours(3)\n"
			"t0 : utctime\n"
			"\tspecifies the starttime of the pattern\n"
			"ta : Timeaxis\n"
			"\tthe time-axis for which the pattern is repeated\n"
			"\t e.g. your pattern might be 8 3h values,and you could supply\n"
			"\ta time-axis 'ta' at hourly resolution\n"
			;
		def("create_periodic_pattern_ts", shyft::api::create_periodic_pattern_ts, args("pattern","dt","t0","ta"), docstr);

	}
    void timeseries() {
        enum_<timeseries::point_interpretation_policy>("point_interpretation_policy")
            .value("POINT_INSTANT_VALUE",timeseries::POINT_INSTANT_VALUE)
            .value("POINT_AVERAGE_VALUE",timeseries::POINT_AVERAGE_VALUE)
            .export_values()
            ;
        enum_<timeseries::convolve_policy>(
            "convolve_policy",
            "Ref Timeseries.convolve_w function, this policy determinte how to handle initial conditions\n"
            "USE_FIRST: value(0) is used for all values before value(0), 'mass preserving'\n"
            "USE_ZERO : fill in zero for all values before value(0):shape preserving\n"
            "USE_NAN  : nan filled in for the first length-1 values of the filter\n"
            )
            .value("USE_FIRST",timeseries::convolve_policy::USE_FIRST)
            .value("USE_ZERO",timeseries::convolve_policy::USE_ZERO)
            .value("USE_NAN",timeseries::convolve_policy::USE_NAN)
            .export_values()
            ;
        class_<timeseries::point> ("Point", "A timeseries point specifying utctime t and value v")
            .def(init<utctime,double>(args("t","v")))
            .def_readwrite("t",&timeseries::point::t)
            .def_readwrite("v",&timeseries::point::v)
            ;
        point_ts<time_axis::fixed_dt>("TsFixed","A time-series with a fixed delta t time-axis");
        point_ts<time_axis::point_dt>("TsPoint","A time-series with a variable delta time-axis");
        TsFactory();
        expose_apoint_ts();
		expose_periodic_ts();
		expose_correlation_functions();
    }
}
