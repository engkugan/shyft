
#include "test_pch.h"
#include "hbv_soil_test.h"
#include "core/hbv_soil.h"


using namespace shyft::core;

void hbv_soil_test::test_regression() {
	hbv_soil::parameter p;
	hbv_soil::calculator<hbv_soil::parameter> calc(p);
	hbv_soil::state s;
	hbv_soil::response r;
	calendar utc;
	utctime t0 = utc.time(2015, 1, 1);
	auto previous_value = s.sm;
	calc.step(s, r, t0, t0 + deltahours(1), 0, 0);
	TS_ASSERT_DELTA(r.outflow, 0.0, 0.0); //TODO: verify some more numbers
	TS_ASSERT_DELTA(s.sm, previous_value, 0.0);
	calc.step(s, r, t0, t0 + deltahours(1), 50.0, 0);
	TS_ASSERT_DELTA(s.sm, 97.22222222, 0.0002);
	TS_ASSERT_DELTA(r.outflow, 1.38888000, 0.05);
}
