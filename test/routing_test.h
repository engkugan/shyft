#pragma once

#include <cxxtest/TestSuite.h>

class routing_test: public CxxTest::TestSuite {
  public:
    void test_hydrograph();
    void test_routing_model();
};
