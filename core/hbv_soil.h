#pragma once
#include "utctime_utilities.h"
namespace shyft {
	namespace core {
		namespace hbv_soil {
			using namespace std;

			struct parameter {
				parameter(double fc = 300.0, double beta = 2.0, double lp = 150) 
					:fc(fc), beta(beta), lp(lp) {
					if (fc < .0) 
						throw runtime_error("fc should be > 0.0");
				}
				double fc=300; // mm
				double beta=2.0;// unit-less
				double lp = 150;// mm
			};

			struct state {
				state(double sm = 0.0) :sm(sm) {}
				double sm=50.0; // mm
				bool operator==(const state&x) const {
					const double eps = 1e-6;
					return fabs(sm - x.sm)<eps;
				}
			};

			struct response {
				double outflow = 0.0; // mm
			};


			/** \brief Hbv soil
			*
			*  reference somewhere..
			*
			* \tparam P Parameter type, implementing the interface:
			* \tparam S State type, implementing the interface:
			* \tparam R Respone type, implementing the interface:
			*/
			template<class P>
			struct calculator {
				P param;
				calculator(const P& p):param(p) {}
				template <class R,class S> 
				void step(S& s, R& r, shyft::core::utctime t0, shyft::core::utctime t1, double insoil, double pot_evap, double act_evap) {
					double fraction = pow(s.sm / param.fc, param.beta);
					r.outflow = fraction*insoil;
					s.sm = s.sm - r.outflow + ((1-fraction)*insoil  - (s.sm < param.lp? act_evap:pot_evap ) );
				}
			};
		}
	} // core
} // shyft