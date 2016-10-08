#pragma once
#include "geo_point.h"
namespace shyft {
    namespace core {

        /** \brief LandTypeFractions are used to describe 'type of land'
         *   like glacier, lake, reservoir and forest.
         *
         *   It is designed as a part of GeoCellData (could be nested, but
         *   we take python/swig into consideration). It's of course
         *   questionable, since we could have individual models for each type
         *   of land, - but current approach is to use a cell-area, and then
         *   describe fractional properties.
         */
        struct land_type_fractions {
            land_type_fractions()
             : glacier_(0), lake_(0), reservoir_(0), forest_(0) {}

            /** \brief LandTypeFraction constructor
             *
             *  Construct the land type fractions based on the sizes of the
             *  different land types. Each size should be greater or equal to
             *  zero.
             */
            land_type_fractions(double glacier_size, double lake_size,
                               double reservoir_size, double forest_size,
                               double unspecified_size) {
                // TODO: Is this needed?
                glacier_size = std::max(0.0, glacier_size);
                lake_size = std::max(0.0, lake_size);
                reservoir_size = std::max(0.0, reservoir_size);
                forest_size = std::max(0.0, forest_size);
                unspecified_size = std::max(0.0, unspecified_size);
                const double sum = glacier_size + lake_size + reservoir_size
                                   + forest_size + unspecified_size;
                if (sum > 0) {
                    glacier_ = glacier_size/sum;
                    lake_ = lake_size/sum;
                    reservoir_ = reservoir_size/sum;
                    forest_ = forest_size/sum;
                }
                else
                    glacier_ = lake_ = reservoir_ = forest_ = 0.0;
            }
            double glacier() const { return glacier_; }
			double lake() const { return lake_; }         // not regulated, assume time-delay until discharge
			double reservoir()const{ return reservoir_; } // regulated, assume zero time-delay to discharge
			double forest() const { return forest_; }
			double unspecified() const { return 1.0 - glacier_ - lake_ - reservoir_ - forest_; }
            void set_fractions(double glacier, double lake, double reservoir, double forest) {
                const double tol = 1.0e-3; // Allow small round off errors
                const double sum = glacier + lake + reservoir + forest;
                if (sum > 1.0 && sum < 1.0 + tol) {
                    glacier /= sum;
                    lake /= sum;
                    reservoir /= sum;
                    forest /= sum;
                } else if (sum > 1.0  ||
                          (glacier < 0.0 || lake < 0.0 || reservoir < 0.0 || forest < 0.0))
                   throw std::invalid_argument("LandTypeFractions:: must be >=0.0 and sum <= 1.0");
                glacier_ = glacier; lake_ = lake; reservoir_ = reservoir; forest_ = forest;
            }
          private:
			double glacier_;
			double lake_;   // not regulated, assume time-delay until discharge
			double reservoir_;// regulated, assume zero time-delay to discharge
			double forest_;
        };

        /** The routing_info contains the geo-static parts of the relation between
         * the cell and the routing sink point */
        struct routing_info {
            int id=0; ///< target routing input identifier (similar to catchment_id), 0 means nil,none
            double distance=1.0; ///< static routing distance[m] to the routing point
            double velocity=1.0; ///< static 'velocity' parameter[m/s], s= v*t, t=s/v, t= avg. travel time[s]
        };

        const double default_radiation_slope_factor=0.9;

        /** \brief geo_cell_data represents common constant geo_cell properties across several possible models and cell assemblies.
         *  The idea is that most of our algorithms uses one or more of these properties,
         *  so we provide a common aspect that keeps this together.
         *  Currently it keep the
         *   - mid-point geo_point, (x,y,z) (default 0)
         *     the area in m^2, (default 1000 x 1000 m2)
         *     land_type_fractions (unspecified=1)
         *     catchment_id   def (-1)
         *     radiation_slope_factor def 0.9
		 */
		struct geo_cell_data {
		    static const int default_area_m2=1000000;
			geo_cell_data() :catchment_ix(0),area_m2(default_area_m2),catchment_id_(-1),radiation_slope_factor_(default_radiation_slope_factor){}

			geo_cell_data(geo_point mid_point,double area=default_area_m2,
                int catchment_id = -1, double radiation_slope_factor=default_radiation_slope_factor,const land_type_fractions& land_type_fractions=land_type_fractions()):
				mid_point_(mid_point), area_m2(area), catchment_id_(catchment_id),radiation_slope_factor_(radiation_slope_factor),fractions(land_type_fractions)
			{}
			const geo_point& mid_point() const { return mid_point_; }
			size_t catchment_id() const { return catchment_id_; }
			void set_catchment_id(size_t catchmentid) {catchment_id_=catchmentid;}
			double radiation_slope_factor() const { return radiation_slope_factor_; }
			const land_type_fractions& land_type_fractions_info() const { return fractions; }
			void set_land_type_fractions(const land_type_fractions& ltf) { fractions = ltf; }
			double area() const { return area_m2; }
			size_t catchment_ix; // internally generated zero-based catchment index, used to correlate to calc-filter, ref. region_model
            routing_info routing;///< keeps the geo-static routing info, where it routes to, and routing distance.
		  private:

			geo_point mid_point_; // midpoint
			double area_m2; //m2
			size_t catchment_id_;
			double radiation_slope_factor_;
			land_type_fractions fractions;
			// geo-type  parts, interesting for some/most response routines, sum fractions should be <=1.0
		};
    }
} // shyft
