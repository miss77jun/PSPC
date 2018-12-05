/*********************************************************************************
  *FileName: roulette_wheel_select.hh
  *Author:  Tengyu
  *Date:  2018/11/19 
  *Description:  Roulette Wheel Selection
  							 Calculate S = the sum of a finesses.
  							 Generate a random number R between 0 and S.
  							 Starting from the top of the population, keep adding the finesses to the partial sum P, till P<R.
  							 The individual for which P exceeds R is the chosen individual.
**********************************************************************************/

#ifndef ROULETTE_WHEEL_SELECT_HH
#define ROULETTE_WHEEL_SELECT_HH

#include <string>
#include <vector>


namespace protocols {
	namespace abinitio {
		class roulette_wheel_select {
			public:
				roulette_wheel_select(std::vector<double> fitness_);
				int selection();
			private:
				double sum_items_in;
				std::vector<double> fitness_in;
		};
	} // abinitio
} // protocols
#endif

