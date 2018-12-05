/*********************************************************************************
  *FileName: roulette_wheel_select.cc
  *Author:  Tengyu
  *Date:  2018/11/19 
  *Description:  Roulette Wheel Selection
**********************************************************************************/
#include "roulette_wheel_select.hh"
#include <numeric>
#include <numeric/random/random.hh>
#include <string>
namespace protocols {
	namespace abinitio {
		roulette_wheel_select::roulette_wheel_select(std::vector<double> fitness_)
			{
				fitness_in = fitness_;
				sum_items_in = 0; 
				for(std::vector<double>::iterator it = fitness_in.begin(); it != fitness_in.end(); it ++)
					{
						sum_items_in = sum_items_in + *it;
					}
			}
		int roulette_wheel_select::selection()
			{
				double rand_num = ( numeric::random::rg().uniform() * sum_items_in );
				double P = 0;
				unsigned int index = 0;
				do{
					if( index <= fitness_in.size() )
						{
							P = P + fitness_in[index];
							index = index + 1;
						}
					else
						{
							std::cout << "roulette_wheel_select:access violation\n";
							return false;
						}
				}while( P < rand_num);
				return index - 1;
			}	
		} // abinitio
	} // protocols