//Tengyu
//


#ifndef CONTACT_PREDICTION_hh
#define CONTACT_PREDICTION_hh

#include <string>
#include <vector>
#include <core/types.hh>
#include <core/pose/Pose.fwd.hh>
#include <core/chemical/ResidueTypeSet.fwd.hh>

using namespace std;
namespace protocols {
	namespace abinitio {
	
		class contact_predicted {
			public:
				contact_predicted( int c_length );
				bool read_contact();
				vector<vector<int> > get_contact_index();
				vector<int>  get_contact_index(unsigned int index_);
				double score_C( core::pose::Pose &pose);
				vector<double> get_prob();
			private:
				double score_one_C( core::pose::Pose &pose , int i);
				core::Real get_distance_paired( core::pose::Pose &pose , core::Size pos1, core::Size pos2 );
				
				
				vector<vector<int> > contact_index_in;
				vector<double> probability_in;
				int length_size;//length of contacts
		};
	} // abinitio
} // protocols
#endif