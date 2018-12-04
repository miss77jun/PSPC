//miss77
//contact
#include <iostream>
#include <fstream>
#include <iostream>
#include <string>
#include <sstream>
#include "contact_predicted.hh"
#include<cmath>
#include <numeric>
#include <numeric/xyzVector.hh>
#include <core/chemical/ResidueTypeSet.fwd.hh>
#include <core/chemical/ChemicalManager.hh>
#include <core/util/SwitchResidueTypeSet.hh>
#include <core/pose/Pose.hh>
#include <core/conformation/Residue.hh>

using namespace std;
namespace protocols {
namespace abinitio {
	contact_predicted::contact_predicted( int c_length )
	{
		length_size = c_length;
	}
	bool contact_predicted::read_contact()
	{
		ifstream contactmap_file("../inputFolder/contactmap.txt");
		string input;
		bool error ; 
		
		for( int i = 0; i < length_size; )
		{
			error = false;
			
			getline(contactmap_file,input);
			if (!contactmap_file) break; //check for eof
			
			int first_col = (int)input[0];
			if (first_col < 48 || first_col > 57 ) 
			{
				continue;
			}
			istringstream buffer(input);
			int index1,index2,zero,eight;
			double probability;
			
			buffer >> index1 >> index2 >> zero >> eight >> probability;
			//cout << "set: " << index1 << "\t" << index2 << "\t" << probability << "\n";
			
			//set index of contact pairs
			vector<int> contact_index;
			contact_index.push_back( index1 );
			contact_index.push_back( index2 );
			contact_index_in.push_back( contact_index );
			
			//set probability
			probability_in.push_back( probability );
			
			string str;
			while(buffer >> str){
				cout << str<<"\n";
			}
			//check for non numerical input
			//and less/more input than needed
//			if ( !buffer || !buffer.eof()) 
//			{
//				error = true;
//				break;
//			}
			
			i ++;
		}
		if (error)
				cout << "contactmap file  is corrupted..." << endl;
		else
			  cout << "Read contactmap successfully!!!" << endl; 
		contactmap_file.close();
		
		return 1;
	}
	vector<vector<int> > contact_predicted::get_contact_index()
	{
		return contact_index_in;
	}
	vector<int>  contact_predicted::get_contact_index(unsigned int index_)
	{
		if( index_ <= contact_index_in.size()  )
			{
				return contact_index_in[index_];
			}
		else
			{
				std::cout << "contact_predicted.cc: contact_index_in access violation!!!\n";
				return *( contact_index_in.end() );
			}
	}
	//score with the probability of contact
	double contact_predicted::score_C( core::pose::Pose &pose)
	{
		double score = 0;
		for(int i = 0; i < length_size ; i ++)
		{
			score = score + score_one_C( pose , i);
		}
		return score;
	}
	//score one contact i
	double contact_predicted::score_one_C( core::pose::Pose &pose , int i)
	{
		//position of contact i
		int pos1 = contact_index_in[i][0];
		int pos2 = contact_index_in[i][1];
		double prob = probability_in[i];
		
		core::Real dis = get_distance_paired( pose, pos1, pos2 );
		//double score = -(( 1 - dis) / ( 1 + exp(-(dis - 8))) + dis / ( 1 + exp( -( 8 - dis))));
		double score;
		if(dis < 20)
			{
				score = 1 - prob + (2 * prob - 1) / (1 + exp( 8 -dis ));
			}
		else
			{
//				score = 1 - prob + (2 * prob - 1) / (1 + exp(8 - 20)) + 5 * log(dis / 20);
					score = 1 - prob + (2 * prob - 1) / (1 + exp(8 - 20)) + exp(dis / 20) - exp(1);
			}
		return score;
	}
	core::Real contact_predicted::get_distance_paired( core::pose::Pose &pose , core::Size pos1, core::Size pos2 )
	{
			std::string atom_name1 = pose.residue_type( pos1 ).name1() == 'G' ? "CA": "CB";
			std::string atom_name2 = pose.residue_type( pos2 ).name1() == 'G' ? "CA": "CB";
			numeric::xyzVector<core::Real> v1 = pose.residue( pos1 ).atom( atom_name1 ).xyz();
			numeric::xyzVector<core::Real> v2 = pose.residue( pos2 ).atom( atom_name2 ).xyz();
			core::Real distance = v1.distance(v2);
			return distance;
	}
	vector<double> contact_predicted::get_prob()
		{
			return probability_in;
		}
} //abinitio
} //protocols