// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file ClassicAbinitio.cc
/// @brief ab-initio fragment assembly protocol for proteins
/// @details
///   Contains currently: Classic Abinitio
///
///
/// @author Oliver Lange
/// @author James Thompson
/// @author Mike Tyka

// Unit Headers
#include <protocols/abinitio/ClassicAbinitio.hh>
#include <protocols/simple_moves/SymmetricFragmentMover.hh>

// Package Headers
#include <protocols/simple_moves/GunnCost.hh>

// Project Headers
#include <core/pose/Pose.hh>
#include <core/kinematics/MoveMap.hh>
#include <core/types.hh>
#include <core/scoring/ScoreFunction.fwd.hh>
#include <core/scoring/ScoreType.hh>
#include <core/scoring/ScoreFunctionFactory.hh>
#include <core/scoring/rms_util.hh>
#include <core/scoring/constraints/ConstraintSet.hh>
#include <core/scoring/constraints/ConstraintSet.fwd.hh>
#include <core/scoring/dssp/Dssp.hh>

#include <protocols/moves/Mover.hh>
#include <protocols/simple_moves/BackboneMover.hh>
#include <protocols/moves/MoverContainer.hh>
#include <protocols/moves/TrialMover.hh>
#include <protocols/moves/RepeatMover.hh>
#include <protocols/moves/WhileMover.hh>
#include <protocols/abinitio/AllResiduesChanged.hh>

//for cenrot
#include <protocols/moves/CompositionMover.hh>
#include <core/pack/task/PackerTask.hh>
#include <core/pack/task/TaskFactory.hh>
#include <core/pack/task/operation/TaskOperations.hh>
#include <protocols/simple_moves/PackRotamersMover.hh>
#include <basic/options/keys/corrections.OptionKeys.gen.hh>
#include <protocols/simple_moves/SwitchResidueTypeSetMover.hh>
//#include <protocols/simple_moves/BackboneMover.hh>


// ObjexxFCL Headers
#include <ObjexxFCL/string.functions.hh>

// Utility headers
#include <utility/exit.hh>
#include <utility/vector1.fwd.hh>
#include <utility/pointer/ReferenceCount.hh>
#include <utility/io/ozstream.hh>
#include <numeric/numeric.functions.hh>
#include <basic/prof.hh>
#include <basic/Tracer.hh>
#include <basic/options/option.hh>
#include <basic/options/keys/abinitio.OptionKeys.gen.hh>
#include <basic/options/keys/run.OptionKeys.gen.hh>
#include <basic/options/keys/templates.OptionKeys.gen.hh>

//// C++ headers
#include <cstdlib>
#include <string>
#ifdef WIN32
#include <ctime>
#endif

//debug

#include <protocols/moves/MonteCarlo.hh>
#include <utility/vector0.hh>
#include <utility/vector1.hh>

#include<map>
#include<cstdlib>
#include<algorithm>
#include<math.h>

//@miss77

#include<map>
#include<cstdlib>
#include<algorithm>
#include<math.h>

//@xty
#include <core/io/silent/SilentFileData.hh>
#include <core/io/silent/silent.fwd.hh>
#include <core/io/silent/SilentStructFactory.hh>
#include <core/io/silent/SilentStruct.hh>
#include <ObjexxFCL/string.functions.hh>
#include <core/io/silent/SilentFileOptions.hh>
#include <core/import_pose/import_pose.hh>
#include <basic/options/keys/in.OptionKeys.gen.hh>
#include <core/scoring/Energies.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/chemical/ResidueTypeSet.fwd.hh>
#include <core/chemical/ChemicalManager.hh>
#include <core/util/SwitchResidueTypeSet.hh>
#include <ObjexxFCL/FArray1D.hh>
#include <ObjexxFCL/format.hh>
#include <numeric/random/random.hh>
#include <numeric>

#include <basic/options/keys/in.OptionKeys.gen.hh>
#include <basic/options/keys/out.OptionKeys.gen.hh>
#include <basic/options/keys/cluster.OptionKeys.gen.hh>
#include <basic/options/keys/symmetry.OptionKeys.gen.hh>

#include <utility/vector1.hh>
#include <utility/excn/Exceptions.hh>


#if defined(WIN32) || defined(__CYGWIN__)
#include <ctime>
#endif
#include <core/scoring/symmetry/SymmetricScoreFunction.hh>
#include <protocols/cluster/cluster.hh>
#include <core/import_pose/pose_stream/MetaPoseInputStream.hh>
#include <core/import_pose/pose_stream/util.hh>
#include <core/import_pose/pose_stream/PDBPoseInputStream.hh>
#include <core/import_pose/pose_stream/PDBPoseInputStream.fwd.hh>
#include <core/chemical/ResidueTypeSet.fwd.hh>
#include <core/chemical/ChemicalManager.hh>
#include <core/import_pose/pose_stream/MetaPoseInputStream.hh>
#include <core/import_pose/pose_stream/SilentFilePoseInputStream.hh>
#include <core/import_pose/pose_stream/LazySilentFilePoseInputStream.hh>
#include <core/import_pose/pose_stream/SilentFilePoseInputStream.fwd.hh>
#include <utility/io/izstream.hh>

#include <utility/vector1.hh>
#include <core/scoring/ScoringManager.hh>
#include <core/scoring/Ramachandran.hh>
#include <string>
#include <iostream>
#include <stdlib.h>
#include <stdio.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <string.h>
#include <malloc.h>
#include <unistd.h>


#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <iomanip>
#include <map>

//@xty
#include <core/io/silent/SilentFileData.hh>
#include <core/io/silent/silent.fwd.hh>
#include <core/io/silent/SilentStructFactory.hh>
#include <core/io/silent/SilentStruct.hh>
#include <ObjexxFCL/string.functions.hh>
#include <core/io/silent/SilentFileOptions.hh>

#include <numeric/random/random.hh>
#include<numeric>
#include <sstream>

#include <time.h> 
#include <numeric/xyzVector.hh>
#include "roulette_wheel_select.hh"

//miss77
bool flag_stage2to4 , flag_stage5 ;
core::pose::Pose  native_pose_temp;

using namespace std;


static THREAD_LOCAL basic::Tracer tr( "protocols.abinitio" );

using core::Real;
using namespace core;
using namespace basic;
using namespace basic::options;
using namespace basic::options::OptionKeys;

/*!
@detail call this:
ClassicAbinitio::register_options() before devel::init().
Derived classes that overload this function should also call Parent::register_options()
*/

// This method of adding options with macros is a pain in the ass for people
// trying to nest ClassicAbinitio as part of other protocols. If you don't call
// ClassicAbinitio::register_options() in your main function, you get a really
// unintuitive segfault as the options system doesn't know about the options
// listed below. The solution is to call register_options() in your main method
// before devel::init(), which is really ugly as the main method shouldn't need
// to know what protocols are called, and it's prone to error because it's an
// easy thing to forget.
// This should get some more thought before it becomes the standard way to add options.

void protocols::abinitio::ClassicAbinitio::register_options() {
	Parent::register_options();
	using namespace basic::options;
	using namespace basic::options::OptionKeys;
	option.add_relevant( OptionKeys::abinitio::increase_cycles );
	option.add_relevant( OptionKeys::abinitio::smooth_cycles_only );
	option.add_relevant( OptionKeys::abinitio::debug );
	option.add_relevant( OptionKeys::abinitio::skip_convergence_check );
	option.add_relevant( OptionKeys::abinitio::log_frags );
	option.add_relevant( OptionKeys::abinitio::only_stage1 );
	option.add_relevant( OptionKeys::abinitio::end_bias );
	option.add_relevant( OptionKeys::abinitio::symmetry_residue );
	option.add_relevant( OptionKeys::abinitio::vdw_weight_stage1 );
	option.add_relevant( OptionKeys::abinitio::override_vdw_all_stages );
	option.add_relevant( OptionKeys::abinitio::recover_low_in_stages );
	option.add_relevant( OptionKeys::abinitio::close_chbrk );
}


namespace protocols {
namespace abinitio {

//little helper function
bool contains_stageid( utility::vector1< ClassicAbinitio::StageID > vec, ClassicAbinitio::StageID query ) {
	return find( vec.begin(), vec.end(), query) != vec.end();
}

/// @detail  large (stage1/stage2)
/// small(stage2/stage3/stage4)
/// smooth_small ( stage3/stage4)
ClassicAbinitio::ClassicAbinitio(
	simple_moves::FragmentMoverOP brute_move_small,
	simple_moves::FragmentMoverOP brute_move_large,
	simple_moves::FragmentMoverOP smooth_move_small,
	int  /*dummy otherwise the two constructors are ambiguous */
) :
	brute_move_small_( brute_move_small ),
	brute_move_large_( brute_move_large ),
	smooth_move_small_( smooth_move_small )
{
	BaseClass::type( "ClassicAbinitio" );
	get_checkpoints().set_type("ClassicAbinitio");
	// std::cerr << "ClassicAbinitio::constructor has stubbed out...(fatal) see code file";
	// runtime_assert( 0 ); //---> needs more implementation to use this constructor: e.g. read out movemap from FragmentMover...
	movemap_ = brute_move_large->movemap();
	//  set_defaults( pose ); in constructor virtual functions are not called
	bSkipStage1_ = false;
	bSkipStage2_ = false;

	close_chbrk_ = false;

	stage4_cycles_pack_rate_ = 0.25;
}

ClassicAbinitio::ClassicAbinitio(
	core::fragment::FragSetCOP fragset_small,
	core::fragment::FragSetCOP fragset_large,
	core::kinematics::MoveMapCOP movemap
)  :
	movemap_( movemap )
{
	BaseClass::type( "ClassicAbinitio" );
	get_checkpoints().set_type("ClassicAbinitio");
	using namespace basic::options;
	simple_moves::ClassicFragmentMoverOP bms, bml, sms;
	using simple_moves::FragmentCostOP;
	using simple_moves::ClassicFragmentMover;
	using simple_moves::SymmetricFragmentMover;
	using simple_moves::SmoothFragmentMover;
	using simple_moves::SmoothSymmetricFragmentMover;
	using simple_moves::GunnCost;
	if ( option[ OptionKeys::abinitio::log_frags ].user() ) {
		if ( !option[ OptionKeys::abinitio::debug ] ) utility_exit_with_message( "apply option abinitio::log_frags always together with abinitio::debug!!!");
		bms = simple_moves::ClassicFragmentMoverOP( new simple_moves::LoggedFragmentMover( fragset_small, movemap ) );
		bml = simple_moves::ClassicFragmentMoverOP( new simple_moves::LoggedFragmentMover( fragset_large, movemap ) );
		sms = simple_moves::ClassicFragmentMoverOP( new SmoothFragmentMover( fragset_small, movemap, FragmentCostOP( new GunnCost ) ) );
	} else if ( option[ OptionKeys::abinitio::symmetry_residue ].user() ) {
		Size const sr (  option[ OptionKeys::abinitio::symmetry_residue ] );
		bms = simple_moves::ClassicFragmentMoverOP( new SymmetricFragmentMover( fragset_small, movemap, sr ) );
		bml = simple_moves::ClassicFragmentMoverOP( new SymmetricFragmentMover( fragset_large, movemap, sr ) );
		sms = simple_moves::ClassicFragmentMoverOP( new SmoothSymmetricFragmentMover( fragset_small, movemap, FragmentCostOP( new GunnCost ), sr ) );
	} else {
		bms = simple_moves::ClassicFragmentMoverOP( new ClassicFragmentMover( fragset_small, movemap ) );
		bml = simple_moves::ClassicFragmentMoverOP( new ClassicFragmentMover( fragset_large, movemap ) );
		sms = simple_moves::ClassicFragmentMoverOP( new SmoothFragmentMover ( fragset_small, movemap, FragmentCostOP( new GunnCost ) ) );
	}

	bms->set_end_bias( option[ OptionKeys::abinitio::end_bias ] ); //default is 30.0
	bml->set_end_bias( option[ OptionKeys::abinitio::end_bias ] );
	sms->set_end_bias( option[ OptionKeys::abinitio::end_bias ] );

	brute_move_small_ = bms;
	brute_move_large_ = bml;
	smooth_move_small_ = sms;

	using namespace core::pack::task;
	//init the packer
	pack_rotamers_ = simple_moves::PackRotamersMoverOP( new protocols::simple_moves::PackRotamersMover() );
	TaskFactoryOP main_task_factory( new TaskFactory );
	main_task_factory->push_back( operation::TaskOperationCOP( new operation::RestrictToRepacking ) );
	//main_task_factory->push_back( new operation::PreserveCBeta );
	pack_rotamers_->task_factory(main_task_factory);

	bSkipStage1_ = false;
	bSkipStage2_ = false;

	close_chbrk_ = false;

	stage4_cycles_pack_rate_ = 0.25;
}

/// @details Call parent's copy constructor and perform a shallow
/// copy of all the data.  NOTE: Shallow copy is only to preserve
/// behavior pre 9/7/2009 when the compiler-provided copy constructor
/// was being invoked.
ClassicAbinitio::ClassicAbinitio( ClassicAbinitio const & src ) :
	//utility::pointer::ReferenceCount(),
	Parent( src )
{
	stage1_cycles_ = src.stage1_cycles_;
	stage2_cycles_ = src.stage2_cycles_;
	stage3_cycles_ = src.stage3_cycles_;
	stage4_cycles_ = src.stage4_cycles_;
	stage5_cycles_ = src.stage5_cycles_;
	score_stage1_ = src.score_stage1_;
	score_stage2_ = src.score_stage2_;
	score_stage3a_ = src.score_stage3a_;
	score_stage3b_ = src.score_stage3b_;
	score_stage4_ = src.score_stage4_;
	score_stage4rot_ = src.score_stage4rot_;
	score_stage5_ = src.score_stage5_;
	apply_large_frags_ = src.apply_large_frags_;
	short_insert_region_ = src.short_insert_region_;
	just_smooth_cycles_ = src.just_smooth_cycles_;
	bQuickTest_ = src.bQuickTest_;
	close_chbrk_ = src.close_chbrk_;
	temperature_ = src.temperature_;
	movemap_ = src.movemap_;
	mc_ = src.mc_;
	brute_move_small_ = src.brute_move_small_;
	brute_move_large_ = src.brute_move_large_;
	smooth_move_small_ = src.smooth_move_small_;
	trial_large_ = src.trial_large_;
	trial_small_ = src.trial_small_;
	smooth_trial_small_ = src.smooth_trial_small_;
	total_trials_ = src.total_trials_;
	bSkipStage1_ = src.bSkipStage1_;
	bSkipStage2_ = src.bSkipStage2_;
	bSkipStage3_ = src.bSkipStage3_;
	bSkipStage4_ = src.bSkipStage4_;
	bSkipStage5_ = src.bSkipStage5_;
	recover_low_stages_ = src.recover_low_stages_;
}

/// @brief Explicit destructor is needed to destroy all the OPs
/// The compiler does all the work, but it requires that we place
/// the destructor in the .cc file.
ClassicAbinitio::~ClassicAbinitio()
{}

/// @brief setup moves, mc-object, scores
/// @details can't call this from constructor; virtual functions don't operate until construction has completed.

void
ClassicAbinitio::init( core::pose::Pose const& pose ) {
	// Parent::init( pose );
	set_defaults( pose );
	// bInitialized_ = true;
}

/// @brief ClassicAbinitio has virtual functions... use this to obtain a new instance
moves::MoverOP
ClassicAbinitio::clone() const
{
	return moves::MoverOP( new ClassicAbinitio( *this ) );
}

void ClassicAbinitio::apply( pose::Pose & pose ) {
	using namespace moves;
	using namespace scoring;
	using namespace scoring::constraints;

	Parent::apply( pose );
	if ( option[ OptionKeys::run::dry_run ]() ) return;
	
	//miss77 output straight chain
	//pose.dump_pdb( "initialization.pdb" );
	

	//basic::prof_reset();
	//get native strucuture
	core::import_pose::pose_from_file( native_pose_temp,option[ in::file::native ](),core::import_pose::PDB_file);
	core::util::switch_to_residue_type_set( native_pose_temp,chemical::CENTROID,true);

	
	const char *native_ss;
  core::scoring::dssp::Dssp dssp_native( native_pose_temp );
  std::string pose_ss = dssp_native.get_dssp_secstruct();
  native_ss = pose_ss.c_str();
	

	bool success( true );
	total_trials_ = 0;
	//contact from RaptorX: top L
	protocols::abinitio::contact_predicted RaptorXcontact( pose.total_residue() );
	RaptorXcontact.read_contact();
	ofstream max_distance_file("max_distance.csv");
	vector<vector<int> > contact_index = RaptorXcontact.get_contact_index();
	for(unsigned int i = 0; i < contact_index.size(); i ++)
	{
		max_distance_file << contact_index[i][0] << "-" <<contact_index[i][1] << "," ;
		cout << "contact_index: "<< contact_index[i][0] << "\t" << contact_index[i][1] << "\n";
	}
	max_distance_file << "\n";
	//output the distance of native structure
	if( native_pose_temp.total_residue() ==  pose.total_residue() )
	{
			for(unsigned int i = 0; i < contact_index.size(); i ++)
		{
			Real distance_native =  get_distance_paired( native_pose_temp , contact_index[i][0], contact_index[i][1]);
			max_distance_file << distance_native << "," ;
		}
  }
  else max_distance_file << "The lengths of target sequence and  experimental strucutre are different!!!" << "," ;
	max_distance_file << "\n";
	max_distance_file.close();
	
	//get status_info
	//start time
	time_t rawtime;
  struct tm * timeinfo;
  time ( &rawtime );
  timeinfo = localtime ( &rawtime );
  printf ( "The current date/time is: %s", asctime (timeinfo) );
	char record_start_time[50];
	sprintf(record_start_time, "%s", asctime (timeinfo));
	
	//prepare spicker
	count_decoy_spicker = 0;
	file_num = 1;
	spicker_file( pose );
	//initialization population
	const Size NP = 100;//100
	const Size GMAX = 300;//contact score
	const	Size GMAX2 = 100;//contact ditribution 
	//test	
	const float cluster_propotion = 0.8; //( 1 - cluster_propotion ) * GMAX for clustering
	//parameters in population evolution
	const Real CR = 0.5;
	const Real per_contact = 0.7;
	//prepare contact_bin
	int start = 4;
	int end = 20;
	int interval = 1;
	initiaize_contact_bin( start, end, interval);
	//record these parameters
	ofstream para("parameters.txt");
	para << "NP:" << NP << "\n" 
			 << "GMAX:" << GMAX << "\n" 
			 << "GMAX2:" << GMAX2 << "\n" 
			 << "cluster_propotion:" << cluster_propotion << "\n" 
			 << "CR:" << CR << "\n" 
			 << "per_contact:" << per_contact << "\n"
			 <<  "start:" << start << "\n" << "end: " << end << "\ninterval: " << interval << "\n"
			 << "multiplier: 10" << "\n" ;
	para.close();
	std::vector< core::pose::Pose > population;
	for(Size i = 0; i < NP; i ++)
	{
		population.push_back(pose);
	}
	
	for(Size i = 0; i < NP; i ++)
	{
		pose = population[i];
		if ( !bSkipStage1_ ) {
		flag_stage2to4 = false;
		PROF_START( basic::STAGE1 );
		clock_t starttime = clock();

		if ( !prepare_stage1( pose ) ) {
			set_last_move_status( moves::FAIL_RETRY );
			return;
		}
		// part 1 ----------------------------------------
		tr.Info <<  "\n===================================================================\n";
		tr.Info <<  "   Stage 1                                                         \n";
		tr.Info <<  "   Folding with score0 for max of " << stage1_cycles() << std::endl;

		if ( option[ basic::options::OptionKeys::run::profile ] ) prof_show();
		if ( option[ basic::options::OptionKeys::abinitio::debug ]() ) {
			output_debug_structure( pose, "stage0" );
		}
		if ( !get_checkpoints().recover_checkpoint( pose, get_current_tag(), "stage_1", false /* fullatom*/, true /*fold tree */ ) ) {
			ConstraintSetOP orig_constraints(NULL);
			orig_constraints = pose.constraint_set()->clone();
			success = do_stage1_cycles( pose );

			if ( tr.Info.visible() ) current_scorefxn().show( tr, pose );
			recover_low( pose, STAGE_1 );
			mc().show_counters();
			total_trials_+=mc().total_trials();
			mc().reset_counters();

			pose.constraint_set( orig_constraints ); // restore constraints - this is critical for checkpointing to work
			get_checkpoints().checkpoint( pose, get_current_tag(), "stage_1", true /*fold tree */ );
		} //recover checkpoint
		get_checkpoints().debug( get_current_tag(), "stage_1", current_scorefxn()( pose ) );

		clock_t endtime = clock();
		PROF_STOP( basic::STAGE1 );
		if ( option[ basic::options::OptionKeys::run::profile ] ) prof_show();
		if ( option[ basic::options::OptionKeys::abinitio::debug ]() ) {
			tr.Info << "Timeperstep: " << (double(endtime) - starttime )/(CLOCKS_PER_SEC ) << std::endl;
			output_debug_structure( pose, "stage1" );
		}
		} //skipStage1
		population[i] = pose;
	}//end initialization
	cout << "Success: initialization population" << endl;
	//miss77
	//process information: rmsd vs energy
	ofstream accept_rmsd_energy_1("accept_rmsd_energy_1.csv"); 

	ofstream accept_rmsd_energy("accept_rmsd_energy.csv"); 
	ofstream average_rmsd_energy("average_rmsd_energy.csv"); 
	//contact based perturbation
	//population evolution
	vector<double> probability = RaptorXcontact.get_prob();
	roulette_wheel_select contac_roulette( probability );
	for (Size g = 0; g < GMAX2; g ++)
	{	
		for (Size i = 0; i < NP; i ++)
		{
			//prepare for perturbation
			//chose one contact pair according its contact probability
			int contac_indx = contac_roulette.selection();
			vector<int> pertb_pair = RaptorXcontact.get_contact_index( contac_indx );
			int expected_bin = prob_chose_bin( probability[contac_indx] );			
			//set the range of fragment assembly 
			int fragAss_size = 12;
			vector<vector<int> > frag_range = set_frag_range( pertb_pair , fragAss_size, population[i] );
			
			//out put frag_range
			cout << frag_range[0][0]<<"*"<< frag_range[0][1]<<"*"<<frag_range[1][0]<<"*"<<frag_range[1][1]<<"*\n";
			//set_move_map
			core::kinematics::MoveMapOP map_part( new core::kinematics::MoveMap( *movemap() ) );
			map_part->set_bb(false);
			//set_move_map adverse
		    core::kinematics::MoveMapOP map_part_adverse( new core::kinematics::MoveMap( *movemap() ) );
			map_part_adverse->set_bb(true);
			for(vector<vector<int> >::iterator it = frag_range.begin() ; it != frag_range.end() ; it ++)
				{
					for(int pos_b = (*it)[0] ; pos_b <= (*it)[1]; pos_b ++ )
						{
							// cout << pos_b << "&" ;
							map_part->set_bb(pos_b, true);
							map_part_adverse->set_bb(pos_b, false);
						}
				}
			//set fragmentMoverOP
	   		simple_moves::FragmentMoverOP fragment_assembly_part;
			fragment_assembly_part = brute_move_large_;
			fragment_assembly_part->set_movemap( map_part );
			//fragment assembly, length 9

			int count_num = 100;
			int current_count = 0;
			int continu_count = 0;
			int accepted_count = 0;
			//
			double current_distance = get_distance_paired( population[i] , pertb_pair[0], pertb_pair[1]);
     		 //int current_bin = floor( current_distance ) - 4;
			double delta = abs( current_distance - expected_bin - 4.5);//caclulate the difference between the two residue to the expected bin center 


			core::pose::Pose trial_current = population[i];
			core::pose::Pose trial_temp = trial_current;
			
			cout << "residue index: " << pertb_pair[0] << "\t" << pertb_pair[1] << "\tcontact prob: " << probability[contac_indx] <<"\texpected_bin: " << expected_bin << "\n";
			do
			{
				fragment_assembly_part->apply( trial_temp );
				current_distance = get_distance_paired( trial_temp , pertb_pair[0], pertb_pair[1]);
			  
				cout << "current_distance: " << current_distance << "expected center: " << expected_bin + 3.5 << "\n";
				if( (abs( current_distance - expected_bin -4.5)) < delta && (current_distance > 3.8))
					{
						trial_current = trial_temp;//accept this fragment assembly
						delta = abs( current_distance - expected_bin -4.5);//apdate dalta
						
						accept_rmsd_energy << format( core::scoring::CA_rmsd( native_pose_temp, trial_current ) )  << ","
										   << format( (*score_stage4_)( trial_current ) ) << "\n";
						current_count = 0;
						accepted_count ++;    			

			  	}
		  	else
				{
					trial_temp = trial_current;//reject this fragment assembly
					current_count ++;	
				}
				continu_count ++;
			}while(current_count < count_num && abs(current_distance - expected_bin -4.5) > 0.5 );

			
			if(abs(current_distance - expected_bin -4.5) <= 0.5)
				{
					cout << "Success in perturbation of expacted bin!\n Accepted Frequency: " << (float)accepted_count / continu_count << "\n";
				}
			else
				{
					cout << "Fail in perturbation of expacted bin!\n Accepted Frequency: " << (float)accepted_count / continu_count << "\n";
				}
				//mutation: fragment assembly, length 3 in contact region 
				//energy function evolution	
				int trajectory_length_small = 100;
				for(int k = 0; k < trajectory_length_small; k ++)
				{
					//set fragmentMoverOP

			    	simple_moves::FragmentMoverOP fragment_assembly_part_small;
				  	fragment_assembly_part_small = brute_move_small_;
				 	fragment_assembly_part_small->set_movemap( map_part );
		  
					fragment_assembly_part_small->apply( trial_temp );
					bool accept_reject = Boltzmann_pose(trial_current , trial_temp , RaptorXcontact , score_stage4_ );
					if(accept_reject) 
					{
						accept_rmsd_energy << format( core::scoring::CA_rmsd( native_pose_temp, trial_temp ) )  << ","
	    											   << format( (*score_stage4_)( trial_temp ) ) << ","
	    											   << format( RaptorXcontact.score_C( trial_temp ) )<< ","
													   << "fragment_assembly_part_small" << "\n";
	    											   
		    		//record decoys for clustering									 
/* 		    		if( ((float)g / GMAX2) >	cluster_propotion	)
					{
						const int COUNT_MAX = 20000;
						if( count_decoy_spicker < COUNT_MAX )
						{
								count_decoy_spicker ++;
						}
						else
						{
								count_decoy_spicker = 1;
								file_num ++;
						}
						output_spicker( trial_temp );
					} */
	    		}
	    	//energy function evolution
				//mutation: fragment assembly, length 9
				int  trajectory_length_large = 30;
				for(int k = 0; k < trajectory_length_large; k ++)
				{
 					//set fragmentMoverOP
			      simple_moves::FragmentMoverOP fragment_assembly_part_large;
				  fragment_assembly_part_large = brute_move_large_;
				  fragment_assembly_part_large->set_movemap( map_part_adverse );
					
					fragment_assembly_part_large->apply( trial_temp );
					bool accept_reject = Boltzmann_pose(trial_current , trial_temp , RaptorXcontact , score_stage4_ );
					if(accept_reject) 
					{
						accept_rmsd_energy << format( core::scoring::CA_rmsd( native_pose_temp, trial_temp ) )  << ","
	    											   << format( (*score_stage4_)( trial_temp ) ) << ","
	    											   << format( RaptorXcontact.score_C( trial_temp ) )<< ","
													   << "fragment_assembly_part_large" << "\n";
	    											   
		    		//record decoys for clustering									 
/* 		    		if( ((float)g / GMAX2) >	cluster_propotion	)
		    			{
			    			const int COUNT_MAX = 20000;
								if( count_decoy_spicker < COUNT_MAX )
								{
										count_decoy_spicker ++;
								}
								else
								{
										count_decoy_spicker = 1;
										file_num ++;
								}
								output_spicker( trial_temp );
		    			}*/
	    			} 
				}
				}
				population[i] = trial_current;
		}//end once generation
		Real average_rmsd = 0, average_energy = 0;
		get_average_rmsd_energy( population ,  average_rmsd, average_energy);
		//cout << "average rmsd and energy: " << format( average_rmsd ) << "\t"<< format( average_energy) << endl;
		average_rmsd_energy << format( average_rmsd ) << ","
							<< format( average_energy) << "\n";
		if( g % 3 == 0 ) average_rmsd_energy.flush();//clear buffer
		output_status_info( timeinfo , g , GMAX + GMAX2);	
		contact_dis( population, RaptorXcontact );
		
	}//end population evolution
	cout << "End stage1. Start contact based perturbation.\n"; 
	//out put the information of the last generation
	ofstream rmsd_energy_after_stage1("rmsd_energy_after_stage1.csv");
	for (Size i = 0; i < NP; i ++)
	{
		rmsd_energy_after_stage1 << format( core::scoring::CA_rmsd( native_pose_temp, population[i] ) )  << ","
	    						 << format( (*score_stage4_)( population[i] ) ) << ","
	    						 << format( RaptorXcontact.score_C( population[i] ) )<< "n";
	}
	rmsd_energy_after_stage1.close();

	//population evolution:stage2
	for (Size g = 0; g < GMAX; g ++)
	{	
		for (Size i = 0; i < NP; i ++)
		{
			//crossover
			core::pose::Pose trial_cross;
				
			if( numeric::random::uniform() < CR)
			{
				crossover( population , i,  trial_cross);
			}
			else trial_cross = population[i];
			
			Size count_Max = 150, count = 0;
			bool accept_state = false;
			while ((!accept_state) && (count < count_Max)) {
				//mutation: fragment assembly, length 9
				core::pose::Pose trial_mutate = trial_cross;
				moves::MoverOP trial= trial_large();
				trial->apply( trial_mutate );
				//selection: score3 and score_c
				Real pro_c = 1 - per_contact * (float)g / GMAX;
//				pro_c = 1;//only contact score funtion 
				bool score_c_state = false;
				if( numeric::random::rg().uniform() < pro_c)
				{
						accept_state = Boltzmann_pose(population[i] , trial_mutate , RaptorXcontact);
						if( accept_state ) score_c_state = true;
				}
				else
				{
						accept_state = Boltzmann_pose(population[i] , trial_mutate , score_stage4_ );
				}
								
				if( accept_state ) //record process information
				{
					accept_rmsd_energy_1 << format( core::scoring::CA_rmsd( native_pose_temp, trial_mutate ) )  << ","
	    								 << format( (*score_stage4_)( trial_mutate ) ) << ",";
	    		if ( score_c_state ) accept_rmsd_energy_1 << format( RaptorXcontact.score_C( trial_mutate ) ) << ",";
	    		accept_rmsd_energy_1 << "\n";
	    		
	    		//record decoys for clustering									 
	    		if( ((float)g / GMAX) >	cluster_propotion	)
	    			{
		    			const int COUNT_MAX = 20000;
							if( count_decoy_spicker < COUNT_MAX )
							{
									count_decoy_spicker ++;
							}
							else
							{
									count_decoy_spicker = 1;
									file_num ++;
							}
							output_spicker( trial_mutate );
	    			}					 
				}				
				count ++;
			}//mutation and selection
			//calculate 
			//cout <<"iteration-g: " << g << "\t" << "individual: " << i << "\t" << "trail frequency: " << (format)( (float)count / count_Max ) << endl;
		}//end once generation
		//average rmsd vs energy of current population
		Real average_rmsd = 0, average_energy = 0;
		get_average_rmsd_energy( population ,  average_rmsd, average_energy);
		//cout << "average rmsd and energy: " << format( average_rmsd ) << "\t"<< format( average_energy) << endl;
		average_rmsd_energy << format( average_rmsd ) << ","
											 	<< format( average_energy) << "\n";
		if( g % 3 == 0 ) average_rmsd_energy.flush();//clear buffer
		output_status_info( timeinfo , g + GMAX2 , GMAX + GMAX2);	
		//record the max distance for each distance at g-th generation
		//record_max_distance( population , RaptorXcontact);
		//record the distribution of distance of contact pairs:4-5,5-6,6-7,7-8,8-9,9-10,10-11,11-12,12-13,13-14,14-15,15-16,16-17,17-18,18-19,19-20,20- 
		contact_dis( population, RaptorXcontact );
	}//end population evolution

	accept_rmsd_energy.close();
	average_rmsd_energy.close();
	accept_rmsd_energy_1.close();

	
	if (success) 
	{
	  //end time
//		time_t rawtime;
	  struct tm * timeinfo2;
	  time ( &rawtime );
	  timeinfo2 = localtime ( &rawtime );
	  output_status_info( record_start_time , timeinfo2 , GMAX+ GMAX2 , GMAX+ GMAX2);	

		spicker();
		return;
	}
	
	if ( !success ) {
		set_last_move_status( moves::FAIL_RETRY );
		return;
	}


	if ( !bSkipStage2_ ) {
		flag_stage2to4 = true;
		//
		// part 2 ----------------------------------------
		tr.Info <<  "\n===================================================================\n";
		tr.Info <<  "   Stage 2                                                         \n";
		tr.Info <<  "   Folding with score1 for " << stage2_cycles() << std::endl;

		PROF_START( basic::STAGE2 );
		clock_t starttime = clock();


		if ( close_chbrk_ ) {
			Real const setting( 0.25 );
			set_score_weight( scoring::linear_chainbreak, setting, STAGE_2 );
			tr.Info <<  " Chain_break score assigned " << std::endl;
		}


		if ( !prepare_stage2( pose ) )  {
			set_last_move_status( moves::FAIL_RETRY );
			return;
		}

		if ( !get_checkpoints().recover_checkpoint( pose, get_current_tag(), "stage_2", false /* fullatom */, true /*fold tree */ ) ) {
			ConstraintSetOP orig_constraints(NULL);
			orig_constraints = pose.constraint_set()->clone();

			success = do_stage2_cycles( pose );
			recover_low( pose, STAGE_2 );                    //default OFF: seems to be a bad choice after score0

			if  ( tr.visible() ) current_scorefxn().show( tr, pose );
			mc().show_counters();
			total_trials_+=mc().total_trials();
			mc().reset_counters();

			pose.constraint_set( orig_constraints ); // restore constraints - this is critical for checkpointing to work
			get_checkpoints().checkpoint( pose, get_current_tag(), "stage_2", true /*fold tree */ );
		}
		get_checkpoints().debug( get_current_tag(), "stage_2", current_scorefxn()( pose ) );

		clock_t endtime = clock();
		PROF_STOP( basic::STAGE2 );
		if ( option[ basic::options::OptionKeys::run::profile ] ) prof_show();
		if ( option[ basic::options::OptionKeys::abinitio::debug ]() ) {
			output_debug_structure( pose, "stage2" );
			tr << "Timeperstep: " << (double(endtime) - starttime )/(CLOCKS_PER_SEC ) << std::endl;
		}
	} //bSkipStage2

	if ( !success ) {
		set_last_move_status( moves::FAIL_RETRY );
		return;
	}

	if ( !bSkipStage3_ ) {
		flag_stage2to4 = true;
		// moved checkpointing into do_stage3_cycles because of structure store

		// part 3 ----------------------------------------
		tr.Info <<  "\n===================================================================\n";
		tr.Info <<  "   Stage 3                                                         \n";
		tr.Info <<  "   Folding with score2 and score5 for " << stage3_cycles() <<std::endl;

		PROF_START( basic::STAGE3 );
		clock_t starttime = clock();

		if ( !prepare_stage3( pose ) ) {
			set_last_move_status( moves::FAIL_RETRY );
			return;
		}
		// this is not the final score-function.. only known after prepare_loop_in_stage3
		// because this is confusing rather not show.if ( tr.Info.visible() ) current_scorefxn().show( tr, pose );

		success = do_stage3_cycles( pose );
		recover_low( pose, STAGE_3b );

		if ( tr.Info.visible() ) current_scorefxn().show( tr, pose );
		mc().show_counters();
		total_trials_+=mc().total_trials();
		mc().reset_counters();

		clock_t endtime = clock();
		PROF_STOP( basic::STAGE3);
		if ( option[ basic::options::OptionKeys::run::profile ] ) prof_show();
		if ( option[ basic::options::OptionKeys::abinitio::debug ]() ) {
			output_debug_structure( pose, "stage3" );
			tr << "Timeperstep: " << (double(endtime) - starttime )/( CLOCKS_PER_SEC) << std::endl;
		}

		//  pose.dump_pdb("stage3.pdb");

	}

	if ( !success ) {
		set_last_move_status( moves::FAIL_RETRY );
		return;
	}
	if ( !bSkipStage4_ ) {
		flag_stage2to4 = true;
		// part 4 ------------------------------------------
		tr.Info <<  "\n===================================================================\n";
		tr.Info <<  "   Stage 4                                                         \n";
		tr.Info <<  "   Folding with score3 for " << stage4_cycles() <<std::endl;

		PROF_START( basic::STAGE4 );
		clock_t starttime = clock();

		if ( !prepare_stage4( pose ) ) {
			set_last_move_status( moves::FAIL_RETRY );
			return;
		}

		//score-fxn may be changed in do_stage4_cycles...
		// confusing if shown here already... if ( tr.Info.visible() ) current_scorefxn().show( tr, pose);
		success = do_stage4_cycles( pose );
		recover_low( pose, STAGE_4  );

		if ( tr.Info.visible() ) current_scorefxn().show( tr, pose);
		mc().show_counters();
		total_trials_+=mc().total_trials();
		mc().reset_counters();

		clock_t endtime = clock();
		PROF_STOP( basic::STAGE4 );
		if ( option[ basic::options::OptionKeys::run::profile ] ) prof_show();
		if ( option[ basic::options::OptionKeys::abinitio::debug ]() ) {
			output_debug_structure( pose, "stage4" );
			tr << "Timeperstep: " << (double(endtime) - starttime )/( CLOCKS_PER_SEC ) << std::endl;
		}

		tr.Info <<  "\n===================================================================\n";
		tr.Info <<  "   Finished Abinitio                                                 \n";
		tr.Info <<  std::endl;
		//  pose.dump_pdb("stage4.pdb");
	}
	
	flag_stage5 = false;
	if ( !bSkipStage5_ ) {
		flag_stage2to4 = false;
		flag_stage5 = true;

		// part 5 ------------------------------------------
		tr.Info <<  "\n===================================================================\n";
		tr.Info <<  "   Stage 5                                                         \n";
		tr.Info <<  "   Folding with score3 for " << stage5_cycles() <<std::endl;

		PROF_START( basic::STAGE5 );
		clock_t starttime = clock();

		if ( !prepare_stage5( pose ) ) {
			set_last_move_status( moves::FAIL_RETRY );
			return;
		}

		success = do_stage5_cycles( pose );
		recover_low( pose, STAGE_5 );

		if ( tr.Info.visible() ) current_scorefxn().show( tr, pose);
		//  current_scorefxn().show(tr, pose);
		mc().show_counters();
		total_trials_+=mc().total_trials();
		mc().reset_counters();

		clock_t endtime = clock();
		PROF_STOP( basic::STAGE5 );
		if ( option[ basic::options::OptionKeys::run::profile ] ) prof_show();
		if ( option[ basic::options::OptionKeys::abinitio::debug ]() ) {
			output_debug_structure( pose, "stage5" );
			tr << "Timeperstep: " << (double(endtime) - starttime )/( CLOCKS_PER_SEC ) << std::endl;
		}

		tr.Info <<  "\n===================================================================\n";
		tr.Info <<  "   Now really finished Abinitio                                                 \n";
		tr.Info <<  std::endl;
		//  pose.dump_pdb("stage5.pdb");

	}
	
	//@miss77
	char fname[256];

	sprintf(fname , "spicker_decoy");
	string savepath(fname);
	ofstream fp(savepath.c_str() , std::ios::app);

	fp << setw(8)  << pose.total_residue()
		 << format( (*score_stage4_)(pose) )
		 << setw(8) << "&%*" << "\n";

	for(Size i = 1; i <= pose.total_residue() ; i ++ )
	{
			numeric::xyzVector<core::Real> v = pose.residue( i ).xyz(" CA ");
			fp << format( v[0] )
				 << format( v[1] )
				 << format( v[2] )<< "\n";
	}
	//comopare rmsd&&energy
	char f_energy[256];
	sprintf( f_energy , "rmsd_energy%.3f.csv", Real( option[ OptionKeys::abinitio::increase_cycles  ]) );
	string f_energy_path( f_energy );
	ofstream fout( f_energy_path.c_str() , std::ios::app );
	fout <<  format( core::scoring::CA_rmsd( native_pose_temp, pose ) )  << ","
	     <<  format( (*score_stage4_)(pose) ) << "\n";
	
	//compare second structure
	char f_name_sec[256];				
	sprintf( f_name_sec , "second_structure_%.3f.csv" , Real( option[ OptionKeys::abinitio::increase_cycles ]) );
	string f_name_sec_path( f_name_sec );
	ofstream fout_sec( f_name_sec_path.c_str() , std::ios::app );
	for(Size i = 1; i <= pose.total_residue() ; i ++ )
	{				
		fout_sec << pose.secstruct(i) << native_ss[i - 1] << "," ;
	}
	fout_sec << "\n";

	get_checkpoints().flush_checkpoints();

	if ( !success ) set_last_move_status( moves::FAIL_RETRY );

	//basic::prof_show();

	return;
}// ClassicAbinitio::apply( pose::Pose & pose )

char* ClassicAbinitio::format(double value)
{
		char *temp = new char[50];
		sprintf(temp, "%10.3f", value);
		return temp;

}
std::string
ClassicAbinitio::get_name() const {
	return "ClassicAbinitio";
}

//@brief return FramgentMover for smooth_small fragment insertions (i.e., stage4 moves)
simple_moves::FragmentMoverOP
ClassicAbinitio::smooth_move_small() {
	return smooth_move_small_;
}

//@brief return FragmentMover for small fragment insertions ( i.e., stage3/4 moves )
simple_moves::FragmentMoverOP
ClassicAbinitio::brute_move_small() {
	return brute_move_small_;
}

//@brief return FragmentMover for large fragment insertions (i.e., stage1/2 moves )
simple_moves::FragmentMoverOP
ClassicAbinitio::brute_move_large() {
	return brute_move_large_;
}

//@brief change the movemap ( is propagated to mover-objects )
//@detail overload if your extension stores additional moves as member variables
void
ClassicAbinitio::set_movemap( core::kinematics::MoveMapCOP mm )
{
	movemap_ = mm;
	if ( smooth_move_small_ ) smooth_move_small_->set_movemap( mm );
	if ( brute_move_small_  ) brute_move_small_ ->set_movemap( mm );
	if ( brute_move_large_  ) brute_move_large_ ->set_movemap( mm );
}

//@brief set new instances of FragmentMovers
void
ClassicAbinitio::set_moves(
	simple_moves::FragmentMoverOP brute_move_small,
	simple_moves::FragmentMoverOP brute_move_large,
	simple_moves::FragmentMoverOP smooth_move_small
)
{
	smooth_move_small_ = smooth_move_small;
	brute_move_small_  = brute_move_small;
	brute_move_large_  = brute_move_large;
	update_moves();
}

//@brief returns current movemap
core::kinematics::MoveMapCOP
ClassicAbinitio::movemap() {
	return movemap_;
}

//@detail read cmd_line options and set default versions for many protocol members: trials/moves, score-functions, Monte-Carlo
void ClassicAbinitio::set_defaults( pose::Pose const& pose ) {
	temperature_ = 2.0;
	bSkipStage1_ = false;
	bSkipStage2_ = false;
	bSkipStage3_ = false;
	bSkipStage4_ = false;
	bSkipStage5_ = true; //vats is turned off by default
	set_default_scores();
	set_default_options();
	set_default_mc( pose, *score_stage1_ );
	update_moves();
}

//@detail called to notify about changes in Movers: new movemap or Moverclass
void ClassicAbinitio::update_moves() {
	/* set apply_large_frags_ and
	short_insert_region_
	*/
	/* what about move-map ? It can be set manually for all Fragment_Moves .. */
	// set_move_map();
	set_trials();
}

//@detail create instances of TrialMover for our FragmentMover objects
void ClassicAbinitio::set_trials() {
	// setup loop1
	runtime_assert( brute_move_large_ != 0 );
	trial_large_ = moves::TrialMoverOP( new moves::TrialMover( brute_move_large_, mc_ ) );
	//trial_large_->set_keep_stats( true );
	trial_large_->keep_stats_type( moves::accept_reject );

	runtime_assert( brute_move_small_ != 0 );
	trial_small_ = moves::TrialMoverOP( new moves::TrialMover( brute_move_small_, mc_ ) );
	//trial_small_->set_keep_stats( true );
	trial_small_->keep_stats_type( moves::accept_reject );

	runtime_assert( smooth_move_small_ != 0 );
	smooth_trial_small_ = moves::TrialMoverOP( new moves::TrialMover( smooth_move_small_, mc_ ) );
	//smooth_trial_small_->set_keep_stats( true );
	smooth_trial_small_->keep_stats_type( moves::accept_reject );

	//build trial_pack mover
	moves::SequenceMoverOP combo_small( new moves::SequenceMover() );
	combo_small->add_mover(brute_move_small_);
	combo_small->add_mover(pack_rotamers_);
	trial_small_pack_ = moves::TrialMoverOP( new moves::TrialMover(combo_small, mc_) );
	moves::SequenceMoverOP combo_smooth( new moves::SequenceMover() );
	combo_smooth->add_mover(smooth_move_small_);
	combo_smooth->add_mover(pack_rotamers_);
	smooth_trial_small_pack_ = moves::TrialMoverOP( new moves::TrialMover(combo_smooth, mc_) );
}

//@detail sets Monto-Carlo object to default
void ClassicAbinitio::set_default_mc(
	pose::Pose const & pose,
	scoring::ScoreFunction const & scorefxn
) {
	set_mc( moves::MonteCarloOP( new moves::MonteCarlo( pose, scorefxn, temperature_ ) ) );
}

//@detail sets Monto-Carlo object
void ClassicAbinitio::set_mc( moves::MonteCarloOP mc_in ) {
	mc_ = mc_in;
	if ( trial_large_ ) trial_large_->set_mc( mc_ );
	if ( trial_small_ ) trial_small_->set_mc( mc_ );
	if ( smooth_trial_small_ ) smooth_trial_small_->set_mc( mc_ );
}

//@detail override cmd-line setting for "increase_cycling"
void ClassicAbinitio::set_cycles( Real increase_cycles ) {
	stage1_cycles_ = static_cast< int > (2000 * increase_cycles);
	stage2_cycles_ = static_cast< int > (2000 * increase_cycles);
	stage3_cycles_ = static_cast< int > (2000 * increase_cycles);
	stage4_cycles_ = static_cast< int > (4000 * increase_cycles);
	stage5_cycles_ = static_cast< int > (50000* increase_cycles);//vats

	using namespace basic::options;
	if ( option[ OptionKeys::abinitio::only_stage1 ]() ) {
		stage2_cycles_ = 0;
		stage3_cycles_ = 0;
		stage4_cycles_ = 0;
		bSkipStage2_ = bSkipStage3_ = /*bSkipStage3_ =*/ true;  // Was bSkipStage4_ meant? ~Labonte
	}
}

void ClassicAbinitio::set_default_scores() {
	using namespace scoring;
	using namespace basic::options;
	using namespace basic::options::OptionKeys;

	tr.Debug << "creating standard scoring functions" << std::endl;

	if ( option[ OptionKeys::abinitio::stage1_patch ].user() ) {
		score_stage1_  = ScoreFunctionFactory::create_score_function( "score0", option[ OptionKeys::abinitio::stage1_patch ]() );
	} else {
		score_stage1_  = ScoreFunctionFactory::create_score_function( "score0" );
	}

	if ( option[ OptionKeys::abinitio::stage2_patch ].user() ) {
		score_stage2_  = ScoreFunctionFactory::create_score_function( "score1", option[ OptionKeys::abinitio::stage2_patch ]() );
	} else {
		score_stage2_  = ScoreFunctionFactory::create_score_function( "score1" );
	}

	if ( option[ OptionKeys::abinitio::stage3a_patch ].user() ) {
		score_stage3a_ = ScoreFunctionFactory::create_score_function( "score2", option[ OptionKeys::abinitio::stage3a_patch ]() );
	} else {
		score_stage3a_ = ScoreFunctionFactory::create_score_function( "score2" );
	}

	if ( option[ OptionKeys::abinitio::stage3b_patch ].user() ) {
		score_stage3b_ = ScoreFunctionFactory::create_score_function( "score5", option[ OptionKeys::abinitio::stage3b_patch ]() );
	} else {
		score_stage3b_ = ScoreFunctionFactory::create_score_function( "score5" );
	}

	if ( option[ OptionKeys::abinitio::stage4_patch ].user() ) {
		score_stage4_  = ScoreFunctionFactory::create_score_function( "score3", option[ OptionKeys::abinitio::stage4_patch ]() );
	} else {
		score_stage4_  = ScoreFunctionFactory::create_score_function( "score3" );
	}

	//loading the cenrot score
	score_stage4rot_ = ScoreFunctionFactory::create_score_function( "score4_cenrot_relax" );
	//score_stage4rot_->set_weight(core::scoring::cen_rot_dun, 0.0);
	score_stage4rot_sc_ = ScoreFunctionFactory::create_score_function( "score4_cenrot_repack" );
	//score_stage4rot_sc_->set_weight(core::scoring::cen_rot_dun, 1.0);

	if ( option[ OptionKeys::abinitio::stage5_patch ].user() ) { //vats
		score_stage5_  = ScoreFunctionFactory::create_score_function( "score3", option[ OptionKeys::abinitio::stage5_patch ]() );
	} else {
		score_stage5_  = ScoreFunctionFactory::create_score_function( "score3" );
	}


	if ( option[ OptionKeys::abinitio::override_vdw_all_stages ] ) {
		set_score_weight( scoring::vdw, option[ OptionKeys::abinitio::vdw_weight_stage1 ], ALL_STAGES );
	}
}


/// @brief sets a score weight for all stages of abinitio
void ClassicAbinitio::set_score_weight( scoring::ScoreType type, Real setting, StageID stage ) {
	tr.Debug << "set score weights for ";
	if ( stage == ALL_STAGES ) tr.Debug << "all stages ";
	else tr.Debug << "stage " << (stage <= STAGE_3a ? stage : ( stage-1 ) ) << ( stage == STAGE_3b ? "b " : " " );
	tr.Debug << scoring::name_from_score_type(type) << " " << setting << std::endl;
	if ( score_stage1_  && ( stage == STAGE_1  || stage == ALL_STAGES ) ) score_stage1_ ->set_weight(type, setting);
	if ( score_stage2_  && ( stage == STAGE_2  || stage == ALL_STAGES ) ) score_stage2_ ->set_weight(type, setting);
	if ( score_stage3a_ && ( stage == STAGE_3a || stage == ALL_STAGES ) ) score_stage3a_->set_weight(type, setting);
	if ( score_stage3b_ && ( stage == STAGE_3b || stage == ALL_STAGES ) ) score_stage3b_->set_weight(type, setting);
	if ( score_stage4_  && ( stage == STAGE_4  || stage == ALL_STAGES ) ) score_stage4_ ->set_weight(type, setting);
	if ( score_stage4rot_  && ( stage == STAGE_4  || stage == ALL_STAGES ) ) score_stage4rot_ ->set_weight(type, setting);
	if ( score_stage5_  && ( stage == STAGE_5  || stage == ALL_STAGES ) ) score_stage5_ ->set_weight(type, setting);//vats
}

//@brief currently used score function ( depends on stage )
scoring::ScoreFunction const& ClassicAbinitio::current_scorefxn() const {
	return mc().score_function();
}

//@brief set current scorefunction
void ClassicAbinitio::current_scorefxn( scoring::ScoreFunction const& scorefxn ) {
	mc().score_function( scorefxn );
}

//@brief set individual weight of current scorefunction --- does not change the predefined scores: score_stageX_
void ClassicAbinitio::set_current_weight( core::scoring::ScoreType type, core::Real setting ) {
	scoring::ScoreFunctionOP scorefxn ( mc().score_function().clone() );
	scorefxn->set_weight( type, setting );
	mc().score_function( *scorefxn ); //trigger rescore
}

void ClassicAbinitio::set_default_options() {
	bSkipStage1_ = bSkipStage2_ = bSkipStage3_ = bSkipStage4_ = false;
	bSkipStage5_ = true; //vats turned off by default
	using namespace basic::options;
	just_smooth_cycles_ = option[ OptionKeys::abinitio::smooth_cycles_only ]; // defaults to false
	bQuickTest_ = basic::options::option[ basic::options::OptionKeys::run::test_cycles ]();

	if ( bQuickTest() ) {
		set_cycles( 0.001 );
	} else {
		set_cycles( option[ OptionKeys::abinitio::increase_cycles ] ); // defaults to factor of 1.0
	}

	if ( just_smooth_cycles_ ) {
		bSkipStage1_ = bSkipStage2_ = bSkipStage3_ = bSkipStage5_ = true;
	}
	if ( option[ OptionKeys::abinitio::only_stage1 ] ) {
		bSkipStage2_ = bSkipStage3_ = bSkipStage4_ = bSkipStage5_= true;
	}

	if ( option[ OptionKeys::abinitio::include_stage5 ] ) {
		bSkipStage5_ = false;
	}

	apply_large_frags_   = true;  // apply large frags in phase 2!

	// in rosetta++ switched on in fold_abinitio if contig_size < 30 in pose_abinitio never
	short_insert_region_ = false;  // apply small fragments in phase 2!

	if ( option[ OptionKeys::abinitio::recover_low_in_stages ].user() ) {
		for ( IntegerVectorOption::const_iterator it = option[ OptionKeys::abinitio::recover_low_in_stages ]().begin(),
				eit = option[ OptionKeys::abinitio::recover_low_in_stages ]().end(); it!=eit; ++it ) {
			if ( *it == 1 ) recover_low_stages_.push_back( STAGE_1 );
			else if ( *it == 2 ) recover_low_stages_.push_back( STAGE_2 );
			else if ( *it == 3 ) {
				recover_low_stages_.push_back( STAGE_3a );
				recover_low_stages_.push_back( STAGE_3b );
			} else if ( *it == 4 ) recover_low_stages_.push_back( STAGE_4 );
		}
	} else {
		recover_low_stages_.clear();
		recover_low_stages_.push_back( STAGE_1 );
		recover_low_stages_.push_back( STAGE_2 );
		recover_low_stages_.push_back( STAGE_3a );
		recover_low_stages_.push_back( STAGE_3b );
		recover_low_stages_.push_back( STAGE_4 );
		recover_low_stages_.push_back( STAGE_5 );
	}

	close_chbrk_ = option[ OptionKeys::abinitio::close_chbrk ];

}


/// @brief (helper) functor class which keeps track of old pose for the
/// convergence check in stage3 cycles
/// @detail
/// calls of operator ( pose ) compare the
class hConvergenceCheck;
typedef  utility::pointer::shared_ptr< hConvergenceCheck >  hConvergenceCheckOP;

class hConvergenceCheck : public moves::PoseCondition {
public:
	hConvergenceCheck() : bInit_( false ), ct_( 0 ) {}
	void reset() { ct_ = 0; bInit_ = false; }
	void set_trials( moves::TrialMoverOP trin ) {
		trials_ = trin;
		runtime_assert( trials_->keep_stats_type() < moves::no_stats );
		last_move_ = 0;
	}
	virtual bool operator() ( const core::pose::Pose & pose );
private:
	pose::Pose very_old_pose_;
	bool bInit_;
	Size ct_;
	moves::TrialMoverOP trials_;
	Size last_move_;
};

// keep going --> return true
bool hConvergenceCheck::operator() ( const core::pose::Pose & pose ) {
	if ( !bInit_ ) {
		bInit_ = true;
		very_old_pose_ = pose;
		return true;
	}
	runtime_assert( trials_ != 0 );
	tr.Trace << "TrialCounter in hConvergenceCheck: " << trials_->num_accepts() << std::endl;
	if ( numeric::mod(trials_->num_accepts(),100) != 0 ) return true;
	if ( (Size) trials_->num_accepts() <= last_move_ ) return true;
	last_move_ = trials_->num_accepts();
	// change this later to this: (after we compared with rosetta++ and are happy)
	// if ( numeric::mod(++ct_, 1000) != 0 ) return false; //assumes an approx acceptance rate of 0.1

	// still here? do the check:

	core::Real converge_rms = core::scoring::CA_rmsd( very_old_pose_, pose );
	very_old_pose_ = pose;
	if ( converge_rms >= 3.0 ) {
		return true;
	}
	// if we get here thing is converged stop the While-Loop
	tr.Info << " stop cycles in stage3 due to convergence " << std::endl;
	return false;
}


bool ClassicAbinitio::do_stage1_cycles( pose::Pose &pose ) {
	AllResiduesChanged done( pose, brute_move_large()->insert_map(), *movemap() );
	moves::MoverOP trial( stage1_mover( pose, trial_large() ) );

	// FragmentMoverOP frag_mover = brute_move_large_;
	// fragment::FragmentIO().write("stage1_frags_classic.dat",*frag_mover->fragments());

	Size j;
	for ( j = 1; j <= stage1_cycles(); ++j ) {
		trial->apply( pose ); // apply a large fragment insertion, accept with MC boltzmann probability
		if ( done(pose) ) {
			tr.Info << "Replaced extended chain after " << j << " cycles." << std::endl;
			mc().reset( pose ); // make sure that we keep the final structure
			return true;
		}
	}
	tr.Warning << "Warning: extended chain may still remain after " << stage1_cycles() << " cycles!" << std::endl;
	done.show_unmoved( pose, tr.Warning );
	mc().reset( pose ); // make sure that we keep the final structure
	return true;
}

bool ClassicAbinitio::do_stage2_cycles( pose::Pose &pose ) {

	//setup cycle
	moves::SequenceMoverOP cycle( new moves::SequenceMover() );
	if ( apply_large_frags_   ) cycle->add_mover( trial_large_->mover() );
	if ( short_insert_region_ ) cycle->add_mover( trial_small_->mover() );

	Size nr_cycles = stage2_cycles() / ( short_insert_region_ ? 2 : 1 );
	moves::TrialMoverOP trials( new moves::TrialMover( cycle, mc_ptr() ) );
	moves::RepeatMover( stage2_mover( pose, trials ), nr_cycles ).apply(pose);

	//is there a better way to find out how many steps ? for instance how many calls to scoring?
	return true; // as best guess
}

/*! @detail stage3 cycles:
nloop1 : outer iterations
nloop2 : inner iterations
stage3_cycle : trials per inner iteration
every inner iteration we switch between score_stage3a ( default: score2 ) and score_stage3b ( default: score 5 )

prepare_loop_in_stage3() is called before the stage3_cycles() of trials are started.

first outer loop-iteration is done with TrialMover trial_large()
all following iterations with trial_small()

start each iteration with the lowest_score_pose. ( mc->recover_low() -- called in prepare_loop_in_stage3() )

*/
bool ClassicAbinitio::do_stage3_cycles( pose::Pose &pose ) {
	using namespace ObjexxFCL;

	// interlaced score2 / score 5 loops
	// nloops1 and nloops2 could become member-variables and thus changeable from the outside
	int nloop1 = 1;
	int nloop2 = 10; //careful: if you change these the number of structures in the structure store changes.. problem with checkpointing
	// individual checkpoints for each stage3 iteration would be a remedy. ...

	if ( short_insert_region_ ) {
		nloop1 = 2;
		nloop2 = 5;
	}

	hConvergenceCheckOP convergence_checker ( NULL );
	if ( !option[ basic::options::OptionKeys::abinitio::skip_convergence_check ] ) {
		convergence_checker = hConvergenceCheckOP( new hConvergenceCheck );
	}

	moves::TrialMoverOP trials = trial_large();
	int iteration = 1;
	for ( int lct1 = 1; lct1 <= nloop1; lct1++ ) {
		if ( lct1 > 1 ) trials = trial_small(); //only with short_insert_region!
		for ( int lct2 = 1; lct2 <= nloop2; lct2++, iteration++  ) {
			tr.Debug << "Loop: " << lct1 << "   " << lct2 << std::endl;

			if ( !prepare_loop_in_stage3( pose, iteration, nloop1*nloop2 ) ) return false;

			if ( !get_checkpoints().recover_checkpoint( pose, get_current_tag(), "stage_3_iter"+string_of( lct1)+"_"+string_of(lct2),
					false /*fullatom */, true /*fold tree */ ) ) {


				tr.Debug << "  Score stage3 loop iteration " << lct1 << " " << lct2 << std::endl;
				if ( convergence_checker ) {
					moves::TrialMoverOP stage3_trials = stage3_mover( pose, lct1, lct2, trials );
					convergence_checker->set_trials( stage3_trials ); //can be removed late
					moves::WhileMover( stage3_trials, stage3_cycles(), convergence_checker ).apply( pose );
				} else {    //no convergence check -> no WhileMover
					moves::RepeatMover( stage3_mover( pose, lct1, lct2, trials ), stage3_cycles() ).apply( pose );
				}

				if ( numeric::mod( (int)iteration, 2 ) == 0 || iteration > 7 ) recover_low( pose, STAGE_3a );
				recover_low( pose, STAGE_3b );

				get_checkpoints().checkpoint( pose, get_current_tag(), "stage_3_iter"+string_of( lct1)+"_"+string_of(lct2), true /*fold tree */ );
			}//recover_checkpoint
			get_checkpoints().debug( get_current_tag(), "stage_3_iter"+string_of( lct1)+"_"+string_of(lct2), current_scorefxn()( pose ) );

			//   structure_store().push_back( mc_->lowest_score_pose() );
		} // loop 2
	} // loop 1
	return true;
}


// interlaced score2 / score 5 loops
/*! @detail stage4 cycles:
nloop_stage4: iterations
stage4_cycle : trials per  iteration

first iteration: use trial_small()
following iterations: use trial_smooth()
only trial_smooth() if just_smooth_cycles==true

prepare_loop_in_stage4() is called each time before the stage4_cycles_ of trials are started.

start each iteration with the lowest_score_pose. ( mc->recover_low()  in prepare_loop_in_stage4()  )

*/
bool ClassicAbinitio::do_stage4_cycles( pose::Pose &pose ) {
	Size nloop_stage4 = 3;

	using namespace basic::options;
	using namespace basic::options::OptionKeys;

	if ( option[corrections::score::cenrot]() ) nloop_stage4=2;

	for ( Size kk = 1; kk <= nloop_stage4; ++kk ) {
		tr.Debug << "prepare ..." << std::endl ;
		if ( !prepare_loop_in_stage4( pose, kk, nloop_stage4 ) ) return false;

		if ( !get_checkpoints().recover_checkpoint( pose, get_current_tag(), "stage4_kk_" + ObjexxFCL::string_of(kk), false /* fullatom */, true /* fold_tree */ ) ) {
			moves::TrialMoverOP trials;
			if ( kk == 1 && !just_smooth_cycles_ ) {
				trials = trial_small();
			} else {
				tr.Debug << "switch to smooth moves" << std::endl;
				trials = trial_smooth();
			}

			tr.Debug << "start " << stage4_cycles() << " cycles" << std::endl;
			moves::RepeatMover( stage4_mover( pose, kk, trials ), stage4_cycles() ).apply(pose);
			tr.Debug << "finished" << std::endl;
			recover_low( pose, STAGE_4 );

			get_checkpoints().checkpoint( pose, get_current_tag(), "stage4_kk_" + ObjexxFCL::string_of(kk), true /*fold tree */ );
		}
		get_checkpoints().debug( get_current_tag(), "stage4_kk_" + ObjexxFCL::string_of(kk),  current_scorefxn()( pose ) );

		//don't store last structure since it will be exactly the same as the final structure delivered back via apply
		//  if( kk < nloop_stage4 ) // <-- this line was missing although the comment above was existant.
		//   structure_store().push_back( mc_->lowest_score_pose() );
	}  // loop kk

	if ( option[corrections::score::cenrot] ) {
		//switch to cenrot model
		tr.Debug << "switching to cenrot model ..." << std::endl;
		protocols::simple_moves::SwitchResidueTypeSetMover to_cenrot(chemical::CENTROID_ROT);
		to_cenrot.apply(pose);

		//init pose
		(*score_stage4rot_)( pose );
		pack_rotamers_->score_function(score_stage4rot_sc_);
		pack_rotamers_->apply(pose);

		mc_->reset(pose);
		replace_scorefxn( pose, STAGE_4rot, 0 );
		//mc_->set_temperature(1.0);
		//mc_->set_autotemp(true, 1.0);

		//debug
		//tr.Debug << "starting_energy: " << (*score_stage4rot_)( pose ) << std::endl;
		//tr.Debug << "starting_temperature: " << mc_->temperature() << std::endl;

		for ( Size rloop=1; rloop<=3; rloop++ ) {
			//change vdw weight
			switch (rloop) {
			case 1 :
				score_stage4rot_->set_weight(core::scoring::vdw, score_stage4rot_->get_weight(core::scoring::vdw)/9.0);
				break;
			case 2 :
				score_stage4rot_->set_weight(core::scoring::vdw, score_stage4rot_->get_weight(core::scoring::vdw)*3.0);
				break;
			case 3 :
				score_stage4rot_->set_weight(core::scoring::vdw, score_stage4rot_->get_weight(core::scoring::vdw)*3.0);
				break;
			}

			//stage4rot
			//for (Size iii=1; iii<=100; iii++){
			//pose::Pose startP = pose;
			//tr << "temperature: " << mc_->temperature() << std::endl;
			moves::RepeatMover( stage4rot_mover( pose, rloop, trial_smooth() ), stage4_cycles()/100 ).apply(pose);
			//tr << "delta_rms: " << core::scoring::CA_rmsd( startP, pose ) << std::endl;
			//}
		}
	}

	return true;
}

bool ClassicAbinitio::do_stage5_cycles( pose::Pose &pose ) {//vats

	Size nmoves = 1;
	core::kinematics::MoveMapOP mm_temp( new core::kinematics::MoveMap( *movemap() ) );
	simple_moves::SmallMoverOP small_mover( new simple_moves::SmallMover( mm_temp, temperature_, nmoves) );
	small_mover->angle_max( 'H', 2.0 );
	small_mover->angle_max( 'E', 2.0 );
	small_mover->angle_max( 'L', 5.0 );

	moves::TrialMoverOP trials( new moves::TrialMover( small_mover, mc_ptr() ) );
	moves::RepeatMover( stage5_mover( pose, trials ), stage5_cycles() ).apply( pose );

	// moves::MoverOP trial( stage5_mover( pose, small_mover ) );
	// Size j;
	// for( j = 1; j <= stage5_cycles(); ++j ) {
	//  trial->apply( pose );
	// }
	mc().reset( pose );
	return true;

}


moves::TrialMoverOP
ClassicAbinitio::stage1_mover( pose::Pose &, moves::TrialMoverOP trials ) {
	return trials;
}


moves::TrialMoverOP
ClassicAbinitio::stage2_mover( pose::Pose &, moves::TrialMoverOP trials ) {
	return trials;
}


moves::TrialMoverOP
ClassicAbinitio::stage3_mover( pose::Pose &, int, int, moves::TrialMoverOP trials ) {
	return trials;
}

moves::TrialMoverOP
ClassicAbinitio::stage4_mover( pose::Pose &, int, moves::TrialMoverOP trials ) {
	return trials;
}

moves::TrialMoverOP
ClassicAbinitio::stage4rot_mover( pose::Pose &, int, moves::TrialMoverOP trials ) {
	if ( trials == trial_small_ ) {
		return trial_small_pack_;
	} else {
		return smooth_trial_small_pack_;
	}
}

moves::TrialMoverOP //vats
ClassicAbinitio::stage5_mover( pose::Pose &, moves::TrialMoverOP trials ) {
	return trials;
}

void ClassicAbinitio::recover_low( core::pose::Pose& pose, StageID stage ){
	if ( contains_stageid( recover_low_stages_, stage ) ) {
		mc_->recover_low( pose );
	}
}

// anything you want to have done before the stages ?
void ClassicAbinitio::replace_scorefxn( core::pose::Pose& pose, StageID stage, core::Real /*intra_stage_progress */ ) {
	// must assume that the current pose is the one to be accepted into the next stage! (this change was necessary for
	// checkpointing to work correctly.

	//intra_stage_progress = intra_stage_progress;
	if ( score_stage1_  && ( stage == STAGE_1 ) ) current_scorefxn( *score_stage1_ );
	if ( score_stage2_  && ( stage == STAGE_2 ) ) current_scorefxn( *score_stage2_ );
	if ( score_stage3a_ && ( stage == STAGE_3a) ) current_scorefxn( *score_stage3a_ );
	if ( score_stage3b_ && ( stage == STAGE_3b) ) current_scorefxn( *score_stage3b_ );
	if ( score_stage4_  && ( stage == STAGE_4 ) ) current_scorefxn( *score_stage4_ );
	if ( score_stage4rot_  && ( stage == STAGE_4rot ) ) current_scorefxn( *score_stage4rot_ );
	if ( score_stage5_  && ( stage == STAGE_5 ) ) current_scorefxn( *score_stage5_ );//vats
	Real temperature( temperature_ );
	if ( stage == STAGE_5 ) temperature = 0.5;
	mc_->set_autotemp( true, temperature );
	mc_->set_temperature( temperature ); // temperature might have changed due to autotemp..
	mc_->reset( pose );
}


moves::TrialMoverOP ClassicAbinitio::trial_large() {
	return ( apply_large_frags_ ? trial_large_ : trial_small_ );
}

moves::TrialMoverOP ClassicAbinitio::trial_small() {
	return trial_small_;
}

moves::TrialMoverOP ClassicAbinitio::trial_smooth() {
	return smooth_trial_small_;
}

// prepare stage1 sampling
bool ClassicAbinitio::prepare_stage1( core::pose::Pose &pose ) {
	replace_scorefxn( pose, STAGE_1, 0.5 );
	mc_->set_autotemp( false, temperature_ );
	// mc_->set_temperature( temperature_ ); already done in replace_scorefxn
	// mc_->reset( pose );
	(*score_stage1_)( pose );
	/// Now handled automatically.  score_stage1_->accumulate_residue_total_energies( pose ); // fix this
	return true;
}

bool ClassicAbinitio::prepare_stage2( core::pose::Pose &pose ) {
	replace_scorefxn( pose, STAGE_2, 0.5 );

	(*score_stage2_)(pose);
	/// Now handled automatically.  score_stage2_->accumulate_residue_total_energies( pose );
	return true;
}


bool ClassicAbinitio::prepare_stage3( core::pose::Pose &pose ) {
	replace_scorefxn( pose, STAGE_3a, 0 );
	//score for this stage is changed in the do_stage3_cycles explicitly
	if ( option[ templates::change_movemap ].user() && option[ templates::change_movemap ] == 3 ) {
		kinematics::MoveMapOP new_mm( new kinematics::MoveMap( *movemap() ) );
		new_mm->set_bb( true );
		set_movemap( new_mm ); // --> store it in movemap_ --> original will be reinstated at end of apply()
	}
	return true;
}


bool ClassicAbinitio::prepare_stage4( core::pose::Pose &pose ) {
	replace_scorefxn( pose, STAGE_4, 0 );
	(*score_stage4_)( pose );
	/// Now handled automatically.  score_stage4_->accumulate_residue_total_energies( pose ); // fix this

	if ( option[ templates::change_movemap ].user() && option[ templates::change_movemap ] == 4 ) {
		kinematics::MoveMapOP new_mm( new kinematics::MoveMap( *movemap() ) );
		new_mm->set_bb( true );
		tr.Debug << "option: templates::change_movemap ACTIVE: set_movemap" << std::endl;
		set_movemap( new_mm ); // --> store it in movemap_ --> original will be reinstated at end of apply()
	}
	return true;
}

bool ClassicAbinitio::prepare_stage5( core::pose::Pose &pose ) {//vats
	// temperature_ = 0.5; //this has to be reset to original temperature!!!
	// no special if-statement in replace_scorefxn...OL
	replace_scorefxn( pose, STAGE_5, 0 );
	(*score_stage5_)( pose );
	return true;
}


bool ClassicAbinitio::prepare_loop_in_stage3( core::pose::Pose &pose/*pose*/, Size iteration, Size total ){
	// interlace score2/score5

	Real chbrk_weight_stage_3a = 0;
	Real chbrk_weight_stage_3b = 0;

	if ( numeric::mod( (int)iteration, 2 ) == 0 || iteration > 7 ) {
		Real progress( iteration );
		chbrk_weight_stage_3a = 0.25 * progress;
		tr.Debug << "select score_stage3a..." << std::endl;
		recover_low( pose, STAGE_3a );
		replace_scorefxn( pose, STAGE_3a, 1.0* iteration/total );
	} else {
		Real progress( iteration );
		chbrk_weight_stage_3b = 0.05 * progress;
		tr.Debug << "select score_stage3b..." << std::endl;
		recover_low( pose, STAGE_3b );
		replace_scorefxn( pose, STAGE_3b, 1.0* iteration/total );
	}

	if ( close_chbrk_ ) {

		set_score_weight( scoring::linear_chainbreak, chbrk_weight_stage_3a , STAGE_3a );
		set_score_weight( scoring::linear_chainbreak, chbrk_weight_stage_3b , STAGE_3b );

	}


	return true;
}

bool ClassicAbinitio::prepare_loop_in_stage4( core::pose::Pose &pose, Size iteration, Size total ){
	replace_scorefxn( pose, STAGE_4, 1.0* iteration/total );

	Real chbrk_weight_stage_4 (iteration*0.5+2.5);

	if ( close_chbrk_ ) {
		set_current_weight( scoring::linear_chainbreak, chbrk_weight_stage_4 );
	}

	return true;
}

//@brief obtain currently used monte-carlo object --> use to obtain current score-func: mc().score_function()
moves::MonteCarloOP
ClassicAbinitio::mc_ptr() {
	return mc_;
}


void ClassicAbinitio::output_debug_structure( core::pose::Pose & pose, std::string prefix ) {
	using namespace core::io::silent;

	mc().score_function()( pose );
	Parent::output_debug_structure( pose, prefix );

	if ( option[ basic::options::OptionKeys::abinitio::explicit_pdb_debug ]() ) {
		pose.dump_pdb( prefix + get_current_tag() + ".pdb" );
	}

	if ( option[ basic::options::OptionKeys::abinitio::log_frags ].user() ) {
		std::string filename = prefix + "_" + get_current_tag() + "_" + std::string( option[ basic::options::OptionKeys::abinitio::log_frags ]() );
		utility::io::ozstream output( filename );
		simple_moves::LoggedFragmentMover& log_frag = dynamic_cast< simple_moves::LoggedFragmentMover& > (*brute_move_large_);
		log_frag.show( output );
		log_frag.clear();
	}

} // ClassicAbinitio::output_debug_structure( core::pose::Pose & pose, std::string prefix )

bool ClassicAbinitio::crossover( std::vector< core::pose::Pose > & population_, const Size index, core::pose::Pose &trial_){
		Size NP_ = population_.size();
		Size seq_len = population_[0].total_residue();
		//select a random conformati                 on
		Size rand_index = index;
		do
		{
			rand_index = numeric::random::random_range( 0, NP_ - 1);
		}while( rand_index == index );
		//random region
		Size cross_len = 3;
		Size cross_region_start = numeric::random::random_range( 1, seq_len - cross_len + 1);
		//exchange phi and psi 
		exchange_dihedral_angle( population_[index], population_[rand_index], trial_, cross_region_start, cross_len);
		return 1;
}
//exchange dihedral angle 
bool ClassicAbinitio::exchange_dihedral_angle( core::pose::Pose &target, core::pose::Pose &rand_conform, 			core::pose::Pose &trial_result, const Size region_start, const Size length)
{
	trial_result = target;
	for(Size i = 0; i < length; i ++)
	{
		Size position = i + region_start;
		trial_result.set_phi( position, rand_conform.phi( position ) );
	}
	return 1;
}
//calculate the score of two poses; if accepted, the pose is substituted.
bool ClassicAbinitio::Boltzmann_pose( core::pose::Pose &pose_before , core::pose::Pose &pose_late , core::scoring::ScoreFunctionOP score_function)
{
		double score_before = (*score_function)(pose_before);
		double score_late = (*score_function)(pose_late);
		bool accept_state = boltzmann( score_before , score_late);
		if(accept_state)
		{
				pose_before = pose_late;
		}
		return accept_state;
}
///calculate the score_contact of two poses; if accepted, the pose is substituted.
bool ClassicAbinitio::Boltzmann_pose( core::pose::Pose &pose_before , core::pose::Pose &pose_late , contact_predicted& contact_)
{
		double score_before = contact_.score_C( pose_before );
		double score_late = contact_.score_C( pose_late );
		//set the KT, change the scale of score_C
		int multiplier = 10;
		bool accept_state = boltzmann( multiplier * score_before , multiplier * score_late);
		if(accept_state)
		{
				pose_before = pose_late;
		}
		return accept_state;
}
bool ClassicAbinitio::Boltzmann_pose( core::pose::Pose &pose_before , core::pose::Pose &pose_late , contact_predicted& contact_, core::scoring::ScoreFunctionOP score_function)
	{
		double score_before1 = (*score_function)(pose_before);
		double score_late1 = (*score_function)(pose_late);
		bool accept_state1 = boltzmann( score_before1 , score_late1);
		
		double score_before2 = contact_.score_C( pose_before );
		double score_late2 = contact_.score_C( pose_late );
		//set the KT, change the scale of score_C
		int multiplier = 10;
		bool accept_state2 = boltzmann( multiplier * score_before2 , multiplier * score_late2);
		
		if( accept_state1 && accept_state2)
			{
				pose_before = pose_late;
				return true;
			}
			return false;
		
	}
//calculate the possibility of acceptance according to the Metropolis criterion
bool ClassicAbinitio::boltzmann( double score_before , double score_late )
{
		Real const score_delta( score_late - score_before );
		Real temperature_ = 2.0;
		Real const boltz_factor = -score_delta / temperature_ ;
		Real const probability = std::exp( std::min (40.0, std::max(-40.0,boltz_factor)) ) ;
		if ( probability < 1 ) {
				if ( numeric::random::rg().uniform() >= probability ) {
						return false; // rejected
				}
		}
		return true;
}

bool ClassicAbinitio::get_average_rmsd_energy( std::vector< core::pose::Pose > & population_, Real &rmsd, Real &energy)
{
	Size NP_ = population_.size();
	rmsd = 0;
	energy = 0;
	for(Size i = 0 ; i < NP_ ; i ++)
		{
			rmsd = rmsd + core::scoring::CA_rmsd( native_pose_temp, population_[i] ) ;
			energy = energy + (*score_stage4_)( population_[i] );
		}
		rmsd = rmsd / NP_;
		energy = energy / NP_;
		return 1;
}			
					
void ClassicAbinitio::spicker_file(core::pose::Pose &pose )
{
		Size residue_len = pose.total_residue();
		ofstream rmsinp("rmsinp");
		rmsinp <<"1" << setw(5) <<  residue_len << "\n"
					 << residue_len << "\n";
		rmsinp.close();

		ofstream seq;
		seq.open("seq.dat");
		for(Size i = 1; i <= residue_len ; i ++)
		{
				seq << setw(5) << i
						<< setw(6) << pose.residue(i).name3() << "\n";
		}
		seq.close();
}				
	
void ClassicAbinitio::output_spicker( core::pose::Pose &pose )
{
		char fname[256];

		sprintf(fname , "spicker_decoy%d",file_num);
		string savepath(fname);
		ofstream fp(savepath.c_str() , std::ios::app);

		fp << setw(8)  << pose.total_residue()
			 << format( (*score_stage4_)(pose) )
			 << setw( 8) << count_decoy_spicker << "\n";

		for(Size i = 1; i <= pose.total_residue() ; i ++ )
		{
				numeric::xyzVector<core::Real> v = pose.residue( i ).xyz(" CA ");
				fp << format( v[0] )
					 << format( v[1] )
					 << format( v[2] )<< "\n";
		}
}
bool ClassicAbinitio::spicker( )
		{
				ofstream tra_in;
				tra_in.open("tra.in");
				tra_in<< setw(2) << file_num << " 1 1" << "\n";
				for(int i = 1; i <= file_num ;i ++)
				{
						tra_in << "spicker_decoy" << i <<"\n";
				}
				tra_in.close();
				cout << "cmd:spicker\n";
				

				return system( "~/spicker/spicker");

		}
bool ClassicAbinitio::output_status_info( struct tm * timeinfo_ , const int g_ , const int G_)
{
		ofstream status_info("status_info.txt");
		status_info << asctime ( timeinfo_  ) << g_ + 1 << "\n" << G_ ;
		status_info.close();
		return 1;
}
bool ClassicAbinitio::output_status_info( char *start_time,  struct tm * timeinfo_2 , const int g_ , const int G_)
{
		ofstream status_info("status_info.txt");
		status_info << start_time << g_  << "\n" << G_  << "\n" << asctime ( timeinfo_2 );
		status_info.close();
		return 1;
}

void ClassicAbinitio::record_max_distance( const std::vector< core::pose::Pose > population_ , protocols::abinitio::contact_predicted &RaptorXcontact_)
{
		vector<vector<int> > contact_index_ = RaptorXcontact_.get_contact_index();
		vector<double> max_distance;
		for( unsigned int i = 0; i < contact_index_.size(); i ++)
		{
				vector<int> current_pair = contact_index_[i];
				double current_max_distance = get_max_distance_P( population_ , current_pair );
				max_distance.push_back( current_max_distance );				
				
		}
		//output max distance to file
		ofstream max_distance_file("max_distance.csv", ios::app);
		for(unsigned int i = 0; i < max_distance.size() ; i ++)
		{
				max_distance_file << max_distance[i] << ",";
				//cout << "max_distance" << max_distance[i] << "\t" ;
		}
		max_distance_file << '\n';
		max_distance_file.close();
		cout << "\n" ;
}		
		
		
		//calculate the max distance of a contact pair in population
double ClassicAbinitio::get_max_distance_P(std::vector< core::pose::Pose > population_ , const std::vector<int> current_pair_ )
{
	double max_dis = 0;
	for(unsigned int i = 0; i < population_.size() ; i ++)
	{
		// Get coordinates of both contact and calculate distance
		core::Real distance = get_distance_paired( population_[i], current_pair_[0], current_pair_[1]);	 
		if( distance > max_dis ) max_dis = distance;
	} 
	return max_dis;
}

core::Real ClassicAbinitio::get_distance_paired( core::pose::Pose &pose , Size pos1, Size pos2 )
{
	  if( pos1 <= pose.total_residue() && pos2 <= pose.total_residue() )
		{
			std::string atom_name1 = pose.residue_type( pos1 ).name1() == 'G' ? "CA": "CB";
			std::string atom_name2 = pose.residue_type( pos2 ).name1() == 'G' ? "CA": "CB";
			numeric::xyzVector<core::Real> v1 = pose.residue( pos1 ).atom( atom_name1 ).xyz();
			numeric::xyzVector<core::Real> v2 = pose.residue( pos2 ).atom( atom_name2 ).xyz();
			core::Real distance = v1.distance(v2);
			return distance;
		}
		else return 0;
}
void ClassicAbinitio::contact_dis( const std::vector< core::pose::Pose > population_ , protocols::abinitio::contact_predicted &RaptorXcontact_ )
{
	vector<vector<int> > contact_index_ = RaptorXcontact_.get_contact_index();
	for( unsigned int i = 0; i < contact_index_.size(); i ++)
	{
		vector<int> current_pair = contact_index_[i];
		//distance of native structure
		Real distance_native =  get_distance_paired( native_pose_temp , current_pair[0], current_pair[1]);
		
		
		count_num( population_,  current_pair);//calculate the distribution of the population
		double entropy = entropy_contact_dis( population_.size() );//calculate the entropy of the contact distribution
		char fname[256];
		//contact_Num_index1_index2_Bin
		sprintf(fname , "./contact_dis/contact%d_%d-%d_%.2f.csv", i , current_pair[0], current_pair[1], distance_native);
		string savepath(fname);
		ofstream fp(savepath.c_str() , std::ios::app);
		for(unsigned int i = 0; i < contact_bin_count.size() ; i ++)
		{
			fp << contact_bin_count[i] << "," ;
		}
		fp << entropy << ",\n";
	}
}
void ClassicAbinitio::initiaize_contact_bin( int start_, int end_, int interval_)
	{
		//4-5,5-6,6-7...
		int length_ = (end_ - start_ + 1) / interval_;
		for(int i = 0; i < length_; i ++)
		{
			contact_bin_count.push_back( 0 );
		}
	}
	
void ClassicAbinitio::reset_contact_bin( )
	{
		for(unsigned int i = 0; i < contact_bin_count.size() ; i ++)
		{
			contact_bin_count[i] = 0;
		}
	}
	
bool ClassicAbinitio::count_num( std::vector< core::pose::Pose > population_, vector<int> current_pair_ )
	{
		reset_contact_bin( );		
//count the distance distribution of the contact pair for the population
		for(unsigned int i = 0; i < population_.size() ; i ++)
		{
			double distance_flo = get_distance_paired( population_[i] , current_pair_[0], current_pair_[1]);
			int dis_bin = floor( distance_flo ) - 4;
			//count and record to contact_bin_count
			if(dis_bin < 0) dis_bin = 0;
			int bin_length = contact_bin_count.size();
			if(dis_bin < bin_length && dis_bin >= 0)
				{
					contact_bin_count[ dis_bin ] ++;
				}
			else if( dis_bin >= bin_length )
				{
					contact_bin_count[ bin_length - 1] ++;
				}
			else
				{
					return false;
				}		
		}
		return true;
	}
	int ClassicAbinitio::prob_chose_bin( double contact_prob )
		{
			//calculate the number of bins less than 8 angstrom
			double prob_less_8 = contact_prob / 4;
			double prob_greater_8 = ( 1 - contact_prob) / (contact_bin_count.size() - 4);
			vector<double> bin_prob;
			for(unsigned int i = 1;i <= contact_bin_count.size(); i ++)
			{
				if(i <= 4)
					{
						bin_prob.push_back( prob_less_8 );
//						cout<< prob_less_8 << "\t";
					}
				else
					{
						bin_prob.push_back( prob_greater_8 );
//						cout<< prob_greater_8 << "\t";
					}
//					cout << "\n";
			}
			roulette_wheel_select bin_roulette( bin_prob );
			int bin_indx = bin_roulette.selection();
			return bin_indx;
		}
	vector<vector<int> > ClassicAbinitio::set_frag_range( std::vector<int> pertb_pair_ , int fragAss_size_ , core::pose::Pose pose)
	{
		int pose_length = pose.total_residue();
		
		vector<vector<int> > perturb_range;

		for(unsigned int i = 0; i < pertb_pair_.size(); i ++)
		{
			vector<int> purturb_range_part;
			int pos = pertb_pair_[i];
			if( (pos - fragAss_size_ / 2) >= 1 &&  (pos + fragAss_size_ / 2) <= pose_length)
				{
					purturb_range_part.push_back( pos - fragAss_size_ / 2 );
					purturb_range_part.push_back( pos + fragAss_size_ / 2 );
				}
			else if( (pos - fragAss_size_ / 2) < 1 &&(pos + fragAss_size_ / 2) <= pose_length )
				{
					purturb_range_part.push_back( 1 );
					purturb_range_part.push_back( fragAss_size_ + 1 );
				}
			else if( (pos - fragAss_size_ / 2) >= 1 && (pos + fragAss_size_ / 2) > pose_length )
				{
					purturb_range_part.push_back( pos - fragAss_size_ / 2 );
					purturb_range_part.push_back( pose_length );
				}
			else
				{
					std::cout << "ClassicAbinitio::set_frag_range: the value of fragAss_size_ is greater than length of length!";
				}
				perturb_range.push_back( purturb_range_part );
		}
		return perturb_range;
	}
	double ClassicAbinitio::entropy_contact_dis(Size population_size)
	{
		double entropy_ = 0;
		for(Size i = 0; i < contact_bin_count.size(); i ++)
		{
			double proba = (double) contact_bin_count[i] / population_size ;
			if(contact_bin_count[i] ==0) proba = 0.0000001;
			entropy_ =  entropy_ + proba* log(proba);
		}
		return -entropy_;
	}
} //abinitio
} //protocols
