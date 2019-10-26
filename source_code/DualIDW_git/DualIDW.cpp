//!-----------------------------------------------------------------------
//!
//! This is a SGeMS (the Stanford Geostatistical Modeling Software) plug-in, which implements 
//! a dual inverse distance weighting (DIDW) interpolation. 

//! This file may be distributed and/or modified under the terms of the
//! GNU General Public License version 2 as published by the Free Software
//! Foundation and appearing in the file LICENSE.GPL included in the
//! packaging of this file.

//! Title of Paper Submitted: 
//	Dual inverse distance weighting interpolation with locally varying exponents
//
//
//
//!-----------------------------------------------------------------------
//!
//!	DUAL INVERSE DISTANCE WEIGHTING
//!	************************************************************************
//!
//!	This plug-in implements a systematical enhancement of IDW to account both for the 
//!	accounts both for the clustering of nearby samples and for the distance to the unknown location.
//!
//!	MODIFIED FROM:
//!	The Stanford Geostatistical Modeling Software (SGeMS, url: http://sgems.sourceforge.net)
//!
//!-----------------------------------------------------------------------


#include "basicStatistics.h"

#include "DualIDW.h"

int DualIDW::nInitialized_ = 0;

DIDW_V1_API std::ofstream &fcout = (std::ofstream()); 

#ifdef OUT_PUT_TO_SCREEN
std::ostream& DualIDW::_cout_file(std::cout);
#else
std::ofstream& DualIDW::_cout_file(fcout);
#endif


OptParameters DualIDW::optPars_ = 
{
	1,						// IDW or Kriging, 0- IDW; 1- Kriging

	2.0,					// default value of sample-to-estimate Power
	2.0,					// default value of sample-to-sample Power
	0,						// Is Anisotropic Distance used? 0- NO; 1- YES;

	std::vector<int>(0),	// Global OPT, 1- yes; 0- no; for [D2E_power; D2D_power;]
	std::vector<int>(0),	// Local  OPT, 1- yes; 0- no; for [D2E_power; D2D_power;]

	// 3 pars: the number of possible values; Min Weight; Max Weight
	std::vector<double>(0), // Sample-to-estimation powers
	std::vector<double>(0), // Sample-to-sample powers 

	std::vector<int>(0),	//	Optimization Goal (only for global OPT)

	false,					// Need output test details? ? 0- NO; 1- YES;
	false,					// Use the D2E power value as D2D power? 0- NO; 1- YES;
};

DualIDW::DualIDW()
{		
	Kconstraints_ = 0;
	combiner_ = 0;
	simul_grid_ = 0;
	neighborhood_ = 0;
	min_neigh_ = 0;
	rhs_covar_blk_ = 0;
	rhs_covar_ = 0;

	_ni_KrgPars = 0;

	var_prop_PowerS2U_ = 0;
	var_prop_PowerS2S_ = 0;

	var_prop_TrueError_ = 0;
	var_prop_AbsError_  = 0;
	var_prop_S2sDis_  = 0;
	mean_prop_S2sDis_  = 0;
	var_prop_SampNum_  = 0;

	dPowerSample2Est_OptResult = 2;
	dPowerSample2Sap_OptResult = 2;


#ifdef OUT_PUT_TO_SCREEN
	//_cout_file = std::cout;
#else
	std::string sAppPath = QCoreApplication::applicationDirPath().toStdString();
	std::string sPath = sAppPath + "\\CalProcess\\ResultLogs.dat";

	_cout_file.open(sPath.c_str(), ios::out|ios::app);

	if(_cout_file.fail())
	{
		std::cout<<" Open files error: ResultLogs.dat"<<std::endl;		
	}
#endif
}

DualIDW::~DualIDW()
{
#ifdef OUT_PUT_TO_SCREEN
	//
#else
	if(_cout_file.fail())
	{
		std::cout<<" Open files error: ResultLogs.dat"<<std::endl;		
	}
	else if ((_cout_file).is_open())
	{
		(_cout_file).close();
	}
#endif

	if( Kconstraints_ )
		delete Kconstraints_;

	if( combiner_ )
		delete combiner_;
}

Named_interface* DualIDW::create_new_interface( std::string& ) 
{

	return  (new DualIDW);
}

void DualIDW::clean( const std::string& prop ) 
{
	simul_grid_->remove_property( prop );
}

// get parameters form the interface 
bool DualIDW::initialize( const Parameters_handler* parameters, Error_messages_handler* errors )
{	
	//	return Kriging::initialize(parameters, errors);
	_basic_KrgPars = (Parameters_handler*)parameters;

	std::string sAppPath = QCoreApplication::applicationDirPath().toStdString();
	_sModelingPath = sAppPath + "\\CalProcess\\DualIDW.par"; 

	if (-1 == access(_sModelingPath.c_str(), 0))
	{
		return false;
	}

	//-----------------
	// Extract the parameters input by the user from the parameter handler
	std::string simul_grid_name = parameters->value( "Grid_Name.value" );
	errors->report( simul_grid_name.empty(), 
		"Grid_Name", "No grid selected" );
	property_name_ = parameters->value( "Property_Name.value" );
	errors->report( property_name_.empty(), 
		"Property_Name", "No property name specified" );


	// Get the simulation grid from the grid manager
	if( !simul_grid_name.empty() ) {
		bool ok = geostat_utils::create( simul_grid_, simul_grid_name,
			"Grid_Name", errors );
		if( !ok ) return false;
	}
	else 
		return false;

	gridTempRegionSelector_.set_temporary_region(
		parameters->value( "Grid_Name.region" ), simul_grid_);


	std::string harddata_grid_name = parameters->value( "Hard_Data.grid" );
	errors->report( harddata_grid_name.empty(), 
		"Hard_Data", "No hard data specified" );
	harddata_property_name_ = parameters->value( "Hard_Data.property" );
	errors->report( harddata_property_name_.empty(), 
		"Hard_Data", "No property name specified" );

	// Get the harddata grid from the grid manager
	if( !harddata_grid_name.empty() ) {
		bool ok = geostat_utils::create( harddata_grid_, harddata_grid_name,
			"Hard_Data", errors );
		if( !ok ) return false;
	}
	else 
		return false;

	// If the hard data is on the same grid than the
	// estimation grid, than it cannot be set to a region others than the one
	// already set on that grid
	// The is only used if the neighborhood is to consider only data within a region
	if(harddata_grid_ !=  simul_grid_)
		hdgridTempRegionSelector_.set_temporary_region(
		parameters->value( "Hard_Data.region" ),harddata_grid_ );

	int max_neigh = 
		String_Op::to_number<int>( parameters->value( "Max_Conditioning_Data.value" ) );

	min_neigh_ = 
		String_Op::to_number<int>( parameters->value( "Min_Conditioning_Data.value" ) );
	errors->report( min_neigh_ >= max_neigh, 
		"Min_Conditioning_Data", "Min must be less than Max" );

	// set-up the covariance
	bool init_cov_ok =
		geostat_utils::initialize_covariance( &covar_, "Variogram", 
		parameters, errors );
	if( !init_cov_ok ) return false;
	do_block_kriging_ = parameters->value("do_block_kriging.value") == "1";
	if( do_block_kriging_ ) {

		RGrid* block_grid = dynamic_cast<RGrid*>(simul_grid_);
		if(!block_grid) {
			errors->report("Grid_Name","Must be a Cartesian grid to use the Block Kriging Option");
			return false;
		}

		nblock_pts_[0] = String_Op::to_number<int>(parameters->value( "npoints_x.value" ) );
		nblock_pts_[1] = String_Op::to_number<int>(parameters->value( "npoints_y.value" ) );
		nblock_pts_[2] = String_Op::to_number<int>(parameters->value( "npoints_z.value" ) );

		errors->report(nblock_pts_[0] <= 0,"npoints_x","At least one point is necessary");
		errors->report(nblock_pts_[1] <= 0,"npoints_y","At least one point is necessary");
		errors->report(nblock_pts_[2] <= 0,"npoints_z","At least one point is necessary");
		if(!errors->empty()) return false;

		rhs_covar_blk_ = new Block_covariance<Location>(covar_,nblock_pts_,block_grid->geometry()->cell_dims());
	}
	else rhs_covar_ = new Covariance<Location>(covar_);

	/*  bool init_blk_cok_ok =  initialize_blk_covariance( &blk_covar_, 
	nblock_pts_,covar_, simul_grid_->geometry->cell_dims());
	if( !init_blk_cok_ok ) {
	errors->errors("Variogram",
	"Bad intialization of the block variogram, check the discretization and/or variogram");
	return false;
	}
	*/

	// More complicated stuff related to the cokriging covaraince setup 

	hdata_prop_ =
		harddata_grid_->select_property( harddata_property_name_ );
	if( !hdata_prop_ ) {
		std::ostringstream error_stream;
		error_stream << harddata_grid_name <<  " does not have a property called " 
			<< harddata_property_name_;
		errors->report( "Hard_Data", error_stream.str() );
	}


	//-------------
	// Set up the search neighborhood

	GsTLTriplet ellips_ranges;
	GsTLTriplet ellips_angles;
	bool extract_ok = 
		geostat_utils::extract_ellipsoid_definition( ellips_ranges, ellips_angles,
		"Search_Ellipsoid.value",
		parameters, errors );
	if( !extract_ok ) return false;

	extract_ok = geostat_utils::is_valid_range_triplet( ellips_ranges );
	errors->report( !extract_ok,
		"Search_Ellipsoid",
		"Ranges must verify: major range >= " 
		"medium range >= minor range >= 0" );
	if( !extract_ok ) return false;


	harddata_grid_->select_property(harddata_property_name_);
	if( dynamic_cast<Point_set*>(harddata_grid_) ) {
		neighborhood_ = SmartPtr<Neighborhood>(
			harddata_grid_->neighborhood( ellips_ranges, ellips_angles, &covar_, true ) );
	} 
	else {
		neighborhood_ =  SmartPtr<Neighborhood>(
			harddata_grid_->neighborhood( ellips_ranges, ellips_angles, &covar_ ));
	}
	neighborhood_->select_property( harddata_property_name_ );
	neighborhood_->max_size( max_neigh );

	geostat_utils::set_advanced_search(neighborhood_.raw_ptr(), 
		"AdvancedSearch", parameters, errors);

	estimation_weights_.reserve( 2 * max_neigh );



	//-----------------
	// The kriging constraints and combiner

	//  std::string kriging_type = parameters->value( "Kriging_Type.type" );
	//  set_kriging_parameters( kriging_type, parameters, errors );

	geostat_utils::KrigTagMap tags_map;
	tags_map[ geostat_utils::SK  ] = "Kriging_Type/parameters.mean";
	tags_map[ geostat_utils::KT  ] = "Kriging_Type/parameters.trend";
	tags_map[ geostat_utils::LVM ] = "Kriging_Type/parameters.property";
	//  tags_map[ geostat_utils::LVM ] = "Kriging_Type/prop_mean_grid;harddata_grid_name;Kriging_Type/prop_mean_hdata";

	//geostat_utils::Kriging_type ktype = 
	ktype_ = 
		geostat_utils::kriging_type( "Kriging_Type.type", parameters, errors );
	geostat_utils::initialize( ktype_, combiner_, Kconstraints_,
		tags_map,
		parameters, errors,
		simul_grid_ );

	if( !errors->empty() )
		return false;

	return true;
}

// // calculate the kriging error variance and the variance of the estimete
template<
	class SymmetricMatrix,
	class MatVector,
	class location,
	class Covariance,
	class Vector
>
double DualIDW::cal_kring_variance_by_weights0(
	SymmetricMatrix& A, 
	MatVector& b,
	Vector& weights,
	const location& center,
	Covariance& covar,
	double* pdVarianceOfEstimate 
	)  
	{
		DEBUG_PRINT_KRIGING_SYSTEM( A,b);
		int nb_conditioning_data = b.size() - 1;

		double dKrgVariance = 0.0;
		double dPart1 = 0.0, dPart2 = 0.0, dPart3 = 0.0;

		dKrgVariance += covar(Location(0, 0, 0), Location(0, 0, 0));
		dPart1 = covar(Location(0, 0, 0), Location(0, 0, 0));

		int dim = A.num_rows();
		// Since the LU algorithm of TNT does not use the fact that A is symetric,
		// we have to copy the upper part of A to its lower part.

		int i=1;
		int j;

		typedef typename Neighborhood::const_iterator InputIterator;

		while(i < dim)
		{
			j=1;

			while(j < dim)
			{
				dPart2 += ( weights[i-1]*weights[j-1]*A(i,j) );

				j++;
			}

			i++;
		}

		i = 1;

		while(i < dim)
		{
			dPart3 -= (weights[i-1]*2.0*b(i));

			i++;
		}		

		dKrgVariance = dPart1 + dPart2 + dPart3;

		if(pdVarianceOfEstimate)
		{
			*pdVarianceOfEstimate = dPart2;
		}

		//if (dPart2 < 0 && dKrgVariance>=0)
		//{
		//	std::cout<< "dPart2 < 0 && dKrgVariance>=0 " <<std::endl;
		//}

		return dKrgVariance;
	}

	int DualIDW::EstimationDIDW( GsTL_project* proj)
	{
		// output basic information
		{
			_cout_file<<"Currently used method: ";

			Parameters_handler* parameters = (Parameters_handler*)_basic_KrgPars.raw_ptr();

			std::string harddata_grid_name = parameters->value( "Hard_Data.grid" );
			_cout_file<<"harddata_data = "<<harddata_grid_name<<"."<<harddata_property_name_<<std::endl;

			std::string sim_grid_name = parameters->value( "Grid_Name.value" );
			_cout_file<<"sim_data = "<<sim_grid_name<<"."<<property_name_<<std::endl;

			this->_cout_file 
				<< "\nbUsePower_D2E_AS_D2D = \t"<< DualIDW::optPars_.bUsePower_D2U_AS_D2D
				<< "\nOptType_G = \t"<< DualIDW::optPars_.nOptType_G[0]<<","<< DualIDW::optPars_.nOptType_G[1]
			<< "\tOptType_L = \t"<< DualIDW::optPars_.nOptType_L[0]<<","<< DualIDW::optPars_.nOptType_L[1]
			<< "\ndPowerSample2Est = \t"<< DualIDW::optPars_.dPowerSample2Est
				<< "\tdPowerSample2Sap = \t"<< DualIDW::optPars_.dPowerSample2Sap

				<< "\nbUseAnisotropicDistance = \t"<< DualIDW::optPars_.bUseAnisotropicDistance
				<< std::endl;
		}

		// get the simulated grid 
		Geostat_grid* simul_grid = simul_grid_; 


		GsTLGridProperty* prop = 0;
		if (prop = simul_grid->property(property_name_))
		{
			simul_grid->select_property(property_name_);
		}
		else
		{
			// create the property
			appli_message("creating new property: " << property_name_ << "..." );
			prop = geostat_utils::add_property_to_grid( simul_grid, property_name_ );
			simul_grid->select_property( prop->name() );
		}

		GsTLGridProperty* var_prop = 0;
		if (var_prop = simul_grid->property(property_name_ + "_krig_var"))
		{

		}
		else
		{
			// create property for kriging variance
			std::string var_prop_name = prop->name() + "_krig_var";
			var_prop = geostat_utils::add_property_to_grid( simul_grid, var_prop_name );
		}

		// the used D-U power
		if (var_prop_PowerS2U_ = simul_grid->property(property_name_ + "_PowerS2U"))
		{

		}
		else
		{
			std::string var_prop_name = prop->name() + "_PowerS2U";
			var_prop_PowerS2U_ = geostat_utils::add_property_to_grid( simul_grid, var_prop_name );
		}

		if (var_prop_PowerS2S_ = simul_grid->property(property_name_ + "_PowerS2S"))
		{

		}
		else
		{
			std::string var_prop_name = prop->name() + "_PowerS2S";
			var_prop_PowerS2S_ = geostat_utils::add_property_to_grid( simul_grid, var_prop_name );
		}

		if (var_prop_TrueError_ = simul_grid->property(property_name_ + "_TrueError"))
		{

		}
		else
		{
			std::string var_prop_name = prop->name() + "_TrueError";
			var_prop_TrueError_ = geostat_utils::add_property_to_grid( simul_grid, var_prop_name );
		}


		if (var_prop_AbsError_ = simul_grid->property(property_name_ + "_AbsError"))
		{

		}
		else
		{
			std::string var_prop_name = prop->name() + "_AbsError";
			var_prop_AbsError_ = geostat_utils::add_property_to_grid( simul_grid, var_prop_name );
		}

		if (var_prop_S2sDis_ = simul_grid->property(property_name_ + "_S2sDisVar"))
		{

		}
		else
		{
			std::string var_prop_name = prop->name() + "_S2sDisVar";
			var_prop_S2sDis_ = geostat_utils::add_property_to_grid( simul_grid, var_prop_name );
		}

		if (mean_prop_S2sDis_ = simul_grid->property(property_name_ + "_S2sDisMean"))
		{

		}
		else
		{
			std::string var_prop_name = prop->name() + "_S2sDisMean";
			mean_prop_S2sDis_ = geostat_utils::add_property_to_grid( simul_grid, var_prop_name );
		}

		if (var_prop_SampNum_ = simul_grid->property(property_name_ + "_SampNum"))
		{

		}
		else
		{
			std::string var_prop_name = prop->name() + "_SampNum";
			var_prop_SampNum_ = geostat_utils::add_property_to_grid( simul_grid, var_prop_name );
		}

		// get the permissible powers 
		// std::vector<double> dPowersS2U, dPowersS2S;
		PowerOptRanges(dPowersS2U_, dPowersS2S_);

		Whole_OptimizationByEnumeration(simul_grid, var_prop);
		this->_cout_file << "Global Opt"<<std::endl;

		return 0;
	}

	double DualIDW::EstimationWithAllPoints(Geostat_grid* simul_grid, GsTLGridProperty* var_prop, 
		interpolation_accuracy_measurement_inx *p_IAM_inx, bool bOutPut /*= false*/)
	{
		// Set up a progress notifier	
		int total_steps = simul_grid->size();
		int frequency = std::max( total_steps / 20, 1 );
		SmartPtr<Progress_notifier> progress_notifier = 0;

		// show the progress 
		if(DualIDW::optPars_.nOptType_G[0] == 0 && DualIDW::optPars_.nOptType_G[1] == 0)
		{
			if(DualIDW::optPars_.nOptType_L[0] == 0 && DualIDW::optPars_.nOptType_L[1] == 0)
			{
				progress_notifier = utils::create_notifier( "Running Estimation", total_steps, frequency );
			}
			else if(DualIDW::optPars_.nOptType_L[0] == 1 && DualIDW::optPars_.nOptType_L[1] == 1)
			{
				progress_notifier = utils::create_notifier( "Running Local Optimization", total_steps, frequency );
			}
			else 
			{
				progress_notifier = utils::create_notifier( "Running Semi-Local Optimization", total_steps, frequency );
			}
		}

		// dvWeightsCR: the correlation coefficient between the OK and DIDW weights 
		std::vector<double> dvMeasuredVals, dvEstimatedVals, dvTrueErrors, dvWeightsCR, dvKVs;
		std::vector<double> dvS2S_DisVariance, dvS2S_DisMean;

		GsTLGridProperty *pReferenceProperty = simul_grid->property(hdata_prop_->name());

		typedef Geostat_grid::iterator iterator;
		iterator begin = simul_grid->begin();
		iterator end = simul_grid->end();

		begin_ = &(*begin);

		for( ; begin != end; ++begin )   
		{
			if( progress_notifier &&  !progress_notifier->notify() ) 
			{
				clean( property_name_ );
				return 1;
			}

			int center_node_id = begin->node_id();

			//_cout_file<<"Point to be Estimated: "<<begin->location()<<"\t"<<begin->property_value()<<std::endl;

			neighborhood_->includes_center(false); // the center of the neighborhood will not be viewed as a sample 

			neighborhood_->find_neighbors( *begin );

			//////////////////////[2018-4-10,16:42 O]//////////////+
			//TODO: calculate the variance of the numbers of neighborhood samples in the four quadrants 
			// SampleNum
			//{
			//	p_IAM_inx = 0;
			//	Search_filter* pSF = neighborhood_->search_neighborhood_filter();

			//	Octant_search_filter *pOSF = (Octant_search_filter*)pSF;

			//	std::vector<int> octant_registrar_ = pOSF->octant_registrar();

			//	in a search plane, only the four quadrants are needed
			//	octant_registrar_.erase(octant_registrar_.begin());
			//	octant_registrar_.erase(octant_registrar_.begin());
			//	octant_registrar_.erase(octant_registrar_.begin());
			//	octant_registrar_.erase(octant_registrar_.begin());


			//	double dVar = 0;
			//	double dSum = 0, dAvg = 0;
			//	{
			//		int nTotalSize = octant_registrar_.size();
			//		dSum = (double)std::accumulate(octant_registrar_.begin(), octant_registrar_.end(), 0);
			//		dAvg = dSum/nTotalSize;
			//		for (int ii = 0; ii<octant_registrar_.size(); ii++)
			//		{
			//			dVar += (octant_registrar_[ii] - dAvg)*(octant_registrar_[ii] - dAvg);
			//		}

			//		dVar = dVar/(nTotalSize - 1);
			//	}

			//	begin->set_property_value( neighborhood_->size() );

			//	var_prop->set_value( dVar, begin->node_id() );

			//	continue;
			//}

			//////////////////////[2018-4-10,16:42 O]//////////////-

			if( neighborhood_->size() < min_neigh_ )  
			{
				begin->set_not_informed();
				var_prop->set_not_informed(center_node_id);

				continue;
			}

			if(!neighborhood_->is_valid()) 
			{
				begin->set_not_informed();
				var_prop->set_not_informed(center_node_id);

				continue;
			};


			double variance = EstimationWithOnePoint();

			if(variance >= 0.0) 
			{
				// the kriging system could be solved
				double estimate = (*combiner_)( estimation_weights_.begin(), 
					estimation_weights_.end(),
					*(neighborhood_.raw_ptr()) );

				begin->set_property_value( estimate );

				var_prop->set_value( variance, begin->node_id() );

				var_prop_SampNum_->set_value( neighborhood_->size(), begin->node_id() );

				//_cout_file<<"Estimated Point Result: "
				// <<"EV:"<<"\t"<<estimate
				// <<"KV:"<<"\t"<<variance<<std::endl;

				dvKVs.push_back(variance);

				if(pReferenceProperty)
				{
					dvMeasuredVals.push_back( pReferenceProperty->get_value(begin->node_id()) );

					dvEstimatedVals.push_back(estimate);

					double dTrueError = estimate - pReferenceProperty->get_value(begin->node_id());
					dvTrueErrors.push_back(dTrueError);

					var_prop_TrueError_->set_value( dTrueError, begin->node_id() );
					var_prop_AbsError_->set_value( fabs(dTrueError), begin->node_id() );

					double dS2S_DisVariance = GetVariance(dvSample2SampleDis_);
					double dS2S_DisMean = GetAverage(dvSample2SampleDis_);
					if(dS2S_DisVariance >= 0)
					{
						dvS2S_DisVariance.push_back(dS2S_DisVariance);
						var_prop_S2sDis_->set_value(dS2S_DisVariance, begin->node_id() );

						dvS2S_DisMean.push_back(dS2S_DisMean);
						mean_prop_S2sDis_->set_value(dS2S_DisMean, begin->node_id() );
					}

					if(estimation_weights_.size() >= 10)	// only for the case that sample size is no smaller than 10
					{
						// calculate the corelation relation between the OK and DIDW weights 
						ori_kriging_weights_.resize(estimation_weights_.size());

						double dTemp = 0;

						if( (dTemp = GetCC(ori_kriging_weights_, estimation_weights_)) < 9.8E12)
						{
							dvWeightsCR.push_back( dTemp);
						}
					}
				}
				else
				{
					var_prop_TrueError_->set_not_informed(center_node_id);
					var_prop_AbsError_->set_not_informed(center_node_id);
					var_prop_S2sDis_->set_not_informed(center_node_id);
					mean_prop_S2sDis_->set_not_informed(center_node_id);
					var_prop_SampNum_->set_not_informed(center_node_id);
				}
			}
			else 
			{
				begin->set_not_informed();
				var_prop->set_not_informed(center_node_id);

				var_prop_TrueError_->set_not_informed(center_node_id);
				var_prop_AbsError_->set_not_informed(center_node_id);
				var_prop_S2sDis_->set_not_informed(center_node_id);
				mean_prop_S2sDis_->set_not_informed(center_node_id);
				var_prop_SampNum_->set_not_informed(center_node_id);
			}
		}

		double dObjValue = 9.9e20;
		if(p_IAM_inx && dvTrueErrors.size() > 1)
		{
			double dErrorMean  = GetAverage(dvTrueErrors);
			double dErrorVariance   = GetVariance(dvTrueErrors);
			dErrorVariance = sqrt(dErrorVariance);// standard dev. 

			double dMTE  = GetMTE (dvMeasuredVals, dvEstimatedVals, 0, false);
			double dMAE  = GetMAE (dvMeasuredVals, dvEstimatedVals, 0, false);
			double dRMSE = GetRMSE(dvMeasuredVals, dvEstimatedVals, 0, false);
			double dCC   = GetCC  (dvMeasuredVals, dvEstimatedVals);
			double dCV   = GetCovariance (dvMeasuredVals, dvEstimatedVals);
			//double dCV   = GetRankCC (dvMeasuredVals, dvEstimatedVals);

			double dMinError = *(std::min_element(dvTrueErrors.begin(), dvTrueErrors.end()));
			double dMaxError = *(std::max_element(dvTrueErrors.begin(), dvTrueErrors.end()));


			double dMeanWeightsCR  = GetAverage(dvWeightsCR);// Mean correlation coefficient

			double dMeanKV  = GetAverage(dvKVs);// Mean Kriging variance


			p_IAM_inx->iam.nNumberOfEstimates = dvTrueErrors.size();
			p_IAM_inx->iam.dMeanS2S_DisVariance = GetAverage(dvS2S_DisVariance);
			p_IAM_inx->iam.dMeanS2S_DisMean = GetAverage(dvS2S_DisMean);


			p_IAM_inx->dPowerSample2Est = DualIDW::optPars_.dPowerSample2Est;
			p_IAM_inx->dPowerSample2Sap = DualIDW::optPars_.dPowerSample2Sap;

			p_IAM_inx->iam.dErrorMean = dErrorMean;
			p_IAM_inx->iam.dErrorVariance = dErrorVariance;

			p_IAM_inx->iam.dMTE = dMTE;
			p_IAM_inx->iam.dMAE = dMAE;
			p_IAM_inx->iam.dRMSE = dRMSE;
			p_IAM_inx->iam.dCC = dCC;

			p_IAM_inx->iam.dMinError = dMinError;
			p_IAM_inx->iam.dMaxError = dMaxError;

			p_IAM_inx->iam.dMeanWeightsCR = dMeanWeightsCR;
			p_IAM_inx->iam.dMeanKV = dMeanKV;
			p_IAM_inx->iam.dCV = dCV;

			std::vector<int>  &nvGlobalOptGoalType = DualIDW::optPars_.nvGlobalOptGoalType;

			dObjValue = 0;

			for (int iTypeInx = 0; iTypeInx<DualIDW::optPars_.nvGlobalOptGoalType.size(); iTypeInx++)
			{
				if(DualIDW::optPars_.nvGlobalOptGoalType[iTypeInx])
				{
					switch(iTypeInx)
					{
					case 0://ErrorMean
						dObjValue += fabs(p_IAM_inx->iam.dErrorMean);
						break;
					case 1://ErrorVariance
						dObjValue += p_IAM_inx->iam.dErrorVariance;
						break;
					case 2://MinError
						dObjValue += fabs(p_IAM_inx->iam.dMinError);
						break;
					case 3://MaxError
						dObjValue += fabs(p_IAM_inx->iam.dMaxError);
						break;

					case 4://MTE
						dObjValue += fabs(p_IAM_inx->iam.dMTE);
						break;
					case 5://MAE
						dObjValue += (p_IAM_inx->iam.dMAE);
						break;
					case 6://RMSE
						dObjValue += (p_IAM_inx->iam.dRMSE);
						break;

					case 7://  -CR
						dObjValue += (-p_IAM_inx->iam.dCC);
						break;
					case 8://  -CV
						dObjValue += (p_IAM_inx->iam.dCV);
						break;

					case 9:// - MeanWeightsCR
						dObjValue += (-p_IAM_inx->iam.dMeanWeightsCR);
						break;
					case 10://MeanKV
						dObjValue += p_IAM_inx->iam.dMeanKV;
						break;
					}
				}						
			}
		}

		//TEST CODE: add by O on [2018-4-15]+
		//FOR: output the interpolation accuracy

		if (p_IAM_inx && bOutPut)
		{
			this->_cout_file 
				//<<"IsAD\t"<<"PowerS2U\t"<<"PowerS2S\t"<<"MeanError\t"<<"Error_Variance\t"<<"MAE\t"<<"RMSE\t"<<"CR\n"
				<< DualIDW::optPars_.bUseAnisotropicDistance<<"\t"
				<< p_IAM_inx->dPowerSample2Est<<"\t" 
				<< p_IAM_inx->dPowerSample2Sap<<"\t"

				<<p_IAM_inx->iam.dErrorMean<<"\t"
				<<p_IAM_inx->iam.dErrorVariance<<"\t"
				<<p_IAM_inx->iam.dMinError<<"\t"
				<<p_IAM_inx->iam.dMaxError<<"\t"

				<<p_IAM_inx->iam.dMTE<<"\t"
				<<p_IAM_inx->iam.dMAE<<"\t"
				<<p_IAM_inx->iam.dRMSE<<"\t"

				<<p_IAM_inx->iam.dCC<<"\t"
				<<p_IAM_inx->iam.dCV<<"\t"

				<<p_IAM_inx->iam.dMeanWeightsCR<<"\t"
				<<p_IAM_inx->iam.dMeanKV<<"\t"

				<< std::endl;
		}
		//TEST CODE: add by O on [2018-4-15]-

		return dObjValue;
	}



	// get all of the possible value ranges for the two powers to be opted
	double DualIDW::PowerOptRanges(std::vector<double> &dPowersS2U, std::vector<double> &dPowersS2S)
	{
		dPowersS2U.clear();

		//if(DualIDW::optPars_.bPowersOptimized[0] == 1)
		{
			double dMinMaxPowerS2U[2] = {DualIDW::optPars_.dPowerS2U_ValueRange[1],DualIDW::optPars_.dPowerS2U_ValueRange[2]};
			int nTotalStepS2U = DualIDW::optPars_.dPowerS2U_ValueRange[0];
			double dStepPowerS2U = (dMinMaxPowerS2U[1] - dMinMaxPowerS2U[0])/nTotalStepS2U;
			for (int i = 0; i<=nTotalStepS2U; i++)
			{
				dPowersS2U.push_back(dMinMaxPowerS2U[0] +  double(i)*dStepPowerS2U);
			}
		}

		//dPowersS2U.push_back(DualIDW::optPars_.dPowerSample2Est);

		dPowersS2S.clear();
		//if(DualIDW::optPars_.bPowersOptimized[1] == 1)
		{
			double dMinMaxPowerS2S[2] = {DualIDW::optPars_.dPowerS2S_ValueRange[1],DualIDW::optPars_.dPowerS2S_ValueRange[2]};

			int nTotalStepS2S = DualIDW::optPars_.dPowerS2S_ValueRange[0];

			double dStepPowerS2S = (dMinMaxPowerS2S[1] - dMinMaxPowerS2S[0])/nTotalStepS2S;

			for (int i = 0; i<=nTotalStepS2S; i++)
			{
				dPowersS2S.push_back(dMinMaxPowerS2S[0] +  double(i)*dStepPowerS2S);
			}
		}

		//dPowersS2S.push_back(DualIDW::optPars_.dPowerSample2Sap);

		return 1.0;
	}

	// get the distances needed for the estimation
	void DualIDW::GetDistance4Estimation()
	{
		int status = -1;
		Geovalue &begin = (*begin_);
		int center_node_id = begin.node_id();

		double variance = -1.0E10;

		std::vector<double>& weights = estimation_weights_;
		double& kriging_variance = variance;
		const Location& center = begin.location();
		const Neighborhood& neighbors = *(neighborhood_.raw_ptr());
		Covariance<Location>& covar= covar_;
		Covariance<Location>& covar_rhs = *rhs_covar_;
		KrigingConstraints& Kconstraints = *Kconstraints_;

		MatrixLib::Symmetric_matrix& A = A_;	// Kriging system: Ax = b
		MatrixLib::Vector &b = b_;

		typedef matrix_lib_traits< GSTL_TNT_lib > MatrixLib;

		// dPowerSample2Est; dPowerSample2Sap
		double dPowerSample2Est = 2.0, dPowerSample2Sap = 2.0;
		dPowerSample2Est = DualIDW::optPars_.dPowerSample2Est;
		dPowerSample2Sap = DualIDW::optPars_.dPowerSample2Sap;

		MatrixLib::Symmetric_matrix sample_to_sample_dis_MatrixA;
		sample_to_sample_dis_MatrixA.resize(neighbors.size(), 0);

		MatrixLib::Vector &sample_to_est_dis_VecB = sample_to_est_dis_VecB_;
		sample_to_est_dis_VecB.resize((int)neighbors.size());
		std::vector<double> &dvSample2SampleDis = dvSample2SampleDis_;  
		dvSample2SampleDis.clear();
		double dDisSumTotal = 0, dDisSum1 = 0; 

		// is the variogram distance used?
		if(DualIDW::optPars_.bUseAnisotropicDistance) 				
		{
			int i=1;
			int j;

			Neighborhood::const_iterator first_neigh = neighbors.begin();
			Neighborhood::const_iterator last_neigh  = neighbors.end();

			double dC0 = covar(first_neigh->location(),first_neigh->location());

			// Only compute the upper triangle of the matrix
			for(Neighborhood::const_iterator row = first_neigh; row != last_neigh; row++ ) 
			{
				j=1;

				//for(Neighborhood::const_iterator col = row ; col != last_neigh ; col++)
				for(Neighborhood::const_iterator col = first_neigh ; col != last_neigh ; col++)
				{
					double dTemp = dC0 - covar( row->location(), col->location() );
					sample_to_sample_dis_MatrixA(i,j++) = dTemp;

					dDisSumTotal += dTemp;
					dDisSum1 += dTemp;
				}

				sample_to_est_dis_VecB(i) = dC0 - covar_rhs( row->location(), center );

				dvSample2SampleDis.push_back( dDisSum1/neighbors.size() ); 
				dDisSum1 = 0;

				i++;
			}
		}
		// Euclidean distance 
		else
		{
			int i=1;
			int j;

			Neighborhood::const_iterator first_neigh = neighbors.begin();
			Neighborhood::const_iterator last_neigh  = neighbors.end();

			// Only compute the upper triangle of the matrix
			for(Neighborhood::const_iterator row = first_neigh; row != last_neigh; row++ ) 
			{
				j=1;

				//for(Neighborhood::const_iterator col = row ; col != last_neigh ; col++)
				for(Neighborhood::const_iterator col = first_neigh ; col != last_neigh ; col++)
				{
					double dTemp = euclidean_norm( row->location() - col->location() );
					sample_to_sample_dis_MatrixA(i,j++) = dTemp;

					dDisSumTotal += dTemp;
					dDisSum1 += dTemp;
				}

				sample_to_est_dis_VecB(i)=euclidean_norm( row->location() - center );

				dvSample2SampleDis.push_back( dDisSum1/neighbors.size() ); 
				dDisSum1 = 0;

				i++;
			}
		}

		return;
	}

	double DualIDW::EstimationWithOnePoint()
	{
		int status = -1;
		Geovalue &begin = (*begin_);
		int center_node_id = begin.node_id();

		double variance = -1.0E10;

		std::vector<double>& weights = estimation_weights_;
		double& kriging_variance = variance;
		const Location& center = begin.location();
		const Neighborhood& neighbors = *(neighborhood_.raw_ptr());
		Covariance<Location>& covar= covar_;
		Covariance<Location>& covar_rhs = *rhs_covar_;
		KrigingConstraints& Kconstraints = *Kconstraints_;

		MatrixLib::Symmetric_matrix& A = A_;	// to build the kiring system: A b = x
		MatrixLib::Vector &b = b_;

		typedef matrix_lib_traits< GSTL_TNT_lib > MatrixLib;

		// dPowerSample2Est; dPowerSample2Sap
		double dPowerSample2Est = DualIDW::optPars_.dPowerSample2Est;
		double dPowerSample2Sap = DualIDW::optPars_.dPowerSample2Sap;

		// If the neighborhood is empty, there is no kriging to be done.
		if( neighbors.is_empty() ) 
		{
			gstl_warning( "Empty neighborhood. No kriging to be done. " );
			return 2;
		}

		// build the kriging system
		{
			int conditioning_data = 
				build_kriging_system(A, b,
				weights, 
				center, neighbors,
				covar, covar_rhs, Kconstraints);
		}

		// get the kriging weights
		if (1)
		{
			MatrixLib::Symmetric_matrix tempA;

			int conditioning_data = 
				build_kriging_system(tempA, b,
				ori_kriging_weights_, 
				center, neighbors,
				covar, covar_rhs, Kconstraints);
			status = 
				kriging_constraints_traits<
				KrigingConstraints,
				GSTL_TNT_lib
				>::const_kriging_solver(tempA, b, ori_kriging_weights_.begin());
		}	

		// to calculate the DIDW weights
		{
			{
				int conditioning_data = 
					build_kriging_system(A, b,
					weights, 
					center, neighbors,
					covar, covar_rhs, Kconstraints);
			}

			MatrixLib::Vector &sample_to_est_dis_VecB = sample_to_est_dis_VecB_;
			sample_to_est_dis_VecB.resize((int)neighbors.size());
			std::vector<double> &dvSample2SampleDis = dvSample2SampleDis_; 

			GetDistance4Estimation();

			if ( 
				(DualIDW::optPars_.nOptType_L[0] == 1 || DualIDW::optPars_.nOptType_L[1] == 1) 
				) 
			{
				Optimization_MinKV(dPowerSample2Est, dPowerSample2Sap);
			}

			CalculateDualIdwWeights(dPowerSample2Est, dPowerSample2Sap, false, false);
		}	

		// for kriging variance
		if(kriging_variance < -1.0E10 + 0.1)
		{
			kriging_variance = cal_kring_variance_by_weights(
				A, 
				b,
				weights,
				center,
				covar);
		}

		var_prop_PowerS2S_->set_value(dPowerSample2Sap, begin_->node_id());
		var_prop_PowerS2U_->set_value(dPowerSample2Est, begin_->node_id());

		return variance;
	}

	void DualIDW::TestOutputInterAccuracy(std::string& str_method_name)// TEST [2019-3-21 lizhanglin]
	{
		Geostat_grid*  simul_grid2 = simul_grid_; 

		std::vector<double> dEsts, dTrus;
		{
			typedef Geostat_grid::iterator iterator;
			iterator begin = simul_grid2->begin();
			iterator end = simul_grid2->end();

			GsTLGridProperty *pReferenceVal  = simul_grid2->property("V"/*harddata_property_name_*/);
			GsTLGridProperty *pEstimatedVal  = simul_grid2->property(property_name_);

			if(pReferenceVal && pEstimatedVal)
				for( ; begin != end; ++begin ) 
				{
					int id = begin->node_id();

					if(pEstimatedVal->is_informed(id) && pReferenceVal->is_informed(id))
					{
						dEsts.push_back(pEstimatedVal->get_value(begin->node_id()));
						dTrus.push_back(pReferenceVal->get_value(begin->node_id()));
					}
				}
		}

		if (dTrus.size() > 0)
		{

			double dMAE = GetMAE(dTrus, dEsts, 0, false);
			double dNRMSE = GetRMSE(dTrus, dEsts, 0, true);
			double dRMSE = GetRMSE(dTrus, dEsts, 0, false);
			double dCR = GetCC(dEsts, dTrus); 

			_cout_file<<str_method_name<<"\t interpolation accuracy (the reference property is V):\n";
			_cout_file<<"\t\t number of estimates:"<<dEsts.size()<<std::endl;
			_cout_file<<"\t\tMAE = "<<dMAE
				<<"\t\tNRMSE = "<<dNRMSE
				<<"\t\tRMSE = "<<dRMSE
				<<"\tCR = "<<dCR<<std::endl;

			_cout_file<<"\t\t"<<dMAE<<"\t\t"<<dNRMSE<<"\t\t"<<dRMSE<<"\t"<<dCR<<std::endl;
		}
	}




	// the main executive function for this plug-in 
	int DualIDW::execute( GsTL_project* proj )
	{
		this->_cout_file<<"\n++++++++++++++++++Optimization Begins!++++++++" ;//<< std::endl;

		char tmp[64]; 
		time_t t = time(0);
		strftime(tmp, sizeof(tmp), "current time: %Y/%m/%d, %H:%M:%S\n", localtime(&t));
		this->_cout_file << tmp<< std::endl;

		if (InitializedOptParsByFile(DualIDW::optPars_) == false)
		{
			_cout_file<<"Error: Cann't find DIDW pars!"<<std::endl;
			return false;
		}

		switch(DualIDW::optPars_.nDIDWOrKrg)
		{
		case 0:
			std::cout<<"Estimation"<<std::endl;
			EstimationDIDW(0);
			break;
		case 1:
			std::cout<<"Kriging"<<std::endl;
			ExecuteKriging(0);
			break;
		}

		//////////////////////[2019-3-18,16:25 lizhanglin]//////////////+
		//TODO: 输出最终图层中的计算结果
		if(DualIDW::optPars_.nDIDWOrKrg == 0)
		{
			TestOutputInterAccuracy(std::string("DIDW"));
			//ExecuteKriging(0);	
			//TestOutputInterAccuracy(std::string("Kriging"));
		}
		else
		{
			TestOutputInterAccuracy(std::string("Kriging"));
		}		

		//////////////////////[2019-3-18,16:25 lizhanglin]//////////////-

		this->_cout_file<<"------------------Optimization Ends!--------" ;

		time_t t2 = time(0);
		strftime(tmp, sizeof(tmp), "current time: %Y/%m/%d, %H:%M:%S,\t",localtime(&t2));
		this->_cout_file << tmp<< std::endl;			
		this->_cout_file << "processed time:"<<(difftime(t2,t)/60.0)<<"(Mins)\n"<< std::endl;


		//QMessageBox message(QMessageBox::Information, "Information","Optimization Successed!",QMessageBox::Ok,NULL);  
		//message.exec();


		return 0;
	} 


	// the classic kriging process, the cross validation here is for all samples. 
	int DualIDW::ExecuteKriging( GsTL_project* proj, bool b4CrossValid /*= false*/)
	{
		Geostat_grid* simul_grid = 0; 

		simul_grid = simul_grid_;

		// those flags will be used to signal if some of the nodes could not be
		// informed
		bool issue_singular_system_warning = false;
		bool issue_no_conditioning_data_warning = false;

		// Set up a progress notifier	
		int total_steps = simul_grid->size();
		int frequency = std::max( total_steps / 20, 1 );

		SmartPtr<Progress_notifier> progress_notifier = 
		  utils::create_notifier( "Running Kriging", 
			    total_steps, frequency );

		GsTLGridProperty* prop = 0;
		if (prop = simul_grid->property(property_name_))
		{
			simul_grid->select_property(property_name_);
		}
		else
		{
			// create the property
			appli_message("creating new property: " << property_name_ << "..." );
			prop = geostat_utils::add_property_to_grid( simul_grid, property_name_ );
			simul_grid->select_property( prop->name() );
		}

		GsTLGridProperty* var_prop = 0;
		if (var_prop = simul_grid->property(property_name_ + "_krig_var"))
		{

		}
		else
		{
			// create property for kriging variance
			std::string var_prop_name = prop->name() + "_krig_var";
			var_prop = geostat_utils::add_property_to_grid( simul_grid, var_prop_name );
		}

		typedef Geostat_grid::iterator iterator;
		iterator begin = simul_grid->begin();
		iterator end = simul_grid->end();

		/*
		Block_covariance<Location>* rhs_covar_blk = 0;
		if(do_block_kriging_)  
		rhs_covar_blk = static_cast<Block_covariance<Location>*>(rhs_covar_);
		*/

		double dReservedVal = 0;
		bool bReserved = false;

		for( ; begin != end; ++begin )   
		{
			if( !progress_notifier->notify() ) {
			  clean( property_name_ );
			  return 1;
			}

			int node_id = begin->node_id();

			bReserved = false;

			//_cout_file<<"Point to be Estimated: "<<begin->location()<<"\t"<<begin->property_value()<<std::endl;

			neighborhood_->includes_center(false); 

			neighborhood_->find_neighbors( *begin );

			if( neighborhood_->size() < min_neigh_ )  
			{
				begin->set_not_informed();
				var_prop->set_not_informed(node_id);
				continue;
			}

			if(!neighborhood_->is_valid()) 
			{
				begin->set_not_informed();
				var_prop->set_not_informed(node_id);
				continue;
			};

			/*
			//   if( neighborhood_->is_empty() ) {
			if( neighborhood_->size() < min_neigh_ ) {
			//if we don't have any conditioning data, skip the node
			issue_no_conditioning_data_warning = true;
			continue;
			}
			*/
			double variance;

			int status;  

			if(rhs_covar_blk_) {
				status  = kriging_weights_2( estimation_weights_, variance,
					begin->location(), *(neighborhood_.raw_ptr()),
					covar_,*rhs_covar_blk_, *Kconstraints_ );
			} 
			else {
				status = kriging_weights_2( estimation_weights_, variance,
					begin->location(), *(neighborhood_.raw_ptr()),
					covar_,*rhs_covar_, *Kconstraints_ );
			}

			if(status == 0) 
			{
				// the kriging system could be solved
				double estimate = (*combiner_)( estimation_weights_.begin(), 
					estimation_weights_.end(),
					*(neighborhood_.raw_ptr()) );
				begin->set_property_value( estimate );
				var_prop->set_value( variance, begin->node_id() );


				//_cout_file<<"Estimated Point Result: "
				// <<"EV:"<<"\t"<<estimate
				// <<"KV:"<<"\t"<<variance<<std::endl;
			}
			else 
			{
				// the kriging system could not be solved, issue a warning and skip the
				// node

				begin->set_not_informed();
				var_prop->set_not_informed(node_id);

				issue_singular_system_warning = true;

			}

		}
		/* This pop-up windows breaks script

		if( issue_singular_system_warning )
		GsTLcerr << "Kriging could not be performed at some locations because\n"
		<< "the kriging system was singular\n" 
		<< gstlIO::end; 
		if( issue_no_conditioning_data_warning )
		GsTLcerr << "Kriging could not be performed at some locations because\n"
		<< "the neighborhood of those locations was empty.\n"
		<< "Try increasing the size of the search ellipsoid.\n"
		<< gstlIO::end; 
		*/

		return 0;
	}

	// Calculate Dual IDW weights
	void DualIDW::CalculateDualIdwWeights(double dPowerSample2Est, double dPowerSample2Sap, bool bStandardized /*= false*/, bool bOutPut /*= false*/)
	{
		MatrixLib::Vector &sample_to_est_dis_VecB = sample_to_est_dis_VecB_;
		std::vector<double> &dvSample2SampleDis = dvSample2SampleDis_; 

		const Neighborhood& neighbors =  *(neighborhood_.raw_ptr());

		int nSamples = neighbors.size();

		std::vector<double>& weights = estimation_weights_;

		// Step 1/3 calculate the D-U weight
		std::vector<double> dWeightSample2Estimate;
		{
			weights.resize(neighbors.size());

			std::copy(sample_to_est_dis_VecB.begin(), sample_to_est_dis_VecB.begin() + weights.size(), weights.begin());

			std::vector<double>::iterator it = weights.begin();

			std::vector<double>::iterator it_MinDis = std::min_element(weights.begin(), weights.end());

			// for the extremely nearby data 
			if (*it_MinDis < 0.1E-6)
			{
				for (it = weights.begin(); it != weights.end(); it++)
				{
					*it = 0;
				}
				(*weights.begin()) = 1.0;
			}
			else
			{
				double dPower = dPowerSample2Est;

				for (it = weights.begin(); it != weights.end(); it++)
				{
					// large D-U distance is related to small weight and thus the negative exponent is applied 
					*it = pow(*it, -dPower);
				}

				double dSumDis = std::accumulate(weights.begin(), weights.end(), 0.0);

				for (it = weights.begin(); it != weights.end(); it++)
				{
					*it = (*it)/dSumDis;
				}

				if (bStandardized)
				{
					double dSumDis = std::accumulate(weights.begin(), weights.end(), 0.0);

					for (it = weights.begin(); it != weights.end(); it++)
					{
						*it = (*it)/dSumDis;
					}
				}
			}

			dWeightSample2Estimate.resize(weights.size()); 
			std::copy(weights.begin(), weights.begin() + weights.size(), dWeightSample2Estimate.begin());
		}

		// Step 2/3 calculate the D-D weight
		std::vector<double> dWeightSample2Sample;
		if( neighbors.size() > 1)
		{
			weights.resize(neighbors.size());

			std::copy(dvSample2SampleDis.begin(), dvSample2SampleDis.begin() + weights.size(), weights.begin());

			std::vector<double>::iterator it = weights.begin();

			std::vector<double>::iterator it_MinDis = std::min_element(weights.begin(), weights.end());

			//if (*it_MinDis < 0.1E-6)
			//{
			//	for (it = weights.begin(); it != weights.end(); it++)
			//	{
			//		*it = 0;
			//	}
			//	(*weights.begin()) = 1.0;
			//}
			//else
			{
				double dPower = dPowerSample2Sap;

				for (it = weights.begin(); it != weights.end(); it++)
				{
					// large D-D distance is related to the large weight, and thus the positive power is applied
					*it = pow(*it, dPower);
				}

				double dSumDis = std::accumulate(weights.begin(), weights.end(), 0.0);

				for (it = weights.begin(); it != weights.end(); it++)
				{
					*it = (*it)/dSumDis;
				}

				if (bStandardized)
				{
					double dSumDis = std::accumulate(weights.begin(), weights.end(), 0.0);

					for (it = weights.begin(); it != weights.end(); it++)
					{
						*it = (*it)/dSumDis;
					}
				}
			}

			dWeightSample2Sample.resize(weights.size()); 
			std::copy(weights.begin(), weights.begin() + weights.size(), dWeightSample2Sample.begin());
		}

		// Step 3/3 combine the D-D and D-U weights
		if( neighbors.size() > 1)
		{
			//std::inner_product(dWeightSample2Estimate.begin(), dWeightSample2Estimate.end(),
			//	dWeightSample2Sample.begin(), dWeightSample2Sample.end(), 0.0);

			for (int iSapInx = 0; iSapInx<dWeightSample2Sample.size(); iSapInx++)
			{
				weights[iSapInx] = dWeightSample2Estimate[iSapInx]*dWeightSample2Sample[iSapInx];						
			}

			double dSumDis = std::accumulate(weights.begin(), weights.end(), 0.0);

			for (int iSapInx = 0; iSapInx<dWeightSample2Sample.size(); iSapInx++)
			{
				weights[iSapInx] = weights[iSapInx]/dSumDis;
			}
		}

		//TEST CODE add by O on [2018-4-14]+
		//FOR: output the D-D& D-U distances and weights
		if(bOutPut)
		{
			this->_cout_file 
				<< "dPowerSample2Est = \t"<< dPowerSample2Est
				<< "dPowerSample2Sap = \t"<< dPowerSample2Sap
				<< std::endl;

			this->_cout_file 
				<< "X\t"<< "Y\t"<< "Z\t"
				<< "samp2sap_dis\t"
				<< "samp2Est_dis\t"

				<< "dWeightSample2Estimate\t"
				<< "dWeightSample2Sample\t"
				<< "weights\t"
				<< std::endl;


			Neighborhood::const_iterator first_neigh = neighbors.begin();
			Neighborhood::const_iterator row = first_neigh; 

			for (int i = 0; i<dvSample2SampleDis.size(); i++)
			{
				if(dWeightSample2Sample.size() > 0)
					this->_cout_file
					<< row->location()[0]<<"\t"
					<< row->location()[1]<<"\t"
					<< row->location()[2]<<"\t"

					<< dvSample2SampleDis[i]<<"\t"
					<< sample_to_est_dis_VecB(i + 1)<<"\t"
					<< dWeightSample2Estimate[i]<<"\t"
					<< dWeightSample2Sample[i]<<"\t"
					<< weights[i]<<"\t"
					<< std::endl;		
				else
					this->_cout_file
					<< row->location()[0]<<"\t"
					<< row->location()[1]<<"\t"
					<< row->location()[2]<<"\t"

					<< dvSample2SampleDis[i]<<"\t"
					<< sample_to_est_dis_VecB(i + 1)<<"\t"
					<< dWeightSample2Estimate[i]<<"\t"
					<< "0.0"<<"\t"
					<< weights[i]<<"\t"
					<< std::endl;	


				row++;
			}
		}
		//TEST CODE add by O on [2018-4-14]-

		return;
	}


	// Calculate Estimation Variance
	double DualIDW::CalculateEstimationVar(double dTempPowerSample2Est, double dTempPowerSample2Sap, bool bStandardized /*= false*/)
	{
		MatrixLib::Vector &sample_to_est_dis_VecB = sample_to_est_dis_VecB_;
		std::vector<double> &dvSample2SampleDis = dvSample2SampleDis_; 

		const Neighborhood& neighbors =  *(neighborhood_.raw_ptr());

		int nSamples = neighbors.size();

		std::vector<double>& weights = estimation_weights_;

		// Step 1/3 calculate the D-U weight 
		std::vector<double> dWeightSample2Estimate;
		{
			weights.resize(neighbors.size());

			std::copy(sample_to_est_dis_VecB.begin(), sample_to_est_dis_VecB.begin() + weights.size(), weights.begin());

			std::vector<double>::iterator it = weights.begin();

			std::vector<double>::iterator it_MinDis = std::min_element(weights.begin(), weights.end());

			// for the extremely nearby data 
			if (*it_MinDis < 0.1E-6)
			{
				for (it = weights.begin(); it != weights.end(); it++)
				{
					*it = 0;
				}
				(*weights.begin()) = 1.0;
			}
			else
			{
				double dPower = dTempPowerSample2Est;

				for (it = weights.begin(); it != weights.end(); it++)
				{
					// large D-U distance is related to small weight and thus the negative exponent is applied 
					*it = pow(*it, -dPower);
				}

				if (bStandardized)
				{
					double dSumDis = std::accumulate(weights.begin(), weights.end(), 0.0);

					for (it = weights.begin(); it != weights.end(); it++)
					{
						*it = (*it)/dSumDis;
					}
				}
			}

			dWeightSample2Estimate.resize(weights.size()); 
			std::copy(weights.begin(), weights.begin() + weights.size(), dWeightSample2Estimate.begin());
		}

		// Step 2/3 calculate the D-D weight
		std::vector<double> dWeightSample2Sample;
		if( neighbors.size() > 1)
		{
			weights.resize(neighbors.size());

			std::copy(dvSample2SampleDis.begin(), dvSample2SampleDis.begin() + weights.size(), weights.begin());

			std::vector<double>::iterator it = weights.begin();

			std::vector<double>::iterator it_MinDis = std::min_element(weights.begin(), weights.end());

			//if (*it_MinDis < 0.1E-6)
			//{
			//	for (it = weights.begin(); it != weights.end(); it++)
			//	{
			//		*it = 0;
			//	}
			//	(*weights.begin()) = 1.0;
			//}
			//else
			{
				double dPower = dTempPowerSample2Sap;

				for (it = weights.begin(); it != weights.end(); it++)
				{
					// large D-D distance is related to the large weight, and thus the positive power is applied
					*it = pow(*it, dPower);
				}


				if (bStandardized)
				{
					double dSumDis = std::accumulate(weights.begin(), weights.end(), 0.0);

					for (it = weights.begin(); it != weights.end(); it++)
					{
						*it = (*it)/dSumDis;
					}
				}
			}

			dWeightSample2Sample.resize(weights.size()); 
			std::copy(weights.begin(), weights.begin() + weights.size(), dWeightSample2Sample.begin());
		}

		// Step 3/3 combine the D-U and D-D weights 
		if( neighbors.size() > 1)
		{
			//std::inner_product(dWeightSample2Estimate.begin(), dWeightSample2Estimate.end(),
			//	dWeightSample2Sample.begin(), dWeightSample2Sample.end(), 0.0);

			for (int iSapInx = 0; iSapInx<dWeightSample2Sample.size(); iSapInx++)
			{
				weights[iSapInx] = dWeightSample2Estimate[iSapInx]*dWeightSample2Sample[iSapInx];						
			}

			double dSumDis = std::accumulate(weights.begin(), weights.end(), 0.0);

			for (int iSapInx = 0; iSapInx<dWeightSample2Sample.size(); iSapInx++)
			{
				weights[iSapInx] = weights[iSapInx]/dSumDis;
			}
		}

		double variance;
		// for kriging variance
		{
			double& kriging_variance = variance;
			const Location& center = (begin_)->location();
			const Neighborhood& neighbors =  *(neighborhood_.raw_ptr());
			Covariance<Location>& covar= covar_;
			Covariance<Location>& covar_rhs = *(rhs_covar_blk_);
			KrigingConstraints& Kconstraints = *(Kconstraints_ );

			MatrixLib::Symmetric_matrix& A = A_;	// Kriging system Ax = b
			MatrixLib::Vector &b = b_;


			kriging_variance = cal_kring_variance_by_weights0(
				A, 
				b,
				weights,
				center,
				covar);
		}

		return variance;
	}



	// read the parameter file and return Parameters_handler
	Parameters_handler* DualIDW::InputPars(const char* filename)
	{
		std::ifstream infile( filename);
		if( !infile ) {
			return 0;
		}

		std::ostringstream file_content;
		char ch;
		while( file_content && infile.get( ch ) )
			file_content.put( ch );

		Parameters_handler *par_handler = (Parameters_handler *)Parameters_handler_xml::create_new_interface(file_content.str());

		appli_assert( par_handler );

		if( !par_handler->is_ready() ) {
			return 0;
		}
		std::string algo_name = par_handler->value( "algorithm.name" );
		if( algo_name.empty() ) {
			return 0;
		}

		return par_handler;
	}

	// Initialized Opt parameters by _ni_KrgPars
	bool DualIDW::InitializedOptParsByFile(OptParameters &optPars)
	{
		// get the value of _ni_KrgPars 
		{
			std::string str = _sModelingPath;

			std::ifstream infile(str.c_str());
			if( infile ) 
			{
				_ni_KrgPars = InputPars(str.c_str());
			}
			else
			{
				_cout_file<<"Error when open:"<<str<<std::endl;
				return false;
			}

			if (_ni_KrgPars.raw_ptr() == 0)
			{
				_cout_file<<"Error when read:"<<str<<std::endl;
				return false;
			}
		}

		Parameters_handler* par_handler = (Parameters_handler*)(_ni_KrgPars.raw_ptr());

		// initialization
		{
			optPars.nDIDWOrKrg = String_Op::to_number<int>( par_handler->value( "nDIDWOrKrg.value" ));

			// for DURAL IDW 
			optPars.dPowerSample2Est = String_Op::to_number<double>( par_handler->value( "dPowerSample2Est.value" ));
			optPars.dPowerSample2Sap = String_Op::to_number<double>( par_handler->value( "dPowerSample2Sap.value" ));
			optPars.bUseAnisotropicDistance = String_Op::to_number<bool>( par_handler->value( "bUseAnisotropicDistance.value" ));

			optPars.dPowerS2U_ValueRange = String_Op::to_numbers<double>( par_handler->value( "dPowerS2U_ValueRange.value" ));
			if ( optPars.dPowerS2U_ValueRange.size() != 3)
			{
				optPars.dPowerS2U_ValueRange.clear();
				optPars.dPowerS2U_ValueRange.push_back(50);
				optPars.dPowerS2U_ValueRange.push_back(0);
				optPars.dPowerS2U_ValueRange.push_back(20);
			}

			optPars.dPowerS2S_ValueRange = String_Op::to_numbers<double>( par_handler->value( "dPowerS2S_ValueRange.value" ));
			if ( optPars.dPowerS2S_ValueRange.size() != 3)
			{
				optPars.dPowerS2S_ValueRange.clear();
				optPars.dPowerS2S_ValueRange.push_back(50);
				optPars.dPowerS2S_ValueRange.push_back(0);
				optPars.dPowerS2S_ValueRange.push_back(20);
			}

			optPars.nOptType_G = String_Op::to_numbers<int>( par_handler->value( "nOptType_G.value" ));
			if ( optPars.nOptType_G.size() != 2)
			{
				optPars.nOptType_G.clear();
				optPars.nOptType_G.push_back(0);
				optPars.nOptType_G.push_back(0);
			}

			optPars.nOptType_L = String_Op::to_numbers<int>( par_handler->value( "nOptType_L.value" ));
			if ( optPars.nOptType_L.size() != 2)
			{
				optPars.nOptType_L.clear();
				optPars.nOptType_L.push_back(0);
				optPars.nOptType_L.push_back(0);
			}

			optPars.nvGlobalOptGoalType = String_Op::to_numbers<int>( par_handler->value( "nvGlobalOptGoalType.value" ));
			if ( optPars.nvGlobalOptGoalType.size() == 0)
			{
				optPars.nvGlobalOptGoalType.push_back(1);
			}

			optPars.bOutputTestDetails = String_Op::to_number<bool>( par_handler->value( "bOutputTestDetails.value" ));
			optPars.bUsePower_D2U_AS_D2D = String_Op::to_number<bool>( par_handler->value( "bUsePower_D2U_AS_D2D.value" ));

			if (optPars.bUsePower_D2U_AS_D2D)
			{
				DualIDW::optPars_.nOptType_G[1] = 0; // now 1 is invalid for this parameter
				DualIDW::optPars_.nOptType_L[1] = 0; // now 1 is invalid for this parameter
			}
		}

		return true;
	}

	int DualIDW::Whole_OptimizationByEnumeration(Geostat_grid* simul_grid, GsTLGridProperty* var_prop)
	{
		std::map<double, interpolation_accuracy_measurement_inx> temp_interpolation_result_index;
		std::map<double, interpolation_accuracy_measurement_inx> interpolation_result_index;
		interpolation_accuracy_measurement_inx temp_iam_inx;

		// get the permissible powers 
		std::vector<double> &dPowersS2U = dPowersS2U_, &dPowersS2S = dPowersS2S_;	

		SmartPtr<Progress_notifier> progress_notifier = 0;

		if (DualIDW::optPars_.nOptType_G[0]==1 && DualIDW::optPars_.nOptType_G[1]==1 )
		{
			// Set up a progress notifier	
			int total_steps = dPowersS2U.size()*dPowersS2S.size();
			int frequency = std::max( total_steps / 20, 1 );
			progress_notifier = utils::create_notifier( "Running Whole Optimization", total_steps, frequency );
		}
		else if (DualIDW::optPars_.nOptType_G[0]==1)
		{
			// Set up a progress notifier	
			int total_steps = dPowersS2U.size();
			int frequency = std::max( total_steps / 20, 1 );
			progress_notifier = utils::create_notifier( "Running Semi-Whole Optimization", total_steps, frequency );
		}
		else if (DualIDW::optPars_.nOptType_G[1]==1)
		{
			// Set up a progress notifier	
			int total_steps = dPowersS2S.size();
			int frequency = std::max( total_steps / 20, 1 );
			progress_notifier = utils::create_notifier( "Running Semi-Whole Optimization", total_steps, frequency );
		}

		// output the field names
		//if (DualIDW::optPars_.bOutputTestDetails)
		{

			this->_cout_file 
				<< "bAnisotropicDistance\t d-u_Power\t d-d_Power\t ErrorMean\t ErrorStdVar\t MinError\t MaxError\t MTE\t MAE\t RMSE\t CC\t CV\t MeanWeightsCR\t MeanKV\t"
				<< std::endl;

		}

		// use dPowersS2U 
		std::vector<double>::iterator it_power_s2u = dPowersS2U.begin();	

		do
		{
			DualIDW::optPars_.dPowerSample2Est = (DualIDW::optPars_.nOptType_G[0] == 1)? (*it_power_s2u):(DualIDW::optPars_.dPowerSample2Est);

			std::vector<double>::iterator it_power_s2s = dPowersS2S.begin();

			do
			{
				if( progress_notifier ) {
					if( !progress_notifier->notify() ) {
						clean( property_name_ );
						return 0.0;
					}
				}

				DualIDW::optPars_.dPowerSample2Sap = (DualIDW::optPars_.nOptType_G[1] == 1)?(*it_power_s2s):(DualIDW::optPars_.dPowerSample2Sap);

				if (DualIDW::optPars_.bUsePower_D2U_AS_D2D)  
				{
					DualIDW::optPars_.dPowerSample2Sap = DualIDW::optPars_.dPowerSample2Est;
					DualIDW::optPars_.nOptType_G[1] = 0; // now 1 is invalid for this parameter
				}


				//// use dPowersS2S as the first loop 
				//std::vector<double>::iterator it_power_s2s = dPowersS2S.begin();

				//do
				//{
				//	DualIDW::optPars_.dPowerSample2Sap = (DualIDW::optPars_.nOptType_G[1] == 1)?(*it_power_s2s):(DualIDW::optPars_.dPowerSample2Sap);
				//	std::vector<double>::iterator it_power_s2u = dPowersS2U.begin();	

				//	do
				//	{
				//		if( progress_notifier ) {
				//			if( !progress_notifier->notify() ) {
				//				clean( property_name_ );
				//				return 0.0;
				//			}
				//		}

				//		DualIDW::optPars_.dPowerSample2Est = (DualIDW::optPars_.nOptType_G[0] == 1)? (*it_power_s2u):(DualIDW::optPars_.dPowerSample2Est);

				//		if (DualIDW::optPars_.bUsePower_D2U_AS_D2D)  
				//		{
				//			DualIDW::optPars_.dPowerSample2Est = DualIDW::optPars_.dPowerSample2Sap;
				//			DualIDW::optPars_.nOptType_G[1] = 0; // now 1 is invalid for this parameter
				//		}


				double dMainObj = EstimationWithAllPoints(simul_grid, var_prop, &temp_iam_inx, DualIDW::optPars_.bOutputTestDetails);

				interpolation_result_index.insert(std::pair<double, interpolation_accuracy_measurement_inx> (
					dMainObj, temp_iam_inx 
					));

				temp_interpolation_result_index.insert(std::pair<double, interpolation_accuracy_measurement_inx> (
					dMainObj, temp_iam_inx 
					));

				if (DualIDW::optPars_.nOptType_G[1] == 0)	break;
			}
			while(++it_power_s2s != dPowersS2S.end());  // use dPowersS2U as the first loop 
			//while (++it_power_s2u != dPowersS2U.end());   // use dPowersS2S as the first loop 

			// if the details are not needed
			// output the best results corresponding to the curent D-D (or D-U, decided by the loop struct) and D-U exponent
			if(!DualIDW::optPars_.bOutputTestDetails)
			{
				std::map<double, interpolation_accuracy_measurement_inx>::iterator map_it = temp_interpolation_result_index.begin();

				interpolation_accuracy_measurement_inx *p_IAM_inx = &(map_it->second);

				this->_cout_file 
					//<<"IsAD\t"<<"PowerS2U\t"<<"PowerS2S\t"<<"MeanError\t"<<"Error_Variance\t"<<"MAE\t"<<"RMSE\t"<<"CR\n"
					<< DualIDW::optPars_.bUseAnisotropicDistance<<"\t"
					<< p_IAM_inx->dPowerSample2Est<<"\t" 
					<< p_IAM_inx->dPowerSample2Sap<<"\t"

					<< p_IAM_inx->iam.dErrorMean<<"\t"
					<< p_IAM_inx->iam.dErrorVariance<<"\t"
					<< p_IAM_inx->iam.dMinError<<"\t"
					<< p_IAM_inx->iam.dMaxError<<"\t"

					<< p_IAM_inx->iam.dMTE<<"\t"
					<< p_IAM_inx->iam.dMAE<<"\t"
					<< p_IAM_inx->iam.dRMSE<<"\t"

					<< p_IAM_inx->iam.dCC<<"\t"
					<< p_IAM_inx->iam.dCV<<"\t"

					<< p_IAM_inx->iam.dMeanWeightsCR<<"\t"
					<< p_IAM_inx->iam.dMeanKV<<"\t"

					<< std::endl;
			}

			temp_interpolation_result_index.clear(); // clear the current result 

			if (DualIDW::optPars_.nOptType_G[0] == 0)	break;
		}
		while (++it_power_s2u != dPowersS2U.end()); // use dPowersS2U as the first loop 
		//while( ++it_power_s2s != dPowersS2S.end() );  // use dPowersS2S as the first loop


		int iii = 0;
		for (std::map<double, interpolation_accuracy_measurement_inx>::iterator map_it = interpolation_result_index.begin();
			map_it != interpolation_result_index.end();
			map_it++)
		{
			if (iii++ == 1) // output the first result 
			{
				break;
			}

			temp_iam_inx = (map_it->second);

			this->_cout_file 
				<<"\nNumber of estimates:\t"<<temp_iam_inx.iam.nNumberOfEstimates
				<<"\tMeanS2S_DisVariance:\t"<<temp_iam_inx.iam.dMeanS2S_DisVariance
				<<"\tMeanS2S_DisMean:\t"<<temp_iam_inx.iam.dMeanS2S_DisMean
				<<"\nscore:\t"<<map_it->first
				<<"\tdPowerS-U:\t"<< temp_iam_inx.dPowerSample2Est<<"\tdPowerS-S:\t"<< temp_iam_inx.dPowerSample2Sap
				//<< "dErrorMean\t"<<dErrorMean
				<< "\nError_Variance:\t"<< temp_iam_inx.iam.dErrorVariance << "\tGetMTE:"<<temp_iam_inx.iam.dMTE
				<< "\tError_Min:\t"<< temp_iam_inx.iam.dMinError
				<< "\tError_Max:\t"<< temp_iam_inx.iam.dMaxError

				<< "\nGetMAE:\t"<<temp_iam_inx.iam.dMAE 
				<< "\tGetRMSE:\t"<<temp_iam_inx.iam.dRMSE
				<< "\tGetCC\t:"<<temp_iam_inx.iam.dCC
				<< "\tGetCV\t:"<<temp_iam_inx.iam.dCV
				<< std::endl;
		}

		// if global OPT is applied, output the result and re-estimate based on the opted result.
		if(interpolation_result_index.size() > 2)
		{
			// only the best POWER is adopted 
			temp_iam_inx = interpolation_result_index.begin()->second;

			this->_cout_file <<"best powers:" << std::endl;

			this->_cout_file 
				<<"dPowerSample2Est:"<< temp_iam_inx.dPowerSample2Est
				<<"dPowerSample2Sap:"<< temp_iam_inx.dPowerSample2Sap
				<< std::endl;

			DualIDW::optPars_.dPowerSample2Est = temp_iam_inx.dPowerSample2Est;
			DualIDW::optPars_.dPowerSample2Sap = temp_iam_inx.dPowerSample2Sap;

			// execute the estimation with the optimized results 
			EstimationWithAllPoints(simul_grid, var_prop, &(interpolation_accuracy_measurement_inx()), true);
		}

		return 0.0;
	}

	int DualIDW::Optimization_MinKV(double& dPowerSample2Est, double& dPowerSample2Sap)
	{
		{
			// get the permissible powers 
			std::vector<double> &dPowersS2U = dPowersS2U_, &dPowersS2S = dPowersS2S_;	

			double dKV = 1.0E99;

			int i = 1;
			do
			{
				double dTempPowerSample2Est = (DualIDW::optPars_.nOptType_L[0] == 1)?(dPowersS2U[i-1]):(DualIDW::optPars_.dPowerSample2Est);

				int j = 1;
				do
				{
					double dTempPowerSample2Sap = (DualIDW::optPars_.nOptType_L[1] == 1)?(dPowersS2S[j-1]):(DualIDW::optPars_.dPowerSample2Sap);

					if (DualIDW::optPars_.bUsePower_D2U_AS_D2D)  
					{
						dTempPowerSample2Sap = dTempPowerSample2Est;
						DualIDW::optPars_.nOptType_L[1] = 0; // now 1 is invalid for this parameter
					}

					double kriging_variance = CalculateEstimationVar(dTempPowerSample2Est, dTempPowerSample2Sap);								

					if (kriging_variance < dKV)
					{
						dKV = kriging_variance;

						dPowerSample2Est = dTempPowerSample2Est;
						dPowerSample2Sap = dTempPowerSample2Sap;
					}

					if(DualIDW::optPars_.nOptType_L[1] == 0) break;
					j++;
				}
				while(j<= dPowersS2S.size());

				if(DualIDW::optPars_.nOptType_L[0] == 0) break;
				i++;
			}
			while(i<= dPowersS2U.size());
		}

		return 0;
	}


	//GEOSTAT_PLUGIN(DIDW);
	extern "C" __declspec(dllexport) int plugin_init() 
	{ 
		const std::string geostatAlgo_manager = "/GeostatAlgo";

		SmartPtr<Named_interface> ni = Root::instance()->interface( geostatAlgo_manager ); 

		Manager* dir = dynamic_cast<Manager*>( ni.raw_ptr() ); 

		if( !dir ) 
		{ 
			//GsTLlog << "Directory " << geostatAlgo_manager << " does not exist \n"; 
			return 1; 
		} 

		DualIDW toto2; 
		dir->factory( toto2.name(), DualIDW::create_new_interface ); 

		return 0; 
	}