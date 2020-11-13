#ifndef DualIDW_H
#define DualIDW_H


// for GSTL
#include <GsTL/geometry/covariance.h> 
#include <GsTL/geometry/Block_covariance.h>
#include <GsTL/kriging/kriging_weights.h>
#include <GsTL/geometry/Block_covariance.h>
#include <GsTL/utils/debug_tools.h>

#include <GsTLAppli/geostat/common.h>
#include <GsTLAppli/geostat/geostat_algo.h> 
#include <GsTLAppli/geostat/utilities.h> 
#include <GsTLAppli/geostat/kriging.h>
#include <GsTLAppli/geostat/parameters_handler.h>
#include <GsTLAppli/geostat/parameters_handler_impl.h>

#include <GsTLAppli/grid/grid_model/property_copier.h>
#include <GsTLAppli/grid/grid_model/grid_region_temp_selector.h>
#include <GsTLAppli/grid/grid_model/geostat_grid.h>
#include <GsTLAppli/grid/grid_model/point_set.h>
#include <GsTLAppli/grid/grid_model/combined_neighborhood.h>
#include <GsTLAppli/grid/grid_model/gval_iterator.h> 
#include <GsTLAppli/grid/grid_model/grid_initializer.h>
#include <GsTLAppli/grid/grid_model/cartesian_grid.h>
#include <GsTLAppli/grid/grid_model/point_set.h>

#include <GsTLAppli/utils/gstl_plugins.h>
#include <GsTLAppli/utils/gstl_messages.h>
#include <GsTLAppli/utils/string_manipulation.h>
#include <GsTLAppli/utils/error_messages_handler.h>

#include <GsTLAppli/appli/manager_repository.h>
#include <GsTLAppli/appli/utilities.h>

#include <GsTLAppli\gui\variogram2\variogram_modeler_gui.h>

#include <string> 
#include <map>

#include "common.h"
#include <QFileDialog>
#include <QMessageBox>
#include <io.h>

#include <time.h>


typedef Geostat_grid::location_type Location; 

class Neighborhood; 
class RGrid;
class DualIDW;

// parameters used in optimization
struct OptParameters 
{
	int nDIDWOrKrg;							// IDW or Kriging, 0- IDW; 1- Kriging

	// for DURAL IDW 
	double dPowerSample2Est;				// default value of sample-to-estimate Power
	double dPowerSample2Sap;				// default value of sample-to-sample Power
	bool bUseAnisotropicDistance;			// Is Anisotropic Distance used? 0- NO; 1- YES;

	std::vector<int> 	nOptType_G;			// Global OPT [2], 1- yes; 0- no; for [D2E_power; D2D_power;]
	std::vector<int> 	nOptType_L;			// Local  OPT [2], 1- yes; 0- no; for [D2E_power; D2D_power;]

	// 3 pars: the number of possible values; Min Weight; Max Weight
	std::vector<double> dPowerS2U_ValueRange; // Sample-to-estimation powers
	std::vector<double> dPowerS2S_ValueRange; // Sample-to-sample powers 
	

	//	Optimization Goal (only for global OPT)
	// 0 ErrorMean;
	// 1 ErrorVariance;
	// 2 MinError;
	// 3 MaxError;

	// 4 MTE;
	// 5 MAE;
	// 6 RMSE;

	// 7 (-CC);// Correlation Coefficient 
	// 8 (-CV);// Covariance

	// 9 (-MeanWeightsCR);
	// 10 MeanKV;
	std::vector<int>  nvGlobalOptGoalType;

	bool bOutputTestDetails;				// Need output test details? ? 0- NO; 1- YES;
	bool bUsePower_D2U_AS_D2D;				// Use the D2E power value as D2D power? 0- NO; 1- YES;
};

struct interpolation_accuracy_measurement
{
	int nNumberOfEstimates;
	double dMeanS2S_DisVariance;
	double dMeanS2S_DisMean;

	double dErrorMean;
	double dErrorVariance;
	double dMinError;
	double dMaxError;

	double dMTE;
	double dMAE;
	double dRMSE;

	double dCC;// Correlation Coefficient 
	double dCV;// Covariance

	double dMeanWeightsCR;
	double dMeanKV;
};

struct  interpolation_accuracy_measurement_inx
{
public:
	interpolation_accuracy_measurement iam;

	double dPowerSample2Est;			// default value of sample-to-estimate Power
	double dPowerSample2Sap;			// default value of sample-to-sample Power
};

// Main plug-in class 
class DIDW_V1_API DualIDW: public Geostat_algo 
{
public:
	// calculate the kriging error variance and the variance of the estimete
	template<class SymmetricMatrix,class MatVector,class location,class Covariance,class Vector>
	static double cal_kring_variance_by_weights0(
		SymmetricMatrix& A, 
		MatVector& b,
		Vector& weights,
		const location& center,
		Covariance& covar,
		double* pdVarianceOfEstimate =0
		);

public: 
	DualIDW(); 
	~DualIDW(); 

	virtual bool initialize( const Parameters_handler* parameters, Error_messages_handler* errors ); 
	virtual int execute( GsTL_project* proj=0 ); 
	virtual std::string name() const { return "DualIDW"; } 

public: 
	static Named_interface* create_new_interface( std::string& ); 

protected:
	void clean( const std::string& prop ); 

public:
	std::vector<double> ini_kriging_weights1_;				// ini_kriging_weights1_ 
	std::vector<std::vector<double>> ini_kriging_weightss_;	// ini_kriging_weights1_ 

	double* _dReservedVal;			

	static int nInitialized_;

	Geovalue *begin_;

	typedef matrix_lib_traits< GSTL_TNT_lib > MatrixLib; 
	MatrixLib::Symmetric_matrix A_;	// A, b matrices of the kriging system
	MatrixLib::Vector b_;

public:
	std::vector<double> dPowersS2U_, dPowersS2S_;	// get all of the permissible value ranges for the two powers to be opted
	GsTLGridProperty *var_prop_PowerS2U_, *var_prop_PowerS2S_;
	GsTLGridProperty *var_prop_TrueError_, *var_prop_AbsError_;
	GsTLGridProperty *var_prop_S2sDis_,*mean_prop_S2sDis_, *var_prop_SampNum_;

	MatrixLib::Vector sample_to_est_dis_VecB_;
	std::vector<double> dvSample2SampleDis_; 

	double	dPowerSample2Est_OptResult, dPowerSample2Sap_OptResult; // the optimized power results

public:
	// Calculate Estimation Variance
	double CalculateEstimationVar(double dTempPowerSample2Est, double dTempPowerSample2Sap, bool bStandardized = false);

	// get the distances needed for the estimation
	void GetDistance4Estimation();

	void CalculateDualIdwWeights(double dPowerSample2Est, double dPowerSample2Sap, bool bStandardized = false, bool bOutPut = false);

public:
	double EstimationWithOnePoint();

	// get all of the possible value ranges for the two powers to be opted
	double PowerOptRanges(std::vector<double> &dPowersS2U, std::vector<double> &dPowersS2S); 
	
	double EstimationWithAllPoints(Geostat_grid* simul_grid, GsTLGridProperty* var_prop,
		interpolation_accuracy_measurement_inx *p_IAM_inx , bool bOutPut = false);

	int EstimationDIDW( GsTL_project* proj);

	int Whole_OptimizationByEnumeration(Geostat_grid* simul_grid, GsTLGridProperty* var_prop);

	int Optimization_MinKV(double& dPowerSample2Est, double& dPowerSample2Sap);

public:
	int ExecuteKriging( GsTL_project* proj=0 , bool b4CrossValid = false);
	Parameters_handler* InputPars(const char* filename);
	bool InitializedOptParsByFile(OptParameters &optPars);

public:
	void TestOutputInterAccuracy(std::string& str_method_name);// TEST [2019-3-21 lizhanglin]

public:
	static OptParameters optPars_;

#ifdef OUT_PUT_TO_SCREEN
	static std::ostream& _cout_file;
#else
	static std::ofstream& _cout_file;
#endif

public:
	SmartPtr<Named_interface>  _ni_KrgPars;

	SmartPtr<Named_interface>  _basic_KrgPars;

	std::string _sModelingPath;

protected: 
	typedef Geostat_grid::location_type Location; 
	typedef std::vector<double>::const_iterator weight_iterator; 
	typedef Kriging_combiner< weight_iterator, Neighborhood > KrigingCombiner; 
	typedef Kriging_constraints< Neighborhood, Location > KrigingConstraints;  

public:   
	Geostat_grid* simul_grid_; 
	std::string property_name_; 
	Geostat_grid* harddata_grid_; 
	GsTLGridProperty* hdata_prop_;
	GsTLGridProperty* blk_hdata_prop_;
	std::string harddata_property_name_;

	geostat_utils::Kriging_type ktype_;

protected: 
	SmartPtr<Neighborhood> neighborhood_; 

	Block_covariance<Location>* rhs_covar_blk_;

	KrigingCombiner* combiner_; 
	KrigingConstraints* Kconstraints_; 

	std::vector<double> estimation_weights_;
	std::vector<double> ori_kriging_weights_;

	int min_neigh_;
	GsTLVector<int> nblock_pts_;

	Temporary_gridRegion_Selector gridTempRegionSelector_;
	Temporary_gridRegion_Selector hdgridTempRegionSelector_;

	bool do_block_kriging_;
public:
	Covariance<Location> covar_;		// Covariance between samples 
	Covariance<Location>*  rhs_covar_;	// Covariance between estimated point and samples  
};



extern "C" DIDW_V1_API int plugin_init(); 

#endif