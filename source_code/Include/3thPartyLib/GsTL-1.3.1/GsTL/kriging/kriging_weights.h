/* GsTL: the Geostatistics Template Library
 * 
 * Author: Nicolas Remy
 * Copyright (c) 2000 The Board of Trustees of the Leland Stanford Junior University
 * 
 * All rights reserved.
 * 
 * Redistribution and use in source and binary forms, with or without modification, 
 * are permitted provided that the following conditions are met:
 * 
 *   1.Redistributions of source code must retain the above copyright notice, this 
 *     list of conditions and the following disclaimer. 
 *   2.Redistributions in binary form must reproduce the above copyright notice, this 
 *     list of conditions and the following disclaimer in the documentation and/or other
 *     materials provided with the distribution. 
 *   3.The name of the author may not be used to endorse or promote products derived 
 *     from this software without specific prior written permission. 
 * 
 * THIS SOFTWARE IS PROVIDED BY THE AUTHOR ``AS IS'' AND ANY EXPRESS OR IMPLIED 
 * WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY 
 * AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE AUTHOR BE 
 * LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
 * DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
 * LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
 * THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING 
 * NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF
 * ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 */

 

#ifndef __GSTL_KRIGING_WEIGHTS_H__
#define __GSTL_KRIGING_WEIGHTS_H__
#ifdef __GNUC__
#pragma interface
#endif


/** This class defines a generic kriging Weighting System.
 * It needs a linear algebra library for solving linear systems.
 * The library is up to the user: MatrixLibrary is a template
 * argument (see documentation for requirements on the matrix library).
 * The default matrix library is a slightly modified version of 
 * TNT ( http://math.nist.gov/tnt/ )
 */

#include <GsTL/utils/gstl_error_messages.h>
#include <GsTL/utils/debug_tools.h>
#include <GsTL/matrix_library/matrix_lib_traits.h>
#include <GsTL/matrix_library/gstl_tnt_lib.h>
#include <GsTL/kriging/kriging_constraints_traits.h>

#include <GsTL/kriging/helper_functions.h>


#include <iostream>


/** Compute the kriging weights.
 * @return 0 if no problem was encountered, 1 if the kriging system could not be
 * solved, and 2 if no neighbor were provided (empty neighborhood).
 * 
 * @param weights is the vector where the kriging weights will be output
 * @param kriging_variance is where the kriging variance will be output
 * @param center is the location being estimated
 * @param neighbors points to the neighborhood of conditioning data used to 
 * estimate center. If there are no neighbors, no kriging is performed and  
 * the return code is 2.
 * @param covar is the covariance function
 * @param Kconstraints are the kriging constraints (eg SK, OK, KT, ...).
 */

template<
         class MatrixLibrary,
         class Location,
         class Neighborhood,
         class Covariance,
         class KrigingConstraints,
         class Vector
        >
int kriging_weights(
		    Vector& weights,
		    double& kriging_variance,
		    const Location& center,
		    const Neighborhood& neighbors,
		    Covariance& covar,
		    KrigingConstraints& Kconstraints
		    ) {

  typedef matrix_lib_traits< MatrixLibrary > MatrixLib;
  
  // If the neighborhood is empty, there is no kriging to be done.
  if( neighbors.is_empty() ) {
    gstl_warning( "Empty neighborhood. No kriging to be done. " );
    return 2;
  }
  
  typename MatrixLib::Symmetric_matrix A;
  typename MatrixLib::Vector b;
  
  int conditioning_data = 
    build_kriging_system(A,b,
		         weights, 
		         center, neighbors,
		         covar, Kconstraints);

  // solve the system
  int status = 
         kriging_constraints_traits<
                                    KrigingConstraints,
                                    MatrixLibrary
                                   >::const_kriging_solver(A, b, weights.begin());

  // Compute the kriging variance
  if(status == 0) {
    double C0=covar(center,center);
    kriging_variance = 
      compute_kriging_variance(weights.begin(), 
                               weights.begin()+conditioning_data, weights.end(),
                               b, Kconstraints, center, C0);
  }
  else
    kriging_variance = -99;

  return status;
}



/** This overloaded function does not compute the kriging variance.
 * There are cases where the kriging variance is of no use (eg indicator
 * kriging). This function saves some time by not computing the kriging
 * variance.
 */ 
template<
         class MatrixLibrary,
         class Location,
         class Neighborhood,
         class Covariance,
         class KrigingConstraints,
         class Vector
        >
int kriging_weights(
		    Vector& weights,
		    const Location& center,
		    const Neighborhood& neighbors,
		    Covariance& covar,
		    KrigingConstraints& Kconstraints
		    ) {

  typedef matrix_lib_traits<MatrixLibrary> MatrixLib;
  
  // If the neighborhood is empty, there is no kriging to be done.
  if( neighbors.is_empty() ) {
    gstl_warning( "Empty neighborhood. No kriging to be done. " );
    return 2;
  }


  typename MatrixLib::Symmetric_matrix A;
  typename MatrixLib::Vector b;

  build_kriging_system(A,b,
		       weights, 
		       center, neighbors,
		       covar, Kconstraints);

  // solve the system
  int status = 
       kriging_constraints_traits<
                           KrigingConstraints,
                           MatrixLibrary
                                 >::kriging_solver(A, b, weights.begin());
  return status;
}




/** This overloaded function uses the default matrix library: TNT.
 */ 
template<
         class Location,
         class Neighborhood,
         class Covariance,
         class KrigingConstraints,
         class Vector
        >
inline int
kriging_weights(
		Vector& weights,
		double& kriging_variance,
		const Location& center,
		const Neighborhood& neighbors,
		Covariance& covar,
		KrigingConstraints& Kconstraints
		) {

  return kriging_weights< GSTL_TNT_lib >( weights, kriging_variance,
					  center,neighbors,
					  covar,Kconstraints );
}



/** This overloaded function uses the default matrix library: TNT.
 */ 
template<
         class Location,
         class Neighborhood,
         class Covariance,
         class KrigingConstraints,
         class Vector
        >
inline int
kriging_weights(
		Vector& weights,
		const Location& center,
		const Neighborhood& neighbors,
		Covariance& covar,
		KrigingConstraints& Kconstraints
		) {

  return kriging_weights< GSTL_TNT_lib >( weights,
					  center,neighbors,
					  covar,Kconstraints );
}



//===========================================


/** Compute the kriging weights, using a different covariance function for
 * the kriging matrix and the right-hand side of the kriging system. This
 * is useful for block kriging for example.
 * @return 0 if no problem was encountered, 1 if the kriging system could not be
 * solved, and 2 if no neighbors were provided (empty neighborhood).
 * 
 * @param weights is the vector where the kriging weights will be output
 * @param kriging_variance is where the kriging variance will be output
 * @param center is the location being estimated
 * @param neighbors points to the neighborhood of conditioning data used to 
 * estimate \c center. If there are no neighbors, no kriging is performed and  
 * the return code is 2.
 * @param covar is the covariance function used for the kriging matrix
 * @param covar_rhs is the covariance function used to compute the right-hand
 * side of the kriging system, ie the covariance between the data and the 
 * unknown.
 * @param Kconstraints are the kriging constraints (eg SK, OK, KT, ...).
 */

template<
         class MatrixLibrary,
         class Location,
         class Neighborhood,
         class Covariance,
         class Covariance2,
         class KrigingConstraints,
         class Vector
        >
int kriging_weights_2(
		    Vector& weights,
		    double& kriging_variance,
		    const Location& center,
		    const Neighborhood& neighbors,
		    Covariance& covar, Covariance2& covar_rhs,
		    KrigingConstraints& Kconstraints
		    ) {

  typedef matrix_lib_traits< MatrixLibrary > MatrixLib;
  
  // If the neighborhood is empty, there is no kriging to be done.
  if( neighbors.is_empty() ) {
    gstl_warning( "Empty neighborhood. No kriging to be done. " );
    return 2;
  }
  
  typename MatrixLib::Symmetric_matrix A;
  typename MatrixLib::Vector b;
  
  int conditioning_data = 
    build_kriging_system(A,b,
     		         weights, 
		         center, neighbors,
		         covar, covar_rhs, Kconstraints);



 //gstl_debug::print_kriging_system(A, b, std::cout);

  // solve the system
  int status = 
         kriging_constraints_traits<
                                    KrigingConstraints,
                                    MatrixLibrary
                                   >::const_kriging_solver(A, b, weights.begin());

  //gstl_debug::print_kriging_system(A, b, std::cout);
  //std::cout<<"Kriging Weights"<<std::endl;
  //gstl_debug::print_range(weights.begin(), weights.end(), std::cout);

  // Compute the kriging variance
  if(status == 0) {
    double C0=covar(center,center);
    kriging_variance =
      compute_kriging_variance(weights.begin(), 
                               weights.begin()+conditioning_data, weights.end(),
                               b, Kconstraints, center, C0);
  }
  else
    kriging_variance = -99;

  return status;
}



/** This overloaded function does not compute the kriging variance.
 * There are cases where the kriging variance is of no use (eg indicator
 * kriging). This function saves some time by not computing the kriging
 * variance.
 */ 
template<
         class MatrixLibrary,
         class Location,
         class Neighborhood,
         class Covariance,
         class Covariance2,
         class KrigingConstraints,
         class Vector
        >
int kriging_weights_2(
		    Vector& weights,
		    const Location& center,
		    const Neighborhood& neighbors,
		    Covariance& covar, Covariance2& covar_rhs,
		    KrigingConstraints& Kconstraints
		    ) {

  typedef matrix_lib_traits<MatrixLibrary> MatrixLib;
  
  // If the neighborhood is empty, there is no kriging to be done.
  if( neighbors.is_empty() ) {
    gstl_warning( "Empty neighborhood. No kriging to be done. " );
    return 2;
  }


  typename MatrixLib::Symmetric_matrix A;
  typename MatrixLib::Vector b;

  build_kriging_system(A,b,
		       weights, 
		       center, neighbors,
		       covar, covar_rhs, Kconstraints);


  // solve the system
  int status = 
       kriging_constraints_traits<
                           KrigingConstraints,
                           MatrixLibrary
                                 >::kriging_solver(A, b, weights.begin());


  return status;
}




/** This overloaded function uses the default matrix library: TNT.
 */ 
template<
         class Location,
         class Neighborhood,
         class Covariance,
         class Covariance2,
         class KrigingConstraints,
         class Vector
        >
inline int
kriging_weights_2(
		Vector& weights,
		double& kriging_variance,
		const Location& center,
		const Neighborhood& neighbors,
		Covariance& covar, Covariance2& covar_rhs,
		KrigingConstraints& Kconstraints
		) {

  return kriging_weights_2< GSTL_TNT_lib >( weights, kriging_variance,
					  center,neighbors,
					  covar, covar_rhs, Kconstraints );
}



/** This overloaded function uses the default matrix library: TNT.
 */ 
template<
         class Location,
         class Neighborhood,
         class Covariance,
         class Covariance2,
         class KrigingConstraints,
         class Vector
        >
inline int
kriging_weights_2(
		Vector& weights,
		const Location& center,
		const Neighborhood& neighbors,
		Covariance& covar, Covariance2& covar_rhs,
		KrigingConstraints& Kconstraints
		) {

  return kriging_weights_2< GSTL_TNT_lib >( weights,
					  center,neighbors,
					  covar, covar_rhs, Kconstraints );
}





//===========================================
// The function that does the job
template<
         class SymmetricMatrix,
         class MatVector,
         class Location,
         class Neighborhood,
         class Covariance,
         class KrigingConstraints,
         class Vector
        >
int build_kriging_system(
			  SymmetricMatrix& A, 
			  MatVector& b,
			  Vector& weights,
			  const Location& center,
			  const Neighborhood& neighbors,
			  Covariance& covar,
			  KrigingConstraints& Kconstraints
			  ) {

  int nb_conditioning_data = Kconstraints(A, b,
	                                  center, neighbors);
  
  build_invariant(A,b,
		  center,
		  neighbors.begin(), neighbors.end(),
		  covar);
  
  DEBUG_PRINT_KRIGING_SYSTEM( A,b);

  // Resize the output vector if necessary
  if( static_cast<int>(weights.size()) != static_cast<int>(A.num_cols()) )
    weights.resize(A.num_cols());

  return nb_conditioning_data;
}



template<
         class SymmetricMatrix,
         class MatVector,
         class Location,
         class Neighborhood,
         class Covariance,
         class Covariance2,
         class KrigingConstraints,
         class Vector
        >
int build_kriging_system(
			  SymmetricMatrix& A, 
			  MatVector& b,
			  Vector& weights,
			  const Location& center,
			  const Neighborhood& neighbors,
			  Covariance& covar, Covariance2& covar_rhs,
			  KrigingConstraints& Kconstraints
			  ) {

  int nb_conditioning_data = Kconstraints(A, b,
	                                   center, neighbors);
  
  build_invariant(A, b,
		              center,
		              neighbors.begin(), neighbors.end(),
		              covar, covar_rhs);
  
  DEBUG_PRINT_KRIGING_SYSTEM( A,b);

  // Resize the output vector if necessary
  if( static_cast<int>(weights.size()) != static_cast<int>(A.num_cols()) )
    weights.resize(A.num_cols());

  return nb_conditioning_data;
}

		 template<
			 class SymmetricMatrix,
			 class MatVector,
			 class Location,
			 class Neighborhood,
			 class Covariance,
			 class Covariance2,
			 class KrigingConstraints,
			 class Vector
		 >
		 double build_idw_system(
		 SymmetricMatrix& A, 
		 MatVector& b,
		 Vector& weights,
		 const Location& center,
		 const Neighborhood& neighbors,
		 Covariance& covar, Covariance2& covar_rhs,
		 KrigingConstraints& Kconstraints
		 ) {		 


			 DEBUG_PRINT_KRIGING_SYSTEM( A,b);
			 int nb_conditioning_data = b.size();

			 // 参考源自文献: statistical apprach to inverse distance interpolation 
			 // 利用 A b 及weights计算KRG方差
			 double dKrgVariance = 0.0;
			 //double dPart1 = 0.0, dPart2 = 0.0, dPart3 = 0.0;

			 dKrgVariance += covar(Location(0, 0, 0), Location(0, 0, 0));
			 //dPart1 = covar(Location(0, 0), Location(0, 0));

			 int dim = A.num_rows();
			 // Since the LU algorithm of TNT does not use the fact that A is symetric,
			 // we have to copy the upper part of A to its lower part.

			 int i=1;
			 int j;

			 typedef typename Neighborhood::const_iterator InputIterator;

			 for(InputIterator row = neighbors.begin(); row != neighbors.end(); row++ ) 
			 {
				 j=1;

				 for(InputIterator col = neighbors.begin() ; col != neighbors.end() ; col++)	  
				 {
					 //double dA_IJ = covar( row->location(), col->location() );
					 //double dA_IJ_temp = A(i,j);
					 //dPart2 += ( weights[i-1]*weights[j-1]*A(i,j) );

					 dKrgVariance += ( weights[i-1]*weights[j-1]*A(i,j) );

					 j++;
				 }

				 //dKrgVariance -= (weights[i-1]*2.0*b(i));

				 i++;
			 }

			 // 下面的过程, 可以通过STD的内积 inner_product来计算得到
			 i = 1;
			 for(InputIterator row = neighbors.begin(); row != neighbors.end(); row++ ) 
			 {
				 //double dB_I = covar_rhs( row->location(), center );
				 //double dB_I_temp = b(i);
				 //dPart3 += (weights[i-1]*2.0*b(i));

				 dKrgVariance -= (weights[i-1]*2.0*b(i));

				 i++;
			 }		

			 return dKrgVariance;
			 }


			 template<
				 class SymmetricMatrix,
				 class MatVector,
				 class Location,
				 class Covariance,
				 class Vector
			 >
			 double cal_kring_variance_by_weights(
			 SymmetricMatrix& A, 
			 MatVector& b,
			 Vector& weights,
			 const Location& center,
			 Covariance& covar
			 )  
				 {

					 DEBUG_PRINT_KRIGING_SYSTEM( A,b);
					 int nb_conditioning_data = b.size() - 1;

					 // 参考源自文献: statistical apprach to inverse distance interpolation 
					 // 利用 A b 及weights计算KRG方差
					 double dKrgVariance = 0.0;
					 //double dPart1 = 0.0, dPart2 = 0.0, dPart3 = 0.0;

					 dKrgVariance += covar(Location(0, 0, 0), Location(0, 0, 0));
					 //dPart1 = covar(Location(0, 0), Location(0, 0));

					 int dim = A.num_rows();
					 // Since the LU algorithm of TNT does not use the fact that A is symetric,
					 // we have to copy the upper part of A to its lower part.

					 int i=1;
					 int j;

					 typedef typename Neighborhood::const_iterator InputIterator;

					 //for(InputIterator row = neighbors.begin(); row != neighbors.end(); row++ ) 
					 while(i < dim)
					 {
						 j=1;

						 //for(InputIterator col = neighbors.begin() ; col != neighbors.end() ; col++)	  
						 while(j < dim)
						 {
							 //double dA_IJ = covar( row->location(), col->location() );
							 //double dA_IJ_temp = A(i,j);
							 //dPart2 += ( weights[i-1]*weights[j-1]*A(i,j) );

							 dKrgVariance += ( weights[i-1]*weights[j-1]*A(i,j) );

							 j++;
						 }

						 //dKrgVariance -= (weights[i-1]*2.0*b(i));

						 i++;
					 }

					 // 下面的过程, 可以通过STD的内积 inner_product来计算得到
					 i = 1;
					 //for(InputIterator row = neighbors.begin(); row != neighbors.end(); row++ ) 
					 while(i < dim)
					 {
						 //double dB_I = covar_rhs( row->location(), center );
						 //double dB_I_temp = b(i);
						 //dPart3 += (weights[i-1]*2.0*b(i));

						 dKrgVariance -= (weights[i-1]*2.0*b(i));

						 i++;
					 }		

					 return dKrgVariance;
				 }		
	
	// 使用常数幂指数 double dPower = 2.0
	// IDW 相关功能的实现, 相对于同一待估点, 不同的样点使用不同的幂指数
	template<
		class MatrixLibrary,
		class Location,
		class Neighborhood,
		class Covariance,
		class KrigingConstraints,
		class Vector
	>
	int idw_weights(
		Vector& weights,
		const Location& center,
		const Neighborhood& neighbors,
		Covariance& covar,
		KrigingConstraints& Kconstraints,
		double dPower = 2.0,
		double* dEstimatedValue = 0,		//	估计值
		double* dKriging_variance = 0,
		double* dInterpolationVariance = 0	//	插值方差
	) 
	{
		typedef matrix_lib_traits< MatrixLibrary > MatrixLib;

		// If the neighborhood is empty, there is no idw to be done.
		if( neighbors.is_empty() ) {
			gstl_warning( "Empty neighborhood. No kriging to be done. " );
			return 2;
		}  	


		std::vector<double> vDduDistance;
		double dTempEduDistance = 0.0, dSumEduDistance = 0.0;


		//typedef typename Neighborhood::location_type location_type;
		for(Neighborhood::const_iterator neigh_it=neighbors.begin(); 
			neigh_it!=neighbors.end(); 
			neigh_it++)
		{
			dTempEduDistance = sqrt( double(			
				(neigh_it->location()[0] - center[0])*(neigh_it->location()[0] - center[0]) +
				(neigh_it->location()[1] - center[1])*(neigh_it->location()[1] - center[1]) +
				(neigh_it->location()[2] - center[2])*(neigh_it->location()[2] - center[2])
				)
				);

			if (dTempEduDistance < 1.0E-10)
			{
				dTempEduDistance = 1.0E-10;
			}

			//dSumEduDistance_ids += dTempEduDistance;

			dTempEduDistance += 0;

			//dTempEduDistance = sqrt( covar(center, center) - covar(neigh_it->get_location(), center) );
			//dTempEduDistance = covar(center, center) - covar(neigh_it->get_location(), center);

			//if (dKriging_variance && *dKriging_variance > -0.1)
			//{
			//	//dTempEduDistance = pow(dTempEduDistance*(dKriging_variance[iiiii++]), -dPower);
			//	dTempEduDistance = pow(dTempEduDistance, -(dKriging_variance[iiiii++]));
			//}
			//else
			{
				dTempEduDistance = pow(dTempEduDistance, -dPower);


				//dTempEduDistance_ids = pow(dTempEduDistance_ids, -1.0);
			}

			dSumEduDistance += dTempEduDistance;
			vDduDistance.push_back(dTempEduDistance);


			//dSumEduDistance_ids += dTempEduDistance_ids;
			//vDduDistance_ids.push_back(dTempEduDistance_ids);
		}

		//if (dKriging_variance && *dKriging_variance > -0.1)
		//{
		//	delete []dKriging_variance;
		//	dKriging_variance = 0;
		//}

		weights.resize(neighbors.size());
		Vector::iterator it = weights.begin();
		for (std::vector<double>::iterator it_D = vDduDistance.begin(); 
			it != weights.end(); 
			it++, it_D++)
		{
			*it = (*it_D)/dSumEduDistance;
		}

		if (dEstimatedValue)
		{
			*dEstimatedValue = linear_combination(weights.begin(), weights.end(), neighbors);
		}


		if (dEstimatedValue && dInterpolationVariance)
		{
			*dInterpolationVariance = 0.0;
			Vector::iterator it = weights.begin();
			for(Neighborhood::const_iterator neigh_it=neighbors.begin(); 
				neigh_it!=neighbors.end(); 
				neigh_it++, it++)
			{
				double dTemp = neigh_it->property_value();
				double dWei = *it;
				*dInterpolationVariance += ( dWei*(*dEstimatedValue - dTemp)*(*dEstimatedValue - dTemp) );
			}
		}

		// 如果需要计算 Kriging_variance 的话,则通过协方差矩阵进行计算
		if (dKriging_variance)
		{
			typename MatrixLib::Symmetric_matrix A;
			typename MatrixLib::Vector b;

			// 之前已经调用过 build_kriging_system 进行计算了
			// 建立克里格方程系统, 写 A b 矩阵的值
			//int conditioning_data = 
			//	build_kriging_system(A,b,
			//	weights, 
			//	center, neighbors,
			//	covar, Kconstraints);

			// 设置 A b 的维数
			int nb_conditioning_data = Kconstraints(A, b,
				center, neighbors);
			// 设置A b中元素的值
			build_invariant(A,b,
				center,
				neighbors.begin(), neighbors.end(),
				covar);

			// 计算估计方差(即KRG方差)
			// 使用通用公式, 由 样品之间\样品与待估点之间的协方差 样品权值, 即矩阵A和B, 计算IDW对应的KRG方差
			*dKriging_variance = cal_kring_variance_by_weights(A, b,
				weights, 
				center,
				covar);
				//,
				//covar, covar, Kconstraints
		}

		return 0;
	}

	// 增加了一个 幂次参数 Vector& dvPowers,
	// IDW 相关功能的实现, 相对于同一待估点, 不同的样点使用不同的幂指数
	template<
		class MatrixLibrary,
		class Location,
		class Neighborhood,
		class Covariance,
		class KrigingConstraints,
		class Vector
	>
	int idw_weights(
	Vector& weights,
	const Location& center,
	const Neighborhood& neighbors,
	Covariance& covar,
	KrigingConstraints& Kconstraints,
	Vector& dvPowers,
	double* dEstimatedValue = 0,		//	估计值
	double* dKriging_variance = 0,
	double* dInterpolationVariance = 0	//	插值方差
	) 
	{
		typedef matrix_lib_traits< MatrixLibrary > MatrixLib;

		// If the neighborhood is empty, there is no idw to be done.
		if( neighbors.is_empty() ) {
			gstl_warning( "Empty neighborhood. No kriging to be done. " );
			return 2;
		}  	

		// 记录待估点与样品点之间的欧氏距离
		std::vector<double> vDduDistance;
		double dTempEduDistance = 0.0, dSumEduDistance = 0.0;

		// 遍历 Neighborhood 中的所有样品数据, 计算其与等估点之间的 距离

		Vector::iterator PowerIt = dvPowers.begin();
		double dPower = *PowerIt;

		for(Neighborhood::const_iterator neigh_it=neighbors.begin(); 
			neigh_it!=neighbors.end(); 
			neigh_it++)
		{
			dTempEduDistance = sqrt( double(			
				(neigh_it->get_location()[0] - center[0])*(neigh_it->get_location()[0] - center[0]) +
				(neigh_it->get_location()[1] - center[1])*(neigh_it->get_location()[1] - center[1]) +
				(neigh_it->get_location()[2] - center[2])*(neigh_it->get_location()[2] - center[2])
				)
				);

			dTempEduDistance += 0;

			if (PowerIt != dvPowers.end())
			{
				dPower = *PowerIt;
				PowerIt++;
			}
			dTempEduDistance = pow(dTempEduDistance, -dPower);
			//dTempEduDistance = dPower*pow(dTempEduDistance, -2.0);

			dSumEduDistance += dTempEduDistance;

			vDduDistance.push_back(dTempEduDistance);
		}

		// 计算IDW权值
		weights.resize(neighbors.size());
		Vector::iterator it = weights.begin();
		for (std::vector<double>::iterator it_D = vDduDistance.begin(); 
			it != weights.end(); 
			it++, it_D++)
		{
			*it = (*it_D)/dSumEduDistance;
		}

		// 记录IDW的估计值
		if (dEstimatedValue)
		{
			*dEstimatedValue = linear_combination(weights.begin(), weights.end(), neighbors);
		}

		// 计算插值方差
		if (dEstimatedValue && dInterpolationVariance)
		{
			*dInterpolationVariance = 0.0;
			Vector::iterator it = weights.begin();
			for(Neighborhood::const_iterator neigh_it=neighbors.begin(); 
				neigh_it!=neighbors.end(); 
				neigh_it++, it++)
			{
				double dTemp = neigh_it->get_property_value();
				double dWei = *it;
				*dInterpolationVariance += ( dWei*(*dEstimatedValue - dTemp)*(*dEstimatedValue - dTemp) );
			}
		}

		// 如果需要计算 Kriging_variance 的话,则通过协方差矩阵进行计算
		if (dKriging_variance)
		{
			typename MatrixLib::Symmetric_matrix A;
			typename MatrixLib::Vector b;

			// 之前已经调用过 build_kriging_system 进行计算了
			// 建立克里格方程系统, 写 A b 矩阵的值
			//int conditioning_data = 
			//	build_kriging_system(A,b,
			//	weights, 
			//	center, neighbors,
			//	covar, Kconstraints);

			// 设置 A b 的维数
			int nb_conditioning_data = Kconstraints(A, b,
				center, neighbors);
			// 设置A b中元素的值
			build_invariant(A,b,
				center,
				neighbors.begin(), neighbors.end(),
				covar);

			// 计算估计方差(即KRG方差)
			// 使用通用公式, 由 样品之间\样品与待估点之间的协方差 样品权值, 即矩阵A和B, 计算IDW对应的KRG方差
			*dKriging_variance= cal_kring_variance_by_weights(A,b,
				weights, 
				center,
				covar);
				//,
				//covar, covar, Kconstraints
		}


		// solve the system
		int status = 0;
		return status;
	}

#endif
