

#ifndef _BASIC_STATISTICS_SUPPORT_			
#define _BASIC_STATISTICS_SUPPORT_
						

#include <stdio.h>
#include <tchar.h>
#include <math.h>

#include <vector>
#include <string>
#include <algorithm>

using namespace std;

#define PAI				3.1415926535897932384626433832795	// ¦Ð
#define PiOver180		1.74532925199433E-002				 
#define PiUnder180		5.72957795130823E+001				 
#define SMALL_NUMBER_POSITIVE	0.000001					 
#define LARGE_NUMBER	1E20	


double GetAverage(vector<double>& daZ, vector<double>* paAlpha = NULL); 

double GetVariance(vector<double>& daX);

double GetVariance(vector<double>& daX, vector<double>& daExceptVals);

double GetMedian(vector<double>& daVals);

double GetAverage(vector<double>& daZ, vector<double>& daExceptVals);

double GetSum(vector<double>& daVals);

double GetMTE(vector<double>& daMeasuredVals, vector<double>& daEvalutedVals, vector<double>* pWeights = 0, bool bStandard = true); 


// mean absolute error 
double GetMAE(vector<double>& daMeasuredVals, vector<double>& daEvalutedVals, vector<double>* pWeights = 0, bool bStandard = true); 

// mean relative error
double GetMRE(vector<double>& daMeasuredVals, vector<double>& daEvalutedVals, vector<double>* pWeights = 0);  

// root mean squared error
double GetRMSE(vector<double>& daMeasuredVals, vector<double>& daEvalutedVals, vector<double>* pWeights = 0, bool bStandard = true);

// correlation coefficient
double GetCC(vector<double>& daMeasuredVals, vector<double>& daEvalutedVals);

// rank correlation coefficient
double GetRankCC(vector<double>& daMeasuredVals, vector<double>& daEvalutedVals);

void GetRank(vector<double>& dvOriVals, vector<double>& viRanks);

// Covariance
double GetCovariance(vector<double>& daMeasuredVals, vector<double>& daEvalutedVals);

// Standardized squared deviation
double GetSQD(vector<double>& daMeasuredVals, vector<double>& daEvalutedVals, vector<double>& daKrigingVariance);

// Coefficient of variation
double GetCV(vector<double>& daVals); 

#endif