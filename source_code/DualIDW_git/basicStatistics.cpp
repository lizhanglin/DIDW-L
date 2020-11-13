#include "basicStatistics.h"

 double GetAverage(vector<double>& daZ, vector<double>* paAlpha/* = NULL*/)
 {
	 double dAverage = 0.0;
	 int i=0, 
		 n = daZ.size();
	 if(n>0)
	 {
		 if(paAlpha && (*paAlpha).size()==n)
		 {
			 for(i=0; i<n; i++)
			 {
				 dAverage += (*paAlpha)[i]*daZ[i];
			 }
		 }
		 else
		 {
			 for(i=0; i<n; i++)
			 {
				 dAverage += daZ[i];
			 }
			 dAverage /= n;
		 }
	 }
	 return dAverage;
 }


double GetVariance(vector<double>& daX)
{
	double S2 = 0.0, dAver = 0.0;
	int i=0, n = daX.size();
	//if(n>0)
	//{
	//	dAver = GetAverage(daX);
	//	for(i=0; i<n; i++)
	//	{
	//		S2 += daX[i]*daX[i];
	//	}
	//	S2 = S2/n - dAver*dAver;
	//}

	if(n>1)
	{
		dAver = GetAverage(daX);
		for(i=0; i<n; i++)
		{
			S2 += ( (daX[i] - dAver)*(daX[i] - dAver) ) ;
		}
		S2 = S2/(n-1);
	}
	else
	{
		S2 = -1;
	}
	
	return S2;
}

double GetVariance(vector<double>& daX, vector<double>& daExceptVals)
{
	vector<double> daXCopy(daX);

	for (int i = 0; i<daExceptVals.size(); i++)
	{
		vector<double>::iterator result = find( daXCopy.begin( ), daXCopy.end( ), daExceptVals[i] );

		if  ( result != daXCopy.end() ) 
		{
			daXCopy.erase(result);
		}		
	}
	return GetVariance(daXCopy);	
}


double GetMedian(vector<double>& daVals)
{
	vector<double> daZ(daVals);
	daZ.resize(daVals.size());

	std::copy(daVals.begin(), daVals.end(), daZ.begin());

	double  dMedianV = 0.0;

	std::sort(daZ.begin(), daZ.end());

	if (daZ.size() %2 == 1)
	{
		dMedianV = *(daZ.begin() + daZ.size()/2);
	}
	else
	{
		dMedianV = *(daZ.begin() + daZ.size()/2) + *(daZ.begin() + daZ.size()/2 - 1);
		dMedianV *= 0.5;
	}

	return dMedianV;
}


double GetAverage(vector<double>& daZ, vector<double>& daExceptVals)
{
	vector<double> daXCopy(daZ);

	for (int i = 0; i<daExceptVals.size(); i++)
	{
		vector<double>::iterator result = find( daXCopy.begin( ), daXCopy.end( ), daExceptVals[i] );

		if  ( result != daXCopy.end() ) 
		{
			daXCopy.erase(result);
		}		
	}
	return GetAverage(daXCopy);
}

double GetSum(vector<double>& daVals)
{
	double dSum = 0.0;

	vector<double>::iterator vdIterator = daVals.begin();

	for (; vdIterator != daVals.end(); vdIterator++)
	{
		dSum += (*vdIterator);
	}

	return dSum;
}


double GetMTE(vector<double>& daMeasuredVals, vector<double>& daEvalutedVals, vector<double>* pWeights /*= 0*/, bool bStandard /*= true*/)
{
	double dStdVariance = 1.0;
	if (bStandard)
	{
		dStdVariance = sqrt(GetVariance(daMeasuredVals));
	} 

	vector<double> daTemp;
	daTemp.resize(daMeasuredVals.size());


	if (pWeights && pWeights->size() == daMeasuredVals.size())
	{
		vector<double>::iterator daTempIterator, daMeasuredValsIterator, daEvalutedValsIterator,
			daWeightIterator = pWeights->begin();

		daTempIterator = daTemp.begin();
		daMeasuredValsIterator = daMeasuredVals.begin();
		daEvalutedValsIterator = daEvalutedVals.begin();

		for (; daTempIterator != daTemp.end(); daTempIterator++)
		{
			*daTempIterator = (*daWeightIterator)*(*daEvalutedValsIterator - *daMeasuredValsIterator);

			daMeasuredValsIterator++;
			daEvalutedValsIterator++;
			daWeightIterator++;
		}

		return GetSum(daTemp)/dStdVariance;
	}
	else
	{
		vector<double>::iterator daTempIterator, daMeasuredValsIterator, daEvalutedValsIterator;

		daTempIterator = daTemp.begin();
		daMeasuredValsIterator = daMeasuredVals.begin();
		daEvalutedValsIterator = daEvalutedVals.begin();

		for (; daTempIterator != daTemp.end(); daTempIterator++)
		{
			*daTempIterator = (*daEvalutedValsIterator - *daMeasuredValsIterator);

			daMeasuredValsIterator++;
			daEvalutedValsIterator++;
		}

		return GetAverage(daTemp)/dStdVariance;
	}

}

// mean absolute error 
double GetMAE(vector<double>& daMeasuredVals, vector<double>& daEvalutedVals, vector<double>* pWeights /*= 0*/, bool bStandard /*= true*/)
{
	double dStdVariance = 1.0;
	if (bStandard)
	{
		dStdVariance = sqrt(GetVariance(daMeasuredVals));
	} 

	vector<double> daTemp;
	daTemp.resize(daMeasuredVals.size());


	if (pWeights && pWeights->size() == daMeasuredVals.size())
	{
		vector<double>::iterator daTempIterator, daMeasuredValsIterator, daEvalutedValsIterator,
			daWeightIterator = pWeights->begin();

		daTempIterator = daTemp.begin();
		daMeasuredValsIterator = daMeasuredVals.begin();
		daEvalutedValsIterator = daEvalutedVals.begin();

		for (; daTempIterator != daTemp.end(); daTempIterator++)
		{
			*daTempIterator = (*daWeightIterator)*fabs(*daEvalutedValsIterator - *daMeasuredValsIterator);

			daMeasuredValsIterator++;
			daEvalutedValsIterator++;
			daWeightIterator++;
		}

		return GetSum(daTemp)/dStdVariance;
	}
	else
	{
		vector<double>::iterator daTempIterator, daMeasuredValsIterator, daEvalutedValsIterator;

		daTempIterator = daTemp.begin();
		daMeasuredValsIterator = daMeasuredVals.begin();
		daEvalutedValsIterator = daEvalutedVals.begin();

		for (; daTempIterator != daTemp.end(); daTempIterator++)
		{
			*daTempIterator = fabs(*daEvalutedValsIterator - *daMeasuredValsIterator);

			daMeasuredValsIterator++;
			daEvalutedValsIterator++;
		}

		return GetAverage(daTemp)/dStdVariance;
	}

} 

// mean relative error
double GetMRE(vector<double>& daMeasuredVals, vector<double>& daEvalutedVals, vector<double>* pWeights /*= 0*/)
{
	vector<double> daTemp;
	daTemp.resize(daMeasuredVals.size());

	if (pWeights && pWeights->size() == daMeasuredVals.size())
	{
		vector<double>::iterator daTempIterator, daMeasuredValsIterator, daEvalutedValsIterator,
			daWeightIterator = pWeights->begin();

		daTempIterator = daTemp.begin();
		daMeasuredValsIterator = daMeasuredVals.begin();
		daEvalutedValsIterator = daEvalutedVals.begin();

		for (; daTempIterator != daTemp.end(); daTempIterator++)
		{
			*daTempIterator = (*daWeightIterator)*fabs(*daEvalutedValsIterator - *daMeasuredValsIterator);

			if (fabs(*daMeasuredValsIterator) < 0.000001)
			{
				*daTempIterator = (*daTempIterator)/0.000001;	
			}
			else
			{
				*daTempIterator = (*daTempIterator)/(*daMeasuredValsIterator);
			}

			daMeasuredValsIterator++;
			daEvalutedValsIterator++;
			daWeightIterator++;
		}
	}
	else
	{
		vector<double>::iterator daTempIterator, daMeasuredValsIterator, daEvalutedValsIterator;

		daTempIterator = daTemp.begin();
		daMeasuredValsIterator = daMeasuredVals.begin();
		daEvalutedValsIterator = daEvalutedVals.begin();

		for (; daTempIterator != daTemp.end(); daTempIterator++)
		{
			*daTempIterator = fabs(*daEvalutedValsIterator - *daMeasuredValsIterator);

			if (fabs(*daMeasuredValsIterator) < 0.000001)
			{
				*daTempIterator = (*daTempIterator)/0.000001;
			}
			else
			{
				*daTempIterator = (*daTempIterator)/(*daMeasuredValsIterator);
			}
			
			daMeasuredValsIterator++;
			daEvalutedValsIterator++;
		}
	}

	return GetAverage(daTemp);
}  

// root mean squared error
double GetRMSE(vector<double>& daMeasuredVals, vector<double>& daEvalutedVals, vector<double>* pWeights /*= 0*/, bool bStandard /*= true*/)
{
	double dStdVariance = 1.0;
	if (bStandard)
	{
		dStdVariance = sqrt(GetVariance(daMeasuredVals));
	} 

	double dTemp = 0.0;

	if (pWeights && pWeights->size() == daMeasuredVals.size())
	{
		vector<double>::iterator daMeasuredValsIterator, daEvalutedValsIterator,
			daWeightIterator = pWeights->begin();

		daMeasuredValsIterator = daMeasuredVals.begin();
		daEvalutedValsIterator = daEvalutedVals.begin();

		for (; daMeasuredValsIterator != daMeasuredVals.end(); 
			daMeasuredValsIterator++, daEvalutedValsIterator++, daWeightIterator++)
		{
			dTemp += /*(*daWeightIterator)**/(*daWeightIterator)*(*daEvalutedValsIterator - *daMeasuredValsIterator)*(*daEvalutedValsIterator - *daMeasuredValsIterator);			
		}

		dTemp = sqrt(dTemp)/dStdVariance;		
	}
	else
	{
		vector<double>::iterator daMeasuredValsIterator, daEvalutedValsIterator;

		daMeasuredValsIterator = daMeasuredVals.begin();
		daEvalutedValsIterator = daEvalutedVals.begin();

		for (; daMeasuredValsIterator != daMeasuredVals.end(); 
			daMeasuredValsIterator++, daEvalutedValsIterator++)
		{
			dTemp += (*daEvalutedValsIterator - *daMeasuredValsIterator)*(*daEvalutedValsIterator - *daMeasuredValsIterator);
		}

		dTemp = sqrt(dTemp/daMeasuredVals.size())/dStdVariance;
	}

	return dTemp;
} 

// correlation coefficient
double GetCC(vector<double>& daMeasuredVals, vector<double>& daEvalutedVals)
{
	double dMeasuredMeanVals = GetAverage(daMeasuredVals);
	double dEvalutedMeanVals = GetAverage(daEvalutedVals);

	vector<double> daTemp1, daTemp2, daTemp3;
	daTemp1.resize(daMeasuredVals.size());
	daTemp2.resize(daMeasuredVals.size());
	daTemp3.resize(daMeasuredVals.size());

	vector<double>::iterator daTempIterator1 = daTemp1.begin();
	vector<double>::iterator daTempIterator2 = daTemp2.begin();
	vector<double>::iterator daTempIterator3 = daTemp3.begin();

	vector<double>::iterator daMeasuredValsIterator, daEvalutedValsIterator;

	daMeasuredValsIterator = daMeasuredVals.begin();
	daEvalutedValsIterator = daEvalutedVals.begin();

	for (; daMeasuredValsIterator != daMeasuredVals.end(); 
				daMeasuredValsIterator++, daEvalutedValsIterator++,
				daTempIterator1++,daTempIterator2++,daTempIterator3++)
	{
		*daTempIterator1 = (*daEvalutedValsIterator - dEvalutedMeanVals)*(*daMeasuredValsIterator - dMeasuredMeanVals);

		*daTempIterator2 = (*daEvalutedValsIterator - dEvalutedMeanVals)*(*daEvalutedValsIterator - dEvalutedMeanVals);

		*daTempIterator3 = (*daMeasuredValsIterator - dMeasuredMeanVals)*(*daMeasuredValsIterator - dMeasuredMeanVals);
	}

	double dNumerator = GetSum(daTemp1);

	double dDenominator = sqrt( GetSum(daTemp2)*GetSum(daTemp3) );

	if(fabs(dDenominator) < 0.1e-12)
	{
		return 9.9E12;
	}
	
	return dNumerator/dDenominator;
} 

void GetRank(vector<double>& dvOriVals, vector<double>& viRanks)
{
	std::vector<double> dvVals(dvOriVals);
	std::sort(dvVals.begin(), dvVals.end());

	for (int i = 0; i<dvOriVals.size(); i++)
	{
		std::vector<double>::iterator iter = std::lower_bound(dvVals.begin(), dvVals.end(), dvOriVals[i]);
		//std::cout << "index " << int(iter - v.begin()) << std::endl;
		viRanks.push_back(int( 1 + iter - dvVals.begin()));		
	}

	return;
}

// rank correlation coefficient
double GetRankCC(vector<double>& daMeasuredVals, vector<double>& daEvalutedVals)
{
	vector<double> viRanks1, viRanks2;

	GetRank(daMeasuredVals, viRanks1);
	GetRank(daEvalutedVals, viRanks2);

	return GetCC(viRanks1, viRanks2);
} 

// Covariance
double GetCovariance(vector<double>& daMeasuredVals, vector<double>& daEvalutedVals)
{
	double dMeasuredMeanVals = GetAverage(daMeasuredVals);
	double dEvalutedMeanVals = GetAverage(daEvalutedVals);

	vector<double> daTemp1;
	daTemp1.resize(daMeasuredVals.size());

	vector<double>::iterator daTempIterator1 = daTemp1.begin();

	vector<double>::iterator daMeasuredValsIterator = daMeasuredVals.begin();
	vector<double>::iterator daEvalutedValsIterator = daEvalutedVals.begin();

	for (; daMeasuredValsIterator != daMeasuredVals.end(); 
		daMeasuredValsIterator++, daEvalutedValsIterator++, daTempIterator1++)
	{
		*daTempIterator1 = (*daEvalutedValsIterator - dEvalutedMeanVals)*(*daMeasuredValsIterator - dMeasuredMeanVals);
	}

	double dNumerator = GetAverage(daTemp1);

	return dNumerator;
} 


// Standardized squared deviation
double GetSQD(vector<double>& daMeasuredVals, vector<double>& daEvalutedVals, vector<double>& daKrigingVariance)
{
	double dMeasuredMeanVals = GetAverage(daMeasuredVals);
	double dEvalutedMeanVals = GetAverage(daEvalutedVals);

	vector<double> daTemp1;
	daTemp1.resize(daMeasuredVals.size());
	vector<double>::iterator daTemp1It = daTemp1.begin();

	vector<double>::iterator daMeasuredValsIt = daMeasuredVals.begin();
	vector<double>::iterator daEvalutedValsIt = daEvalutedVals.begin();
	vector<double>::iterator daKrigingVarianceIt = daKrigingVariance.begin();

	double dTempV = 0;

	for (; daTemp1It != daTemp1.end(); daTemp1It++,
		daMeasuredValsIt++,daEvalutedValsIt++,daKrigingVarianceIt++)
	{
		double dTempV = (*daEvalutedValsIt - *daMeasuredValsIt)*(*daEvalutedValsIt - *daMeasuredValsIt);
		
		(*daTemp1It) = dTempV/(*daKrigingVarianceIt);
	}

	double dAverageV = GetAverage(daTemp1);
	double  dMedianV = GetMedian(daTemp1);
	

	double dScore1 = fabs(1.0 - dAverageV);

	double dScore2 = fabs(0.445 - dMedianV);
	
	return (dScore1 + dScore2);
} 


double GetCV(vector<double>& daVals)
{
	if ( daVals.size() < 2 || fabs(GetAverage(daVals) ) < 0.000001)
	{
		return 0.0;
	}
	else
	{
		return ( sqrt(GetVariance(daVals)) / GetAverage(daVals));
	}
} 