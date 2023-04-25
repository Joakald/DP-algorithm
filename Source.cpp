#include <vector>
#include <cassert>
#include <fstream>
#include <iostream>
#include <iomanip>
#include <string>
#include <random>
#include <algorithm>
#include <limits>

#include "Matrix.h"
#include "Settings.h"

// interpolation
double Interpolate(double w, int time, MatrixInt& CMP)
{
	int IntPart = (int)w;
	double alpha = w - (double)IntPart;
	while (IntPart >= CMP.rows)
		IntPart--;
	
	if ( IntPart + 1 < CMP.rows )
		return std::max(0.0,(1.0 - alpha) * CMP.Get(IntPart,time) + alpha * CMP.Get(IntPart + 1,time));
	
	return CMP.Get(IntPart, time);
}

// get objective given time, wealth and specific consumptions, given that future consumptions are indirectly determined by current consumptions
double GetObjectiveCooperationStochastic(int time, int wealth, int cp1, int cp2, MatrixInt& CMP1, MatrixInt& CMP2, int MaxWealth)
{
	// assertions
	assert(time >= 0);
	assert(time < nPeriods);
	assert(wealth >= 0);
	assert(wealth <= MaxWealth);
	assert(cp1 >= 0);
	assert(cp2 >= 0);
	assert(cp1 + cp2 <= wealth);
	assert(CMP1.rows >= MaxWealth + 1);
	assert(CMP1.columns == nPeriods);
	assert(CMP2.rows >= MaxWealth + 1);
	assert(CMP2.columns == nPeriods);

	int nShocks = nPeriods - time - 1;
	int nOutcomes = (int)pow(2, nShocks);
	std::vector<int> ShockVector(nShocks, 1);

	double Objective = 0.0;
	
	// holders
	double hcp1;
	double hcp2;
	double hV;
	double hw;

	// loop over all the outcomes and add the probability weighted value function to the objective
	for (int i = 1; i <= nOutcomes; i++)
	{	
		// reset holders at the start of every outcome
		hcp1 = cp1;
		hcp2 = cp2;
		hV = 0.0;
		hw = wealth;

		// loop from time t to nPeriods - 1
		for (int t = time; t < nPeriods; t++)
		{
			// get consumptions for the period t by looking in the consumption matrices	
			if (t > time)
			{
				hcp1 = Interpolate(hw, t, CMP1);
				hcp2 = Interpolate(hw, t, CMP2);
			}
			// check if the terminal condition is satisfied
			if (t == nPeriods - 1)
			{
				if (!TerminalCondition(hw - hcp1 - hcp2))
				{
					return -std::numeric_limits<double>::max();
				}
			}

			if (hcp1 < 0.0)
			{
				int x = hw;
			}
			assert(hcp1 >= 0.0);
			assert(hcp2 >= 0.0);
			assert(hcp1 + hcp2 <= hw);

			// get the instantaneuos utility of period t, discounted to period time, and add it to the objective holder
			hV += Weight(wealth, GammaP1, LambdaP1) * DiscountFunction(t - time, DeltaP1, BetaP1, wealth) * UtilityFunction(hcp1, hcp2, hw, GammaP1, alphaP1);
			hV += Weight(wealth, GammaP2, LambdaP2) * DiscountFunction(t - time, DeltaP2, BetaP2, wealth) * UtilityFunction(hcp2, hcp1, hw, GammaP2, alphaP2);

			// set the wealth holder to the wealth of the next period
			if (t < nPeriods - 1)
			{
				hw = Dynamics(hw, hcp1, hcp2, ShockVector[t - time]);
				assert(hw >= 0);
			}
		}
		// add the value of the value function for this partiular shock to the objetive function
		Objective =+ hV;

		// modify shocks to get the next outcome
		int j = nShocks;
		while (j > 0)
		{
			if (i % (int)pow(2,nShocks - j) == 0)
				ShockVector[j - 1] *= -1;
			j--;
		}
	}
	// divide the objective by the number of potential shock outcomes to get the probability weighted average
	Objective /= nOutcomes;
	return Objective;
}

// get objective given time, wealth and specific consumptions, given that future consumptions are indirectly determined by current consumptions
double GetObjectiveCooperation(int time, int wealth, int cp1, int cp2, MatrixInt& CMP1, MatrixInt& CMP2, int MaxWealth)
{
	// assertions
	assert(time >= 0);
	assert(time < nPeriods);
	assert(wealth >= 0);
	assert(wealth <= MaxWealth);
	assert(cp1 >= 0);
	assert(cp2 >= 0);
	assert(cp1 + cp2 <= (int)(wealth * (1.0 + Growth)));
	assert(CMP1.rows >= MaxWealth + 1);
	assert(CMP1.columns == nPeriods);
	assert(CMP2.rows >= MaxWealth + 1);
	assert(CMP2.columns == nPeriods);

	double Objective = 0.0;

	// holders
	double hcp1 = cp1;
	double hcp2 = cp2;
	double hV = 0.0;
	double hw = wealth;

	// loop from time to nPeriods - 1
	for (int t = time; t < nPeriods; t++)
	{
		// get consumptions for the period t by looking in the consumption matrices	
		if (t > time)
		{
			hcp1 = Interpolate(hw, t,CMP1);
			hcp2 = Interpolate(hw, t, CMP2);
		}

		assert(hcp1 >= 0.0);
		assert(hcp2 >= 0.0);
		assert(hcp1 + hcp2 <= hw);

		// get the instantaneous utility of period t, discounted to period time, and add it to the objective holder
		hV += Weight(wealth, GammaP1, LambdaP1) * DiscountFunction(t - time, DeltaP1, BetaP1, wealth) * UtilityFunction(hcp1, hcp2, hw, GammaP1, alphaP1);
		hV += Weight(wealth, GammaP2, LambdaP2) * DiscountFunction(t - time, DeltaP2, BetaP2, wealth) * UtilityFunction(hcp2, hcp1, hw, GammaP2, alphaP2);

		// set the wealth holder to the wealth of the next period
		if (t < nPeriods - 1)
		{
			hw = Dynamics(hw, hcp1, hcp2, 0);
			assert(hw >= 0.0);
		}
		// get the value of wealth after consumption but before growth
		if (t == nPeriods - 1)
			hw = hw - hcp1 - hcp2;
	}
	Objective = hV;
	// check that the terminal condition is satisfied for the strategy (c1,c2)
	if (TerminalCondition(hw))
		return Objective;
	return -std::numeric_limits<double>::max();
}

// get objective for player i given, time, wealth and specific consumptions, given that future consumptions are indirectly determined by current consumptions
double GetObjectiveNoncooperation(int time, int wealth, int cpi, int cpj, 
	MatrixInt& CMPi, MatrixInt& CMPj, int MaxWealth,
	double GammaPi, double LambdaPi, double DeltaPi, double BetaPi, double alphaPi)
{
	// assertions
	assert(time >= 0);
	assert(time < nPeriods);
	assert(wealth >= 0);
	assert(wealth <= MaxWealth);
	assert(cpi >= 0);
	assert(cpj >= 0);
	assert(cpi + cpj <= (int)(wealth * (1.0 + Growth)));
	assert(CMPi.rows >= MaxWealth + 1);
	assert(CMPi.columns == nPeriods);
	assert(CMPj.rows >= MaxWealth + 1);
	assert(CMPj.columns == nPeriods);

	double Objective = 0.0;

	// holders
	double hcpi = cpi;
	double hcpj = cpj;
	double hV = 0.0;
	double hw = wealth;

	// loop from time to nPeriods - 1
	for (int t = time; t < nPeriods; t++)
	{
		// get consumptions for the period t by looking in the consumption matrices	
		if (t > time)
		{
			hcpi = Interpolate(hw, t, CMPi);
			hcpj = Interpolate(hw, t, CMPj);
		}

		assert(hcpi >= 0.0);
		assert(hcpj >= 0.0);
		assert(hcpi + hcpj <= hw);

		// get the instantaneous utility of period t for player i, discounted to period time, and add it to the objective holder
		hV += Weight(wealth, GammaPi, LambdaPi) * DiscountFunction(t - time, DeltaPi, BetaPi, wealth) 
			* UtilityFunction(hcpi, hcpj, hw, GammaPi, alphaPi);

		// set the wealth holder to the wealth of the next period
		if (t < nPeriods - 1)
		{
			hw = Dynamics(hw, hcpi, hcpj, 0);
			assert(hw >= 0.0);
		}
		// get the value of wealth after consumption but before growth
		if (t == nPeriods - 1)
			hw = hw - hcpi - hcpj;
	}
	Objective = hV;
	// check that the terminal condition is satisfied for the strategy (c1,c2)
	if (TerminalCondition(hw))
		return Objective;
	return -std::numeric_limits<double>::max();
}

// get optimal consumption vector given time and wealth
std::vector<int> GetOptimalConsumptionCooperation(int time, int wealth, MatrixInt& CMP1, MatrixInt& CMP2, int MaxWealth)
{
	// assertions
	assert(time >= 0);
	assert(time < nPeriods);
	assert(wealth >= 0);
	assert(wealth <= MaxWealth);
	assert(CMP1.rows >= MaxWealth + 1);
	assert(CMP1.columns == nPeriods);
	assert(CMP2.rows >= MaxWealth + 1);
	assert(CMP2.columns == nPeriods);

	// initialize holders
	double V_old = -std::numeric_limits<double>::max();
	double V_new = -std::numeric_limits<double>::max();
	std::vector<int> Consumption(2,0);

	// loop over all consumptions from 0 to current wealth and retrieve the one
	// that yields the highest objective
	for (int c1 = 0; c1 <= wealth; c1++)
	{
		for (int c2 = 0; c2 <= wealth; c2++)
		{
			int wealthlimit = wealth;
			if (c1 + c2 <= wealthlimit)
			{
				// get objective for that time, wealth and consumptions, given that the consumptions
				// of later periods are indirectly determined by current consumptions
				V_new = GetObjectiveCooperationStochastic(time, wealth, c1, c2, CMP1, CMP2, MaxWealth);
				// if the objective is higher then the current holder value, we update the holder and the consumptions
				if (V_new > V_old)
				{
					V_old = V_new;
					Consumption[0] = c1;
					Consumption[1] = c2;
				}
			}
		}
	}
	std::cout << "W: " << wealth << ", T: " << time << "\n";
	return Consumption;
}

// get optimal consumption under noncooperation
std::vector<int> GetOptimalConsumptionNonCooperation(int time, int wealth, MatrixInt& CMP1, MatrixInt& CMP2, int MaxWealth)
{
	// assertions
	assert(time >= 0);
	assert(time < nPeriods);
	assert(wealth >= 0);
	assert(wealth <= MaxWealth);
	assert(CMP1.rows >= MaxWealth + 1);
	assert(CMP1.columns == nPeriods);
	assert(CMP2.rows >= MaxWealth + 1);
	assert(CMP2.columns == nPeriods);

	// initialize holders
	double V_old = -std::numeric_limits<double>::max();
	double V_new = -std::numeric_limits<double>::max();
	
	// initialize matrices to hold the response strategies for each player
	MatrixInt ResponseMatrixP1(2, wealth + 1, -1);
	MatrixInt ResponseMatrixP2(2, wealth + 1, -1);

	// find response strategies for player 1
	for (int c2 = 0; c2 <= wealth; c2++)
	{
		for (int c1 = 0; c1 <= wealth; c1++)
		{
			int wealthlimit = wealth;
			if (c1 + c2 <= wealthlimit)
			{
				// get objective for player 1 at that time and response, given wealth and the consumption of player 2
				// and that the consumption of later periods are indirectly determined by current consumptions
				V_new = GetObjectiveNoncooperation(time, wealth, c1, c2, CMP1, CMP2, MaxWealth,
					GammaP1, LambdaP1, DeltaP1, BetaP1, alphaP1);
				// if the objective is higher then the current holder value, we update the holder and the consumptions
				if (V_new > V_old)
				{
					V_old = V_new;
					ResponseMatrixP1.Set(c2, 0, c1);
					ResponseMatrixP1.Set(c2, 1, c2);
				}
			}
		}
	}

	// find response strategies for player 2
	for (int c1 = 0; c1 <= wealth; c1++)
	{
		for (int c2 = 0; c2 <= wealth; c2++)
		{
			int wealthlimit = wealth;
			if (c1 + c2 <= wealthlimit)
			{
				// get objective for player 2 at that time and response, given wealth and the consumption of player 1
				// and that the consumption of later periods are indirectly determined by current consumptions
				V_new = GetObjectiveNoncooperation(time, wealth, c2, c1, CMP1, CMP2, MaxWealth,
					GammaP2, LambdaP2, DeltaP2, BetaP2, alphaP2);
				// if the objective is higher then the current holder value, we update the holder and the consumptions
				if (V_new > V_old)
				{
					V_old = V_new;
					ResponseMatrixP2.Set(c1, 0, c2);
					ResponseMatrixP2.Set(c1, 1, c1);
				}
			}
		}
	}

	// create holder for solution and distance
	std::vector<int> SolutionIndex(2, 0);
	double Diff_old = std::numeric_limits<double>::max();
	double Diff_new = std::numeric_limits<double>::max();

	// find the intersection of the response functions
	for (int i = 0; i <= wealth; i++)
	{
		for (int j = 0; j <= wealth; j++)
		{
			int wealthlimit = wealth;
			double c1_diff = std::pow(ResponseMatrixP1.Get(i, 0) - ResponseMatrixP2.Get(j, 1), 2);
			double c2_diff = std::pow(ResponseMatrixP1.Get(i, 1) - ResponseMatrixP2.Get(j, 0), 2);
			Diff_new = c1_diff + c2_diff;
			if (Diff_new < Diff_old)
			{
				SolutionIndex[0] = i;
				SolutionIndex[1] = j;
			}
		}
	}	
	std::vector<int> Solution(2, 0);
	Solution[0] = ResponseMatrixP1.Get(SolutionIndex[0], 0);
	Solution[1] = ResponseMatrixP1.Get(SolutionIndex[0], 1);
	// return solution
	return Solution;
}

// get the maximum attainable wealth for each period
std::vector<int> GetMaxByPeriod()
{
	std::vector<int> MaxWealth;
	double wealth = InitialWealth;
	for (int i = 0; i < nPeriods; i++)
	{
		MaxWealth.push_back((int)(wealth+1.0));
		wealth = Dynamics(wealth, 0, 0, 1); // change in stochastic case
	}
	assert(MaxWealth.size() == nPeriods);
	return MaxWealth;
}

int main()
{
	assert(1.0 + Growth - Sigma >= 0.0);

	/// CREATE SHOCK VECTOR ///

	std::random_device rd;
	std::mt19937 rng(rd());

	std::vector<int> Shocks;

	for (int i = 0; i < nPeriods; i++)
	{
		std::uniform_int_distribution<int> uni(-1, 0);
		int shock = uni(rng);
		if (shock == 0)
			shock = 1;
		Shocks.push_back(shock);
	}

	/// END CREATE SHOCK VECTOR ///



	/// CREATE CONTAINERS ///

	// get the maximum attainable wealth for each period
	std::vector<int> MaxWealthVector = GetMaxByPeriod();
	int FinalMaxWealth = MaxWealthVector[nPeriods - 1];
	int RoundingError = 0;
	int rounding_error_temp = 0;

	// create consumption matrices
	MatrixInt ConsumptionMatrixP1(FinalMaxWealth + 1, nPeriods, -1);
	MatrixInt ConsumptionMatrixP2(FinalMaxWealth + 1, nPeriods, -1);

	// create consumption matrices with fraction of wealth
	MatrixInt ConsumptionMatrixFractionsP1(FinalMaxWealth + 1, nPeriods, -1);
	MatrixInt ConsumptionMatrixFractionsP2(FinalMaxWealth + 1, nPeriods, -1);
	
	// create optimal path vectors
	std::vector<int> OptimalStatePath(1, InitialWealth);
	std::vector<int> OptimalConsumptionPathP1;
	std::vector<int> OptimalConsumptionPathP2;

	// create holders
	std::vector<int> Consumption(2,0);
	int index(0);
	int TotalConsumptionP1(0);
	int TotalConsumptionP2(0);
	int DiscountedConsumptionP1(0);
	int DiscountedConsumptionP2(0);

	/// END CREATE CONTAINERS ///



	/// FILL CONTAINERS ///

	// loop over the periods and wealth levels and get the optimal consumption for each
	// time period, wealth level and player
	for (int i = nPeriods - 1; i >= 0; i--) // nPeriods time periods
	{
		for (int j = 0; j <= MaxWealthVector[i]; j++)
		{
			// get consumptions for the two players for the given time period and wealth level
			Consumption = GetOptimalConsumptionCooperation(i, j, ConsumptionMatrixP1, ConsumptionMatrixP2, MaxWealthVector[i]);
			
			assert(Consumption[0] >= 0);
			assert(Consumption[1] >= 0);
			assert(Consumption[0] + Consumption[1] <= j);

			// add to the matrices
			ConsumptionMatrixP1.Set(j, i, Consumption[0]);
			ConsumptionMatrixP2.Set(j, i, Consumption[1]);

			if (j >= 1)
			{
				rounding_error_temp = ConsumptionMatrixP1.Get(j - 1, i) - ConsumptionMatrixP1.Get(j, i);
				if (rounding_error_temp > RoundingError)
					RoundingError = rounding_error_temp;
			}
			if (j >= 1)
			{
				rounding_error_temp = ConsumptionMatrixP2.Get(j - 1, i) - ConsumptionMatrixP2.Get(j, i);
				if (rounding_error_temp > RoundingError)
					RoundingError = rounding_error_temp;
			}
		}
	}

	// create the optimal paths by pulling data from the consumption matrices
	for (int i = 0; i < nPeriods; i++)
	{
		index = OptimalStatePath[i]; // get the index for the consumption matrices

		assert(index <= MaxWealthVector[i]);

		Consumption[0] = ConsumptionMatrixP1.Get(index, i); 
		Consumption[1] = ConsumptionMatrixP2.Get(index, i);
		
		// get the consumptions

		assert(Consumption[0] >= 0);
		assert(Consumption[1] >= 0);
		assert(Consumption[0] + Consumption[1] <= OptimalStatePath[i]);
		
		// place consumptions in the vectors
		OptimalConsumptionPathP1.push_back(Consumption[0]);
		OptimalConsumptionPathP2.push_back(Consumption[1]);

		// add the next period's wealth to the wealth vector (set Z = 0 because we want the excepted change)
		OptimalStatePath.push_back((int)(Dynamics(OptimalStatePath[i], Consumption[0], Consumption[1], 0)+0.5));
	}
	
	assert(OptimalConsumptionPathP1.size() == nPeriods);
	assert(OptimalConsumptionPathP2.size() == nPeriods);
	assert(OptimalStatePath.size() == nPeriods + 1);

	// get total consumptions
	for (int i = 0; i < nPeriods; i++)
	{
		TotalConsumptionP1 += OptimalConsumptionPathP1[i];
		TotalConsumptionP2 += OptimalConsumptionPathP2[i];

		DiscountedConsumptionP1 += (int)(DiscountFunction(i, DeltaP1, BetaP1, InitialWealth) * OptimalConsumptionPathP1[i]);
		DiscountedConsumptionP2 += (int)(DiscountFunction(i, DeltaP2, BetaP2, InitialWealth) * OptimalConsumptionPathP2[i]);
	}

	/// END FILL CONTAINERS ///



	/// OUTPUT DATA ///

	// output consumption matrix for player 1
	std::ofstream ostP1{ "Output/ConsumptionMatrixP1.txt" };
	for (int j = 0; j < nPeriods; j++)
	{
		ostP1 << "T" + std::to_string(j) << ",\t";
	}
	ostP1 << "Stock" << "\n";
	for (int i = 0; i <= FinalMaxWealth; i++)
	{
		for (int j = 0; j < nPeriods; j++)
		{
			ostP1 << ConsumptionMatrixP1.Get(i, j) << ",\t";
		}
		ostP1 << i << "\n";
	}
	ostP1.close();

	// output consumption matrix for player 2
	std::ofstream ostP2{ "Output/ConsumptionMatrixP2.txt" };
	for (int j = 0; j < nPeriods; j++)
	{
		ostP2 << "T" + std::to_string(j) << ",\t";
	}
	ostP2 << "Stock" << "\n";
	for (int i = 0; i <= FinalMaxWealth; i++)
	{
		for (int j = 0; j < nPeriods; j++)
		{
			ostP2 << ConsumptionMatrixP2.Get(i, j) << ",\t";
		}
		ostP2 << i << "\n";
	}
	ostP2.close();

	// output consumption matrix for fractions for player 1
	std::ofstream ostFP1{ "Output/ConsumptionMatrixFractionsP1.txt" };
	for (int j = 0; j < nPeriods; j++)
	{
		ostFP1 << "T" + std::to_string(j) << ",\t";
	}
	ostFP1 << "Stock" << std::fixed << std::setprecision(3) << "\n";
	for (int i = 1; i <= FinalMaxWealth; i++)
	{
		for (int j = 0; j < nPeriods; j++)
		{
			ostFP1 << (double)ConsumptionMatrixP1.Get(i, j) / (double)i << ",\t";
		}
		ostFP1 << i << "\n";
	}
	ostFP1.close();

	// output consumption matrix for fractions for player 2
	std::ofstream ostFP2{ "Output/ConsumptionMatrixFractionsP2.txt" };
	for (int j = 0; j < nPeriods; j++)
	{
		ostFP2 << "T" + std::to_string(j) << ",\t";
	}
	ostFP2 << "Stock" << std::fixed << std::setprecision(3) << "\n";
	for (int i = 1; i <= FinalMaxWealth; i++)
	{
		for (int j = 0; j < nPeriods; j++)
		{
			ostFP2 << (double)ConsumptionMatrixP2.Get(i, j) / (double)i << ",\t";
		}
		ostFP2 << i << "\n";
	}
	ostFP2.close();

	// latex table data for P1
	std::ofstream ostP1Latex{ "Output/ConsumptionMatrixP1Latex.txt" };
	for (int j = 0; j < nPeriods; j++)
	{
		ostP1Latex << "s = " + std::to_string(j) << " & ";
	}
	ostP1Latex << "Stock" << " \\" << "\\" << "\n";
	for (int i = 0; i <= FinalMaxWealth; i+= 20)
	{
		for (int j = 0; j < nPeriods; j++)
		{
			ostP1Latex << ConsumptionMatrixP1.Get(i, j) << " & ";
		}
		ostP1Latex << i << " \\" << "\\" << "\n";;
	}
	ostP1Latex.close();

	// latex table data for P2
	std::ofstream ostP2Latex{ "Output/ConsumptionMatrixP2Latex.txt" };
	for (int j = 0; j < nPeriods; j++)
	{
		ostP2Latex << "s = " + std::to_string(j) << " & ";
	}
	ostP2Latex << "Stock" << " \\" << "\\" << "\n";;
	for (int i = 0; i <= FinalMaxWealth; i+= 20)
	{
		for (int j = 0; j < nPeriods; j++)
		{
			ostP2Latex << ConsumptionMatrixP2.Get(i, j) << " & ";
		}
		ostP2Latex << i << " \\" << "\\" << "\n";;
	}
	ostP2Latex.close();

	// combined latex table
	std::ofstream ostCombinedLatex{ "Output/ConsumptionMatrixCombinedLatex.txt" };
	for (int j = 0; j < nPeriods; j++)
	{
		ostCombinedLatex << "$c_1(" << std::to_string(j) << ")$" << " & " << "$c_2(" << std::to_string(j) << ")$" << " & ";
	}
	ostCombinedLatex << "x(s)" << " \\" << "\\" << "\n";;
	for (int i = 0; i <= FinalMaxWealth; i += 10)
	{
		for (int j = 0; j < nPeriods; j++)
		{
			ostCombinedLatex << ConsumptionMatrixP1.Get(i, j) << " & " << ConsumptionMatrixP2.Get(i, j) << " & ";
		}
		ostCombinedLatex << i << " \\" << "\\" << "\n";;
	}
	ostCombinedLatex.close();

	// latex graph data for P1
	std::ofstream ostP1LatexGraph{ "Output/ConsumptionMatrixP1LatexGraph.txt" };
	for (int j = 0; j < nPeriods; j++)
	{
		for (int i = 0; i <= FinalMaxWealth; i += 20)
		{
			ostP1LatexGraph << "( " << i << " , " << ConsumptionMatrixP1.Get(i, j) << " )\n";
		}
	}
	ostP1LatexGraph.close();

	// latex graph data for P2
	std::ofstream ostP2LatexGraph{ "Output/ConsumptionMatrixP2LatexGraph.txt" };
	for (int j = 0; j < nPeriods; j++)
	{
		for (int i = 0; i <= FinalMaxWealth; i += 20)
		{
			ostP2LatexGraph << "( " << i << " , " << ConsumptionMatrixP2.Get(i, j) << " )\n";
		}
	}
	ostP2LatexGraph.close();


	// output matrix with optimal paths
	std::ofstream ostOpt{ "Output/OptimalPaths.txt" };
	ostOpt << "TC-equilibrium in a 2-player game" << "\n\n";
	ostOpt << "Periods: " << nPeriods << "\n";
	ostOpt << "Initial Stock: " << InitialWealth << "\n";
	ostOpt << "Gamma1: " << GammaP1 << "\n";
	ostOpt << "Gamma2: " << GammaP2 << "\n";
	ostOpt << "Delta1: " << DeltaP1 << "\n";
	ostOpt << "Delta2: " << DeltaP2 << "\n";
	ostOpt << "Beta1: " << BetaP1 << "\n";
	ostOpt << "Beta2: " << BetaP2 << "\n";
	ostOpt << "Lambda1: " << LambdaP1 << "\n";
	ostOpt << "Lambda2: " << LambdaP2 << "\n";
	ostOpt << "Alpha: " << Alpha << "\n";
	ostOpt << "bb: " << bb << "\n";
	ostOpt << "Sigma: " << Sigma << "\n";
	ostOpt << "Quadratic Max P1: " << (a * alphaP2) / (2 * (b + alphaP1) * (b + alphaP2) - 2 * pow(b, 2)) << "\n";
	ostOpt << "Quadratic Max P2: " << (a * alphaP1) / (2 * (b + alphaP1) * (b + alphaP2) - 2 * pow(b, 2)) << "\n";
	ostOpt << "Total Consumption P1: " << TotalConsumptionP1 << "\n";
	ostOpt << "Total Consumption P2: " << TotalConsumptionP2 << "\n";
	ostOpt << "Discounted Consumption P1: " << DiscountedConsumptionP1 << "\n";
	ostOpt << "Discounted Consumption P2: " << DiscountedConsumptionP2 << "\n";
	ostOpt << "Max Rounding Error: " << RoundingError << "\n\n";

	ostOpt << "Time" << "\t\t";
	for (int j = 0; j < nPeriods; j++)
	{
		ostOpt << "T" + std::to_string(j) << ",\t";
	}
	ostOpt << "\nConsumptionP1" << "\t";
	for (int j = 0; j < nPeriods; j++)
	{
		ostOpt << OptimalConsumptionPathP1[j] << ",\t";
	}
	ostOpt << "\nConsumptionP2" << "\t";
	for (int j = 0; j < nPeriods; j++)
	{
		ostOpt << OptimalConsumptionPathP2[j] << ",\t";
	}
	ostOpt << "\nConsumptionT" << "\t";
	for (int j = 0; j < nPeriods; j++)
	{
		ostOpt << OptimalConsumptionPathP1[j] + OptimalConsumptionPathP2[j] << ",\t";
	}
	ostOpt << "\nConsumptionF" << std::fixed << std::setprecision(2) << "\t";
	for (int j = 0; j < nPeriods; j++)
	{
		if (OptimalStatePath[j] == 0)
		{
			ostOpt << (double)0 << "\t";
		}
		else
		{
			ostOpt << ((double)OptimalConsumptionPathP1[j] + (double)OptimalConsumptionPathP2[j]) 
				/ (double)OptimalStatePath[j] << ",\t";
		}
	}
	ostOpt << "\nStock" << "\t\t";
	for (int j = 0; j <= nPeriods; j++)
	{
		ostOpt << OptimalStatePath[j] << ",\t";
	}

	/// END OUTPUT DATA ///



	/// OUTOUT LATEX TABLE

	ostOpt << "\n\nTime";
	for (int j = 0; j < nPeriods; j++)
	{
		ostOpt << " & " << "t" + std::to_string(j);
	}
	ostOpt << "\\" << "\nConsumptionP1";
	for (int j = 0; j < nPeriods; j++)
	{
		ostOpt << " & " << OptimalConsumptionPathP1[j];
	}
	ostOpt << "\\" << "\nConsumptionP2";
	for (int j = 0; j < nPeriods; j++)
	{
		ostOpt << " & " << OptimalConsumptionPathP2[j];
	}
	ostOpt << "\\" << "\nFractionConsumptionP1";
	for (int j = 0; j < nPeriods; j++)
	{
		ostOpt << " & " << (double)OptimalConsumptionPathP1[j] / (double)OptimalStatePath[j];
	}
	ostOpt << "\\" << "\nFractionConsumptionP2";
	for (int j = 0; j < nPeriods; j++)
	{
		ostOpt << " & " << (double)OptimalConsumptionPathP2[j] / (double)OptimalStatePath[j];
	}
	ostOpt << "\\" << "\nConsumptionT";
	for (int j = 0; j < nPeriods; j++)
	{
		ostOpt << " & " << OptimalConsumptionPathP1[j] + OptimalConsumptionPathP2[j];
	}
	ostOpt << "\\" << "\nConsumptionF" << std::fixed << std::setprecision(2);
	for (int j = 0; j < nPeriods; j++)
	{
		if (OptimalStatePath[j] == 0)
		{
			ostOpt << " & " << (double)0;
		}
		else
		{
			ostOpt << " & " << ((double)OptimalConsumptionPathP1[j] + (double)OptimalConsumptionPathP2[j])
				/ (double)OptimalStatePath[j];
		}
	}
	ostOpt << "\\" << "\nStock" << " & ";
	for (int j = 0; j <= nPeriods; j++)
	{
		ostOpt << " & " << OptimalStatePath[j];
	}
	ostOpt << "\\";
	ostOpt.close();
	return 0;
}
