#pragma once
#include <cmath>

/// PARAMETERS ///

// time and wealth
static int nPeriods = 5;
static int InitialWealth = 200;
static int TerminalWealthCondition = 0;

// preferences
static double DeltaP1 = 0.4;
static double DeltaP2 = 0.5;
static double BetaP1 = 0.9;
static double BetaP2 = 0.95;
static double GammaP1 = 1.8;
static double GammaP2 = 2.0;

static double a = 30;
static double b = 0.1;
static double alphaP1 = 0.2;
static double alphaP2 = 0.1;


// planner parameters
static double LambdaP1 = 1.0;
static double LambdaP2 = 1.0;

// dynamics parameters
static double Sigma = 0.6;
static double bb = 5.0;
static double Alpha = 0.7;

static double A = 100.0;
static double B = 1.0;
static double Growth = 0.00;
static double Theta = 0.5;

/// END PARAMETERS ///



/// FUNCTIONS ///

double UtilityFunction(double c1, double c2, double w, double gamma, double alpha)
{
	// power
	if (gamma == 1.0)
		return log(c1);
	else
	{
		if (c1 == 0.0)
		{
			if (gamma < 1.0)
				return 0.0;
			else
				return -std::numeric_limits<double>::max();
		}
		else
			return (1.0 / (1.0 - gamma)) * pow(c1, 1.0 - gamma);
	}

	// quadratic
	//return a * c1 - b * c1 * c2 - (b + alpha) * pow(c1, 2);
}
double TerminalFunction(double w, double gamma, double alpha)
{
	// power
	if (w == 0)
		return 0.0;
	return (1.0 / (1.0 - gamma)) * pow(w, 1.0 - gamma);
}
bool TerminalCondition(double w)
{
	// return true if wealth above the minumum value
	if (w >= TerminalWealthCondition)
		return true;
	return false;
}
double DiscountFunction(int time, double Delta, double Beta, int w)
{
	// exponential
	//return std::pow(Delta, time);

	// beta-delta
	if (time == 0)
		return 1.0;
	else
		return Beta * std::pow(Delta, time);
}
double Dynamics(double w, double c1, double c2, int Z)
{
	
	// fish war
	assert(c1 >= 0);
	assert(c2 >= 0);
	if (Z == 0)
	{
		// deterministic case
		// return bb * std::pow(w - c1 - c2, Alpha);
		// stochastic case with expected value
		double Expected = 0.5 * std::pow(1.0 + Sigma, Alpha)
			+ 0.5 * std::pow(1.0 - Sigma, Alpha);
		return bb * Expected * std::pow(w - c1 - c2, Alpha);
	}
	return bb * std::pow((w - c1 - c2) * (1.0 + Sigma * Z), Alpha);
	
	/*
	// cake eating problem
	assert(c1 >= 0);
	assert(c2 >= 0);
	if (Z == 0)
		return w - c1 - c2;
	return (w - c1 - c2) * (1.0 + Sigma * Z);
	*/
}

double Weight(int w, double gamma, double lambda)
{
	// constant weight
	return lambda;
}

/// END FUNCTIONS ///