#include <iostream>
#include <cmath>
#include <vector>

using namespace std;

const int n1 = 41;
const int n2 = 41;
const int n3 = 41;
const int tmax = 100;

const double h1 = 0.01;
const double h2 = 0.01;
const double h3 = 0.01;
const double tau = 0.01;

vector<vector<vector<vector<double> > > > y =
	vector<vector<vector<vector<double> > > >(n1,
		vector<vector<vector<double> > >(n2,
			vector<vector<double> >(n3,
				vector<double>(tmax, 0.))));

vector<vector<vector<double> > > a0 =
	vector<vector<vector<double> > >(n2,
		vector<vector<double> >(n3,
			vector<double>(tmax)));

vector<vector<vector<double> > > a1 =
	vector<vector<vector<double> > >(n2,
		vector<vector<double> >(n3,
			vector<double>(tmax)));

vector<vector<vector<double> > > b0 =
	vector<vector<vector<double> > >(n1,
		vector<vector<double> >(n3,
			vector<double>(tmax)));

vector<vector<vector<double> > > b1 =
	vector<vector<vector<double> > >(n1,
		vector<vector<double> >(n3,
			vector<double>(tmax)));

vector<vector<vector<double> > > c0 =
	vector<vector<vector<double> > >(n1,
		vector<vector<double> >(n2,
			vector<double>(tmax)));

vector<vector<vector<double> > > c1 =
	vector<vector<vector<double> > >(n1,
		vector<vector<double> >(n2,
			vector<double>(tmax)));

const double lambda1 = 1;
const double lambda2 = 3;

int i, i1, i2, i3, j;

// u = e ^ (lambda1 * (x1 + x2 + x3) + lambda2 * t)

void setBorderConditions()
{
	for (i2 = 0; i2 < n2; ++i2)
	{
		for (i3 = 0; i3 < n3; ++i3)
		{
			for (j = 0; j < tmax; ++j)
			{
				a0[i2][i3][j] =
					exp(lambda1 * (i2 * h2 + i3 * h3) + lambda2 * (j + (1. / 3.)) * tau);
				a1[i2][i3][j] =
					exp(lambda1 * ((n1 - 1) * h1 + i2 * h2 + i3 * h3) + lambda2	* (j + (1. / 3.)) * tau);
			}
		}
	}

	for (i1 = 0; i1 < n1; ++i1)
	{
		for (i3 = 0; i3 < n3; ++i3)
		{
			for (j = 0; j < tmax; ++j)
			{
				b0[i1][i3][j] =
					exp(lambda1 * (i1 * h1 + i3 * h3) + lambda2 * (j + (2. / 3.)) * tau);
				b1[i1][i3][j] =
					exp(lambda1 * (i1 * h1 + (n2 - 1) * h2 + i3 * h3) + lambda2 * (j + (2. / 3.)) * tau);
			}
		}
	}

	for (i1 = 0; i1 < n1; ++i1)
	{
		for (i2 = 0; i2 < n2; ++i2)
		{
			for (j = 0; j < tmax; ++j)
			{
				c0[i1][i2][j] =
					exp(lambda1 * (i1 * h1 + i2 * h2) + lambda2 * (j + 1) * tau);
				c1[i1][i2][j] =
					exp(lambda1 * (i1 * h1 + i2 * h2 + (n3 - 1) * h3) + lambda2 * (j + 1) * tau);
			}
		}
	}
}

void setInitialApproximation()
{
	for (i1 = 0; i1 < n1; ++i1)
	{
		for (i2 = 0; i2 < n2; ++i2)
		{
			for (i3 = 0; i3 < n3; ++i3)
			{
				y[i1][i2][i3][0] = exp(lambda1 * (i1 * h1 + i2 * h2 + i3 * h3));
			}
		}
	}
}

void main()
{
	setBorderConditions();

	setInitialApproximation();

	double epsilon1 = 2 * h1 * h1 / tau;
	double epsilon2 = 2 * h2 * h2 / tau;
	double epsilon3 = 2 * h3 * h3 / tau;

	double *alpha;
	double *beta;

	vector<vector<vector<vector<double> > > > tempY =
		vector<vector<vector<vector<double> > > >(n1,
			vector<vector<vector<double> > >(n2,
				vector<vector<double> >(n3,
					vector<double>(2))));

	for (j = 0; j < tmax - 1; ++j)
	{
		alpha = new double[n1];
		beta = new double[n1];
		for (i2 = 0; i2 < n2; ++i2)
		{
			for (i3 = 0; i3 < n3; ++i3)
			{
				alpha[0] = 0;
				beta[0] = a0[i2][i3][j];
				for (i = 1; i < n1 - 1; ++i)
				{
					alpha[i] = 1 / (2 + epsilon1 - alpha[i - 1]);
					beta[i] =
						((y[i + 1][i2][i3][j] + y[i - 1][i2][i3][j] + beta[i - 1]) +
						(epsilon1 - 2) * y[i][i2][i3][j]) /
							(2 + epsilon1 - alpha[i - 1]);
				}
				tempY[n1 - 1][i2][i3][0] = a1[i2][i3][j];
				for (i = n1 - 2; i >= 0; --i)
				{
					tempY[i][i2][i3][0] =
						alpha[i] * tempY[i + 1][i2][i3][0] + beta[i];
				}
			}
		}

		alpha = new double[n2];
		beta = new double[n2];
		for (i1 = 0; i1 < n1; ++i1)
		{
			for (i3 = 0; i3 < n3; ++i3)
			{
				alpha[0] = 0;
				beta[0] = b0[i1][i3][j];
				for (i = 1; i < n2 - 1; ++i)
				{
					alpha[i] = 1 / (2 + epsilon2 - alpha[i - 1]);
					beta[i] =
						((tempY[i1][i + 1][i3][0] + tempY[i1][i - 1][i3][0] + beta[i - 1]) +
						(epsilon2 - 2) * tempY[i1][i][i3][0]) /
							(2 + epsilon2 - alpha[i - 1]);
				}
				tempY[i1][n2 - 1][i3][1] = b1[i1][i3][j];
				for (i = n2 - 2; i >= 0; --i)
				{
					tempY[i1][i][i3][1] =
						alpha[i] * tempY[i1][i + 1][i3][1] + beta[i];
				}
			}
		}

		alpha = new double[n3];
		beta = new double[n3];
		for (i1 = 0; i1 < n1; ++i1)
		{
			for (i2 = 0; i2 < n2; ++i2)
			{
				alpha[0] = 0;
				beta[0] = c0[i1][i2][j];
				for (i = 1; i < n3 - 1; ++i)
				{
					alpha[i] = 1 / (2 + epsilon3 - alpha[i - 1]);
					beta[i] =
						((tempY[i1][i2][i + 1][1] + tempY[i1][i2][i - 1][1] + beta[i - 1]) +
						(epsilon3 - 2) * tempY[i1][i2][i][1]) /
							(2 + epsilon3 - alpha[i - 1]);
				}
				y[i1][i2][n3 - 1][j + 1] = c1[i1][i2][j];
				for (i = n3 - 2; i >= 0; --i)
				{
					y[i1][i2][i][j + 1] =
						alpha[i] * y[i1][i2][i + 1][j + 1] + beta[i];
				}
			}
		}

		double maxDifference = 0;
		double koord[3] = {0, 0, 0};

		for (i1 = 0; i1 < n1; ++i1)
		{
			for (i2 = 0; i2 < n2; ++i2)
			{
				for (i3 = 0; i3 < n3; ++i3)
				{
					if (fabs(exp(lambda1 * (i1 * h1 + i2 * h2 + i3 * h3) + lambda2 * (j + 1) * tau) -
						y[i1][i2][i3][j + 1]) > maxDifference)
					{
						maxDifference = 
							fabs(exp(lambda1 * (i1 * h1 + i2 * h2 + i3 * h3) + lambda2 * (j + 1) * tau) -
							y[i1][i2][i3][j + 1]);
						koord[0] = i1;
						koord[1] = i2;
						koord[2] = i3;
					}
				}
			}
		}

		std::cout << maxDifference << std::endl;
	}

	system("pause");
}
