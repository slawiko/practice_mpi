#include <iostream>
#include <cmath>
#include <vector>
#include <mpi.h>

using namespace std;

const int OFFSET = 3;

int n1 = 41;
int n2 = 41;
int n3 = 41;
int tmax = 100;
int r = 1;
int blockSize1, blockSize2, buffSize;

const double h1 = 0.01;
const double h2 = 0.01;
const double h3 = 0.01;
const double tau = 0.01;

const double lambda1 = 1;
const double lambda2 = 3;

void setBorderConditions(
	vector<vector<vector<double>>> &a0,
	vector<vector<vector<double>>> &a1,
	vector<vector<vector<double>>> &b0,
	vector<vector<vector<double>>> &b1,
	vector<vector<vector<double>>> &c0,
	vector<vector<vector<double>>> &c1
) {
	for (int i2 = 0; i2 < n2; ++i2) {
		for (int i3 = 0; i3 < n3; ++i3) {
			for (int j = 0; j < tmax; ++j) {
				a0[i2][i3][j] = exp(lambda1 * (i2 * h2 + i3 * h3) + lambda2 * (j + (1. / 3.)) * tau);
				a1[i2][i3][j] = exp(lambda1 * ((n1 - 1) * h1 + i2 * h2 + i3 * h3) + lambda2 * (j + (1. / 3.)) * tau);
			}
		}
	}

	for (int i1 = 0; i1 < n1; ++i1) {
		for (int i3 = 0; i3 < n3; ++i3) {
			for (int j = 0; j < tmax; ++j) {
				b0[i1][i3][j] = exp(lambda1 * (i1 * h1 + i3 * h3) + lambda2 * (j + (2. / 3.)) * tau);
				b1[i1][i3][j] = exp(lambda1 * (i1 * h1 + (n2 - 1) * h2 + i3 * h3) + lambda2 * (j + (2. / 3.)) * tau);
			}
		}
	}

	for (int i1 = 0; i1 < n1; ++i1) {
		for (int i2 = 0; i2 < n2; ++i2) {
			for (int j = 0; j < tmax; ++j) {
				c0[i1][i2][j] = exp(lambda1 * (i1 * h1 + i2 * h2) + lambda2 * (j + 1) * tau);
				c1[i1][i2][j] = exp(lambda1 * (i1 * h1 + i2 * h2 + (n3 - 1) * h3) + lambda2 * (j + 1) * tau);
			}
		}
	}
}

void setInitialApproximation(vector<vector<vector<vector<double>>>> &y) {
	for (int i1 = 0; i1 < n1; ++i1) {
		for (int i2 = 0; i2 < n2; ++i2) {
			for (int i3 = 0; i3 < n3; ++i3) {
				y[i1][i2][i3][0] = exp(lambda1 * (i1 * h1 + i2 * h2 + i3 * h3));
			}
		}
	}
}

double get(double *arr, int m, int l, int i, int j, int k) {
	return arr[i * m * l + j * l + k];
}

void set(double *arr, int m, int l, int i, int j, int k, double val) {
	arr[i * m * l + j * l + k] = val;
}

void swapPointers(double*& a, double*& b) {
	double* tmp;
	tmp = a;
	a = b;
	b = tmp;
}

void executeMaster() {
	MPI_Status status;
	MPI_Request sendRequest;
	MPI_Request receiveRequest;

	vector<vector<vector<vector<double>>>> y =
		vector<vector<vector<vector<double>>>>(n1, vector<vector<vector<double>>>(n2, vector<vector<double>>(n3, vector<double>(tmax, 0.))));

	vector<vector<vector<double>>> tempY = vector<vector<vector<double>>>(n1, vector<vector<double>>(n2, vector<double>(n3)));

	setInitialApproximation(y);

	double *input1 = new double[buffSize];
	double *input2 = new double[buffSize];
	double *output1 = new double[buffSize];
	double *output2 = new double[buffSize];
	double *inArr, *outArr;
	double block, difference, maxDifference;

	for (int t = 0; t < tmax - 1; t++) {
		output1[0] = t;
		output2[0] = t;
		output1[1] = 1;
		output2[1] = 1;
		MPI_Irecv(NULL, 0, MPI_DOUBLE, MPI_PROC_NULL, 0, MPI_COMM_WORLD, &receiveRequest);
		MPI_Isend(NULL, 0, MPI_DOUBLE, MPI_PROC_NULL, 0, MPI_COMM_WORLD, &sendRequest);
		for (int i2 = 0; i2 < r + 1; i2++) {
			if (i2 < r) {
				output2[2] = i2;
				outArr = output2 + OFFSET;
				for (int i = 0; i < n1; i++) {
					for (int j = 0; j < blockSize2 && i2 * blockSize2 + j < n2; j++) {
						for (int k = 0; k < n3; k++) {
							set(outArr, blockSize2, n3, i, j, k, y[i][i2 * blockSize2 + j][k][t]);
						}
					}
				}
				MPI_Wait(&sendRequest, &status);
				swapPointers(output1, output2);
				MPI_Recv(NULL, 0, MPI_DOUBLE, MPI_ANY_SOURCE, 1, MPI_COMM_WORLD, &status);
				MPI_Isend(output1, OFFSET + n1 * blockSize2 * n3, MPI_DOUBLE, status.MPI_SOURCE, 0, MPI_COMM_WORLD, &sendRequest);
			}
			MPI_Wait(&receiveRequest, &status);
			swapPointers(input1, input2);
			if (i2 < r) {
				MPI_Irecv(input1, OFFSET + n1 * blockSize2 * n3, MPI_DOUBLE, MPI_ANY_SOURCE, 0, MPI_COMM_WORLD, &receiveRequest);
			}
			if (i2 > 0) {
				block = input2[2];
				inArr = input2 + OFFSET;
				for (int i = 0; i < n1; i++) {
					for (int j = 0; j < blockSize2 && block * blockSize2 + j < n2; j++) {
						for (int k = 0; k < n3; k++) {
							tempY[i][block * blockSize2 + j][k] = get(inArr, blockSize2, n3, i, j, k);
						}
					}
				}
			}
		}

		output1[1] = 2;
		output2[1] = 2;
		MPI_Irecv(NULL, 0, MPI_DOUBLE, MPI_PROC_NULL, 0, MPI_COMM_WORLD, &receiveRequest);
		MPI_Isend(NULL, 0, MPI_DOUBLE, MPI_PROC_NULL, 0, MPI_COMM_WORLD, &sendRequest);
		for (int i1 = 0; i1 < r + 1; i1++) {
			if (i1 < r) {
				output2[2] = i1;
				outArr = output2 + OFFSET;
				for (int i = 0; i < blockSize1 && i1 * blockSize1 + i < n1; i++) {
					for (int j = 0; j < n2; j++) {
						for (int k = 0; k < n3; k++) {
							set(outArr, n2, n3, i, j, k, tempY[i1 * blockSize1 + i][j][k]);
						}
					}
				}
				MPI_Wait(&sendRequest, &status);
				swapPointers(output1, output2);
				MPI_Recv(NULL, 0, MPI_DOUBLE, MPI_ANY_SOURCE, 1, MPI_COMM_WORLD, &status);
				MPI_Isend(output1, OFFSET + blockSize1 * n2 * n3, MPI_DOUBLE, status.MPI_SOURCE, 0, MPI_COMM_WORLD, &sendRequest);
			}
			MPI_Wait(&receiveRequest, &status);
			swapPointers(input1, input2);
			if (i1 < r) {
				MPI_Irecv(input1, OFFSET + blockSize1 * n2 * n3, MPI_DOUBLE, MPI_ANY_SOURCE, 0, MPI_COMM_WORLD, &receiveRequest);
			}
			if (i1 > 0) {
				block = input2[2];
				inArr = input2 + OFFSET;
				for (int i = 0; i < blockSize1 && block * blockSize1 + i < n1; i++) {
					for (int j = 0; j < n2; j++) {
						for (int k = 0; k < n3; k++) {
							y[block * blockSize1 + i][j][k][t + 1] = get(inArr, n2, n3, i, j, k);
						}
					}
				}
			}
		}

		maxDifference = 0;
		for (int i1 = 0; i1 < n1; i1++) {
			for (int i2 = 0; i2 < n2; i2++) {
				for (int i3 = 0; i3 < n3; i3++) {
					difference = fabs(exp(lambda1 * (i1 * h1 + i2 * h2 + i3 * h3) + lambda2 * (t + 1) * tau) - y[i1][i2][i3][t + 1]);
					if (difference > maxDifference) {
						maxDifference = difference;
					}
				}
			}
		}

		std::cout << maxDifference << std::endl;
	}

	delete[] input1;
	delete[] input2;
	delete[] output1;
	delete[] output2;
}

void executeWorker() {
	MPI_Status status;
	MPI_Request sendRequest;
	MPI_Request receiveRequest;

	vector<vector<vector<double>>> a0 = vector<vector<vector<double>>>(n2, vector<vector<double>>(n3, vector<double>(tmax)));
	vector<vector<vector<double>>> a1 = vector<vector<vector<double>>>(n2, vector<vector<double>>(n3, vector<double>(tmax)));
	vector<vector<vector<double>>> b0 = vector<vector<vector<double>>>(n1, vector<vector<double>>(n3, vector<double>(tmax)));
	vector<vector<vector<double>>> b1 = vector<vector<vector<double>>>(n1, vector<vector<double>>(n3, vector<double>(tmax)));
	vector<vector<vector<double>>> c0 = vector<vector<vector<double>>>(n1, vector<vector<double>>(n2, vector<double>(tmax)));
	vector<vector<vector<double>>> c1 = vector<vector<vector<double>>>(n1, vector<vector<double>>(n2, vector<double>(tmax)));

	setBorderConditions(a0, a1, b0, b1, c0, c1);

	double epsilon1 = 2 * h1 * h1 / tau;
	double epsilon2 = 2 * h2 * h2 / tau;
	double epsilon3 = 2 * h3 * h3 / tau;

	double *alpha1 = new double[n1];
	double *beta1 = new double[n1];
	double *alpha2 = new double[n2];
	double *beta2 = new double[n2];
	double *alpha3 = new double[n3];
	double *beta3 = new double[n3];

	double *input1 = new double[buffSize];
	double *input2 = new double[buffSize];
	double *output1 = new double[buffSize];
	double *output2 = new double[buffSize];
	double *inArr, *outArr;
	double t, task, block, messageSize;

	vector<vector<vector<double>>> tempY = vector<vector<vector<double>>>(blockSize1, vector<vector<double>>(n2, vector<double>(n3)));

	MPI_Send(NULL, 0, MPI_DOUBLE, 0, 1, MPI_COMM_WORLD);
	MPI_Irecv(input1, buffSize, MPI_DOUBLE, 0, MPI_ANY_TAG, MPI_COMM_WORLD, &receiveRequest);
	MPI_Isend(NULL, 0, MPI_DOUBLE, MPI_PROC_NULL, 0, MPI_COMM_WORLD, &sendRequest);

	while (true) {
		MPI_Wait(&receiveRequest, &status);
		if (status.MPI_TAG == 1) {
			break;
		}
		swapPointers(input1, input2);
		MPI_Send(NULL, 0, MPI_DOUBLE, 0, 1, MPI_COMM_WORLD);
		MPI_Irecv(input1, buffSize, MPI_DOUBLE, 0, MPI_ANY_TAG, MPI_COMM_WORLD, &receiveRequest);

		t = input2[0];
		task = input2[1];
		block = input2[2];
		inArr = input2 + OFFSET;
		output2[0] = t;
		output2[1] = task;
		output2[2] = block;
		outArr = output2 + OFFSET;

		if (task == 1) {
			for (int i2 = 0; i2 < blockSize2 && block * blockSize2 + i2 < n2; i2++) {
				for (int i3 = 0; i3 < n3; i3++) {
					alpha1[0] = 0;
					beta1[0] = a0[block * blockSize2 + i2][i3][t];
					for (int i = 1; i < n1 - 1; i++) {
						alpha1[i] = 1 / (2 + epsilon1 - alpha1[i - 1]);
						beta1[i] =
							((get(inArr, blockSize2, n3, i + 1, i2, i3) + get(inArr, blockSize2, n3, i - 1, i2, i3) + beta1[i - 1]) +
							(epsilon1 - 2) * get(inArr, blockSize2, n3, i, i2, i3)) /
								(2 + epsilon1 - alpha1[i - 1]);
					}
					set(outArr, blockSize2, n3, n1 - 1, i2, i3, a1[block * blockSize2 + i2][i3][t]);
					for (int i = n1 - 2; i >= 0; i--) {
						set(outArr, blockSize2, n3, i, i2, i3, alpha1[i] * get(outArr, blockSize2, n3, i + 1, i2, i3) + beta1[i]);
					}
				}
			}
		} else {
			for (int i1 = 0; i1 < blockSize1 && block * blockSize1 + i1 < n1; i1++) {
				for (int i3 = 0; i3 < n3; i3++) {
					alpha2[0] = 0;
					beta2[0] = b0[block * blockSize1 + i1][i3][t];
					for (int i = 1; i < n2 - 1; i++) {
						alpha2[i] = 1 / (2 + epsilon2 - alpha2[i - 1]);
						beta2[i] =
							((get(inArr, n2, n3, i1, i + 1, i3) + get(inArr, n2, n3, i1, i - 1, i3) + beta2[i - 1]) +
							(epsilon2 - 2) * get(inArr, n2, n3, i1, i, i3)) /
								(2 + epsilon2 - alpha2[i - 1]);
					}
					tempY[i1][n2 - 1][i3] = b1[block * blockSize1 + i1][i3][t];
					for (int i = n2 - 2; i >= 0; i--) {
						tempY[i1][i][i3] =
							alpha2[i] * tempY[i1][i + 1][i3] + beta2[i];
					}
				}
			}

			for (int i1 = 0; i1 < blockSize1 && block * blockSize1 + i1 < n1; i1++) {
				for (int i2 = 0; i2 < n2; ++i2) {
					alpha3[0] = 0;
					beta3[0] = c0[block * blockSize1 + i1][i2][t];
					for (int i = 1; i < n3 - 1; i++) {
						alpha3[i] = 1 / (2 + epsilon3 - alpha3[i - 1]);
						beta3[i] =
							((tempY[i1][i2][i + 1] + tempY[i1][i2][i - 1] + beta3[i - 1]) +
							(epsilon3 - 2) * tempY[i1][i2][i]) /
								(2 + epsilon3 - alpha3[i - 1]);
					}
					set(outArr, n2, n3, i1, i2, n3 - 1, c1[block * blockSize1 + i1][i2][t]);
					for (int i = n3 - 2; i >= 0; i--) {
						set(outArr, n2, n3, i1, i2, i, alpha3[i] * get(outArr, n2, n3, i1, i2, i + 1) + beta3[i]);
					}
				}
			}
		}

		MPI_Wait(&sendRequest, &status);
		swapPointers(output1, output2);
		messageSize = (task == 1) ? (OFFSET + n1 * blockSize2 * n3) : (OFFSET + blockSize1 * n2 * n3);
		MPI_Isend(output1, messageSize, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD, &sendRequest);
	}

	delete[] input1;
	delete[] input2;
	delete[] output1;
	delete[] output2;

	delete[] alpha1;
	delete[] beta1;
	delete[] alpha2;
	delete[] beta2;
	delete[] alpha3;
	delete[] beta3;
}

int main(int argc, char* argv[]) {
	for (int i = 1; i < argc - 1; i++) {
		if (strcmp(argv[i], "-n1") == 0) {
			n1 = atoi(argv[i + 1]);
		}
		if (strcmp(argv[i], "-n2") == 0) {
			n2 = atoi(argv[i + 1]);
		}
		if (strcmp(argv[i], "-n3") == 0) {
			n3 = atoi(argv[i + 1]);
		}
		if (strcmp(argv[i], "-t") == 0) {
			tmax = atoi(argv[i + 1]) + 1;
		}
		if (strcmp(argv[i], "-r") == 0) {
			r = atoi(argv[i + 1]);
		}
	}
	
	blockSize1 = n1 / r + (n1 % r ? 1 : 0);
	blockSize2 = n2 / r + (n2 % r ? 1 : 0);
	buffSize = OFFSET + n1 * n2 * n3;

	MPI_Init(&argc, &argv);
	int rank, processCount;
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	MPI_Comm_size(MPI_COMM_WORLD, &processCount);

	if (rank == 0) {
		executeMaster();
		for (int i = 1; i < processCount; i++) {
			MPI_Send(NULL, 0, MPI_DOUBLE, i, 1, MPI_COMM_WORLD);
		}
	} else {
		executeWorker();
	}

	MPI_Finalize();
	return 0;
}
