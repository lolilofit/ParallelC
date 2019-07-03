#include<stdio.h>
#include<malloc.h>
#include<omp.h>
#include<stdlib.h>


double sqrt(const double x) {

	if (x <= 0) {
		return 0;
	}

	int i = 0;
	while ((i * i) <= x) {
		i++;
	}
	i--;
	double z = x - i * i;
	double y = z / (2 * i);
	double w = i + y;
	return w - (y * y) / (2 * w);
}


int main(int argc, char* argv[]) {

	omp_set_num_threads(4);

	int i, j, k, N = 1000;
	double rN = 0.0, bN = 0.0, is_end;

	double procs_s1 = 0.0, procs_s2 = 0.0;
	double* A = (double*)malloc(sizeof(double)*N*N);
	double* b = (double*)malloc(sizeof(double)*N);
	double* x = (double*)malloc(sizeof(double)*N);
	double* r = (double*)malloc(sizeof(double)*N);
	double epsilon = 0.00001;

	double alpha, beta;
	double *buf_z = (double*)malloc(sizeof(double)*N);
	double* buf_res = (double*)malloc(sizeof(double)*N);
	double* buf_prev_r = (double*)malloc(sizeof(double)*N);

	for (i = 0; i < N; ++i) {
		for (j = 0; j < N; ++j) {
			if (i == j)
				A[i*N + j] = 2.0;
			else
				A[i*N + j] = 1.0;
		}
	}
	for (i = 0; i<N; i++) {
		x[i] = 0;
		b[i] = N + 1;
		r[i] = 0;

	}

	for (i = 0; i < N; ++i) {
		for (j = 0; j < N; ++j) {
			r[i] += A[i*N + j] * x[j];
		}
	}

	for (i = 0; i < N; ++i)
		r[i] = b[i] + (-1)*r[i];

	for (i = 0; i<N; ++i)
		buf_z[i] = r[i];

	double sum_one = 0.0, sum_two = 0.0;

	for (i = 0; i<N; ++i) {
		sum_one += r[i] * r[i];
	}
	for (i = 0; i<N; ++i) {
		sum_two += b[i] * b[i];
	}

	procs_s1 += sum_one;
	procs_s2 += sum_two;

	is_end = sqrt(procs_s1) / sqrt(procs_s2);
	procs_s1 = 0.0;
	procs_s2 = 0.0;

	sum_one = 0.0;
	sum_two = 0.0;

	double start = omp_get_wtime();

#pragma omp parallel private(i, j,sum_one, sum_two)
	{
		int for_k, sum_num;
		for (for_k = 0; for_k<3000; for_k++) {
			//while (epsilon <= is_end) {

#pragma omp for
			for (i = 0; i < N; ++i) {
				sum_num = 0;
#pragma omp parallel for reduction(+: sum_num)
				for (j = 0; j < N; ++j) {
					sum_num += A[i*N + j] * buf_z[j];
				}
				buf_res[i] = sum_num;
			}

#pragma omp for
			for (i = 0; i < N; ++i) {
				sum_one += r[i] * r[i];
			}

#pragma omp for
			for (i = 0; i < N; ++i) {
				sum_two += buf_res[i] * buf_z[i];
			}


#pragma omp barrier

#pragma omp reduction(+: procs_s1)
			procs_s1 += sum_one;
#pragma omp reduction(+: procs_s2)
			procs_s2 += sum_two;


#pragma omp barrier
#pragma omp master
			{
				alpha = procs_s1 / procs_s2;
				procs_s1 = 0.0;
				procs_s2 = 0.0;
			}
			sum_one = 0.0;
			sum_two = 0.0;

#pragma omp barrier

#pragma omp for
			for (i = 0; i < N; ++i)
				x[i] = x[i] + buf_z[i] * alpha;

#pragma omp for
			for (i = 0; i< N; ++i)
				buf_prev_r[i] = r[i];

#pragma omp for
			for (i = 0; i < N; ++i)
				r[i] = r[i] + (-1)*alpha*buf_res[i];

#pragma omp barrier

#pragma omp for
			for (i = 0; i < N; ++i)
				sum_one += r[i] * r[i];
#pragma omp for
			for (i = 0; i < N; ++i)
				sum_two += buf_prev_r[i] * buf_prev_r[i];


#pragma omp barrier
#pragma omp reduction(+: procs_s1)
			procs_s1 += sum_one;
#pragma omp reduction(+: procs_s2)
			procs_s2 += sum_two;


#pragma omp barrier
#pragma omp master
			{
				beta = procs_s1 / procs_s2;
				procs_s1 = 0.0;
				procs_s2 = 0.0;
			}
			sum_one = 0.0;
			sum_two = 0.0;

#pragma omp barrier
#pragma omp for
			for (i = 0; i < N; ++i)
				buf_z[i] = r[i] + buf_z[i] * beta;

#pragma omp master
			alpha = 0.0;
#pragma omp master
			beta = 0.0;

#pragma omp for
			for (i = 0; i<N; ++i)
				sum_one += r[i] * r[i];
#pragma omp for
			for (i = 0; i<N; ++i)
				sum_two += b[i] * b[i];


#pragma omp barrier
#pragma omp reduction(+: procs_s1)
			procs_s1 += sum_one;
#pragma omp reduction(+: procs_s2)
			procs_s2 += sum_two;

#pragma omp barrier
#pragma omp master
			{
				is_end = sqrt(procs_s1) / sqrt(procs_s2);
				procs_s1 = 0.0;
				procs_s2 = 0.0;
			}

			sum_one = 0.0;
			sum_two = 0.0;

#pragma omp for
			for (i = 0; i<N; ++i)
				buf_res[i] = 0;

#pragma omp barrier

		}

	}

	double end = omp_get_wtime();

	for (i = 0; i< N; ++i) {
		printf("%f, \n", x[i]);
	}

	printf("Time taken = %.16g\n", end - start);

	free(A);
	free(r);
	free(buf_prev_r);
	free(buf_z);
	free(x);
	free(b);
	free(buf_res);

	return 0;

}
