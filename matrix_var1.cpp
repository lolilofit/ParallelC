#include<mpi.h>
#include<stdio.h>
#include<malloc.h>
#include <stdio.h>
#include <time.h>
#include<mpe.h>


void plus(double* a1, double* a2, double k2, double* res, int size) {
	int i;
	for (i = 0; i < size; i++)
		res[i] = a1[i] + a2[i] * k2;
}


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


double modul(double* a, int size) {
	int i;
	double sum = 0.0;
	for (i = 0; i < size; i++)
		sum += a[i] * a[i];

	return sqrt(sum); //sqrt!
}

double scal_mul(double* a1, double* a2, int size) {
	int i;
	double sum = 0.0;
	for (i = 0; i < size; i++) {
		sum += a1[i] * a2[i];
	}
	return sum;
}

void copy(double *vec, double* vec_cpy, int N) {
	int i;
	for (i = 0; i< N; i++) {
		vec_cpy[i] = vec[i];
	}
}



int main(int argc, char* argv[]) {

	int N, i, j;

	N = 16;
	double* x;
	double* A;
	double* b;

	b = (double*)malloc(sizeof(double)*N);
	x = (double*)malloc(sizeof(double)*N);
	double* r = (double*)malloc(sizeof(double)*N);

	for (i = 0; i < N; i++) {

		x[i] = 0;
		b[i] = N + 1;
		r[i] = 0;
	}


	double epsilon = 1.0 / (10.0*10.0*10.0*10.0*10.0);
	double alpha, beta;
	int proc_num, proc_rank, flag = 0, rest, row_in_tread;
	int evtid_beginPhase1, evtid_endPhase1;

	MPI_Init(&argc, &argv);
	MPI_Comm_size(MPI_COMM_WORLD, &proc_num);
	MPI_Comm_rank(MPI_COMM_WORLD, &proc_rank);

	struct timespec start, end;
	clock_gettime(CLOCK_MONOTONIC_RAW, &start);



	rest = N;

	int *send_num;
	send_num = (int*)malloc(sizeof(int)*proc_num);
	int *send_data_begin;
	send_data_begin = (int*)malloc(sizeof(int)*proc_num);
	int *send = (int*)malloc(sizeof(int)*proc_num);

	//Заполнение массивов с количеством элементов и индексами начала блоков, которые //отдаются каждому процессу

	for (i = 0; i < proc_num; i++)
		send_num[i] = N / proc_num;

	for (i = 0; i<N%proc_num; i++)
		send_num[i] = N / proc_num + 1;
	send_data_begin[0] = 0;

	for (i = 1; i<proc_num; i++) {
		send_data_begin[i] = send_data_begin[i - 1] + send_num[i - 1];

	}
	double *res = (double*)malloc(sizeof(double)*N);
	double *proc_row = (double*)malloc(sizeof(double)*send_num[proc_rank] * N);
	double *proc_res = (double*)malloc(sizeof(double)*send_num[proc_rank]);
	double* prev_r = (double*)malloc(sizeof(double)*N);
	double* z = (double*)malloc(sizeof(double)*N);

	//Умножение Ax

	for (i = 0; i < send_num[proc_rank]; i++) {

		for (j = 0; j < N; j++) {
			if ((send_data_begin[proc_rank] + i) == j) {
				proc_row[i*N + j] = 2.0;


			}
			else
				proc_row[i*N + j] = 1.0;
		}
	}



	for (i = 0; i<send_num[proc_rank]; i++)
		proc_res[i] = 0;

	for (i = 0; i < (send_num[proc_rank] / N); i++) {

		for (j = 0; j<N; j++) {
			proc_res[i] += proc_row[i*N + j] * b[j];
		}
	}

	int *get_num = (int*)malloc(sizeof(int)*proc_num);
	int *get_data_begin = (int*)malloc(sizeof(int)*proc_num);

	get_data_begin[0] = 0;
	for (i = 0; i< proc_num; i++)
		get_num[i] = N / proc_num;
	for (i = 0; i< N%proc_num; i++)
		get_num[i] = N / proc_num + 1;
	for (i = 1; i<proc_num; i++)
		get_data_begin[i] = get_num[i - 1] + get_data_begin[i - 1];

	MPI_Allgatherv(proc_res, get_num[proc_rank], MPI_DOUBLE, r, get_num, get_data_begin, MPI_DOUBLE, MPI_COMM_WORLD);

	plus(b, res, -1, r, N);
	for (i = 0; i<N; i++)
		z[i] = r[i];




	while (epsilon <= modul(r, N) / modul(b, N)) {
		MPE_Init_log();
		evtid_beginPhase1 = MPE_Log_get_event_number();
		evtid_endPhase1 = MPE_Log_get_event_number();
		MPE_Describe_state(evtid_beginPhase1, evtid_endPhase1, "cycle_begin", "red");

		MPE_Log_event(evtid_beginPhase1, proc_rank, (char*)0);

		//Умножение Az

		for (i = 0; i< send_num[proc_rank]; i++)
			proc_res[i] = 0;
		for (i = 0; i < send_num[proc_rank]; i++) {
			for (j = 0; j<N; j++)
				proc_res[i] += proc_row[i*N + j] * z[j];

		}

		MPI_Allgatherv(proc_res, get_num[proc_rank], MPI_DOUBLE, res, get_num, get_data_begin, MPI_DOUBLE, MPI_COMM_WORLD);

		//Вычисление alpha, r(n+1), x(n+1), beta, z(n+1)

		alpha = scal_mul(r, r, N) / scal_mul(res, z, N);
		plus(x, z, alpha, x, N);

		copy(r, prev_r, N);

		plus(r, res, -1 * alpha, r, N);
		beta = scal_mul(r, r, N) / scal_mul(prev_r, prev_r, N);

		plus(r, z, beta, z, N);


		MPE_Log_event(evtid_endPhase1, proc_rank, (char*)0);
		MPE_Finish_log("lab1.clog2");

	}

	if (proc_rank == 0) {

		for (i = 0; i < N; i++)
			printf("%.01f\n", x[i]);
	}


	clock_gettime(CLOCK_MONOTONIC_RAW, &end);
	printf("Time taken: %lf sec.\n", end.tv_sec - start.tv_sec + 0.000000001*(end.tv_nsec - start.tv_nsec));

	free(A);
	free(r);
	free(prev_r);
	free(z);
	free(x);
	free(b);
	free(res);
	free(proc_row);
	free(proc_res);
	free(get_num);
	free(send_num);
	free(get_data_begin);
	free(send_data_begin);

	MPI_Finalize();

	return 0;
}

