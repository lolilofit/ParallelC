#include<stdio.h>
#include<malloc.h>
#include<mpi.h>
//#include<mpe.h>
#include <time.h>

const int N = 32768;

//линейна€ операци€ a+k*b
void plus(double* a1, double* a2, double k2, double* res, int size) {
	int i;
	for (i = 0; i < size; i++)
		res[i] = a1[i] + a2[i] * k2;
}

//квадратный корень
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

//скал€рное умножение
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
double modul(double* a, int size) {

	int i;
	double sum = 0.0;
	for (i = 0; i<size; i++)
		sum += a[i] * a[i];
	return sqrt(sum);

}

int main(int argc, char* argv[]) {

	int i, j;
	double* b = (double*)malloc(sizeof(double)*N);
	double* x = (double*)malloc(sizeof(double)*N);
	double* r = (double*)malloc(sizeof(double)*N);

	double epsilon = 1.0 / (10.0*10.0*10.0*10.0*10.0);
	double alpha, beta;
	double *bufA, *bufB, *bufC, *buf_r, *buf_x, *buf_z, *buf_prev_r, *buf_res, *z_copy;
	int block_size, proc_num, proc_rank, rest;
	int evtid_beginPhase1, evtid_endPhase1;


	MPI_Init(&argc, &argv);

	MPI_Comm_size(MPI_COMM_WORLD, &proc_num);
	MPI_Comm_rank(MPI_COMM_WORLD, &proc_rank);

	struct timespec start, end;
	clock_gettime(CLOCK_MONOTONIC_RAW, &start);

	int *send_num, *send_numA;
	send_num = (int*)malloc(sizeof(int)*proc_num);
	send_numA = (int*)malloc(sizeof(int)*proc_num);
	int *send_data_begin, *send_data_beginA;
	send_data_begin = (int*)malloc(sizeof(int)*proc_num);
	send_data_beginA = (int*)malloc(sizeof(int)*proc_num);

	for (i = 0; i < proc_num; i++) {
		send_num[i] = N / proc_num;
		send_numA[i] = N / proc_num * N;
	}
	for (i = 0; i<N%proc_num; i++) {
		send_num[i] = N / proc_num + 1;
		send_numA[i] = N / proc_num * N + 1;
	}
	send_data_begin[0] = 0;
	send_data_beginA[0] = 0;
	for (i = 1; i<proc_num; i++) {
		send_data_begin[i] = send_data_begin[i - 1] + send_num[i - 1];
		send_data_beginA[i] = send_data_beginA[i - 1] + send_numA[i - 1];
	}



	block_size = send_num[proc_rank];
	bufA = (double*)malloc(sizeof(double)*N*block_size);
	bufB = (double*)malloc(sizeof(double)*block_size);
	bufC = (double*)malloc(sizeof(double)*block_size);
	buf_x = (double*)malloc(sizeof(double)*block_size);
	buf_r = (double*)malloc(sizeof(double)*block_size);
	buf_z = (double*)malloc(sizeof(double)*block_size);
	buf_res = (double*)malloc(sizeof(double)*block_size);
	buf_prev_r = (double*)malloc(sizeof(double)*block_size);
	z_copy = (double*)malloc(sizeof(double)*block_size);

	for (i = 0; i< block_size; i++) {
		bufC[i] = 0;
		buf_res[i] = 0;
	}


	//инициализаци€ частей A, b, r, x
	for (i = 0; i < send_num[proc_rank]; i++) {
		for (j = 0; j < N; j++) {
			if ((send_data_begin[proc_rank] + i) == j)
				bufA[i*N + j] = 2.0;
			else
				bufA[i*N + j] = 1.0;
		}

		buf_x[i] = 0;
		buf_r[i] = 0;
		bufB[i] = N + 1;
	}

	//умножение Ax
	double sum = 0.0;

	for (i = 0; i< send_num[proc_rank]; i++) {
		sum = 0.0;
		for (j = 0; j<send_num[proc_rank]; j++)
			sum += bufA[N*i + send_data_begin[proc_rank] + j] * buf_x[j];
		buf_r[i] += sum;

	}

	int next_proc, prev_proc, k;
	for (i = 1; i< proc_num; i++) {
		MPI_Sendrecv_replace(buf_x, N / proc_num + 1, MPI_DOUBLE, (proc_rank + 1) % proc_num, 123, (proc_rank + proc_num - 1) % proc_num, 123, MPI_COMM_WORLD, MPI_STATUS_IGNORE);



		for (k = 0; k< send_num[proc_rank]; k++) {
			sum = 0;
			for (j = 0; j<send_num[proc_rank]; j++)
				sum += bufA[k*N + send_data_begin[proc_rank] + j] * buf_x[j];
			buf_r[k] += sum;

		}
	}

	//вычисление частей r(0), z(0) дл€ каждого потока
	plus(bufB, buf_r, -1, buf_r, send_num[proc_rank]);
	for (i = 0; i<send_num[proc_rank]; i++)
		buf_z[i] = buf_r[i];


	int *get_num;
	get_num = (int*)malloc(sizeof(int)*proc_num);
	int *get_data_begin;
	get_data_begin = (int*)malloc(sizeof(int)*proc_num);

	get_data_begin[0] = 0;
	for (i = 0; i< proc_num; i++)
		get_num[i] = N / proc_num;
	for (i = 0; i< N%proc_num; i++)
		get_num[i] = N / proc_num + 1;
	for (i = 1; i<proc_num; i++)
		get_data_begin[i] = get_num[i - 1] + get_data_begin[i - 1];

	//собираютс€ матрицы r, b, x с каждого потока
	MPI_Allgatherv(buf_r, get_num[proc_rank], MPI_DOUBLE, r, get_num, get_data_begin, MPI_DOUBLE, MPI_COMM_WORLD);
	MPI_Allgatherv(bufB, get_num[proc_rank], MPI_DOUBLE, b, get_num, get_data_begin, MPI_DOUBLE, MPI_COMM_WORLD);
	MPI_Allgatherv(buf_x, get_num[proc_rank], MPI_DOUBLE, x, get_num, get_data_begin, MPI_DOUBLE, MPI_COMM_WORLD);



	MPI_Comm comm_gr;
	int *index, *edges;
	int v, reord = 1;
	index = (int*)malloc(proc_num * sizeof(int));
	edges = (int*)malloc(proc_num*(proc_num - 1) * sizeof(int));

	for (i = 0; i < proc_num; i++)
	{
		index[i] = (proc_num - 1)*(i + 1);
		v = 0;
		for (j = 0; j < proc_num; j++)
		{
			if (i != j)
				edges[i * (proc_num - 1) + v++] = j;
		}
	}
	MPI_Graph_create(MPI_COMM_WORLD, proc_num, index, edges, reord, &comm_gr);



	while (epsilon <= modul(r, N) / modul(b, N)) {

		MPE_Init_log();
		evtid_beginPhase1 = MPE_Log_get_event_number();
		evtid_endPhase1 = MPE_Log_get_event_number();
		MPE_Describe_state(evtid_beginPhase1, evtid_endPhase1, "Phase1", "cycle_begin");
		MPE_Log_event(evtid_beginPhase1, proc_rank, (char*)0);

		//рассылка частей матриц каждому процессу
		MPI_Scatterv(b, send_num, send_data_begin, MPI_DOUBLE, bufB, send_num[proc_rank], MPI_DOUBLE, 0, MPI_COMM_WORLD);
		MPI_Scatterv(r, send_num, send_data_begin, MPI_DOUBLE, buf_r, send_num[proc_rank], MPI_DOUBLE, 0, MPI_COMM_WORLD);
		MPI_Scatterv(x, send_num, send_data_begin, MPI_DOUBLE, buf_x, send_num[proc_rank], MPI_DOUBLE, 0, MPI_COMM_WORLD);

		double scalar = 0.0, scal_res = 0.0;
		for (i = 0; i<send_num[proc_rank]; i++)
			scalar += buf_r[i] * buf_r[i];
		MPI_Allreduce(&scalar, &scal_res, 1, MPI_DOUBLE, MPI_SUM, comm_gr);
		alpha = scal_res;

		//”множение јz

		for (i = 0; i<send_num[proc_rank]; i++) {
			buf_res[i] = 0;
			z_copy[i] = buf_z[i];
		}

		for (i = 0; i< send_num[proc_rank]; i++) {
			sum = 0.0;
			for (j = 0; j<send_num[proc_rank]; j++)
				sum += bufA[i*N + send_data_begin[proc_rank] + j] * z_copy[j];

			buf_res[i] += sum;
		}


		for (i = 1; i< proc_num; i++) {

			MPI_Sendrecv_replace(z_copy, N / proc_num + 1, MPI_DOUBLE, (proc_rank - 1 + proc_num) % proc_num, 111, (proc_rank + 1) % proc_num, 111, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

			for (k = 0; k< send_num[proc_rank]; k++) {
				sum = 0.0;
				for (j = 0; j<send_num[(proc_rank + i) % proc_num]; j++)
					sum += bufA[k*N + send_data_begin[(proc_rank + i) % proc_num] + j] * z_copy[j];

				buf_res[k] += sum;

			}
		}

		scalar = 0.0; scal_res = 0.0;
		for (i = 0; i<send_num[proc_rank]; i++)
			scalar += buf_res[i] * buf_z[i];
		MPI_Allreduce(&scalar, &scal_res, 1, MPI_DOUBLE, MPI_SUM, comm_gr);

		//¬ычисление alpha, r(n+1), x(n+1), beta, z(n+1) 
		alpha = alpha / scal_res;
		plus(buf_x, buf_z, alpha, buf_x, send_num[proc_rank]);

		for (i = 0; i< send_num[proc_rank]; i++)
			buf_prev_r[i] = buf_r[i];
		plus(buf_r, buf_res, -1 * alpha, buf_r, send_num[proc_rank]);


		scalar = 0.0; scal_res = 0.0;
		for (i = 0; i<send_num[proc_rank]; i++)
			scalar += buf_r[i] * buf_r[i];
		MPI_Allreduce(&scalar, &scal_res, 1, MPI_DOUBLE, MPI_SUM, comm_gr);
		beta = scal_res;
		scalar = 0.0; scal_res = 0.0;
		for (i = 0; i<send_num[proc_rank]; i++)
			scalar += buf_prev_r[i] * buf_prev_r[i];
		MPI_Allreduce(&scalar, &scal_res, 1, MPI_DOUBLE, MPI_SUM, comm_gr);
		beta = beta / scal_res;

		plus(buf_r, buf_z, beta, buf_z, send_num[proc_rank]);

		MPI_Allgatherv(buf_x, get_num[proc_rank], MPI_DOUBLE, x, get_num, get_data_begin, MPI_DOUBLE, MPI_COMM_WORLD);
		MPI_Allgatherv(buf_r, get_num[proc_rank], MPI_DOUBLE, r, get_num, get_data_begin, MPI_DOUBLE, MPI_COMM_WORLD);
		MPI_Allgatherv(bufB, get_num[proc_rank], MPI_DOUBLE, b, get_num, get_data_begin, MPI_DOUBLE, MPI_COMM_WORLD);

		MPE_Log_event(evtid_endPhase1, proc_num, (char*)0);
		MPE_Finish_log("full_par.clog2");


	}


	if (proc_rank == 0) {
		for (i = 0; i< N; i++) {
			printf("%f, \n", x[i]);
		}

	}

	clock_gettime(CLOCK_MONOTONIC_RAW, &end);
	printf("Time taken: %lf sec.\n", end.tv_sec - start.tv_sec + 0.000000001*(end.tv_nsec - start.tv_nsec));

	free(r);
	free(buf_r);
	free(buf_x);
	free(buf_prev_r);
	free(buf_z);
	free(x);
	free(b);
	free(bufB);
	free(bufA);
	free(buf_res);
	free(bufC);
	free(z_copy);
	MPI_Finalize();

	return 0;

}
