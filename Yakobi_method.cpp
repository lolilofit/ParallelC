#include<stdio.h>
#include<malloc.h>
#include<mpi.h>
#include<time.h>

#define n1 3
#define n2 3
#define n3 3

double Xmax = 1.0, Xmin = -1.0, Ymax = 1.0, Ymin = -1.0, Zmax = 1.0, Zmin = -1.0;
const double a = 10000.0;
const double epsilon = 0.0000001;
const double x0 = -1.0;
const double y0 = -1.0;
const double z0 = -1.0;
double hx = -1.0, hy = -1.0, hz = -1, max = 0.0, global_max;

double f(double x, double y, double z) {
	return x * x + y * y + z * z;
}

double p(double x, double y, double z) {
	return 6.0 - a * f(x, y, z);
}

double coord_x(int i) {
	return x0 + i * hx;
}

double coord_y(int j) {
	return y0 + j * hy;
}

double coord_z(int k) {
	return z0 + k * hz;
}

double modul(double val) {
	if (val > 0.0)
		return val;
	else
		return (-1.0)*val;
}

double next_iter(double i_next, double i_prev, double j_next, double j_prev, double k_next, double k_prev, double current, double x, double y, double z) {

	return ((1 / (2 / (hx*hx) + 2 / (hy*hy) + 2 / (hz*hz)) + a)*((i_next + i_prev) / (hx*hx) + (j_next + j_prev) / (hy*hy) + (k_next + k_prev) / (hz*hz) - p(x, y, z)));
}


void find_local_max(double *val, double *prev_val, int high) {
	int i, j, k;

	for (i = 0; i < high; i++) {
		for (j = 0; j < n2; j++) {
			for (k = 0; k < n3; k++) {
				if ((val[i*n2*n3 + j * n3 + k] - prev_val[i*n2*n3 + j * n3 + k]) > max)
					max = val[i*n2*n3 + j * n3 + k] - prev_val[i*n2*n3 + j * n3 + k];
			}
		}
	}
}


void calculate_value(double* val, double prev, double next, int i, int j, int k, int beg) {
	double j_next = 0.0, j_prev = 0.0, k_next = 0.0, k_prev = 0.0;
	int i1, j1, k1;


	if (j != n2 - 1)
		j_next = val[i*n2*n3 + (j + 1)*n3 + k];
	if (j != 0)
		j_prev = val[i*n2*n3 + ((j - 1)*n3) + k];
	if (k != n3 - 1)
		k_next = val[i*n2*n3 + j * n3 + k + 1];
	if (k != 0)
		k_prev = val[i*n2*n3 + j * n3 + k - 1];

	val[i*n2*n3 + j * n3 + k] = next_iter(next, prev, j_next, j_prev, k_next, k_prev, val[i*n2*n3 + j * n3 + k], coord_x(beg + i), coord_y(j), coord_z(k));

}



int main(int argc, char* argv[]) {

	int i, j, k, proc_num, proc_rank, real_proc_num;
	double Dx, Dy, Dz;
	double *next = (double*)malloc(sizeof(double) * (n2*n3 + 1));
	double *prev = (double*)malloc(sizeof(double) * (n2*n3 + 1));

	Dx = Xmax - Xmin;
	Dy = Ymax - Ymin;
	Dz = Zmax - Zmin;

	hx = Dx / (n1 - 1);
	hy = Dy / (n2 - 1);
	hz = Dz / (n3 - 1);




	MPI_Init(&argc, &argv);

	MPI_Comm_size(MPI_COMM_WORLD, &proc_num);
	MPI_Comm_rank(MPI_COMM_WORLD, &proc_rank);


	double* val = (double*)malloc(sizeof(double) * 2 * (n1 + 1)*(n2 + 1)*(n3 + 1));
	double *prev_val = (double*)malloc(sizeof(double) * 2 * (n1 + 1)*(n2 + 1)*(n3 + 1));
	int *send_num = (int*)malloc(sizeof(int)*proc_num);
	int *send_data_begin = (int*)malloc(sizeof(int)*proc_num);

	for (i = 0; i < proc_num; ++i)
		send_num[i] = n1 / proc_num;

	for (i = 0; i<n1%proc_num; ++i)
		send_num[i] = n1 / proc_num + 1;
	send_data_begin[0] = 0;

	for (i = 1; i<proc_num; ++i) {
		send_data_begin[i] = send_data_begin[i - 1] + send_num[i - 1];
	}
	
	
	//если число процессов больше, чем число строк, посчитать сколько процессов //реально задействуется в вычислениях
	real_proc_num = proc_num;
	for (i = proc_num - 1; i >= 0; --i) {
		if (send_num[i] != 0) {
			real_proc_num = i + 1;
			i = -1;
		}
	}

	
	//инициализация начальными значениями в каждом потоке
	for (i = 0; i< send_num[proc_rank]; ++i) {
		for (j = 0; j<n2; ++j) {
			for (k = 0; k<n3; ++k) {
				if (j != 0 && j != (n2 - 1) && k != 0 && k != (n3 - 1) && (send_data_begin[proc_rank] + i) != 0 && (send_data_begin[proc_rank] + i) != (n1 - 1)) {
					val[i*(n2*n3) + j * n3 + k] = 0;
				}
				else
					val[i*(n2*n3) + j * n3 + k] = f(coord_x(send_data_begin[proc_rank] + i), coord_y(j), coord_z(k));

			}
		}
	}


	double global_max, j_next = 0.0, j_prev = 0.0, k_next = 0.0, k_prev = 0.0;

	for (i = 0; i < n2; ++i) {
		for (j = 0; j<n3; ++j) {
			prev[i*n3 + j] = 0.0;
			next[i*n3 + j] = 0.0;
		}
	}


	for (i = 0; i < send_num[proc_rank]; ++i) {
		for (j = 0; j < n2; ++j) {
			for (k = 0; k<n3; ++k)
				prev_val[i*n2*n3 + j * n3 + k] = val[i*n2*n3 + j * n3 + k];
		}
	}

	MPI_Request reqf[2];
	MPI_Request reqs[2];
	MPI_Status statf[2];
	MPI_Status stats[2];

	//обмен граничными значениями
	if (proc_rank < real_proc_num) {
		if (proc_rank != (real_proc_num - 1)) {
			MPI_Isend(val, send_num[proc_rank] * n2*n3, MPI_DOUBLE, proc_rank + 1, 123, MPI_COMM_WORLD, &reqf[0]);
			MPI_Irecv(next, n1*n2*n3, MPI_DOUBLE, proc_rank + 1, 124, MPI_COMM_WORLD, &reqs[1]);
		}

		if (proc_rank != 0) {
			MPI_Isend(val, send_num[proc_rank] * n2*n3, MPI_DOUBLE, proc_rank - 1, 124, MPI_COMM_WORLD, &reqs[0]);
			MPI_Irecv(prev, n1*n3*n2, MPI_DOUBLE, proc_rank - 1, 123, MPI_COMM_WORLD, &reqf[1]);
		}


		if (proc_rank != (real_proc_num - 1)) {
			MPI_Wait(&reqf[0], &statf[0]);
			MPI_Wait(&reqs[1], &stats[1]);
		}
		if (proc_rank != 0) {
			MPI_Wait(&reqf[1], &statf[1]);
			MPI_Wait(&reqs[0], &stats[0]);

		}
	}


	global_max = epsilon + 1;


	while (global_max >= epsilon) {

		max = 0.0;

		if (proc_rank < real_proc_num) {


			j_next = 0.0; j_prev = 0.0; k_next = 0.0; k_prev = 0.0;


			for (i = 0; i< send_num[proc_rank]; ++i) {
				for (j = 0; j < n2; ++j) {
					for (k = 0; k < n3; ++k)
						prev_val[i*n2*n3 + j * n3 + k] = val[i*n2*n3 + j * n3 + k];
				}
			}
			
			//пересчитать точки на границе
			for (j = 1; j < (n2 - 1); ++j) {
				for (k = 1; k< (n3 - 1); ++k) {
					if (proc_rank != 0 && (send_num[proc_rank] > 1 || proc_rank != (real_proc_num - 1)))
						calculate_value(val, prev[j*n3 + k], next[j*n3 + k], 0, j, k, send_data_begin[proc_rank]);


					if (send_num[proc_rank] > 1 && proc_rank != (real_proc_num - 1))
						calculate_value(val, prev[(send_num[proc_rank] - 1)*n2*n3 + j * n3 + k], next[(send_num[proc_num] - 1)*n2*n3 + j * n3 + k], send_num[proc_rank] - 1, j, k, send_data_begin[proc_rank]);

				}
			}


			for (i = 0; i < n2; ++i) {
				for (j = 0; j<n3; ++j) {
					prev[i*n3 + j] = 0.0;
					next[i*n3 + j] = 0.0;
				}
			}

			
			//обмен граничными значениями
			if (proc_rank != (real_proc_num - 1)) {
				MPI_Isend(val, send_num[proc_rank] * n2*n3, MPI_DOUBLE, proc_rank + 1, 123, MPI_COMM_WORLD, &reqf[0]);
				MPI_Irecv(next, n1*n2*n3, MPI_DOUBLE, proc_rank + 1, 124, MPI_COMM_WORLD, &reqs[1]);
			}

			if (proc_rank != 0) {
				MPI_Isend(val, send_num[proc_rank] * n2*n3, MPI_DOUBLE, proc_rank - 1, 124, MPI_COMM_WORLD, &reqs[0]);
				MPI_Irecv(prev, n1*n3*n2, MPI_DOUBLE, proc_rank - 1, 123, MPI_COMM_WORLD, &reqf[1]);
			}


			//пересчет внутренних точек области
			for (i = 1; i<send_num[proc_rank] - 1; ++i) {
				for (j = 1; j < (n2 - 1); ++j) {
					for (k = 1; k < (n3 - 1); ++k) {
						calculate_value(val, val[(i - 1)*n2*n3 + j * n3 + k], val[(i + 1)*n2*n3 + j * n3 + k], i, j, k, send_data_begin[proc_rank]);
					}
				}
			}


			//найти максимум в каждом потоке для своего куска
			find_local_max(val, prev_val, send_num[proc_rank]);


			if (proc_rank != (real_proc_num - 1)) {
				MPI_Wait(&reqf[0], &statf[0]);
				MPI_Wait(&reqs[1], &stats[1]);
			}
			if (proc_rank != 0) {
				MPI_Wait(&reqf[1], &statf[1]);
				MPI_Wait(&reqs[0], &stats[0]);

			}
		}

		
		//разослать максимумы и найти общий максимум
		if (proc_rank != 0) {
			MPI_Send(&max, 1, MPI_DOUBLE, 0, 222, MPI_COMM_WORLD);
			MPI_Recv(&global_max, 1, MPI_DOUBLE, 0, 223, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
		}
		else {
			global_max = max;
			double p_max = -1.0;

			for (i = 1; i < proc_num; ++i) {
				MPI_Recv(&p_max, 1, MPI_DOUBLE, i, 222, MPI_COMM_WORLD, MPI_STATUS_IGNORE);


				if (global_max < p_max)
					global_max = p_max;
			}
			for (i = 1; i<proc_num; ++i)
				MPI_Send(&global_max, 1, MPI_DOUBLE, i, 223, MPI_COMM_WORLD);
		}


		MPI_Barrier(MPI_COMM_WORLD);

	}



	free(prev);
	free(next);
	free(val);
	free(prev_val);
	free(send_num);
	free(send_data_begin);

	MPI_Finalize();


	return 0;
}

