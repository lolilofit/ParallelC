#include<stdio.h>
#include<malloc.h>
#include<mpi.h>
#include<time.h>
#include<mpe.h>

const int n1 = 32, n2 = 32, n3 = 32, p1 = 4, p2 = 2;


int main(int argc, char* argv[]) {


	int i, j, k, proc_num, proc_rank;
	int evtid_beginPhase1, evtid_endPhase1, evtid_beginPhase2, evtid_endPhase2, evtid_beginPhase3, evtid_endPhase3, evtid_beginPhase4, evtid_endPhase4;

	double* A = (double*)malloc(sizeof(double)*n1*n2);
	double* B = (double*)malloc(sizeof(double)*n2*n3);
	double* C = (double*)malloc(sizeof(double)*(n1*n3));

	int periods[] = { 0, 0 };
	int dims[] = { p1, p2 };
	int p[2], n[3];

	MPI_Comm com;


	MPI_Init(&argc, &argv);
	//      MPE_Init_log();

	MPI_Comm_size(MPI_COMM_WORLD, &proc_num);
	MPI_Comm_rank(MPI_COMM_WORLD, &proc_rank);

	struct timespec start, end;
	clock_gettime(CLOCK_MONOTONIC_RAW, &start);

	MPI_Cart_create(MPI_COMM_WORLD, 2, dims, periods, 0, &com);
	MPI_Barrier(MPI_COMM_WORLD);

	/*      evtid_beginPhase1 = MPE_Log_get_event_number();
	evtid_endPhase1 = MPE_Log_get_event_number();
	MPE_Describe_state(evtid_beginPhase1, evtid_endPhase1, "prepearing", "red");

	MPE_Log_event(evtid_beginPhase1, proc_rank, (char*)0);
	*/

	int *send_num = (int*)malloc(sizeof(int)*proc_num);
	int *send_data_begin = (int*)malloc(sizeof(int)*proc_num);

	for (i = 0; i < p1; ++i) {
		send_num[i] = (n1 / p1)*n2;
		send_data_begin[i] = 0;
	}

	for (i = 0; i<n1%p1; ++i)
		send_num[i*p2] = (n1 / p1 + 1)*n2;

	for (i = 1; i<p1; ++i)
		send_data_begin[i] = send_data_begin[i - 1] + send_num[i - 1];


	int *send_num_b = (int*)malloc(sizeof(int)*proc_num);
	int *send_data_begin_b = (int*)malloc(sizeof(int)*proc_num);

	for (i = 0; i<p2; ++i)
		send_num_b[i] = (n3 / p2)*n2;

	for (i = 0; i<n1%p2; ++i)
		send_num_b[i] = (n3 / p2 + 1)*n2;

	send_data_begin_b[0] = 0;
	for (i = 1; i<p2; ++i)
		send_data_begin_b[i] = send_data_begin_b[i - 1] + send_num_b[i - 1];


	double* bufA = (double*)malloc(sizeof(double)*(n1 / p1 + 1)*n2);
	double* bufB = (double*)malloc(sizeof(double)*(n3 / p2 + 1)*n2);
	double* bufC = (double*)malloc(sizeof(double)*(n3 / p2 + 1)*(n1 / p1 + 1));
	double* B_trans = (double*)malloc(sizeof(double)*n2*n3);


	if (proc_rank == 0) {

		for (i = 0; i<n1; ++i) {
			for (j = 0; j<n2; ++j) {
				if (i == j)
					A[i*n2 + j] = 2.0;
				else
					A[i*n2 + j] = 1.0;
			}
		}

		for (i = 0; i<n2; ++i) {
			for (j = 0; j<n3; ++j)
				B[i*n3 + j] = 1.0;
		}

		for (i = 0; i<n2; ++i) {
			for (j = 0; j<n3; ++j)
				B_trans[j*n2 + i] = B[i*n3 + j];
		}
	}

	int remain_dims[] = { 1 , 0 };
	int remain[] = { 0, 1 };
	MPI_Comm groupsA, groupsB;
	MPI_Cart_sub(com, remain_dims, &groupsA);
	MPI_Cart_sub(com, remain, &groupsB);

	/*      MPE_Log_event(evtid_endPhase1, proc_rank, (char*)0);

	evtid_beginPhase2 = MPE_Log_get_event_number();
	evtid_endPhase2 = MPE_Log_get_event_number();
	MPE_Describe_state(evtid_beginPhase2, evtid_endPhase2, "Send A and B", "blue");

	MPE_Log_event(evtid_beginPhase2, proc_rank, (char*)0);
	*/      if (proc_rank %p2 == 0)
	MPI_Scatterv(A, send_num, send_data_begin, MPI_DOUBLE, bufA, (n1 / p1)*n2, MPI_DOUBLE, 0, groupsA);

	if (proc_rank / p2 == 0)
		MPI_Scatterv(B_trans, send_num_b, send_data_begin_b, MPI_DOUBLE, bufB, (n3 / p2)*n2, MPI_DOUBLE, 0, groupsB);


	MPI_Bcast(bufA, send_num[proc_rank / p2], MPI_DOUBLE, 0, groupsB);
	MPI_Bcast(bufB, send_num_b[proc_rank%p2], MPI_DOUBLE, 0, groupsA);
	//      MPE_Log_event(evtid_endPhase2,  proc_rank, (char*)0);

	MPI_Barrier;
	free(B_trans);

	/*      evtid_beginPhase3 = MPE_Log_get_event_number();
	evtid_endPhase3 = MPE_Log_get_event_number();

	MPE_Describe_state(evtid_beginPhase3, evtid_endPhase3, "mult", "green");
	MPE_Log_event(evtid_beginPhase3, proc_rank, (char*)0);
	*/
	for (i = 0; i< send_num[proc_rank / p2] / n2; ++i) {
		for (j = 0; j< send_num_b[proc_rank%p2] / n2; ++j) {
			bufC[i*send_num_b[proc_rank%p2] / n2 + j] = 0;
			for (k = 0; k<n2; ++k)
				bufC[i* send_num_b[proc_rank%p2] / n2 + j] += bufA[i*n2 + k] * bufB[j* n2 + k];
		}
	}

	MPI_Barrier;
	//      MPE_Log_event(evtid_endPhase3, proc_rank, (char*)0);

	int* get_num = (int*)malloc(sizeof(int)*p1);
	int* get_data_begin = (int*)malloc(sizeof(int)*p1);
	double* subC = (double*)malloc(sizeof(double)*n1*n3);

	/*      evtid_beginPhase4 = MPE_Log_get_event_number();
	evtid_endPhase4 = MPE_Log_get_event_number();
	MPE_Describe_state(evtid_beginPhase4, evtid_endPhase4, "get C", "yellow");

	MPE_Log_event(evtid_beginPhase4, proc_rank, (char*)0);
	*/      for (i = 0; i< p2; ++i) {
		for (j = 0; j< p1; ++j)
			get_num[j] = (send_num_b[proc_rank%p2] / n2)*(n1 / p1);
		get_data_begin[0] = 0;
		for (j = 1; j<p1; ++j)
			get_data_begin[j] = get_data_begin[j - 1] + get_num[j - 1];


		MPI_Gatherv(bufC, (n3 / p2)*(n1 / p1), MPI_DOUBLE, subC, get_num, get_data_begin, MPI_DOUBLE, 0, groupsA);



		if (proc_rank / p2 == 0 && proc_rank != 0)
			MPI_Send(subC, n1*(send_num_b[i] / n2), MPI_DOUBLE, 0, 123, MPI_COMM_WORLD);
		if (proc_rank == 0 && i != 0)
			MPI_Recv(subC, n1*n3, MPI_DOUBLE, i, 123, MPI_COMM_WORLD, MPI_STATUS_IGNORE);


		if (proc_rank == 0) {

			for (j = 0; j<n1; ++j) {
				for (k = 0; k<send_num_b[i] / n2; ++k) {
					C[j*n3 + send_data_begin_b[i] / n2 + k] = subC[j*(send_num_b[i] / n2) + k];
				}
			}
		}

		MPI_Barrier;
	}

	//      MPE_Log_event(evtid_endPhase4, proc_rank, (char*)0);

	clock_gettime(CLOCK_MONOTONIC_RAW, &end);

	if (proc_rank == 0) {
		for (i = 0; i<n1; ++i) {
			for (j = 0; j<n3; ++j)
				printf("%f ", C[i*n3 + j]);
			printf("\n");
		}
	}

	if (proc_rank == 0)
		printf("%d : Time taken: %lf sec.\n", n1, end.tv_sec - start.tv_sec + 0.000000001*(end.tv_nsec - start.tv_nsec));


	free(A);
	free(C);
	free(B);
	free(bufA);
	free(bufB);
	free(bufC);
	free(get_num);
	free(get_data_begin);
	free(subC);
	free(send_num);
	free(send_num_b);
	free(send_data_begin_b);
	free(send_data_begin);

	//       MPE_Finish_log("main");
	MPI_Finalize();

	return 0;

}
