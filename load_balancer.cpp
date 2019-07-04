#include<mpi.h>
#include<stdio.h>
#include<math.h>
#include<unistd.h>
#include <pthread.h>
#include<stdlib.h>
#include <stddef.h>

#define ITER_NUM 2
#define TASK_NUM 16

pthread_t* threads;
pthread_mutex_t iter_mutext;
int global_iter, proc_num, proc_rank;
int *ids;
struct List task_list;
int calculating_now, given_num = -1, is_doing = 1;
int f = 0, signal = 0, signal_rec = 0;


struct Task {
	double val;
	int weight;
};

struct List {
	struct Task *tasks[TASK_NUM];
};

void init_list(struct List *l) {
	int i;
	for (i = 0; i < TASK_NUM; ++i)
		l->tasks[i] = NULL;
}

void add_task(struct Task *task, struct List* l, int pos) {
	if (l->tasks[pos] == NULL) {
		l->tasks[pos] = task;
	}
}

void remove_task(struct List* l, int pos) {
	l->tasks[pos] = NULL;
}

double _sqrt(const double x) {
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


double do_task(struct Task task) {
	int i, res = 0;
	res = task.val;
	for (i = 0; i < task.weight; ++i) {
		res = _sqrt(res);
		res += 4;
		res *= 2;
		res = res / 4;
		res = res / 2;
	}
	return res;
}

void create_task(struct Task* task) {
	if (proc_rank == 1) {
		task->weight = proc_rank + 100 + 5000;
		task->val = rand() % 100 + 1;
	}
	else {
		task->weight = proc_rank + 1;
		task->val = rand() % 100 + 1;
	}
}

void* calc(void* args) {
	int iter_num, j, count, i, proc = 0, get = 0, question;
	double res = 0.0;

	for (iter_num = 0; iter_num < ITER_NUM; iter_num++) {

		pthread_mutex_lock(&iter_mutext);
		init_list(&task_list);

		for (j = 0; j < TASK_NUM; ++j) {
			struct Task *t = (struct Task*)malloc(sizeof(struct Task));
			create_task(t);
			add_task(t, &task_list, j);
		}



		for (proc = 0; proc < proc_num; proc++) {

			if (proc != proc_rank) {
				struct Task *given = (struct Task*)malloc(sizeof(struct Task));
				int answer;

				for (j = 0; j < TASK_NUM; ++j) {
					if (task_list.tasks[j] != NULL) {
						struct Task do_it = *(task_list.tasks[j]);

						calculating_now = j;
						given_num = -1;
						if (get == 0) {
							pthread_mutex_unlock(&iter_mutext);
							res = do_task(do_it);
							pthread_mutex_lock(&iter_mutext);
						}
						else {
							res = do_task(do_it);
						}
						task_list.tasks[j] = NULL;
						get = 0;
					}
				}
				calculating_now = -1;
				pthread_mutex_unlock(&iter_mutext);

				question = 1;
				printf("%d : mes_send to %d\n", proc_rank, proc);
				MPI_Send(&question, 1, MPI_INT, proc, 123, MPI_COMM_WORLD);


				MPI_Recv(&get, 1, MPI_INT, proc, 10, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

				if (get == 1) {
					int w; double v;
					pthread_mutex_lock(&res_mutext);
					MPI_Recv(&w, 1, MPI_INT, proc, 125, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
					MPI_Recv(&v, 1, MPI_DOUBLE, proc, 126, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
					pthread_mutex_unlock(&res_mutext);
					given->weight = w;
					given->val = v;
				}

				free(given);
				pthread_mutex_lock(&iter_mutext);
				if (get == 1) {
					for (i = 0; i < TASK_NUM; i++) {
						if (task_list.tasks[i] == NULL) {
							add_task(given, &task_list, i);
							i = TASK_NUM;
						}
					}
				}


			}
		}

		for (j = 0; j < TASK_NUM; ++j) {
			if (task_list.tasks[j] != NULL) {
				struct Task do_it = *(task_list.tasks[j]);

				calculating_now = j;
				given_num = -1;
				if (get == 0) {
					pthread_mutex_unlock(&iter_mutext);
					res = do_task(do_it);
					pthread_mutex_lock(&iter_mutext);
				}
				else {
					res = do_task(do_it);
				}
				task_list.tasks[j] = NULL;
				get = 0;
			}
		}
		calculating_now = -1;
		pthread_mutex_unlock(&iter_mutext);


		MPI_Barrier(MPI_COMM_WORLD);

	}

	return NULL;
}


void* give_task(void* args) {
	int i, j, k, proc_num_copy, flag = 0, proc_rank_copy, ending_count = 0, sending_count = 0;
	int w; double v;
	int question = 0, answer = 0;
	MPI_Status status;

	proc_num_copy = proc_num;
	proc_rank_copy = proc_rank;


	while (1) {

		if ((ending_count % (proc_num_copy - 1) == 0) && (ending_count / (proc_num_copy - 1) == ITER_NUM))

			//  if(ending_count == (proc_num_copy - 1))
			//        break;

			MPI_Recv(&question, 1, MPI_INT, MPI_ANY_SOURCE, 123, MPI_COMM_WORLD, &status);

		ending_count++;


		pthread_mutex_lock(&iter_mutext);


		for (k = 0; k < TASK_NUM; ++k) {
			if (task_list.tasks[k] != NULL && k != calculating_now && k != given_num) {
				printf("%d : 12\n", proc_rank);
				w = task_list.tasks[k]->weight;
				v = task_list.tasks[k]->val;
				task_list.tasks[k] = NULL;
				given_num = k;
				k = TASK_NUM + 1;
				flag = 1;
				printf("%d : 13\n", proc_rank_copy);
			}
		}

		pthread_mutex_unlock(&iter_mutext);

		if (flag == 0) {
			//no tasks
			answer = 0;
			MPI_Send(&answer, 1, MPI_INT, status.MPI_SOURCE, 10, MPI_COMM_WORLD);

		}
		else {
			answer = 1;


			MPI_Send(&answer, 1, MPI_INT, status.MPI_SOURCE, 10, MPI_COMM_WORLD);

			MPI_Send(&w, 1, MPI_INT, status.MPI_SOURCE, 125, MPI_COMM_WORLD);
			MPI_Send(&v, 1, MPI_DOUBLE, status.MPI_SOURCE, 126, MPI_COMM_WORLD);

		}

		flag = 0;


	}
	return NULL;
}


int create_threads(int proc_rank) {
	pthread_attr_t attr;

	if (pthread_attr_init(&attr) != 0) {
		perror("Cannot initialize attributes");
		return 1;
	};

	if (pthread_attr_setdetachstate(&attr, PTHREAD_CREATE_JOINABLE) != 0) {
		perror("Error in setting attributes");
		return 2;
	}

	if (pthread_create(&threads[proc_rank], &attr, give_task, &ids[proc_rank]) != 0) {
		perror("Cannot create calc thread");
		return 3;
	}

	pthread_attr_destroy(&attr);

	calc(NULL);


	if (pthread_join(threads[proc_rank], NULL) != 0) {
		perror("Cannot join a thread");
		return 4;
	}

	return 0;
}


int main(int argc, char* argv[]) {
	int i, j, k, check_gar;

	MPI_Init_thread(&argc, &argv, MPI_THREAD_MULTIPLE, &check_gar);


	MPI_Comm_size(MPI_COMM_WORLD, &proc_num);
	MPI_Comm_rank(MPI_COMM_WORLD, &proc_rank);

	threads = (pthread_t*)malloc(sizeof(pthread_t)*proc_num);
	ids = (int*)malloc(sizeof(int)*proc_num);
	for (i = 0; i < proc_num; i++)
		ids[i] = i;


	pthread_mutex_init(&iter_mutext, NULL);
	srand(time(NULL));


	int count;

	if (create_threads(proc_rank) != 0) {
		pthread_mutex_destroy(&iter_mutext);
		MPI_Finalize();
		return 2;
	}


	pthread_mutex_destroy(&iter_mutext);
	MPI_Finalize();

	return 0;
}

