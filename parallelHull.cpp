#include <iostream>
#include <cstdlib>
#include <cstdio>
#include <vector> //std::vector
#include <cmath> //sqrt
#include <stdio.h>
#include <stdlib.h>
#include <pthread.h>
#include <unistd.h>
#include <chrono>
#include <thread>
#include <algorithm>

using namespace std; //std::


pthread_mutex_t mut;
int TERMINATE = 0;

class Point{
    public:
        double x, y;
};


vector<Point> convex_hull;

typedef struct Task{
	void (*task_function)(vector<Point>,Point, Point);
	vector<Point> subset;
	Point p1;
	Point p2;
} Task;


//-----------------------------PTHREAD UTILITIES:START-----------------------------------------


Task taskQueue[500];
int task_count=0;

pthread_mutex_t queue_mutex;
pthread_cond_t condQueue;


void executeTask(Task* task){
	task->task_function(task->subset,task->p1,task->p2);
}

void submitTask(Task task) {
    pthread_mutex_lock(&queue_mutex);
    taskQueue[task_count] = task;
    task_count++;
	cout<<"\nNew Task Submitted, Current Tasks: "<<task_count<<endl;
    pthread_mutex_unlock(&queue_mutex);
    pthread_cond_signal(&condQueue);
}


void* start_thread(void* args){
	cout <<"\nThread "<<(int *) args<<" starting..." <<endl;

	while (TERMINATE==0) {
        Task task;

        pthread_mutex_lock(&queue_mutex);

        while (task_count == 0) {
            pthread_cond_wait(&condQueue, &queue_mutex);
        }

        task = taskQueue[0];
        int i;
        for (i = 0; i < task_count - 1; i++) {
            taskQueue[i] = taskQueue[i + 1];
        }
        task_count--;

        pthread_mutex_unlock(&queue_mutex);

        executeTask(&task);
		
		this_thread::sleep_for(chrono::milliseconds(2000));
		if(task_count==0) TERMINATE=1;
    }

}


//-----------------------------QUICKHULL UTILITIES:START-----------------------------------------
//----------------- input and output ----------------------*/
void get_input(vector<Point> &p){
    Point aux; //does the pushback
    //fill the vector from a input file (syntax of a line: double double)
    while (fscanf(stdin, "%lf", &aux.x) != EOF) {    // read a double: x coordinate
	    getc(stdin);                             // read a single space character
	    fscanf(stdin, "%lf", &aux.y);            // read a double: y coordinate
	    p.push_back(aux);
    }
}

void print_vector_of_points(vector<Point> &p){
	for(auto i = p.begin(); i != p.end(); ++i){
		cout << (*i).x << " " << (*i).y << "\n";
	}
	cout << "=================================\n";
}


/*--------------------- Point operations ----------------------*/
// sums two vectors
Point vector_sum(Point a, Point b){
	Point sum;
	sum.x = a.x + b.x;
	sum.y = a.y + b.y;
	return sum;
}

// the product of a vector by a scalar
Point scalar_product(double l, Point a){
	Point product;
	product.x = l * a.x;
	product.y = l * a.y;
	return product;
}

// cross product (only Z-axis value)
double cross_product(Point a, Point b){
    return (((a.x)*(b.y)) - ((b.x)*(a.y)));
}

// dot product
double dot_product(Point a, Point b){
    return (((a.x)*(b.x)) + ((a.y)*(b.y)));
}

// euclidean norm: || . ||
double norm(Point a){
    return (sqrt(dot_product(a,a)));
}

// euclidean distance: || b - a ||
double euclidean_distance(Point a, Point b){
	return(norm(vector_sum(b,scalar_product(-1,a))));
}

// the shortest distance between a point and a line (defined by two points)
// see references (1)
double line_distance(Point a, Point p1, Point p2){
	return(abs(((p2.y - p1.y)*a.x)-((p2.x - p1.x)*a.y)+(p2.x * p1.y)-(p2.y * p1.x))/euclidean_distance(p1,p2));
}


/*-------------------- auxiliary functions --------------------*/
// says if a point a is at the left side, right side or ir colinear to the oriented line p1->p2
// see references (2)
int set_local(Point a, Point p1, Point p2){
	//calculates:	p1->p2 x p1->a (the cross product between the vectors)
	double position = cross_product((vector_sum(p2,scalar_product(-1,p1))),(vector_sum(a,scalar_product(-1,p1))));
	if(position > 0){
		return -1; //point c is at the left of the line p2->p1
	}
	else if(position < 0){
		return 1;  //point c is at the left of the line p1->p2
	}
	return 0; 	   //the point c is collinear to the line p1->p2
}

//returns the farthest point from a Set of points and a line (defined by two points)
Point farthest_point(vector<Point> S, Point p1, Point p2){
	double max_value = 0, aux;
	Point max_point;
	//caculates and compares distance of points
	for(auto i = S.begin(); i != S.end(); ++i){
		aux = line_distance((*i),p1,p2);
		if(max_value <= aux){
			max_value = aux;
			max_point = (*i);
		}
	}
	return max_point;
}

//find the set of points at the left side of a line (defined by two points)
void find_left_set(vector<Point> P, vector<Point>& S, Point p1, Point p2){
	int flag;
	for(auto i = P.begin(); i != P.end(); ++i){
		flag = set_local((*i),p1,p2);
		if(flag==1){
			S.push_back(*i);
		}
	}
}

void find_hull(vector<Point>S, Point p1, Point p2){
	//creates two regions
	vector<Point> S1,S2;
	//the farthest point
	Point c;

	//if S is empty there's no hull to be found

	if(S.empty()) return;

    pthread_mutex_lock(&queue_mutex);

	//finds the farthest point from the line p1->p2
	c = farthest_point(S,p1,p2);

    pthread_mutex_unlock(&queue_mutex);

	//separates the set S in three regions
	find_left_set(S,S1,p1,c);
	find_left_set(S,S2,c,p2);   
   

	//find_hull(CH,S1,p1,c); //sequential
	submitTask((Task) { //parallel
		.task_function = &find_hull,
		.subset = S1,
		.p1 = p1,
		.p2 = c,
	}); 

    pthread_mutex_lock(&queue_mutex);

	convex_hull.push_back(c); 
    pthread_mutex_unlock(&queue_mutex);

	//find_hull(convex_hull,S2,c,p2);	//sequential
	submitTask((Task){ //parallel
		.task_function = &find_hull,
		.subset = S2,
		.p1 = c,
		.p2 = p2,
	}); 
}

//-----------------------------QUICKHULL UTILITIES:END-----------------------------------------

//-----------------------------MAIN THREAD-----------------------------------------------------

int main(int argc, char *argv[]){
	auto start = chrono::high_resolution_clock::now();

    pthread_t *tid;
	pthread_mutex_init(&queue_mutex, NULL);
    pthread_cond_init(&condQueue, NULL);
	int i, numOfThreads;

	if (argc != 2) {
		cout<<"Provide the number of threads to create\n";
		exit(0);
	}

	numOfThreads = atoi(argv[1]); 

	tid = (pthread_t *)malloc(numOfThreads * sizeof(pthread_t));
	if (tid == NULL) { 
		cout<<"<error in memory allocation";
    }
    

    //------------FIRST STEPS: START------------------
    vector<Point> points;
    get_input(points);

   // cout << "\n------------- INPUT -------------\n";
	//print_vector_of_points(points);

	// Find points with smallest and biggest value for x (leftmost and rightmost points)
	Point min = *(points.begin());
	Point max = *(points.begin());

	for (auto i = points.begin(); i != points.end(); i++) {
		if ((*i).x < min.x) {
			min = *i;
		}
		if ((*i).x > max.x) {
			max = *i;
		}
	}

	vector<Point> S1,S2;

	//S1: points at the left of line min->max
	find_left_set(points,S1,min,max);
	//S2: points at the left of line max->min 
	find_left_set(points,S2,max,min);

    //--------------FIRST STEPS:END------------------
	

	for (i = 0; i < numOfThreads; i++) {
		int *thread_id;
		thread_id = (int *)malloc(sizeof(int));
		*thread_id= i;
		pthread_create(&tid[i], NULL, &start_thread,(void *)thread_id);
	}

	convex_hull.push_back(min); 
	
	//find_hull(convex_hull,subset1,min,max); //sequential
	submitTask((Task) { //parallel
		.task_function = &find_hull,
		.subset = S1,
		.p1 = min,
		.p2 = max,
	});

	//adds right most point to the convex hull 
	convex_hull.push_back(max); 

	//find_hull(convex_hull,subset2,min,max); //sequential
	submitTask((Task){ //parallel
		.task_function = &find_hull,
		.subset = S2,
		.p1 = min,
		.p2 = max,
	}); 

	convex_hull.push_back(min); 

	for (i = 0; i < numOfThreads; i++) {
	if (pthread_join(tid[i], NULL) != 0) {
		perror("Failed to join the thread");
		}
	}

	cout << "\n========== CONVEX HULL ==========\n";
	print_vector_of_points(convex_hull);

	auto stop = chrono::high_resolution_clock::now();
	auto duration = chrono::duration_cast<chrono::milliseconds>(stop - start);
	cout<< duration.count() << "ms elapsed."<<endl;
	pthread_mutex_destroy(&queue_mutex);
    pthread_cond_destroy(&condQueue);
	//ends program
	return 0;
}

