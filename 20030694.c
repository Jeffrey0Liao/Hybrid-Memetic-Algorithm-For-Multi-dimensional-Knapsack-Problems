#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <time.h>
#include <stdbool.h>
#include <string.h>
#include <sys/time.h>

/* global parameters */
int RAND_SEED[] = {1,20,30,40,50,60,70,80,90,100,110, 120, 130, 140, 150, 160, 170, 180, 190, 200};
int NUM_OF_RUNS = 1;
int MAX_TIME = 300;  //max amount of time permited (in sec)
int num_of_problems;

/* declare parameters for local search here */
int SA_MAX_ITER = 2000; //total number of runs.

/* parameters for evlutionary algorithms*/
static int POP_SIZE = 50;   //population size
static int POOL_SIZE = 40;	//size of mating pool and children pool
static int tourm_size = 5;	// tournament size
int MAX_NUM_OF_GEN = 100000000; //max number of generations
float CROSSOVER_RATE = 0.8;
float MUTATION_RATE = 0.05;
float LS_RATE1 = 0.0;
float LS_RATE2 = 1.0;
float VNS_RATE = 0.1;

/* data structure definition */
struct item_struct{
    int dim; //no. of dimensions
    int* size; //volume of item in all dimensions
    int p;
    double ratio;
    int indx;
};

struct problem_struct{
    int n; //number of items
    int dim; //number of dimensions
    struct item_struct* items;
    int* capacities;  //knapsack capacities
};

struct solution_struct{
    struct problem_struct* prob; //maintain a shallow copy of the problem data
    float objective;
    int feasibility; //indicate the feasiblity of the solution
    int* x; //chromosome vector
    int* cap_left; //capacity left in all dimensions
    int flag;
};

struct solution_struct best_sln;  //global best solution

//return a random number between 0 and 1
float rand_01(){
    float number;
    number = (float) rand();
    number = number/RAND_MAX;
    //printf("rand01=%f\n", number);
    return number;
}

//return a random nunber ranging from min to max (inclusive)
int rand_int(int min, int max){
    int div = max-min+1;
    int val =rand() % div + min;
    //printf("rand_range= %d \n", val);
    return val;
}

// function to create an empty problem instance
void init_problem(int n, int dim, struct problem_struct** my_prob){
    struct problem_struct* new_prob = malloc(sizeof(struct problem_struct));
    new_prob->n=n; new_prob->dim=dim;
    new_prob->items=malloc(sizeof(struct item_struct)*n);
    for(int j=0; j<n; j++)
        new_prob->items[j].size= malloc(sizeof(int)*dim);
    new_prob->capacities = malloc(sizeof(int)*dim);
    *my_prob = new_prob;
}

//example to create problem instances, actual date should come from file
struct problem_struct** load_problems(char* data_file){
    int i,j,k;
    //int num_of_probs;
    FILE* pfile = fopen(data_file, "r");
    if(pfile==NULL)
        {printf("Data file %s does not exist. Please check!\n", data_file); exit(2); }
    fscanf (pfile, "%d", &num_of_problems);
 
    struct problem_struct** my_problems = malloc(sizeof(struct problem_struct*)*num_of_problems);
    for(k=0; k<num_of_problems; k++)
    {
        int n, dim, obj_opt;
        fscanf (pfile, "%d", &n);
        fscanf (pfile, "%d", &dim); fscanf (pfile, "%d", &obj_opt);
        
        init_problem(n, dim, &my_problems[k]);  //allocate data memory
        for(j=0; j<n; j++)
        {
            my_problems[k]->items[j].dim=dim;
            my_problems[k]->items[j].indx=j;
            fscanf(pfile, "%d", &(my_problems[k]->items[j].p)); //read profit data
            //printf("item[j].p=%d\n",my_problems[k]->items[j].p); 
        }
        for(i=0; i<dim; i++)
        {
            for(j=0; j<n; j++)
            {
                fscanf(pfile, "%d", &(my_problems[k]->items[j].size[i])); //read size data
                //printf("my_problems[%i]->items[%i].size[%i]=%d\n",k,j,i,my_problems[k]->items[j].size[i]);
            }
        }
        for(i=0; i<dim; i++){
            fscanf(pfile, "%d", &(my_problems[k]->capacities[i]));
            //printf("capacities[i]=%d\n",my_problems[k]->capacities[i] );
        }
        double sum;
        for(j=0; j<n; j++)
        {
            for (i=0; i<dim; i++){
            	sum += (double)(my_problems[k]->items[j].size[i])/(my_problems[k]->capacities[i]);
            	//sum += (double)(my_problems[k]->items[j].size[i]);
			}
			my_problems[k]->items[j].ratio = ((double)(my_problems[k]->items[j].p))/(sum);
			sum = 0.0;
        }
    }
    fclose(pfile); //close file
    return my_problems;
}

//create an empty solution
struct solution_struct* create_empty_sol(struct problem_struct* my_prob){
    struct solution_struct* my_sln = malloc(sizeof(struct solution_struct));
    my_sln->prob = my_prob;
    my_sln->objective=0; 
	my_sln->cap_left= malloc(sizeof(int)*my_prob->dim);
	for (int i=0; i<my_prob->dim; i++) my_sln->cap_left[i] = my_prob->capacities[i];
    my_sln->x = malloc(sizeof(int)*my_prob->n);
    for(int i=0; i<my_prob->n; i++) my_sln->x[i]=0;
    my_sln->flag = 0;
    return my_sln;
}

// create a random solution based on an empty solution
void random_sol(struct solution_struct* my_sln){
	int i, j, item;
	bool flag = false;
	for (i=0; i<(my_sln->prob->n)*5; i++){
		do {
			item = rand_int(0, my_sln->prob->n-1);
		} while (my_sln->x[item] == 1);
		if (rand_01()<0.50){
			flag = true;
		}
		for (j=0; j< my_sln->prob->dim; j++){
			if (my_sln->prob->items[item].size[j]>my_sln->cap_left[j]){
				flag = false;
				break;
			}
		}
		if (flag){
			my_sln->x[item] = 1;
			for (j=0; j< my_sln->prob->dim; j++){
				my_sln->cap_left[j] = my_sln->cap_left[j] - my_sln->prob->items[item].size[j];				
			}
			my_sln->objective += my_sln->prob->items[item].p;
		}
		flag = false;
	}
}

// copy one solution to another one
bool copy_solution(struct solution_struct* dest_sln, struct solution_struct* source_sln){
    if(source_sln ==NULL) return false;
    if(dest_sln==NULL)
    {
        dest_sln = malloc(sizeof(struct solution_struct));
    }
    else{
        free(dest_sln->cap_left);
        free(dest_sln->x);
    }
    int n = source_sln->prob->n;
    int m =source_sln->prob->dim;
    dest_sln->x = malloc(sizeof(int)*n);
    dest_sln->cap_left=malloc(sizeof(int)*m);
    for(int i=0; i<m; i++)
        dest_sln->cap_left[i]= source_sln->cap_left[i];
    for(int j=0; j<n; j++)
        dest_sln->x[j] = source_sln->x[j];
    dest_sln->prob= source_sln->prob;
    dest_sln->feasibility=source_sln->feasibility;
    dest_sln->objective=source_sln->objective;
    dest_sln->flag = source_sln->flag;
    return true;
}

//update global best solution from sln
void update_best_solution(struct solution_struct* sln){
    if(best_sln.objective < sln->objective)
    copy_solution(&best_sln, sln);
}

//output a given solution to a file
void output_solution(struct solution_struct* sln, char* out_file){
    if(out_file !=NULL){
        FILE* pfile = fopen(out_file, "a"); //append solution data
        fprintf(pfile, "%i\n", (int)sln->objective);
        for(int i=0; i<sln->prob->n; i++)
        {
            fprintf(pfile, "%i ", sln->x[i]);
        }
        fprintf(pfile, "\n");
        /*for(int j=0; j<sln->prob->n; j++)
            fprintf(pfile, "%i ", sln->prob->items[j].p);
        fprintf(pfile, "\n");*/
        fclose(pfile);
    }
    else
        printf("sln.feas=%d, sln.obj=%f\n", sln->feasibility, sln->objective);
}

// free the allocated memories of a solution 
void free_solution(struct solution_struct* sln)
{
    if(sln!=NULL)
    {
        free(sln->x);
        free(sln->cap_left);
        sln->objective=0;
        sln->prob=NULL;
        sln->feasibility=false;
        sln->flag=0;
    }
}

// free allocated memories of a problem instance
void free_problem(struct problem_struct* prob)
{
    if(prob!=NULL)
    {
        if(prob->capacities !=NULL) free(prob->capacities);
        if(prob->items!=NULL)
        {
            for(int j=0; j<prob->n; j++)
            {
                if(prob->items[j].size != NULL)
                    free(prob->items[j].size);
            }
            free(prob->items);
        }
        free(prob);
    }
}

// first-descent local search
void first_descent(struct solution_struct* current_sln){
	int iter =0;
    struct solution_struct* curt_sln = current_sln;
    struct problem_struct* prob = curt_sln->prob;
    struct solution_struct* rand_neighb = create_empty_sol(prob);
    while(iter<SA_MAX_ITER)
    {
        copy_solution(rand_neighb, curt_sln);
        int item1, item2;
        item1 = rand_int(0, prob->n-1);
        if(curt_sln->x[item1] ==1){
            item2 = rand_int(0, prob->n-1);
            while(curt_sln->x[item2] ==1){//careful, this might lead to deadloop
                item2 = rand_int(0, prob->n-1);
            }
        }
        else{
            item2 = rand_int(0, prob->n-1);
            while(curt_sln->x[item2] ==0){//careful, this might lead to deadloop
                item2 = rand_int(0, prob->n-1);
            }
            int temp = item1;
            item1 = item2;
            item2 = temp;
        }
		
        //testing potential constraint violations after swap
        bool flag=true;
        for(int d=0; d<prob->dim; d++){
            if(rand_neighb->cap_left[d] + prob->items[item1].size[d] <
               prob->items[item2].size[d]){
                flag=false;
                break;
            }
        }
        if(flag){//can swap
            float delta = prob->items[item2].p - prob->items[item1].p;
            if(delta>=0){
                rand_neighb->x[item1]=0;
                rand_neighb->x[item2]=1;
                rand_neighb->objective += delta;
                for(int d=0; d<prob->dim; d++){
                    rand_neighb->cap_left[d] +=  prob->items[item1].size[d] - prob->items[item2].size[d];
                }
            }
            copy_solution(curt_sln, rand_neighb);
        }
        iter++;
    }
    free_solution(rand_neighb);
    free(rand_neighb);
}

bool can_swap(struct solution_struct* sln, int out, int in){
    for(int d =0; d<sln->prob->dim; d++)
    {
        if(sln->cap_left[d]+sln->prob->items[out].size[d] < sln->prob->items[in].size[d])
            return false;
    }
    return true;
}

bool can_move(int nb_indx, int* move, struct solution_struct* curt_sln ){
    bool ret=true;
    if(nb_indx==3){//3-item swap
        int i= move[0], j= move[1], k= move[2];
        if(i<0 || j<0 || k<0) return false;
        if(curt_sln->x[j]>0) {//2-1 swap
            for(int d=0; d<curt_sln->prob->dim; d++){
                if(curt_sln->cap_left[d] + curt_sln->prob->items[i].size[d] +
                   curt_sln->prob->items[j].size[d] < curt_sln->prob->items[k].size[d])
                    return false;
            }
        }
        else {//1-2 swap
            for(int d=0; d<curt_sln->prob->dim; d++){
                if(curt_sln->cap_left[d] + curt_sln->prob->items[i].size[d] <
                   curt_sln->prob->items[j].size[d] + curt_sln->prob->items[k].size[d])
                    return false;
            }
        }
        
    }
    else ret=false;
    return ret;
}

bool apply_move(int nb_indx, int* move, struct solution_struct* sln ){
    bool ret=true;
    if(nb_indx==3){//3-item swap
        int i= move[0], j= move[1], k= move[2];
        if(i<0 || j<0 || k<0) return false;
        if(sln->x[j]>0) {//2-1 swap
            for(int d=0; d<sln->prob->dim; d++){
                sln->cap_left[d] = sln->cap_left[d]+sln->prob->items[i].size[d] +
                    sln->prob->items[j].size[d] - sln->prob->items[k].size[d];
            }
            sln->objective += sln->prob->items[k].p-sln->prob->items[i].p-sln->prob->items[j].p;
            sln->x[i]=0; sln->x[j]=0; sln->x[k]=1;
        }
        else {//1-2 swap
            for(int d=0; d<sln->prob->dim; d++){
                sln->cap_left[d] = sln->cap_left[d]+sln->prob->items[i].size[d] -
                    sln->prob->items[j].size[d] - sln->prob->items[k].size[d];
            }
            sln->objective += sln->prob->items[j].p+sln->prob->items[k].p-sln->prob->items[i].p;
            sln->x[i]=0; sln->x[j]=1; sln->x[k]=1;
        }
        
    }
    else ret=false;
    return ret;
}

void unequal_swap(int nb_indx,struct solution_struct * sln){
    //2-1 swap
    int curt_move[]={-1,-1,-1}, best_move []={-1,-1,-1}, delta=0, best_delta=0;
   /* struct solution_struct* best_neighb = malloc(sizeof(struct solution_struct));
       best_neighb->cap_left = malloc(sizeof(int)*sln->prob->dim);
       best_neighb->x = malloc(sizeof(int)*sln->prob->n);
       copy_solution(best_neighb, sln);*/
    
    for(int i=0; i<sln->prob->n; i++){
        if(sln->x[i]==0) continue;
        for(int j=0; j!=i&&j<sln->prob->n; j++){
            if(sln->x[j]==0) continue;
            for(int k=0;k<sln->prob->n;k++){
                if(sln->x[k] == 0){
                    curt_move[0]=i; curt_move[1]=j; curt_move[2]=k;
                    if(can_move(nb_indx, &curt_move[0], sln)){
                        delta = sln->prob->items[k].p -sln->prob->items[i].p-sln->prob->items[j].p;
                        if(delta > best_delta){
                            best_delta = delta; best_move[0] = i; best_move[1] = j; best_move[2]=k;
                        }
                    }
                }
            }
        }
    }
    //1-2 swap
    for(int i=0; i<sln->prob->n; i++){
        if(sln->x[i]==0) continue;
        for(int j=0; j<sln->prob->n; j++){
            if(sln->x[j]>0) continue;
            for(int k=0;k!=j&&k<sln->prob->n;k++){
                if(sln->x[k] == 0)
                {
                    curt_move[0]=i; curt_move[1]=j; curt_move[2]=k;
                    if(can_move(nb_indx, &curt_move[0], sln)){
                        delta = sln->prob->items[k].p +sln->prob->items[j].p-sln->prob->items[i].p;
                            if(delta > best_delta){
                                best_delta = delta; best_move[0] = i; best_move[1] = j; best_move[2]=k;
                            }
                    }
                }
            }
        }
    }
    if(best_delta>0) { apply_move(nb_indx, &best_move[0], sln);}
    update_best_solution(sln);
}

// uniform crossover
void mate(struct solution_struct* c1, struct solution_struct* c2){
	int i, d, temp;
	for (i=0; i<c1->prob->n; i++){
		
		if (rand_01()<CROSSOVER_RATE){
			temp = c1->x[i];
			c1->x[i] = c2->x[i];
			c2->x[i] = temp;		// swap
			
			//update objectives and capacities
			if (c1->x[i] == 0 && c2->x[i] == 1){
				c1->objective -= c1->prob->items[i].p;
				c2->objective += c2->prob->items[i].p;
				
				for(d=0; d<c1->prob->dim; d++){
                    c1->cap_left[d] += c1->prob->items[i].size[d];
                    c2->cap_left[d] -= c2->prob->items[i].size[d];
                }
			} else if (c1->x[i] == 1 && c2->x[i] == 0){
				c1->objective += c1->prob->items[i].p;
				c2->objective -= c2->prob->items[i].p;
				
				for(d=0; d<c1->prob->dim; d++){
                    c1->cap_left[d] -= c1->prob->items[i].size[d];
                    c2->cap_left[d] += c2->prob->items[i].size[d];
                }
			}
		}
	}
	c1->flag = 0;
	c2->flag = 0;
}

// random mutation (based on mutation rate)
void mutate(struct solution_struct* child){
	int i, temp;
	for (i=0; i<child->prob->n; i++){
		if (rand_01()<MUTATION_RATE){
			// mutation occurs
			if (child->x[i] == 0){
				child->x[i] == 1;
			} else {
				child->x[i] == 0;
			}
		}
	}
}

// heuristically repair infeasible solutions 
void repair(struct solution_struct* child){
	int pointer = 0;
	double rate = 1000000.0;
	bool flag = false;
	for(int dim=0; dim<child->prob->dim; dim++){
		if(child->cap_left[dim]<0){
			flag = true;
			break;
		}
	}
	while(flag){
		for(int i=0; i<child->prob->n; i++){
			if((child->x[i] == 1) && (child->prob->items[i].ratio<rate)){
				pointer = i;
				rate = child->prob->items[i].ratio;
			}
		}	// choose the worst item which is currently selected
		child->x[pointer] = 0;
		child->objective -= child->prob->items[pointer].p;
		flag = false;
		for(int dim=0; dim<child->prob->dim; dim++){
			child->cap_left[dim] += child->prob->items[pointer].size[dim];
			if(child->cap_left[dim]<0){
				flag = true;
			}
		}
		pointer = 0;
		rate = 1000000.0;
	}
} 

// add as many items as possible one by one based on their evaluations within capacity contraints
void heuristic_add(struct solution_struct * sln){
    int j = 0;
    for(int i = 0; i < sln->prob->n; i++){
        if(sln->x[i] == 0){
            //check this gene if can be added.
            for( j = 0; j <sln->prob->dim; j++){
                if(sln->cap_left[j]<sln->prob->items[i].size[j]){
                    break;
                }
                
            }
            if(j == sln->prob->dim){
                sln->x[i] = 1;
                for(int k = 0; k <sln->prob->dim; k++){
                    sln->cap_left[k]-= sln->prob->items[i].size[k];
                }
                sln->objective += sln->prob->items[i].p;
            }
        }
    }
}

// replace part of the individuals in the original population with the generated children
void update(struct solution_struct** population, struct solution_struct** children){
		struct solution_struct* rand_sln;
		struct solution_struct* temp;
		int i,j;
		// increasingly sort population
		for (i=0; i<POP_SIZE-1; i++){
			for (j=0; j<POP_SIZE-1-i; j++){
				if (population[j]->objective > population[j+1]->objective){
					temp = population[j];
					population[j] = population[j+1];
					population[j+1] = temp;
				}
			}
		}
		
		// increasingly sort children
		for (i=0; i<POOL_SIZE-1; i++){
			for (j=0; j<POOL_SIZE-1-i; j++){
				if (children[j]->objective > children[j+1]->objective){
					temp = children[j];
					children[j] = children[j+1];
					children[j+1] = temp;
				}
			}
		}
		
		// compare the individuals in population and children correspondingly
		// only the winning children could be pushed into population
		for (i=0; i<POOL_SIZE; i++){
			if (children[i]->objective > population[i]->objective){
				copy_solution(population[i], children[i]);
				children[i]->flag = 0; 
			}
		}
		
		// eliminate duplicated solutions in the updated population
		for(i=0; i<POP_SIZE; i++){
	        for(j=0; j <POP_SIZE; j++){
	            if(population[i]->objective == population[j]->objective && i!=j){
	                rand_sln = create_empty_sol(population[i]->prob);
					random_sol(rand_sln);	// replace one of the parents with a random one
					heuristic_add(rand_sln);
					temp = population[i];
					population[i] = rand_sln;
					rand_sln = temp;
					free_solution(temp); free(temp);
	            }
	        }
	    }	
}

// memetic algorithm
int memetic(struct problem_struct* prob){
	int i, j, k, counter;	
	// data structure defination (pointer array)
	struct solution_struct* population[POP_SIZE];
	struct solution_struct* mating[POOL_SIZE];
	struct solution_struct* children[POOL_SIZE];
	
	// timer defination (gettimeofday)
	struct timeval start;
	struct timeval end;
	int time_spent = 0;
	gettimeofday( &start, NULL );     
	
	// 1. generate initial population randomly
	for (i=0; i<POP_SIZE; i++){
		population[i] = create_empty_sol(prob);
		random_sol(population[i]);
		if (rand_01()<LS_RATE1){
			first_descent(population[i]);	// local search for some of the initial individuals
		}
	}
	
	// initialization for mating and children pool
	for (i=0; i<POOL_SIZE; i++){
		mating[i] = create_empty_sol(prob);
		children[i] = create_empty_sol(prob);
	}
	
	counter = 0; 
	while (counter<MAX_NUM_OF_GEN && time_spent<MAX_TIME){
		j = 0;
		k = 0;
		int winner, max, current;
		max = 0;
		winner =-1;
		
		// 2. select individuals from population pool to mating pool (tournament selection)
		for (i=0; i<POOL_SIZE; i++){
			for (j=0; j<tourm_size; j++){
				current = rand_int(0, POP_SIZE-1);
				if (population[current]->objective > max){
					winner = current;
					max = population[current]->objective;
				}
			}
			copy_solution(mating[k], population[winner]);	// push the winner into mating pool
			k++;			
			max = 0;
			winner = -1;
		}
		
		// 3. pair-wise mate (uniform crossover)
		j = 0;
		k = 0;
		struct solution_struct* rand_sln;
		struct solution_struct* temp;  
		for (i=0; i<(POOL_SIZE/2); i++){
			int p1 = rand_int(0, POOL_SIZE-1);
			int p2;
			do {
				p2 = rand_int(0, POOL_SIZE-1);
			} while (p2 == p1); 	// randomly select two parents from the 
						
			if (fabs(mating[p1]->objective - mating[p2]->objective) < 0.10){	// two selected parents are identical
				rand_sln = create_empty_sol(mating[p2]->prob);
				random_sol(rand_sln);	// replace one of the parents with a random one
				temp = mating[p2];
				mating[p2] = rand_sln;
				rand_sln = temp;
				free_solution(temp); free(temp);
				heuristic_add(mating[p2]);
				copy_solution(children[2*i], mating[p1]);
				copy_solution(children[2*i+1], mating[p2]);	// push children into children pool
				mate(children[2*i], children[2*i+1]);	// mating 
			} else {
				copy_solution(children[2*i], mating[p1]);
				copy_solution(children[2*i+1], mating[p2]);	// push children into children pool
				mate(children[2*i], children[2*i+1]);	// mating
			}
						
			// 4. mutation of children
			mutate(children[2*i]);
			mutate(children[2*i+1]);
			
			// 5. feasibility repair
			repair(children[2*i]);
			repair(children[2*i+1]);

			// 6. local search for each child (a part of children)
			if (rand_01() < LS_RATE2 && counter%100 == 0 && counter>0){
				if (true){
					first_descent(children[2*i]);
			
				}
				if (true){
					first_descent(children[2*i+1]); 
				} 
			}
			if (rand_01() < VNS_RATE && counter%100 == 0 && counter>0){
				if (children[2*i]->flag == 0){
					unequal_swap(3, children[2*i]);
				}
				if (children[2*i+1]->flag == 0){
					unequal_swap(3, children[2*i+1]); 
				}
			}
			
		}
		
		// 7. update population
		update(population, children);
		update_best_solution(population[POP_SIZE-1]);	// update current best solution		
		
		// update timer and counter
		counter++;
		gettimeofday( &end, NULL );  
		time_spent = end.tv_sec - start.tv_sec; 
		
		// ** generation finishes **
		printf("iteration: %d, current optimal objective: %f\n", counter, best_sln.objective);
	}
	
	// free memories (when memetic algorithm terminates)
	for(int k=0; k<POP_SIZE; k++)
    {
       free_solution(population[k]);
       free(population[k]);	// free population data structure 
    }
    
    for(int k=0; k<POOL_SIZE; k++)
    {
       free_solution(mating[k]);
       free_solution(children[k]);
       free(mating[k]);	// free mating pool
       free(children[k]);	// free children pool
    }
}

int main(int argc, const char * argv[]) {

    printf("Starting the run!\n");
    char data_file[50]={"somefile"}, out_file[50]={}, solution_file[50]={};  //max 50 problem instances per run
    if(argc<3)
    {
        printf("Insufficient arguments. Please use the following options:\n   -s data_file (compulsory)\n   -o out_file (default my_solutions.txt)\n   -c solution_file_to_check\n   -t max_time (in sec)\n");
        return 1;
    }
    else if(argc>9)
    {
        printf("Too many arguments.\n");
        return 2;
    }
    else
    {
        for(int i=1; i<argc; i=i+2)
        {
            if(strcmp(argv[i],"-s")==0)
                strcpy(data_file, argv[i+1]);
            else if(strcmp(argv[i],"-o")==0)
                strcpy(out_file, argv[i+1]);
            else if(strcmp(argv[i],"-c")==0)
                strcpy(solution_file, argv[i+1]);
            else if(strcmp(argv[i],"-t")==0)
                MAX_TIME = atoi(argv[i+1]);
        }
        printf("data_file= %s, output_file= %s, sln_file=%s, max_time=%d", data_file, out_file, solution_file, MAX_TIME);
    } 
    struct problem_struct** my_problems = load_problems(data_file); 
    if(strlen(solution_file)<=0)
    {
        if(strcmp(out_file,"")==0) strcpy(out_file, "my_solutions.txt"); //default output
        FILE* pfile = fopen(out_file, "w"); //open a new file
        fprintf(pfile, "%d\n", num_of_problems); fclose(pfile); 
        for(int k=0; k<num_of_problems; k++)
        {
            best_sln.objective=0; best_sln.feasibility=0;
            for(int run=0; run<NUM_OF_RUNS; run++)
            {
                srand(RAND_SEED[run]);
                
                memetic(my_problems[k]); // call memetic method
            }

            output_solution(&best_sln,out_file);
//			for(int i=0; i<my_problems[k]->n; i++){
//				printf("%d ", my_problems[k]->items[i].p);
//			}
//			printf("\n");
        }
    }
    
    for(int k=0; k<num_of_problems; k++)
    {
       free_problem(my_problems[k]); //free problem data memory
    }
    free(my_problems); //free problems array
    if(best_sln.x!=NULL && best_sln.cap_left!=NULL){ free(best_sln.cap_left); free(best_sln.x);} //free global
    
    printf("MEMETIC FINISHES ***\n");
    return 0;
}
