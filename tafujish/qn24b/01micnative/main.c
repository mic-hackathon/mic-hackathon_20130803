/**** XeonPhi NativeMode + OpenMP  by fujish ****/
/************************************************/
/**** N-queens since 2004-09   by Kenji KISE ****/
/**** N-queens OpenMP Version in C           ****/
/************************************************/
#include <stdio.h>
#include <stdlib.h>
#include <sys/time.h>

/************************************************/
#define NAME "qn24b XeonPhiNative"
#define VER  "version 1.0.0 2013-08-03"
#define MAX 29 /** 32 is a real max! **/
#define MIN 2
#define MAX_TASK 100000

/************************************************/
typedef struct array_t{
  unsigned int cdt; /* candidates        */
  unsigned int col; /* column            */
  unsigned int pos; /* positive diagonal */
  unsigned int neg; /* negative diagonal */
} array;

/** N-queens kernel                            **/
/** n: problem size, h: height, r: row         **/
/************************************************/
long long qn(int n, int h, int r, array *a){
  long long answers = 0;

  for(;;){
    if(r){
      int lsb = (-r) & r;
      a[h+1].cdt = (       r & ~lsb);
      a[h+1].col = (a[h].col & ~lsb);
      a[h+1].pos = (a[h].pos |  lsb) << 1;
      a[h+1].neg = (a[h].neg |  lsb) >> 1;
      
      r = a[h+1].col & ~(a[h+1].pos | a[h+1].neg);
      h++;
    }
    else{
      r = a[h].cdt;
      h--;

      if(h==0) break;
      if(h==n) answers++;
    }
  }
  return answers;
}

/** function to return the time                **/
/************************************************/
long long get_time(){
  struct timeval  tp;
  struct timezone tz;
  gettimeofday(&tp, &tz);
  return tp.tv_sec * 1000000ull + tp.tv_usec;
}

/** print the answer                           **/
/************************************************/
void print_result(int n, long long usec,
                  long long solution){
  int i;
  static long long answer[30] = {
    0, 1, 0, 0, 2, 10, 4, 40, 92, 352, 724,
    2680, 14200, 73712, 365596, 2279184,
    14772512, 95815104, 666090624ull,
    4968057848ull, 39029188884ull,
    314666222712ull, 2691008701644ull,
    24233937684440ull,
    0,0,0,0,0,0
  };

  for(i=0;i<45;i++) printf("="); printf("\n");
  printf("%s %s\n", NAME, VER);
  printf("problem size n        : %d\n", n);
  printf("total   solutions     : %lld\n",
         solution);
  printf("correct solutions     : %lld\n",
         answer[n]);
  printf("million solutions/sec : %.3f\n",
         (double)solution/usec);
  printf("elapsed time (sec)    : %.3f\n",
         usec/1000000.0);
  if(solution!=answer[n])
    printf(" ### Wrong answer\n");
  for(i=0;i<45;i++) printf("="); printf("\n");
  fflush(stdout);
}

/** set the problem size                       **/
/************************************************/
int set_problem_size(int argc, char** argv){
  int n; /* problem size */
  printf("%s %s\n", NAME, VER);

  if(argc!=2){
    printf("Usage: %s n\n", argv[0]);
    exit(0);
  }
  n = atoi(argv[1]);
  if(MIN>n || MAX<n){
    printf("board size range: %d-%d\n",MIN, MAX);
    exit(0);
  }
  return n;
}

/** n-queens solver, specify first 4 level     **/
/************************************************/
long long queen_sub(int n, int *in, int mode){
  int i;
  long long solution=0;   /* solutions          */
  unsigned int row;       /* row candidate      */
  unsigned int h;         /* height or level    */
  array a[MAX]; 

  for(i=0; i<MAX; i++){
    a[i].cdt = a[i].col = a[i].pos = a[i].neg = 0;
  }
  row      = (1<<n)-1;
  a[1].col = (1<<n)-1;
  a[1].pos = 0;
  a[1].neg = 0;
  a[1].cdt = 1; /* care */
  a[2].cdt = 0; /* care */
  
  for(h=1; h<=4; h++){  /** Initialize **/
    int lsb = (1<<in[h]);

    if((a[h].col & lsb)==0) return 0;
    if((a[h].pos & lsb)!=0) return 0;
    if((a[h].neg & lsb)!=0) return 0;

    a[h+1].col = (a[h].col & ~lsb);
    a[h+1].pos = (a[h].pos |  lsb) << 1;
    a[h+1].neg = (a[h].neg |  lsb) >> 1;
    row = a[h+1].col & ~(a[h+1].pos | a[h+1].neg);
  }

  if(mode==0) return 1;

  solution = qn(n, h, row, a);
  return solution;
}

/** get the number of sub problems             **/
/************************************************/
int get_sub_problem_num(int n, int *prob){
  int i;
  int problem=0;
  int max = ((n>>1) + (n&1)) * n * n * n;

  for(i=0; i<=max; i++){
    int in[MAX];
    in[1] = (i-1) / n / n / n;
    in[2] = (i-1) / n / n % n;
    in[3] = (i-1) / n % n;
    in[4] = (i-1) % n;
    
    if(queen_sub(n, in, 0)){
      problem++;
      prob[problem] = i;
    }
  }
  printf("There are %d tasks\n", problem);
  if(problem > MAX_TASK){
    printf("The number of tasks is too large.\n");
    exit(1);
  }
  return problem;
}

/************************************************/
long long solver(int n, int id,
                 int problem, int *prob){
  long long ret;
  int p_no = prob[id];
  int in[MAX];

  if(id<1 || id>problem){
    printf("The id %d is out of range.\n", id);
    return 0;
  }

  in[1] = (p_no-1) / n / n / n;
  in[2] = (p_no-1) / n / n % n;
  in[3] = (p_no-1) / n % n;
  in[4] = (p_no-1) % n;
  
  ret = queen_sub(n, in, 1);
  if(in[1]<(n>>1)) ret = ret * 2;

  return ret;
}

/** main function                              **/
/************************************************/
int main(int argc, char *argv[]){
  int prob[MAX_TASK]; /* store problem IDs */
  int rets[MAX_TASK]; /* store each answer */
  int i;
  int n    = set_problem_size(argc, argv);
  int jobs = get_sub_problem_num(n, prob);
  long long usec    = get_time();
  long long answers = 0;

#pragma omp parallel for schedule(dynamic)
  for(i=1; i<=jobs; i++){
    rets[i] = solver(n, i, jobs, prob);
  }

  for(i=1; i<=jobs; i++) answers += rets[i];
  print_result(n, get_time()-usec, answers);
  return 0;
}
/************************************************/
