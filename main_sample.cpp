
#include<limits.h>
#include<time.h>
#include<assert.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>
#include <string.h>
#include <sys/stat.h>
#include <sys/types.h>
#include <algorithm>
#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <map>
#include <vector>
#include <cstdio>
#include <errno.h>  // Para errno
//#include <string.h> // Para strerror

double rand_easy(double a, double b); //Generate a uniformly distribute random number in [a, b] using rand()
double Prob_Connection(double kappa1, double kappa2, double tetha1, double tetha2, double B, double mu, double Rad); //it compute the connection probability of two nodes
double Prob_Connection_C(double kappa1nc, double kappa2c, double mu_C);
double Random_powerlaw(double gamma, double k0);  //It picks a sample from a given power-law distribution
void read_file(const char *file_name , int * Nodes , int * Links); //It reads a file
double get_fun (double r, double k1, double k2, double k1_min, double k2_min , double gamma1, double gamma2, double nu); //This function computes the Eq.14 in the supplementary of "10.1038/NPHYS3812" paper to
                                                                                                                         //generate kappas in S1*S1 model which are correlated with those in the S1 model
double bisect (double r, double x1, double x2, double err,double k1, double k1_min, double k2_min , double gamma1, double gamma2, double nu); //Bisection method to find zeros of a function
double get_fun_truncated (double r, double k1, double k2, double k1_min, double k2_min , double k1_c, double k2_c, double gamma1, double gamma2, double nu);
double bisect_truncated (double r, double x1, double x2, double err,double k1, double k1_min, double k2_min, double k1_c, double k2_c, double gamma1, double gamma2, double nu);
void swap(int* xp, int* yp);
void selectionSort(int arr[], int n);
double Delta_Theta(double Theta1, double Theta2);
void Assign_Labels(double * Theta_s, int NC , int * Label, int NZ, float R, float alpha);
/// Defines ///////////////////////////////////////////////////////////////////////////////////////
#define	sqr(x)		((x)*(x))			    // equis al cuadrado
#define _for(a,b,c)	for(a=(b); a<(c); ++a)
#define ULONGMAXX	4294967296.0			// unsigned long max + 1

long long factorial(long long x) { if(x==0LL) return 1LL; else return x*factorial(x-1LL); }


int L;
#define MM       	sqr(L)		            // numero de nodos (debe ser par)

#define AVERAGES 100

#define kmax        4000


#define prename     "salida"

using namespace std;


/// Tipos /////////////////////////////////////////////////////////////////////////////////////////

struct tipo_nodes
{
		float kappa;
		float theta;
		float clustering;
		int degreec;
		int degree; 				// numero de links (grado)
		int link[kmax];		// lista de nodos a los que esta linkeado

};


struct tipo_features_s1
{
        float kappa;
		float theta;
		float clustering;
		int degreec;
		int degree; 				// numero de links (grado)
		int link[kmax];		// lista de nodos a los que esta linkeado

};

struct tipo_features_c
{
        float kappa;
        float theta;
        float clustering;
		int degree; 				// numero de links (grado)
		int degreec;
		int link[kmax];		// lista de nodos a los que esta linkeado

};




struct tipo_features
{
		float kappa;
		float theta;
		float clustering;
		int degreec;
		int degree; 				// numero de links (grado)
		int link[kmax];		// lista de nodos a los que esta linkeado

};


struct tipo_links
{
		int node;
		int feat;
};


typedef double flt;







/// Variables globales ////////////////////////////////////////////////////////////////////////////

int M;                          // numero total de nodos

			    // cada nodo con su lista de vecinos y numero de vecinos (ahora asigno dinamico)
tipo_nodes *nodes;
tipo_features *features;
tipo_features_c *f_c;
tipo_features_s1 *f_s1;


int *bolsa;			            // bolson de patas (usado para generar la red)
int nbolsa;						// numero actual de patas en bolsa


int *labels;					// vector para re-rotular los labels
int *masas;					    // vector de masas de clusters


unsigned int seed[256];		    // cosas del generador pseudorandom
unsigned int r;
unsigned char irr;



char nombre1[200];
char nombre2[200];
int masa_max,label_masa_max;


int C[kmax];                    // kupla tortuosa (lista de nodos)
int *K;                         // lista de kuplas tortuosas obtenida por busqueda sistematica
int NK,NKMAX;                   // numero de kuplas en la lista y tamaño actual de la lista (puede crecer con realloc)


char buf[200];



#define n_tiempos   1000

flt tiempos[n_tiempos];
flt promedios_tiempos[n_tiempos];
int muestras_tiempos[n_tiempos];



void inicializa_randomm(void)
{
	int i;
// inicializo generador pseudorandom: uso semilla = segundos desde 1970,
// para que genera siempre una secuencia distinta
	srand((unsigned)time(0));
	irr=1;
	_for(i,0,256) seed[i]=rand();
	r=seed[0];
	_for(i,0,70000) r=seed[irr++]+=seed[r>>24];

}


flt randomm(void)
{
// esto genera un nuevo numero pseudoaleatorio
// uniformemente distribuido en el intervalo [0,1)
	r=seed[irr++]+=seed[r>>24];
	return ( (flt)r/ULONGMAXX );
}



int main(int n_args,char *args[])
{
    int av, parameters;
    float Beta_s;
    float gamma_s;
    float gamma_f;
    float gamma_f_c;
    int N_f_c;
    int N_f;
    int N_f_obs;
    double sigma[kmax];

printf("0~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n");


     float BETA[8]   =  {1.5, 3.0, 1.5, 3.0, 1.5, 3.0, 1.5, 3.0};
    float GAMMAC[8] =  {2.5, 2.5, 2.2, 2.2, 3.0, 3.0, 2.2, 2.2};
    float GAMMAS[8] =  {3.0, 3.0, 3.0, 3.0, 2.2, 2.2, 3.0, 3.0};

    int N_NODES[8] = {1000,1000, 1000, 1000, 1000, 1000, 500, 500};
    int N_FF[8] =    {500, 500, 500, 500, 500, 500, 250, 250};
    int N_F_CC[8] =  {500, 500, 500, 500, 500, 500, 250, 250};

printf("1~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n");
 for(parameters = 0; parameters < 6; parameters++)
    {
        Beta_s = BETA[parameters];
        gamma_s = 5.0;
        gamma_f = GAMMAS[parameters];
        gamma_f_c = GAMMAC[parameters];

        N_f_obs = N_FF[parameters];
        N_f = N_FF[parameters];
        N_f_c = N_F_CC[parameters];

        FILE * fp1;
        char file_fp1[200];
        sprintf(file_fp1, "sample_size_percok6_%1d.dat",parameters);
        fp1 = fopen (file_fp1, "a"); if (fp1 == NULL) { perror("Failed: "); return 1; }

        FILE * fp2;
        char file_fp2[200];
        sprintf(file_fp2, "coord_sample_size_perco_k6_%1d_b%.1f_gs1%.1f_gcm%.1f.dat",parameters,BETA[parameters], gamma_f, gamma_f_c);
        fp2 = fopen (file_fp2, "a"); if (fp2 == NULL) { perror("Failed: "); return 1; }

        fprintf(fp2,"type\t label\t kappa\t theta\n");

//~~~~~~~~~~~~~~~~~~~~ Parameters of S1_S1 model
        int Ns_obs = N_NODES[parameters]; // The desired number of observed nodes in the resulting network
        int Ns = N_NODES[parameters];  //The number of nodes in the model
        float kmean_s = 6.0; //Average degree
        float Corr_k0_s = (1-(1/Ns_obs))/(1-pow(Ns_obs, ((2-gamma_s)/(gamma_s-1)))); //Correction of k0_s
        float k0_s = (((gamma_s-2) * kmean_s)/(gamma_s-1))* Corr_k0_s; //Minimum kappa
        float kc_s = k0_s * (pow(Ns_obs, (1/(gamma_s-1))));    //Maximum kappa
        int delta_s = 1;  //Density of nodes
        float R_s = Ns_obs / (2 * M_PI * delta_s); //The Radius of S1 model
        float mu_s = (Beta_s * sin(M_PI/Beta_s)) / (2 * delta_s * kmean_s * M_PI);  //controls the average degree
        float alpha_s = 1.0;
        int NC = 6;

        //~~~~Nodes of Type 1 (n)
        float gamma_n ;
        gamma_n = gamma_s;
        int N_n;
        N_n = Ns_obs;
        int N_n_obs;     // number of nodes
        float kmean_n = 3.0;  // average degree of nodes of type 1 (Nodes)
        float Corr_k0_n; // Correction of k0_n
        float k0_n; // Minimum kappa of nodes of Type 1
        float kc_n; // Maximum kappa of nodes of Type 1
        int delta_n = 1;  //Density of nodes of type 1

        float kmean_f;  //average degree of nodes of type 2 (Features s1) kmeans_n * N_n = K_means_f * N_f
        float Corr_k0_f; //correction of k0_f
        float k0_f;   //Minimum kappa of nodes of Type 2
        float kc_f; //Maximum kappa of nodes of Type 2
        float kmean_n_s = kmean_n; // contribution to the average degree of the nodes from the S1 connections

        float Beta_bi;
        Beta_bi = Beta_s;
        float R_bi ; //The Radius of S1*S1 model
        float mu_bi = (Beta_bi * sin(M_PI/Beta_bi)) / (2.0* M_PI * kmean_n_s);  //controls the average degree

        //~~~~Parameter of correlation between kappas in S1 and kappas in S1*S1
        float Err = 0.001;  //Acceptable error of bisection method
        float nu = 0.0;   //nu in [0, 1] is the correlation strength parameter
        // float nu = atof(argv[10]);   //nu in [0, 1] is the correlation strength parameter


        //~~~~Nodes of Type 3 (features connected via configuration model)

        float kmean_n_c = kmean_s - kmean_n;  // average degree of nodes of type 1 (Nodes)
        float mu_c = 1.0 / (kmean_n_c * (float)(N_n));

        //~~~~expected average degree of the network of nodes
        float kmean_nodes = 1.0 / ( mu_c * (float)(N_n)) +(Beta_bi * sin(M_PI/Beta_bi)) / (2* M_PI * mu_bi);
       // printf("desired mean degree = %f\t kmean_nodes = %f\n",kmean_s, kmean_nodes);


        int N_features = N_f + N_f_c;
        int links, E, links_features[N_features]; //this vector has the number of links/degree for all features

        int coin_feat_c, coin_feat_s1, aa;
        int success_selected, success;
        double  kmean_f_c;
        int i,j,k,nk;
        float KMC[kmax][2],KMS[kmax][2];
        float H[kmax][AVERAGES],L[kmax][AVERAGES];
        int SUP[kmax];



        //~~~Set the seed for random number generator
        time_t t;
        srand((unsigned) time(&t));




printf("2~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n");

for(av = 0; av < 1; av++)
{

        tipo_links* links;
        int N_l = kmean_n * N_NODES[parameters] * 1.15;

        // Asignación de memoria dinámica usando new
        links = new tipo_links[N_l];
        assert(links != NULL);

        nodes=(tipo_nodes *)malloc(N_n*sizeof(tipo_nodes)); assert(nodes!=NULL);
        f_s1=(tipo_features_s1 *)malloc(N_f*sizeof(tipo_features_s1)); assert(f_s1!=NULL);
        f_c=(tipo_features_c *)malloc(N_f_c*sizeof(tipo_features_c)); assert(f_c!=NULL);
        features=(tipo_features *)malloc(N_features*sizeof(tipo_features)); assert(features!=NULL);

        _for(j ,0,N_n){_for(i,0,kmax) nodes[j].link[i] = -1;}
        _for(j ,0,N_f){_for(i,0,kmax) f_s1[j].link[i] = -1;}
        _for(j ,0,N_f_c){_for(i,0,kmax) f_c[j].link[i] = -1;}

        _for(j ,0,N_f) f_s1[j].degree = 0;
        _for(j ,0,N_f_c) f_c[j].degree = 0;
        _for(i ,0,kmax) _for(j,0,2) {KMC[i][j] = 0.0; KMS[i][j] = 2.0;}  //horrible pero quiero ver si lo termino al menos con chapuza
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ Generating kappas
        int iter=0;
        double kappa;


        while(iter < N_n)
        {
            kappa = Random_powerlaw(gamma_s, k0_s);

            if(kappa <= kc_s && kappa >= k0_s)
            {
                nodes[iter].kappa = kappa;
                iter=iter+1;

            }

        }
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~Parameter Tetha/ S1 Model


    for(int i = 0; i < N_n; i++)
    {
       double e = 2* M_PI;

       nodes[i].theta = rand_easy(0, e);
    }


    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ Kappas of bipartite network for nodes of Type 2: features connected via S1
     kmean_f = kmean_n * N_n / N_f_obs;//kmean_f =  (kmean_n * N_n) / N_f_obs; //This is the mean degree of the features connected via S1. 10 is the mean degree of the nodes if the edges were only the ones from the S1 (therefore the mean degree I put here is not the total mean deagree of nodes)
     Corr_k0_f = (1-(1/N_f_obs))/(1-pow(N_f_obs, ((2-gamma_f)/(gamma_f-1))));
     k0_f = (((gamma_f-2) * kmean_f)/(gamma_f-1)) * Corr_k0_f;
     kc_f = k0_f * (pow(N_f_obs, (1/(gamma_f-1))));
     N_f = 0;

     double k_f;

     N_f = N_f_obs;

//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~Kappas of bipartite network for nodes of Type 2: features connected via S1
         int iter2 = 0;
         float kappa2;

         while(iter2 < N_f)
         {
             kappa2 = Random_powerlaw(gamma_f, k0_f);

             if(kappa2 <= kc_f && kappa2 >= k0_f)
             {
                 f_s1[iter2].kappa = kappa2;
                 iter2= iter2+1;

             }
         }

//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~Kappas of bipartite network for nodes of Type 3: features connected via CM
         kmean_f_c = ( kmean_n_c * N_n ) / ( N_f_c );
         double Corr_k0_f_c = (1 - ( 1/N_f_c ) )/( 1 - pow(N_f_c, ( (2-gamma_f_c )/(gamma_f_c -1) ) ) );
         double k0_f_c = ( ( (gamma_f_c - 2) * kmean_f_c ) / ( gamma_f_c -1) ) * Corr_k0_f_c;
         double kc_f_c = k0_f_c * ( pow(N_f_c, ( 1 / ( gamma_f_c - 1) ) ) );



         iter2 = 0;
         float kappa2_c;

         for(i = 0; i < N_f_c; i++)
         {
             kappa2_c = Random_powerlaw(gamma_f_c, k0_f_c);

             if(kappa2_c <= kc_f_c && kappa2_c >= k0_f_c)
             {

                 f_c[i].kappa = kappa2_c;

             }
         }

    //~~~~~~~~~~~~~~~~~~~Thetas of bipartite network for nodes of Type 2

    for(i = 0; i < N_f; i++)
    {
       double e = 2 * M_PI;
       f_s1[i].theta = rand_easy(0, e);
    }

    for(i = 0; i < N_f_c; i++)
    {
       double e = 2 * M_PI;
       f_c[i].theta = rand_easy(0, e);
    }


    printf("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n");

    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~Compute the probability of connection between nodes in S1 * S1 and generate the bipartite network
    printf("Playing with Mu to get a Bipartite network with a desired average degree....!\n");

    R_bi = N_n / (2 * M_PI * delta_n); //The Radius of S1_S1 model
    mu_bi= (Beta_bi * sin(M_PI/Beta_bi)) / (2 * delta_n * kmean_s * M_PI);
    int ii, iii, jj, jjj;
    int it = 0, degree;
    printf("\n Mu is : %f",  kmean_s);
    double AVG_DEG_1, AVG_DEG_2;
    double Num_links, Num_links_c, Num_links_s1;

    while(1)
    {
       // printf("\n itera # %d", it);
        //printf("\n Mu is : %f", mu_bi);



        Num_links_s1 = 0;

        for(ii = 0; ii < N_n ; ii++)
        {
            degree = 0;
            for(jj = 0 ; jj < N_f; jj++)
            {
                double prob = Prob_Connection(nodes[ii].kappa, f_s1[jj].kappa, nodes[ii].theta,f_s1[jj].theta, Beta_bi, mu_bi, R_bi);

                double r = rand_easy(0, 1);
                if (r <= prob)
                {
                    nodes[ii].link[degree] = jj;
                    f_s1[jj].link[f_s1[jj].degree] = ii;
                    degree++;
                    f_s1[jj].degree++;

                    Num_links_s1 = Num_links_s1 + 1;

                }
            }
            nodes[ii].degree = degree;
        }
       // printf("links s1 = %f\n",Num_links_s1);


        AVG_DEG_1 = Num_links_s1/N_n;
        AVG_DEG_2 = Num_links_s1/N_f;
        //printf("\n The av. degree of nodes of Type 1 s1 is: %f \n", AVG_DEG_1);
        //printf(" The av. degree of nodes of Type 2 (features) s1 is: %f \n", AVG_DEG_2);
       // printf("~~~~~~~~~~~");

        if( fabs(AVG_DEG_1-kmean_n_s) < 0.05)
        {
                    break;
        }

        if(it>1000){ break;}

        if(AVG_DEG_1 < kmean_n_s)
        {

            mu_bi = mu_bi + (0.05*mu_bi);

            _for(j,0,N_f) f_s1[j].degree = 0;
            _for(j ,0,N_n){_for(i,0,kmax) nodes[j].link[i] = -1;}
            _for(j ,0,N_f){_for(i,0,kmax) f_s1[j].link[i] = -1;}


        }
        if(AVG_DEG_1 > kmean_n_s)
        {
            // mu_bi = mu_bi - (0.0005*mu_bi);
            mu_bi = mu_bi - (0.005*mu_bi);


            _for(j ,0,N_f) f_s1[j].degree = 0;
            _for(j ,0,N_n){_for(i,0,kmax) nodes[j].link[i] = -1;}
            _for(j ,0,N_f){_for(i,0,kmax) f_s1[j].link[i] = -1;}


        }

        it= it +1;
    }

    degree = 0;
    for(i = 0; i < N_f; i++)
    {
        degree = degree + f_s1[i].degree;

    }


    degree = 0;
    for(i = 0; i < N_n; i++)
    {
        degree = degree + nodes[i].degree;

    }

 printf("degree nodes = %d\n",degree);
 printf("links s1 = %f\n",Num_links_s1);

it=0;

     while(1)
    {
       // printf("\n itera # %d", it);
       // printf("\n Mu is : %f", mu_c);



        Num_links_c = 0;

        for(ii = 0; ii < N_n ; ii++)
        {
            degree = 0;
            for(jj = 0 ; jj < N_f_c; jj++)
            {
                double prob2 = Prob_Connection_C(nodes[ii].kappa, f_c[jj].kappa,mu_c);

                double r = rand_easy(0, 1);
                if (r <= prob2)
                {
                    nodes[ii].link[degree+nodes[ii].degree] = jj + N_f;
                    f_c[jj].link[f_c[jj].degree] = ii;
                    degree++;
                    f_c[jj].degree++;

                    Num_links_c = Num_links_c + 1;
                }
            }
            nodes[ii].degreec = degree;
        }
       // printf("links c = %f\n",Num_links_c);


        AVG_DEG_1 = Num_links_c/N_n;
        AVG_DEG_2 = Num_links_c/N_f_c;
      //  printf("\n The av. degree of nodes of Type 1c s1 is: %f \n", AVG_DEG_1);
      //  printf(" The av. degree of nodes of Type 2c (features) s1 is: %f \n", AVG_DEG_2);
      //  printf("~~~~~~~~~~~");

        if( fabs(AVG_DEG_1-kmean_n_c) < 0.05)
        {
                    break;
        }

        if(it>1000){ break;}

        if(AVG_DEG_1 < kmean_n_c)
        {

            mu_c = mu_c + (0.05*mu_c);

            _for(j,0,N_f_c) f_c[j].degree = 0;
            _for(j ,0,N_n){_for(i,nodes[j].degree,kmax) nodes[j].link[i] = -1;}
            _for(j ,0,N_f_c){_for(i,nodes[j].degree,kmax) f_c[j].link[i] = -1;}


        }
        if(AVG_DEG_1 > kmean_n_c)
        {
            // mu_bi = mu_bi - (0.0005*mu_bi);
            mu_c = mu_c - (0.005*mu_c);


            _for(j,0,N_f_c) f_c[j].degree = 0;
            _for(j ,0,N_n){_for(i,nodes[j].degree,kmax) nodes[j].link[i] = -1;}
            _for(j ,0,N_f_c){_for(i,nodes[j].degree,kmax) f_c[j].link[i] = -1;}

        }

        it= it +1;
    }


     _for(j,0,N_n) nodes[j].degree = nodes[j].degree + nodes[j].degreec;


             degree = 0;
            for(i = 0; i < N_n; i++)
            {
                degree = degree + nodes[i].degree;

            }
        printf(" degree = %f\n",(float)degree/(float)N_n);







//the different features are put together in the struct features


    for(i = 0; i < N_f; i++)
    {
        features[i].kappa = f_s1[i].kappa;
        features[i].theta = f_s1[i].theta;
        features[i].degree = f_s1[i].degree;

    }

    for(i = 0; i < N_f_c; i++)
    {
        features[i + N_f].kappa = f_c[i].kappa;
        features[i + N_f].theta = f_c[i].theta;
        features[i + N_f].degree = f_c[i].degree;

    }


    for(i = 0; i < N_f; i++)
    {
        for(j = 0; j < kmax; j++)
        {features[i].link[j] = f_s1[i].link[j];}
    }

    for(i = 0; i < N_f_c; i++)
    {
        for(j = 0; j < kmax; j++)
        {_for(j,0,kmax) features[i + N_f].link[j] = f_c[i].link[j];}

    }
    for(i = 0; i < N_n; i++)
    {
        for(j = 0; j < nodes[i].degree;  j++)
        {   if(nodes[i].link[j] >= 0 && nodes[i].link[j] <= N_n + N_f + N_f_c){
                                    fprintf(fp1,"%d\t%d\n", i+1, nodes[i].link[j]+N_n+1);
                                    }

            }

    }

    for(i = 0; i < N_n ; i++)
    {
       fprintf(fp2,"node\t%d\t%f\t%f\n",i + 1, nodes[i].kappa, nodes[i].theta);
    }

    for(i = 0; i < N_n / 2 ; i++)
    {
       fprintf(fp2,"F_S1\t%d\t%f\t%f\n",i + N_n + 1, features[i].kappa, features[i].theta);
    }

    for(i = N_n / 2; i < N_n + 1; i++)
    {
       fprintf(fp2,"F_CM\t%d\t%f\t%f\n",i + N_n + 1, features[i].kappa, features[i].theta);
    }


    int a = 0;
    for(i = 0; i < N_n; i++)
    {
        for(j = 0; j < nodes[i].degree; j++)
        { links[a].node = i;
          links[a].feat = N_n + j;}

    }



 printf("av %d\n parameter %d\n",av, parameters);


}

 fclose(fp1);
  fclose(fp2);
}

}
double rand_easy(double a, double b) //Generate a uniformly distribute random number in [a, b] using rand() that generates a number between [0, RAND_Max]
{
    double range= b - a;
    return ((((double)rand() / (double)RAND_MAX) * range) + a);
}

double Prob_Connection(double kappa1, double kappa2, double tetha1, double tetha2, double B, double mu, double Rad) //it compute the connection probability of two nodes
{
    double Delta_T= M_PI - fabs(M_PI- fabs(tetha1- tetha2));
    double X= (Rad * Delta_T)/(mu * kappa1 * kappa2);
    double result = 1 / (1 + pow(X, B));
    return(result);
}

double Prob_Connection_C(double kappa1nc, double kappa2c, double mu_C) //it computes the connection probability of two nodes
{
    double XC = mu_C * kappa1nc * kappa2c;
    double resultC = XC / (1 + XC);
    return(resultC);
}



double Random_powerlaw(double gamma, double k0)
{
    double prob= rand_easy(0, 1); //Select a random number between 0 and 1
    double value = k0 / (pow ((1-prob), (1/(gamma-1))));  //CCDF of a power law degree distribution
    return value;
}

void read_file(const char *file_name , int * Nodes , int * Links) {
    FILE *myfile = fopen(file_name, "r");
    int Num_Nodes = 0;
    int Num_links = 0;
    int Node1;
    int Node2;

    while (1) //Find the number of vertices
    {
        fscanf(myfile, "%d %d", &Node1 , &Node2);

        if (feof(myfile))
            break;

        if (Node1>Num_Nodes)
        {
            Num_Nodes= Node1;
        }

        if (Node2>Num_Nodes)
        {
            Num_Nodes= Node2;
        }
        Num_links = Num_links+1;

     }

    * Nodes = Num_Nodes ;
    * Links = Num_links;
    fclose(myfile);
}

double get_fun (double r, double k1, double k2, double k1_min, double k2_min , double gamma1, double gamma2, double nu)
{
    double phi1, phi2, C, res, temp;
    phi1 =  -log (1-(pow((k1_min/k1), (gamma1-1))));
    phi2 =  -log (1- (pow((k2_min/k2), (gamma2-1))));
    C= (pow(phi1,(nu/(1-nu))) * k1_min * pow(k1, gamma1)) / ((k1_min * pow(k1, gamma1)) - (pow(k1_min, gamma1) * k1));
    temp= pow(phi1, (1/(1-nu))) + pow(phi2, (1/(1-nu)));
    res = (exp(-(pow(temp, (1-nu)))) * pow(temp, -nu) * C) - r ;
    return(res);
}

double get_fun_truncated (double r, double k1, double k2, double k1_min, double k2_min , double k1_c, double k2_c, double gamma1, double gamma2, double nu)
{
    double phi1, phi2, C, res, temp;
    phi1 =  -log (1-((pow(k1, (1-gamma1)) - (pow(k1_c, (1-gamma1)))) / ((pow(k1_min, (1-gamma1)) - (pow(k1_c, (1-gamma1)))))));
    phi2 =  -log (1-((pow(k2, (1-gamma2)) - (pow(k2_c, (1-gamma2)))) / ((pow(k2_min, (1-gamma2)) - (pow(k2_c, (1-gamma2)))))));
    C=  pow(phi1,(nu/(1-nu))) * (1 / (1-((pow(k1, (1-gamma1)) - (pow(k1_c, (1-gamma1)))) / ((pow(k1_min, (1-gamma1)) - (pow(k1_c, (1-gamma1))))))));
    temp= pow(phi1, (1/(1-nu))) + pow(phi2, (1/(1-nu)));
    res = (exp(-(pow(temp, (1-nu)))) * pow(temp, -nu) * C) - r ;
    return(res);
}


double bisect (double r, double x1, double x2, double err,double k1, double k1_min, double k2_min , double gamma1, double gamma2, double nu)
{
   double f1, f2, f;
   double x;
   int iter=0;
   while(1)
  {
      f1= get_fun(r, k1, x1, k1_min,k2_min,  gamma1,  gamma2, nu);
      f2= get_fun(r, k1, x2, k1_min,k2_min,  gamma1,  gamma2, nu);
      if(f1 * f2 >0 || isnan(f1) || isnan(f2))
      {
          printf("\n initial guesses are wrong!\n");
          double T=nan("");
          return T;   //it retruns nan
      }

      x= (x1+x2)/2;
      if(((x2-x1)/x) < err)
      {
          f=get_fun( r, k1, x, k1_min,  k2_min ,  gamma1,  gamma2, nu);
          return x;
      }

      f=get_fun( r, k1, x, k1_min,  k2_min ,  gamma1,  gamma2, nu);

      if( (f*f1)>0)
      {
          x1 = x;
          f1=f;
      }
      else
      {
          x2= x;
          f2 = f;
      }
      iter= iter+1;
   }
}

double bisect_truncated (double r, double x1, double x2, double err,double k1, double k1_min, double k2_min, double k1_c, double k2_c, double gamma1, double gamma2, double nu)
{
   double f1, f2, f;
   double x;
   int iter=0;

   f1= get_fun_truncated(r, k1, x1, k1_min, k2_min, k1_c, k2_c,  gamma1,  gamma2, nu);
   f2= get_fun_truncated(r, k1, x2, k1_min,k2_min, k1_c, k2_c, gamma1,  gamma2, nu);
   if(f1 * f2 >0)
   {
       printf("%f \t %f \n", f1, f2);

       while(1)
       {
           //printf("\nchecking....\n");
           x1= rand_easy(k2_min, k2_c);
           x2= rand_easy(k2_min, k2_c);

           f1= get_fun_truncated(r, k1, x1, k1_min, k2_min, k1_c, k2_c,  gamma1,  gamma2, nu);
           f2= get_fun_truncated(r, k1, x2, k1_min,k2_min, k1_c, k2_c, gamma1,  gamma2, nu);

           if(f1*f2<0)
           {
               printf("******Yes!!********\n");
               printf("%f \t %f \n", f1, f2);
               break;
           }
       }
   }
    while(1)
   {
        f1= get_fun_truncated(r, k1, x1, k1_min, k2_min, k1_c, k2_c,  gamma1,  gamma2, nu);
        f2= get_fun_truncated(r, k1, x2, k1_min,k2_min, k1_c, k2_c, gamma1,  gamma2, nu);
         x= (x1+x2)/2;
         if(((x2-x1)/x) < err)
         {
             f=get_fun_truncated( r, k1, x, k1_min,  k2_min , k1_c, k2_c, gamma1,  gamma2, nu);
              return x;
          }

          f=get_fun_truncated( r, k1, x, k1_min,  k2_min, k1_c, k2_c,  gamma1,  gamma2, nu);

          if( (f*f1)>0)
          {
              x1 = x;
              f1=f;
          }
          else
          {
              x2= x;
              f2 = f;
          }
          iter= iter+1;
   }
}



void swap(int* xp, int* yp)
{
    int temp = *xp;
    *xp = *yp;
    *yp = temp;
}
// Function to perform Selection Sort
void selectionSort(int arr[], int n)
{
    int i, j, min_idx;

    // One by one move boundary of unsorted subarray
    for (i = 0; i < n - 1; i++) {

        // Find the minimum element in unsorted array
        min_idx = i;
        for (j = i + 1; j < n; j++)
            if (arr[j] < arr[min_idx])
                min_idx = j;

        // Swap the found minimum element
        // with the first element
        swap(&arr[min_idx], &arr[i]);
    }
}

double Delta_Theta(double Theta1, double Theta2) //Generate a uniformly distribute random number in [a, b] using rand() that generates a number between [0, RAND_Max]
{
    return (M_PI - fabs(M_PI- fabs(Theta1 - Theta2)));
}


void Assign_Labels(double * Theta_s, int NC , int * Label, int NZ, float R, float alpha)
{
    float c_centers[NC];
    double e= 2* M_PI;

    float r;

    alpha= -1 * alpha;

    printf("\n");
    for (int i=0; i< NC ; i++)     //choose centers at random
    {
        c_centers[i]= rand_easy(0, e);
        printf("Center number %d is located at: %f\n", i, c_centers[i]);

    }

    for (int node=0; node<NZ ; node++)    //for every node
    {
        float d[NC];
        float pr[NC];
        float delta_Th;
        float d_sum;
        d_sum=0;
        for(int l=0; l<NC; l++)        //for every center
        {

            delta_Th= Delta_Theta(Theta_s[node],c_centers[l]);     //compute delta tetha between node and center
            d[l]= R * delta_Th;            //compute distance
            d_sum= d_sum + pow(d[l], alpha);

        }

        for(int l=0; l<NC; l++)     //compute the probability of assigning a node to each of classes (pr keeps the comulative probability)
        {
            if (l==0)
                pr[l]= pow(d[l], alpha) / d_sum;
            else
                pr[l]= (pow(d[l], alpha) / d_sum) + pr[l-1];


        }

        r= rand_easy(0, 1);    //choose a random number and assign label regarding this number

        for(int l=0; l<NC; l++)
        {
            if(l==0)
            {
                if (r>=0 && r<pr[l])
                {
                    Label[node]= l;
                    break;
                }
            }
            else
            {
                if (r>= pr[l-1] && r<pr[l])
                {
                    Label[node]=l;
                    break;
                }

            }
        }


    }



}



