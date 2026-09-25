#include <time.h>
#include <assert.h>
#include <stdio.h>
#include <stdlib.h>
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
#include <ctime>
#include <errno.h>
#include <unordered_set>
#include <cstdint>
#include <cstdlib>


double rand_easy(double a, double b);

/// Defines
#define sqr(x)      ((x)*(x))
#define _for(a,b,c) for(a=(b); a<(c); ++a)
#define ULONGMAXX   4294967296.0

long long factorial(long long x) { if(x==0LL) return 1LL; else return x*factorial(x-1LL); }

#define DATASET_NAME "edge_list_pp"

#define realisations    1000

using namespace std;

/// Types
struct tipo_nodes
{
    float clustering;
    int degree;
    vector<int> links;
};

struct tipo_features
{
    float clustering_rewiring;
    float clustering;
    int degree;
    vector<int> links;
};

struct tipo_links
{
    int node;
    int feat;
};

typedef double flt;

/// Global variables

tipo_nodes *nodes;
tipo_features *features;
tipo_links *links;


unsigned int seed[256];
unsigned int r;
unsigned char irr;

#define n_tiempos   1000

void inicializa_randomm(void)
{
    int i;
    srand((unsigned)time(0));
    irr=1;
    _for(i,0,256) seed[i]=rand();
    r=seed[0];
    _for(i,0,70000) r=seed[irr++]+=seed[r>>24];
}

flt randomm(void)
{
    r=seed[irr++]+=seed[r>>24];
    return ( (flt)r/ULONGMAXX );
}

int main(int n_args,char *args[])
{
    std::srand(static_cast<unsigned int>(std::time(nullptr))); // Seed RNG with current time
    inicializa_randomm();

    // input --------------------

    const char* extensions[] = { ".dat", ".txt", ".edge" };

    FILE* file = NULL;
    char file_nodes_buf[256];

    for (int e = 0; e < 3; e++) {
        sprintf(file_nodes_buf, "%s%s", DATASET_NAME, extensions[e]);
        file = fopen(file_nodes_buf, "r");

        if (file != NULL) {
            printf("Loaded input file: %s\n", file_nodes_buf);
            break;
        }
    }

    if (!file) {
        printf("Error: could not open input file with extensions .dat, .txt or .edge\n");
        printf("Tried base name: %s\n", DATASET_NAME);
        return 1;
    }


    // -------------------- Output file 1 --------------------
    FILE *fp1;
    char file_fp1[256];
    sprintf(file_fp1, "%s_clustering.dat", DATASET_NAME);

    fp1 = fopen(file_fp1, "w");
    if (fp1 == NULL) {
        perror("Failed: ");
        return 1;
    }

    // -------------------- Output file 2 --------------------
    FILE *fp2;
    char file_fp2[256];
    sprintf(file_fp2, "%s_rewired_clustering.dat", DATASET_NAME);

    fp2 = fopen(file_fp2, "w");
    if (fp2 == NULL) {
        perror("Failed: ");
        return 1;
    }


    int N_n = 0;
    int N_features = 0;
    int N_l = 0;

    int node_id, feature_id;
    int a;

    // Pass 1: count edges (a) and detect max IDs
    int max_node = -1;
    int max_feat = -1;
    a = 0;

    while (fscanf(file, "%d %d", &node_id, &feature_id) == 2) {
        if (node_id > max_node) max_node = node_id;
        if (feature_id > max_feat) max_feat = feature_id;
        a++;
    }

    if (a == 0) {
        printf("Error: empty file or wrong format (expected 2 columns)\n");
        fclose(file);
        return 1;
    }

    N_n = max_node + 1;
    N_features = max_feat + 1;
    N_l = a;

    printf("Detected from file: N_n=%d, N_features=%d, N_l=%d\n", N_n, N_features, N_l);

    // Use new[] (NOT malloc) because structs contain std::vector
    nodes    = new tipo_nodes[N_n];
    features = new tipo_features[N_features];
    links    = new tipo_links[N_l];

    void check_multiedges_after_loading(int N_l, tipo_links* links);

    // Initialize degrees
    for (int i = 0; i < N_features; i++) { features[i].degree = 0; }
    for (int i = 0; i < N_n; i++) { nodes[i].degree = 0; }

    // Pass 2: count degrees for reserve()
    rewind(file);
    while (fscanf(file, "%d %d", &node_id, &feature_id) == 2) {
        if (node_id < 0 || node_id >= N_n) continue;
        if (feature_id < 0 || feature_id >= N_features) continue;

        nodes[node_id].degree++;
        features[feature_id].degree++;
    }

    // Reserve exact capacity to avoid reallocations
    for (int i = 0; i < N_n; i++) {
        nodes[i].links.reserve(nodes[i].degree);
    }
    for (int i = 0; i < N_features; i++) {
        features[i].links.reserve(features[i].degree);
    }

    // Pass 3: fill adjacency and links[]
    rewind(file);
    a = 0;
    while (fscanf(file, "%d %d", &node_id, &feature_id) == 2) {
        if (node_id < 0 || node_id >= N_n) continue;
        if (feature_id < 0 || feature_id >= N_features) continue;

        nodes[node_id].links.push_back(feature_id);
        features[feature_id].links.push_back(node_id);

        // Keep the same 'links' array as before
        links[a].node = node_id;
        links[a].feat = feature_id;
        a++;
    }

    fclose(file);


    // Keep 'a' as the actual number of loaded links, because later code uses it.
    N_l = a;

    check_multiedges_after_loading(N_l, links);

    printf("Clustering Coefficient original network\n");

    /***CLUSTERING coefficient***/
    int i, ii, iii, j, jj, triangles;
    float T;


   _for(i,N_n,N_features)  features[i].clustering  = 0.0;

    for(i = N_n; i < N_features; i++){
           // printf("Nodo %d: degree = %d \n ", i,  features[i].degree);
                T = 0.0;
                if(features[i].degree > 1){

                    for(ii = 0; ii < features[i].degree; ii++){

                    for(iii = ii + 1; iii < features[i].degree; iii++){

                    triangles = 0;

                    if(nodes[features[i].links[ii]].degree > 1 && nodes[features[i].links[iii]].degree > 1){
                        for(j = 0; j < nodes[features[i].links[ii]].degree; j++){

                            for(jj = 0; jj < nodes[ features[i].links[iii] ].degree; jj++){
                                if (nodes[features[i].links[ii]].links[j] == nodes[features[i].links[iii]].links[jj] &&  nodes[features[i].links[ii]].links[j]!= i ){
                                    triangles++;}

                                }    }     }


                    if(triangles >0) {
                        if( nodes[features[i].links[ii]].degree >= nodes[features[i].links[iii]].degree){
                                T +=  (float)triangles / (float)(nodes[features[i].links[iii]].degree - 1);
                                }   else {T +=  (float)triangles / (float)(nodes[features[i].links[ii]].degree - 1);}
                                    }

                                } } }

                features[i].clustering = 2.0 * ( T ) / ( (float)features[i].degree * ((float)features[i].degree - 1.0) );


                        }



        for(i=N_n;i<N_features;i++) {fprintf(fp1,"%d\t%d\t%f\n",i, features[i].degree, features[i].clustering);}

fclose(fp1);


printf("Rewiring\n");

int link1, link2, flag, changes = 0;
int node1, node2, feature1, feature2;

for(int re = 0; re < realisations; re++){

    int changes = 0;   // RESET for every random realization

    do{

        do{
            flag = 0;
            link1 = std::rand() % (a);
            link2 = std::rand() % (a);

            for(ii = 0; ii < nodes[links[link1].node].degree ; ii++){
                   if(nodes[links[link1].node].links[ii] ==  links[link2].feat){flag = 1; break;} }

            if(flag == 0){
                        for(ii = 0; ii < nodes[links[link2].node].degree ; ii++){
                           if(nodes[links[link2].node].links[ii] ==  links[link1].feat){flag = 1; break;}  }  }


           } while( link1 == link2 || flag == 1);

       // printf("L1 %d\t  L2 %d\t node1 %d\t  feat1 %d\t node2 %d\t feat2 %d\n", link1, link2, links[link1].node, links[link1].feat, links[link2].node, links[link2].feat   );

        node1 = links[link1].node;
        feature1 = links[link1].feat;
        node2 = links[link2].node;
        feature2 = links[link2].feat;

        links[link1].node = node2;
        links[link2].node = node1;


        auto newEnd = std::remove(nodes[node1].links.begin(), nodes[node1].links.end(), feature1);
        nodes[node1].links.erase(newEnd, nodes[node1].links.end());
        nodes[node1].links.push_back(feature2);

        auto newEnd2 = std::remove(nodes[node2].links.begin(), nodes[node2].links.end(), feature2);
        nodes[node2].links.erase(newEnd2, nodes[node2].links.end());
        nodes[node2].links.push_back(feature1);

        auto newEnd3 = std::remove(features[feature1].links.begin(), features[feature1].links.end(), node1);
        features[feature1].links.erase(newEnd3,features[feature1].links.end());
        features[feature1].links.push_back(node2);

        auto newEnd4 = std::remove(features[feature2].links.begin(), features[feature2].links.end(), node2);
        features[feature2].links.erase(newEnd4,features[feature2].links.end());
        features[feature2].links.push_back(node1);


        changes++;


    }while(changes < a*a);


    _for(i,N_n ,N_features)  features[i].clustering_rewiring  = 0.0;


     _for(i,N_n,N_features) std::sort(features[i].links.begin(), features[i].links.end());


    for(i = N_n; i < N_features; i++){
           //printf("Nodo %d: degree = %d \n ", i,  features[i].degree);
                T = 0.0;
                if(features[i].degree > 1){

                    for(ii = 0; ii < features[i].degree; ii++){

                    for(iii = ii + 1; iii < features[i].degree; iii++){

                    triangles = 0;

                    if(nodes[features[i].links[ii]].degree > 1 && nodes[features[i].links[iii]].degree > 1){
                        for(j = 0; j < nodes[features[i].links[ii]].degree; j++){

                            for(jj = 0; jj < nodes[ features[i].links[iii] ].degree; jj++){
                                if (nodes[features[i].links[ii]].links[j] == nodes[features[i].links[iii]].links[jj] &&  nodes[features[i].links[ii]].links[j]!= i ){
                                    triangles++;}

                                }    }     }


                    if(triangles >0) {
                        if( nodes[features[i].links[ii]].degree >= nodes[features[i].links[iii]].degree){
                                T +=  (float)triangles / (float)(nodes[features[i].links[iii]].degree - 1);
                                }   else {T +=  (float)triangles / (float)(nodes[features[i].links[ii]].degree - 1);}
                                    }

                                } } }

                features[i].clustering_rewiring = 2.0 * ( T ) / ( (float)features[i].degree * ((float)features[i].degree - 1.0) );


                        }


        for(i=N_n;i<N_features;i++) {fprintf(fp2,"%d\t%d\t%f\t%f\t%f\n",i, features[i].degree, features[i].clustering, features[i].clustering_rewiring, features[i].clustering /  features[i].clustering_rewiring);}
        printf("realisation = %d\n",re);

}
fclose(fp2);


}


double rand_easy(double a, double b) //Generate a uniformly distribute random number in [a, b] using rand() that generates a number between [0, RAND_Max]
{
    double range= b - a;
    return ((((double)rand() / (double)RAND_MAX) * range) + a);
}

static inline uint64_t pack_edge_32_32(int u, int v) {
    // Packs two non-negative 32-bit ints into one 64-bit key
    return (uint64_t)(uint32_t)u << 32 | (uint32_t)v;
}

void check_multiedges_after_loading(int N_l, tipo_links* links) {
    std::unordered_set<uint64_t> seen;
    seen.reserve((size_t)N_l * 2);

    bool has_multi = false;

    for (int e = 0; e < N_l; ++e) {
        const int u = links[e].node;
        const int v = links[e].feat;

        // If your IDs could be negative, skip packing (or handle differently)
        if (u < 0 || v < 0) continue;

        const uint64_t key = pack_edge_32_32(u, v);
        if (!seen.insert(key).second) {
            has_multi = true;
            break;
        }
    }

    if (has_multi) {
        printf("warning your edgelist has multiedges, remove them and start again\n");
         exit(EXIT_FAILURE); // stop program
    }else {printf("no multiedges\n");}
}
