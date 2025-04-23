#include <algorithm>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <cmath>
#include <time.h>
#include <random>
#include <chrono>
#include <iostream>
#include <vector>
#include <queue>
#include <utility>
#include <thread>
#include <unistd.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <fcntl.h>
#include <string.h>
#include <stdbool.h>
#include <dirent.h>
#include <signal.h>
#include <errno.h>
#include <wait.h>
#include <sys/wait.h>

using namespace std;
const double pi = 3.14159265358979323846;
const double radicaldindoipi = 2.50662827463100050241;

int nr_apeluri_GA=0;

double Rastrigin(const vector<double> x) {
    double s = 0;
    int length = x.size();
    for (int i = 0; i < length; i++) {
        s += x[i] * x[i] - 10 * cos(2 * pi * x[i]);
    }
    return 10 * length + s;
}

double Schwefel(const vector<double> x) {
    double s = 0;
    int length = x.size();
    for (int i = 0; i < length; i++) {
        s += (x[i] * sin(sqrt(abs(x[i]))));
    }
    return 418.9829 * length - s;
}

double Michalewicz(const vector<double> x) {
    double s = 0;
    double m = 10;
    int length = x.size();
    for (int i = 0; i < length; i++) {
        s += sin(x[i]) * pow(sin((i+1) * x[i] * x[i] / pi), 2 * m);
    }
    return -s;
}

double DeJong(const vector<double> x) {
    int length = x.size();
    double s = 0;
    for (int i = 0; i < length; i++) {
        s += (x[i] * x[i]);
    }
    return s;
}

double Easom(const vector<double> x){
    int length = x.size();
    double aux = 0, p = -1;
    for(int i = 0; i < length; i++){
        p *= cos(x[i]);
        aux += (x[i] - pi) * (x[i] - pi);
    }
    return p * exp(-aux);
}
/////////////////////////////////

/*
//double temperatura, sigma;
double simpli1, simpli2, S;
double Normal_Distribution(const double &x){
    return simpli1 * exp(-x * x / simpli2);
}

double logaritmic(const double &x){
    return S/(1+log2(1+x)) - 100;
}
*/
//(int)(pow(sin(rep + pow(2, rep))
//(int)((15*sin(rep)-3*rep*cos(rep))/pow(2*rep,log2(rep)*0.5-1))

uniform_real_distribution<> unif(0.0, 1.0);
normal_distribution<> norm(0.5, 0.1);

mt19937_64 generator;
void initializeRandom(int rep, mt19937_64& generator) {
    mt19937_64 helper;
    helper.seed(time(NULL) + rep + clock() * 1000 + hash<thread::id>{}(this_thread::get_id()));
    helper.discard(23412 + rep);
    generator.seed(helper());
}

vector<vector<bool>> initializePopulation(mt19937_64& generator,const int popSize, const int& nrbits) {
    
    int random = 0;
    vector<vector<bool>> rez;
    rez.resize(popSize); for (int i = 0; i < popSize; ++i) rez[i].resize(nrbits);

    for(int k = 0; k < popSize; k++){
        for (int i = 0; i < nrbits; i++) {
            if (random <= 3) {
                initializeRandom(i+k, generator);
                random = generator();
            }
            rez[k][i]=(random & 1);
            random = random >> 1;
        }
    }
    return rez;
}

vector<double> Evaluation(double (*FUNCTIE)(vector<double>),
                          const double& a, const double& b, const int& paramBitLength,
                          const vector<vector<bool>>& Population, const int& popSize){

    int length = Population[0].size();
    vector<double> candidatePhenotype, rez;
    //rez.resize(popSize);
    double transform = (b - a) / ((1 << paramBitLength) - 1);
    int x, P, L = length - paramBitLength;
    for(int k = 0; k < popSize; k++){
        for (int i = 0; i <= L; i += paramBitLength) {
            x = 0; P = paramBitLength + i - 1;
            for (int j = i; j < i + paramBitLength; j++) x = x + Population[k][j] * (1 << (P - j));//cout<<x<<" ";
            candidatePhenotype.push_back(a + x * transform);
        }
        rez.push_back(FUNCTIE(candidatePhenotype));
        candidatePhenotype.clear();
    }
    return rez;
}

int randomPosition(const int& length){
    return generator() % length;
}

void Selection(const int& selectionMethod, const int& tournamentSize_OR_nrElites, 
               double (*FUNCTIE)(vector<double>), const int& nrParameters, const double& a, const double& b,
               const int& paramBitLength, vector<vector<bool>>& Population, const int& popSize){
                
    vector<vector<bool>> newPopulation;
    newPopulation.resize(popSize);
    vector<double> fitness;
    double fitness_sum = 0, P_sum = 0;
    double epsilon = 0.0001, R;
    double min = 999999, max = -999999;
    

    //Prima metoda de selectie - Fortune Wheel 1
    if(selectionMethod == 0){

        vector<double> evalPop = Evaluation(FUNCTIE, a, b, paramBitLength, Population, popSize);
    
        for(int i = 0; i < popSize; i++){
            if(evalPop[i] < min) min = evalPop[i];
            if(evalPop[i] > max) max = evalPop[i];
        }

        double M_m = max - min;

        for(int i = 0; i < popSize; i++) fitness.push_back(1 / (evalPop[i] + M_m)), fitness_sum += fitness[i];
        for(int i = 0; i < popSize; i++) fitness[i] = fitness[i] / fitness_sum + P_sum, P_sum = fitness[i];
        for(int i = 0; i < popSize; i++){//fitness devine vectorul de probabilitati cumulate ^^^
            R = unif(generator);
            for(int j = 0; j < popSize; j++)
                if((R <= fitness[j] && j == 0) || (fitness[j-1] < R && R <= fitness[j])){
                    newPopulation[i] = Population[j];
                    break;
                }
        }
    }

    //A doua metoda de selectie - Fortune Wheel 2
    else if(selectionMethod == 1){

        vector<double> evalPop = Evaluation(FUNCTIE, a, b, paramBitLength, Population, popSize);
    
        for(int i = 0; i < popSize; i++){
            if(evalPop[i] < min) min = evalPop[i];
            if(evalPop[i] > max) max = evalPop[i];
        }

        double M_m_e = max - min + epsilon;

        for(int i = 0; i < popSize; i++) fitness.push_back((max - evalPop[i]) / M_m_e), fitness_sum += fitness[i];
        for(int i = 0; i < popSize; i++) fitness[i] = fitness[i] / fitness_sum + P_sum, P_sum = fitness[i];
        for(int i = 0; i < popSize; i++){//fitness devine vectorul de probabilitati cumulate ^^^
            R = unif(generator);
            for(int j = 0; j < popSize; j++)
                if((fitness[j-1] < R && R <= fitness[j]) || (R <= fitness[j] && j == 0)){
                    newPopulation[i] = Population[j];
                    break;
                }
        }
    }

    //A treia metoda de selectie - Tourtament
    else if(selectionMethod == 2){

        vector<double> tournamentResults;
        vector<vector<bool>> tournamentPopulation;
        tournamentPopulation.resize(tournamentSize_OR_nrElites);
        int newPopSize = 0;

        while(newPopSize < popSize) {

            for(int i = 0; i < tournamentSize_OR_nrElites; i++) {
                tournamentPopulation[i] = Population[randomPosition(popSize)];
            }

            tournamentResults = Evaluation(FUNCTIE, a, b, paramBitLength, tournamentPopulation, tournamentSize_OR_nrElites);

            int winnerIndex = 0;
            for(int i = 0; i < tournamentSize_OR_nrElites; i++) {
                if(tournamentResults[i] < tournamentResults[winnerIndex]) {
                    winnerIndex = i;
                }
            }

            newPopulation[newPopSize] = tournamentPopulation[winnerIndex];
            newPopSize++;

            tournamentPopulation.clear();
            tournamentPopulation.resize(tournamentSize_OR_nrElites);
        }
    }

    //A patra metoda de selectie - Elitism
    else{

        vector<double> evalPop = Evaluation(FUNCTIE, a, b, paramBitLength, Population, popSize);

        for(int i = 0; i < tournamentSize_OR_nrElites; i++){

            int minPos = 0;
            for(int j = 1; j < popSize; j++) if(evalPop[j] < evalPop[minPos]) minPos = j;
            
            newPopulation[i] = Population[minPos];
            evalPop[minPos] = 9999;
        }

        //min = 9999;
        for(int i = 0; i < popSize; i++){
            if(evalPop[i] < min) min = evalPop[i];
           //if(evalPop[i] > max) max = evalPop[i];
        }
        //max = 9999;
        double M_m = max - min;

        for(int i = 0; i < popSize; i++) fitness.push_back(1 / (evalPop[i] + M_m)), fitness_sum += fitness[i];
        for(int i = 0; i < popSize; i++) fitness[i] = fitness[i] / fitness_sum + P_sum, P_sum = fitness[i];
        for(int i = tournamentSize_OR_nrElites; i < popSize; i++){//fitness devine vectorul de probabilitati cumulate ^^^
            R = unif(generator);
            for(int j = 0; j < popSize; j++)
                if((R <= fitness[j] && j == 0) || (fitness[j-1] < R && R <= fitness[j])){
                    newPopulation[i] = Population[j];
                    break;
                }
        }
    }
    Population = newPopulation;
}

vector<double> initializeMutationVector(const double& P_mut, const double& sigma, const double& paramBitLength, const int& nrParameters){

    vector<double> mutProbVector;
    int Length = paramBitLength * nrParameters;
    mutProbVector.resize(Length);

    for(int i = 0; i < paramBitLength; i++){
        double x = (double)i / (double)(paramBitLength-1) * sigma;
        mutProbVector[i] = exp(-x * x / (2 * sigma * sigma)) * 1.394 * sin(P_mut * 0.8);
    }

    for(int i = paramBitLength; i < Length; i++){
        mutProbVector[i] = mutProbVector[i - paramBitLength];
    }
    return mutProbVector;
}

void Mutation(const double& P_mut, vector<double>& mutProbVector, const double& sigma, const int& t,
              const int& nrGenerations, const int& candidateBitLength, const int& paramBitLength,
              vector<vector<bool>>& Population, const int& popSize){

    double R;
    for(int i = 0; i < popSize; i++){
        for(int j = 0; j < candidateBitLength; j++){
            R = unif(generator);
            if(R < mutProbVector[j]) Population[i][j] = !Population[i][j];
        }
    }

    //LA FINAL AJUSTARE VECTOR DE PROBABILITATI MUTATIE
    double doisigmapatrat = 2 * sigma * sigma;
    double sin_P_mut_08_1394 =1.394 * sin(P_mut * 0.8);

    for(int i = 0; i < paramBitLength; i++){
        double x_m = (double)i / (double)(paramBitLength-1) - (double)t / (double)nrGenerations; /* t / nrGenerations = mean \in [0,1] */
        mutProbVector[i] = exp(-x_m * x_m / doisigmapatrat) * sin_P_mut_08_1394;
    }

    for(int i = paramBitLength; i < candidateBitLength; i++) mutProbVector[i] = mutProbVector[i - paramBitLength];
}

void Crossover(const double& P_cr, vector<vector<bool>>& Population, const int& popSize){
    
    double R = unif(generator);
    if(R < P_cr){
        int p1 = randomPosition(popSize);
        int p2 = randomPosition(popSize);
        int pos = randomPosition(Population[0].size() - 1);
        for(int i = pos + 1; i < Population[0].size(); i++){
            swap(Population[p1][i], Population[p2][i]);
        }
    }
}

double Print_bestCandidate(double (*FUNCTIE)(vector<double>),
                   const double& a, const double& b, const int& paramBitLength,
                   const vector<vector<bool>>& Population, const int& popSize){

    vector<double> lastGeneration = Evaluation(FUNCTIE, a, b, paramBitLength, Population, popSize);
    double bestResult = 9999;
    int pos = 0;
    for(int  i = 0; i < popSize; i++) if(lastGeneration[i] < bestResult) bestResult = lastGeneration[i], pos = i;

    long long int length = Population[0].size();
    double transform = (b - a) / ((1 << paramBitLength) - 1);
    long long int x, P, L = length - paramBitLength;
    for (long long int i = 0; i <= L; i += paramBitLength) {
        x = 0; P = paramBitLength + i - 1;
        for (long long int j = i; j < i + paramBitLength; j++) x = x + Population[pos][j] * (1 << (P - j));//cout<<x<<" ";
        //cout << (a + x * transform) << ", ";
    }
    //printf(" = %.5f", bestResult);
    return bestResult;
}

double Genetic_Algorithm(double (*FUNCTIE)(vector<double>),
                       const int& nrParameters, const double& a, const double& b,
                       const int& popSize, const int& nrGenerations,
                       int& selectionMethod, const int& tournamentSize_OR_nrElites,
                       const double& P_mut, const double& sigma, //mean = 0 by default at start
                       const double& P_cr){

    initializeRandom(393473, generator);
    int d = 5;// pt ca vrem 5 zecimale dupa virgula
    int paramBitLength = ceil(log2((b - a) * pow(10, d)));
    int Length = nrParameters * paramBitLength;
    vector<vector<bool>> Population = initializePopulation(generator, popSize, Length);
    vector<double> mutProbVector = initializeMutationVector(P_mut, sigma, paramBitLength, nrParameters);
    
    for(int t = 0; t < nrGenerations; t++){
        
        Selection(selectionMethod, tournamentSize_OR_nrElites, FUNCTIE, nrParameters, a, b, paramBitLength, Population, popSize);
        Mutation(P_mut, mutProbVector, sigma, t, nrGenerations, Length, paramBitLength, Population, popSize);
        Crossover(P_cr, Population, popSize);
    }
    return Print_bestCandidate(FUNCTIE, a, b, paramBitLength, Population, popSize);
}


                /*********************** Meta Algoritm Genetic : ***********************/


double transformDouble(const double a, const double b, const vector<vector<bool>>& metaPopulation, const int& metaCandidate_idx, const int& pos, const int& paramBitLength){
    int x = 0, posFinal = pos + paramBitLength - 1;
    for(int i = pos; i <= posFinal; i++) x = x + metaPopulation[metaCandidate_idx][i] * (1 << (posFinal - i));
    return a + x * (b - a) / ((1 << paramBitLength) - 1);
}

vector<double> Meta_Evaluation(double (*FUNCTIE)(vector<double>), const int& nrParameters, const double& a, const double& b,
                               const vector<vector<bool>>& metaPopulation, const int& META_popSize,
                               const int& metaPopSizeBitLength, const int& metaNrGenerationsBitLength, const int& metaSelectionMethodBitLength, 
                               const int& metatournamentSize_OR_nrElites_BitLength, const int& metaP_mutBitLength, const int& metaSigmaBitLength,
                               const int& metaP_crBitLength){

    vector<double> metaEval;
    int popSize, nrGenerations, selectionMethod, tournamentSize_OR_nrElites;
    double P_mut, sigma, P_cr;

    double e1, e2, e3, e4, e5;
    int bestIndex = 0; double bestResult = 999999;
    printf("Generatia %d\n", nr_apeluri_GA); fflush(stdout);
    for(int i = 0; i < META_popSize; i++){

        int pos = 0;
        popSize = floor(transformDouble(10, 200, metaPopulation, i, pos, metaPopSizeBitLength)); pos += metaPopSizeBitLength;
        nrGenerations = floor(transformDouble(1.0, 2000.0, metaPopulation, i, pos, metaNrGenerationsBitLength)); pos += metaNrGenerationsBitLength;
        selectionMethod = metaPopulation[i][pos] * 2 + metaPopulation[i][pos+1] * 1; pos += metaSelectionMethodBitLength;
        tournamentSize_OR_nrElites = ceil(transformDouble(1.0, 0.05 * (double)popSize, metaPopulation, i, pos, metatournamentSize_OR_nrElites_BitLength)); pos += metatournamentSize_OR_nrElites_BitLength;
        P_mut = transformDouble(0.001, 0.4, metaPopulation, i, pos, metaP_mutBitLength); pos += metaP_mutBitLength;
        sigma = transformDouble(0.01, 11.0, metaPopulation, i, pos, metaSigmaBitLength); pos += metaSigmaBitLength;
        P_cr = transformDouble(0.01, 0.9, metaPopulation, i, pos, metaP_crBitLength);

        printf("%d: parametri cu care se apeleaza GA acum sunt: FUNCTIE, %d, %.2f, %.2f, %d, %d, %d, %d, %.5f, %.5f, %.5f;\n", nr_apeluri_GA, nrParameters, a, b, popSize, nrGenerations, selectionMethod, tournamentSize_OR_nrElites, P_mut, sigma, P_cr);fflush(stdout);

        e1 = Genetic_Algorithm(FUNCTIE, nrParameters, a, b, popSize, nrGenerations, selectionMethod, tournamentSize_OR_nrElites, P_mut, sigma, P_cr);
        e2 = Genetic_Algorithm(FUNCTIE, nrParameters, a, b, popSize, nrGenerations, selectionMethod, tournamentSize_OR_nrElites, P_mut, sigma, P_cr);
        e3 = Genetic_Algorithm(FUNCTIE, nrParameters, a, b, popSize, nrGenerations, selectionMethod, tournamentSize_OR_nrElites, P_mut, sigma, P_cr);
        e4 = Genetic_Algorithm(FUNCTIE, nrParameters, a, b, popSize, nrGenerations, selectionMethod, tournamentSize_OR_nrElites, P_mut, sigma, P_cr);
        e5 = Genetic_Algorithm(FUNCTIE, nrParameters, a, b, popSize, nrGenerations, selectionMethod, tournamentSize_OR_nrElites, P_mut, sigma, P_cr);
        metaEval.push_back((e1 + e2 + e3 + e4 + e5) * 0.2);

        if(bestResult > metaEval[i])bestResult = metaEval[i], bestIndex = i;
    }
    nr_apeluri_GA ++;
    printf("%d: Evaluarea populatiei curente: ", nr_apeluri_GA); fflush(stdout);
    for(int i = 0; i < metaEval.size(); i++)printf("%.5f, ", metaEval[i]);
    printf("\n"); fflush(stdout);

    int pos = 0, i = bestIndex;
    popSize = floor(transformDouble(10, 200, metaPopulation, i, pos, metaPopSizeBitLength)); pos += metaPopSizeBitLength;
    nrGenerations = floor(transformDouble(1.0, 2000.0, metaPopulation, i, pos, metaNrGenerationsBitLength)); pos += metaNrGenerationsBitLength;
    selectionMethod = metaPopulation[i][pos] * 2 + metaPopulation[i][pos+1] * 1; pos += metaSelectionMethodBitLength;
    tournamentSize_OR_nrElites = ceil(transformDouble(1.0, 0.05 * (double)popSize, metaPopulation, i, pos, metatournamentSize_OR_nrElites_BitLength)); pos += metatournamentSize_OR_nrElites_BitLength;
    P_mut = transformDouble(0.001, 0.4, metaPopulation, i, pos, metaP_mutBitLength); pos += metaP_mutBitLength;
    sigma = transformDouble(0.01, 11.0, metaPopulation, i, pos, metaSigmaBitLength); pos += metaSigmaBitLength;
    P_cr = transformDouble(0.01, 0.9, metaPopulation, i, pos, metaP_crBitLength);

    printf("Fenotipul celui mai bun candidat al populatiei curente (pS, nrG, sM, t_E, P_m, s, P_c):  %d, %d, %d, %d, %.5f, %.5f, %.5f || GA(candidat): %.5f\n\n", popSize, nrGenerations, selectionMethod, tournamentSize_OR_nrElites, P_mut, sigma, P_cr, bestResult); fflush(stdout);
    return metaEval;
}

void Meta_Selection(double (*FUNCTIE)(vector<double>), const int& nrParameters, const double& a, const double& b,
                    vector<vector<bool>>& metaPopulation, const int& META_popSize,
                    const int& metaPopSizeBitLength, const int& metaNrGenerationsBitLength, const int& metaSelectionMethodBitLength, 
                    const int& metatournamentSize_OR_nrElites_BitLength, const int& metaP_mutBitLength, const int& metaSigmaBitLength,
                    const int& metaP_crBitLength){

    //Selection(selectionMethod, tournamentSize_OR_nrElites, FUNCTIE, nrParameters, a, b, paramBitLength, Population, popSize);

    vector<double> metaEvalPop = Meta_Evaluation(FUNCTIE, nrParameters, a, b,
                                                 metaPopulation, META_popSize,
                                                 metaPopSizeBitLength, metaNrGenerationsBitLength, metaSelectionMethodBitLength,
                                                 metatournamentSize_OR_nrElites_BitLength, metaP_mutBitLength, metaSigmaBitLength,
                                                 metaP_crBitLength);
    vector<vector<bool>> newPopulation;
    newPopulation.resize(META_popSize);
    vector<double> fitness;
    double fitness_sum = 0, P_sum = 0;
    double epsilon = 0.0001, R;
    double min = 9999, max = -9999;

    //Elitism

    int tournamentSize_OR_nrElites = 5;
    for(int i = 0; i < tournamentSize_OR_nrElites; i++){

        int minPos = 0; double minVal = 9999;
        for(int j = 0; j < META_popSize; j++) if(metaEvalPop[j] < minVal) minVal = metaEvalPop[j], minPos = j;
            
        newPopulation[i] = metaPopulation[minPos];
        metaEvalPop[minPos] = 9999;
    }

    /*
    for(int i = 0; i < META_popSize; i++){
        if(metaEvalPop[i] < min) min = metaEvalPop[i];
        if(metaEvalPop[i] > max) max = metaEvalPop[i];
    }

    double M_m_e = max - min + epsilon;

    for(int i = 0; i < META_popSize; i++) fitness.push_back((max - metaEvalPop[i]) / M_m_e), fitness_sum += fitness[i];
    for(int i = 0; i < META_popSize; i++) fitness[i] = fitness[i] / fitness_sum + P_sum, P_sum = fitness[i];
    for(int i = tournamentSize_OR_nrElites; i < META_popSize; i++){//fitness devine vectorul de probabilitati cumulate ^^^
        R = unif(generator);
        for(int j = 0; j < META_popSize; j++)
            if((fitness[j-1] < R && R <= fitness[j]) || (R <= fitness[j] && j == 0)){
                newPopulation[i] = metaPopulation[j];
                break;
            }
    }
    */

        //min = 9999;
    for(int i = 0; i < META_popSize; i++){
        if(metaEvalPop[i] < min) min = metaEvalPop[i];
        //if(evalPop[i] > max) max = evalPop[i];
    }
    max = 9999;
    double M_m = max - min;

    for(int i = 0; i < META_popSize; i++) fitness.push_back(1.0 / (metaEvalPop[i] + M_m)), fitness_sum += fitness[i];
    for(int i = 0; i < META_popSize; i++) fitness[i] = fitness[i] / fitness_sum + P_sum, P_sum = fitness[i];
    for(int i = tournamentSize_OR_nrElites; i < META_popSize; i++){//fitness devine vectorul de probabilitati cumulate ^^^
        R = unif(generator);
        for(int j = 0; j < META_popSize; j++)
            if((R <= fitness[j] && j == 0) || (fitness[j-1] < R && R <= fitness[j])){
                newPopulation[i] = metaPopulation[j];
                break;
            }
    }

    metaPopulation = newPopulation;
}

void Meta_Mutation(const double& META_P_mut, vector<vector<bool>>& metaPopulation, const int& META_popSize, const int META_candidateBitLength){

    double R;
    for(int i = 0; i < META_popSize; i++){
        for(int j = 0; j < META_candidateBitLength; j++){
            R = unif(generator);
            if(R < META_P_mut) metaPopulation[i][j] = !metaPopulation[i][j];
        }
    }
}

void Meta_Crossover(const double& META_P_cr, vector<vector<bool>>& metaPopulation, const int& META_popSize){

    double R = unif(generator);
    if(R < META_P_cr){
        int p1 = randomPosition(META_popSize);
        int p2 = randomPosition(META_popSize);
        int pos = randomPosition(metaPopulation[0].size() - 1);
        for(int i = pos + 1; i < metaPopulation[0].size(); i++){
            swap(metaPopulation[p1][i], metaPopulation[p2][i]);
        }
    }
}

void Print_bestParameters_GA(double (*FUNCTIE)(vector<double>), const int& nrParameters, const double& a, const double& b,
                             const vector<vector<bool>>& metaPopulation, const int& META_popSize,
                             const int& metaPopSizeBitLength, const int& metaNrGenerationsBitLength, const int& metaSelectionMethodBitLength, 
                             const int& metatournamentSize_OR_nrElites_BitLength, const int& metaP_mutBitLength, const int& metaSigmaBitLength,
                             const int& metaP_crBitLength){

    vector<double> metaEvalPop = Meta_Evaluation(FUNCTIE, nrParameters, a, b,
                                                 metaPopulation, META_popSize,
                                                 metaPopSizeBitLength, metaNrGenerationsBitLength, metaSelectionMethodBitLength,
                                                 metatournamentSize_OR_nrElites_BitLength, metaP_mutBitLength, metaSigmaBitLength,
                                                 metaP_crBitLength);

    double bestResult = 999999;
    int bestIndex = 0;
    for(int  i = 0; i < META_popSize; i++) if(metaEvalPop[i] < bestResult) bestResult = metaEvalPop[i], bestIndex = i;

    int popSize, nrGenerations, selectionMethod, tournamentSize_OR_nrElites;
    double P_mut, sigma, P_cr;
    
    int pos = 0, i = bestIndex; 
    popSize = floor(transformDouble(10, 200, metaPopulation, i, pos, metaPopSizeBitLength)); pos += metaPopSizeBitLength;
    nrGenerations = floor(transformDouble(1.0, 2000.0, metaPopulation, i, pos, metaNrGenerationsBitLength)); pos += metaNrGenerationsBitLength;
    selectionMethod = metaPopulation[i][pos] * 2 + metaPopulation[i][pos+1] * 1; pos += metaSelectionMethodBitLength;
    tournamentSize_OR_nrElites = ceil(transformDouble(1.0, 0.05 * (double)popSize, metaPopulation, i, pos, metatournamentSize_OR_nrElites_BitLength)); pos += metatournamentSize_OR_nrElites_BitLength;
    P_mut = transformDouble(0.001, 0.4, metaPopulation, i, pos, metaP_mutBitLength); pos += metaP_mutBitLength;
    sigma = transformDouble(0.01, 11.0, metaPopulation, i, pos, metaSigmaBitLength); pos += metaSigmaBitLength;
    P_cr = transformDouble(0.01, 0.9, metaPopulation, i, pos, metaP_crBitLength);

    printf("Parametri cei mai buni sunt: \npopSize = %d\nnrGenerations = %d\nselectionMethod = %d\ntournamentSize_OR_nrElites = %d\nP_mut = %.5f\nsigma = %.5f\nP_cr = %.5f\nRezultatul cel mai bun este: %.5f\n**************\n",popSize, nrGenerations, selectionMethod, tournamentSize_OR_nrElites, P_mut, sigma, P_cr, bestResult); fflush(stdout);    
}

void Meta_Genetic_Algorithm(double (*FUNCTIE)(vector<double>),
                            const int& nrParameters, const double& a, const double& b,
                            const int& META_popSize, const int& META_nrGenerations,
                            const double& META_P_mut, const double& META_P_cr){

    initializeRandom(12371, generator);
    int metaPopSizeBitLength = ceil(log2(190.0));
    int metaNrGenerationsBitLength = ceil(log2(1999.0));
    int metaSelectionMethodBitLength = 2;
    int metatournamentSize_OR_nrElites_BitLength = ceil(log2(9.0));
    int metaP_mutBitLength = ceil(log2(0.399 * pow(10, 5))); 
    int metaSigmaBitLength = ceil(log2(10.99 * pow(10, 5)));
    int metaP_crBitLength = ceil(log2(0.89 * pow(10, 5)));

    int Length = metaPopSizeBitLength + metaNrGenerationsBitLength +
                 metaSelectionMethodBitLength + metatournamentSize_OR_nrElites_BitLength +
                 metaP_mutBitLength + metaSigmaBitLength + metaP_crBitLength;

/*
    vector<bool> V = {1,0,1,0,0,0,0,1,
1,1,0,1,0,1,1,0,0,0,1,1,1,
0,1,1,0,
0,0,0,0,0,0,0,1,1,0,1,1,0,0,1,1,
0,0,1,1,1,1,0,0,1,1,0,1,0,0,0,1,0,1,1,0,1,
0,1,1,0,1,1,0,0,1,1,1,1,0,1,0,1,1};
vector<vector<bool>> metaPopulation; metaPopulation.resize(META_popSize);
for(int i = 0 ;i < META_popSize; i++)metaPopulation[i]=V;*/
    vector<vector<bool>> metaPopulation = initializePopulation(generator, META_popSize, Length);

    for(int t = 0; t < META_nrGenerations; t++){

        Meta_Selection(FUNCTIE, nrParameters, a, b,
                       metaPopulation, META_popSize,
                       metaPopSizeBitLength, metaNrGenerationsBitLength, metaSelectionMethodBitLength,
                       metatournamentSize_OR_nrElites_BitLength, metaP_mutBitLength, metaSigmaBitLength,
                       metaP_crBitLength);
        Meta_Mutation(META_P_mut, metaPopulation, META_popSize, Length);
        Meta_Crossover(META_P_cr, metaPopulation, META_popSize);
    }
    Print_bestParameters_GA(FUNCTIE, nrParameters, a, b,
                            metaPopulation, META_popSize,
                            metaPopSizeBitLength, metaNrGenerationsBitLength, metaSelectionMethodBitLength,
                            metatournamentSize_OR_nrElites_BitLength, metaP_mutBitLength, metaSigmaBitLength,
                            metaP_crBitLength);
}

int main(){

    //FUNCTIE
    int nrParameters = 30;
    double a, b;
    //a = -100; b = 100; //Easom
    //a = -500, b = 500; //Schwefel
    //a=400;b=422;
    //a = 0, b = pi; //Michalewicz
    a = -5.12, b = 5.12; // DeJong/Rastrigin

    int META_popSize = 15;
    int META_nrGenerations = 100;

    int selectionMethod = 3;
    int tournamentSize_OR_nrElites = 2;   double d = 5;

    double META_P_mut = 3.0/79.0;//0.86728;//1.0/ceil(log2((b - a) * pow(10, d)))/(double)nrParameters;cout<<P_mut<<endl;
    double sigma = 0.18525;

    double META_P_cr = 0.4;//0.028032;

    auto start = std::chrono::high_resolution_clock::now();

    //printf("%.5f",Genetic_Algorithm(Rastrigin, nrParameters, a, b, popSize, nrGenerations, selectionMethod, tournamentSize_OR_nrElites, P_mut, sigma, P_cr));

    Meta_Genetic_Algorithm(Rastrigin, nrParameters, a, b, META_popSize, META_nrGenerations, META_P_mut, META_P_cr);
    
    auto end = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double, std::milli> duration = end - start;
    std::cout << endl << "Execution time: " << duration.count() * 0.001 << " seconds" << std::endl<< "**************"<< endl << endl;

    return 0;
}
