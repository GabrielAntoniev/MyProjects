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
//double temperatura, sigma;
double simpli1, simpli2, S;

///int stop = 0;
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

double Normal_Distribution(const double &x){
    return simpli1 * exp(-x * x / simpli2);
}

double logaritmic(const double &x){
    return S/(1+log2(1+x)) - 100;
}

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
    vector<vector<bool>> rez;//[popSize][nrbits];//(popSize, vector<int>(nrbits));
    rez.resize(popSize); for (int i = 0; i < popSize; ++i) rez[i].resize(nrbits);

    for(int k = 0; k < popSize; k++){
        for (int i = 0; i < nrbits; i++) {
            if (random <= 3) {
                initializeRandom(i+k, generator);
                random = generator();
            }
            rez[k][i]=(random & 1);
            random = random >> 1;
            //rez[k][i]=0;
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
    //max=-999999;
        for(int i = 0; i < popSize; i++){
            if(evalPop[i] < min) min = evalPop[i];
            //if(evalPop[i] > max) max = evalPop[i];
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
            for(int i = 1; i < tournamentSize_OR_nrElites; i++) {
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

        //cout<<"evaluare pop inainte modific. elitism:";for(int i = 0; i < popSize; i++)printf("%.8f, ",evalPop[i]);cout<<endl;

        for(int i = 0; i < tournamentSize_OR_nrElites; i++){

            int minPos = 0, minVal = 999999;
            for(int j = 0; j < popSize; j++) if(evalPop[j] < evalPop[minPos]) minPos = j;//, minVal = evalPop[j];
            
            newPopulation[i] = Population[minPos];
            evalPop[minPos] = 999999;
        }

    /*
        for(int i = 0; i < popSize; i++){
            if(evalPop[i] < min) min = evalPop[i];
            if(evalPop[i] > max) max = evalPop[i];
        }

        double M_m_e = max - min + epsilon;

        for(int i = 0; i < popSize; i++) fitness.push_back((max - evalPop[i]) / M_m_e), fitness_sum += fitness[i];
        for(int i = 0; i < popSize; i++) fitness[i] = fitness[i] / fitness_sum + P_sum, P_sum = fitness[i];
        for(int i = 0; i < tournamentSize_OR_nrElites; i++){//fitness devine vectorul de probabilitati cumulate ^^^
            R = unif(generator);
            for(int j = 0; j < popSize; j++)
                if((fitness[j-1] < R && R <= fitness[j]) || (R <= fitness[j] && j == 0)){
                    newPopulation[i] = Population[j];
                    break;
                }
        }*/

        //cout<<"evaluare pop dupa 'selectare' elite:";for(int i = 0; i < popSize; i++)printf("%.8f, ",evalPop[i]);cout<<endl;

        //min = 999999;
        for(int i = 0; i < popSize; i++){
            if(evalPop[i] < min) min = evalPop[i];
            //if(evalPop[i] > max) max = evalPop[i];
        }
        

        double M_m = max - min;
        
        //printf("maxim = %.8f, minim = %.8f, M_m = %.8f\n",max,min,M_m);

	    max = 999999;
        for(int i = 0; i < popSize; i++) fitness.push_back(1 / (evalPop[i] + M_m)), fitness_sum += fitness[i];

        //cout<<"fitness populatie: ";for(int i = 0; i < popSize; i++)printf("%.8f, ",fitness[i]);cout<<endl;

        for(int i = 0; i < popSize; i++) fitness[i] = fitness[i] / fitness_sum + P_sum, P_sum = fitness[i];

        //cout<<"P_cumul populatie : ";for(int i = 0; i < popSize; i++)printf("%.8f, ",fitness[i]);cout<<endl<<"**************************"<<endl;

        for(int i = tournamentSize_OR_nrElites; i < popSize; i++){//fitness devine vectorul de probabilitati cumulate ^^^
            R = unif(generator);
            for(int j = 0; j < popSize; j++)
                if((R <= fitness[j] && j == 0) || (fitness[j-1] < R && R <= fitness[j])){
                    newPopulation[i] = Population[j];
                    break;
                }
        }

        //stop ++;
        //if(stop == 1)exit(0);
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

void Mutation(const double& P_mut, vector<double>& mutProbVector, const double& sigma, const int& t, const int& nrGenerations, const int& Length, const int& paramBitLength, vector<vector<bool>>& Population, const int& popSize){

    double R;
    for(int i = 0; i < popSize; i++){
        for(int j = 0; j < Length; j++){
            R = unif(generator);
            if(R < mutProbVector[j]) Population[i][j] = !Population[i][j];
        }
    }

    //LA FINAL AJUSTARE VECTOR DE PROBABILITATI MUTATIE
    double doisigmapatrat = 2 * sigma * sigma;
    double sin_P_mut_08_1394 =1.394 * sin(P_mut * 0.8);

    for(int i = 0; i < paramBitLength; i++){
        double x_m = (double)i / (double)(paramBitLength-1) * sigma - (double)t / (double)nrGenerations; /* t / nrGenerations = mean \in [0,1] */
        mutProbVector[i] = exp(-x_m * x_m / doisigmapatrat) * sin_P_mut_08_1394;
    }

    for(int i = paramBitLength; i < Length; i++){
        mutProbVector[i] = mutProbVector[i - paramBitLength];
    }
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
    double bestResult = 999999;
    int pos = 0;
    for(int  i = 0; i < popSize; i++) if(lastGeneration[i] < bestResult) bestResult = lastGeneration[i], pos = i;

    //cout << "Minimul este in: ";

//    if(FUNCTIE == Rastrigin) cout<<"Rastrigin("; else
//    if(FUNCTIE == Schwefel) cout<<"Schwefel("; else
//    if(FUNCTIE == Michalewicz) cout<<"Michalewicz("; else
//    if(FUNCTIE == DeJong) cout<<"DeJong(";

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
    
    for(int t = 0; t <= nrGenerations; t++){
        
        Selection(selectionMethod, tournamentSize_OR_nrElites, FUNCTIE, nrParameters, a, b, paramBitLength, Population, popSize);
        Mutation(P_mut, mutProbVector, sigma, t, nrGenerations, Length, paramBitLength, Population, popSize);
        Crossover(P_cr, Population, popSize);
    }
    return Print_bestCandidate(FUNCTIE, a, b, paramBitLength, Population, popSize);
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

    int popSize = 200;//20;
    int nrGenerations = 1121;//1405;//2000;

    //int selectionMethod = 3;
    int tournamentSize_OR_nrElites = 2;/*8;*/   double d = 5;

    double P_mut = 0.02720;//0.05263;//0.01704;//1.0/ceil(log2((b - a) * pow(10, d))) ;// ((double)(nrParameters+25)); //0.03704
    double sigma = 0.31203;//0.23;

    double P_cr = 0.53137;///0.02;//0.01;

//Genetic_Algorithm(Rastrigin, nrParameters, a, b, popSize, nrGenerations, CURRENT_SELECTION_METHOD, tournamentSize_OR_nrElites, P_mut, sigma, P_cr);
    //cout << "Rastrigin" << endl;
    //cout << "DOMENIUL DE DEFINITIE: [" << a << ", " << b << ']' << endl;
    //cout << "NR. DE PARAMETRI FUNCTIE: " << nrParameters << endl;
    printf("Rastrigin - Genetic Algorithm\nNR. DE PARAMETRI FUNCTIE: %d\nDOMENIU DE DEFINITIE: [%.2f, %.2f]\nPOP. SIZE: %d\nNR. GENERATIONS: %d\nSELECTION METHOD: see below, jus wait :))\nTOURNAMENT SIZE OR NR. ELITES (ONLY FOR SELECTION METHOD = 2 OR 3): %d\nP_MUTATION :%.5f\nSIGMA: %.5f\nP_CROSSOVER: %.5f\n", nrParameters, a, b, popSize, nrGenerations, tournamentSize_OR_nrElites, P_mut, sigma, P_cr);
    fflush(stdout);


///printf("%.5f",Genetic_Algorithm(Rastrigin, nrParameters, a, b, popSize, nrGenerations, selectionMethod, tournamentSize_OR_nrElites, P_mut, sigma, P_cr));
//exit(0);


    pid_t pid; int k;

    for(k = 1; k <= 4; k++){
        pid=fork();
        if(pid<0){ perror("fork failed");exit(1);}

        if(pid == 0){

            if(k == 1){exit(0);
                

                double statistica_rez[40], statistica_timp[40];
                double st_dev = 0, maxim = -999999, minim = 999999;
                double medie_rez = 0, medie_timp = 0, total_timp = 0;


                ////
                int CURRENT_SELECTION_METHOD = 0;
                ////


                for(long long int i = 0; i < 40; i++){

                    //st_dev = 0, minim = -999999, maxim = 999999;
                    //medie_rez = 0, medie_timp = 0;

                    initializeRandom(i, generator);
                    auto start = chrono::high_resolution_clock::now();

                    statistica_rez[i] = Genetic_Algorithm(Rastrigin, nrParameters, a, b, popSize, nrGenerations, CURRENT_SELECTION_METHOD, tournamentSize_OR_nrElites, P_mut, sigma, P_cr);

                    auto end = std::chrono::high_resolution_clock::now();
                    chrono::duration<double, std::milli> duration = end - start;

                    statistica_timp[i] = (double)(duration.count()) * 0.001;
                    total_timp += statistica_timp[i];
                }

                
                for(long long int i = 0; i < 40; i++){
                    medie_rez += statistica_rez[i]; if(statistica_rez[i] < minim)minim = statistica_rez[i]; if(statistica_rez[i] > maxim)maxim = statistica_rez[i];
                    medie_timp += statistica_timp[i];
                }
                medie_rez /= 40; medie_timp /= 40;

                for(long long int i = 0; i < 40; i++){
                    st_dev += (statistica_rez[i]-medie_rez)*(statistica_rez[i]-medie_rez);
                }
                st_dev /= 39; st_dev = sqrt(st_dev);

             
                
                printf("\nREZULTATE GA: selectionMethod = %d:\n\nVALOARE MINIMA: %.5f\nVALOARE MAXIMA: %.5f\nMEDIE VALORI: %.5f\nDEVIATIE STANDARD: %.5f\nTIMP MEDIU: %.2f SECUNDE\nTOTAL TIMP EXECUTIE PROCES: %.2f\n*****************\n",CURRENT_SELECTION_METHOD ,minim, maxim, medie_rez, st_dev, medie_timp, total_timp);


                
            }


            else 
            if(k == 2){
                

                double statistica_rez[40], statistica_timp[40];
                double st_dev = 0, maxim = -999999, minim = 999999;
                double medie_rez = 0, medie_timp = 0, total_timp = 0;


                ////
                int CURRENT_SELECTION_METHOD = 1;
                ////


                for(long long int i = 0; i < 40; i++){

                    //st_dev = 0, minim = -999999, maxim = 999999;
                    //medie_rez = 0, medie_timp = 0;

                    initializeRandom(i, generator);
                    auto start = chrono::high_resolution_clock::now();

                    statistica_rez[i] = Genetic_Algorithm(Rastrigin, nrParameters, a, b, popSize, nrGenerations, CURRENT_SELECTION_METHOD, tournamentSize_OR_nrElites, P_mut, sigma, P_cr);

                    auto end = std::chrono::high_resolution_clock::now();
                    chrono::duration<double, std::milli> duration = end - start;

                    statistica_timp[i] = (double)(duration.count()) * 0.001;
                    total_timp += statistica_timp[i];
                }

                
                for(long long int i = 0; i < 40; i++){
                    medie_rez += statistica_rez[i]; if(statistica_rez[i] < minim)minim = statistica_rez[i]; if(statistica_rez[i] > maxim)maxim = statistica_rez[i];
                    medie_timp += statistica_timp[i];
                }
                medie_rez /= 40; medie_timp /= 40;

                for(long long int i = 0; i < 40; i++){
                    st_dev += (statistica_rez[i]-medie_rez)*(statistica_rez[i]-medie_rez);
                }
                st_dev /= 39; st_dev = sqrt(st_dev);

                printf("\nREZULTATE GA: selectionMethod = %d:\n\nVALOARE MINIMA: %.5f\nVALOARE MAXIMA: %.5f\nMEDIE VALORI: %.5f\nDEVIATIE STANDARD: %.5f\nTIMP MEDIU: %.2f SECUNDE\nTOTAL TIMP EXECUTIE PROCES: %.2f\n*****************\n",CURRENT_SELECTION_METHOD ,minim, maxim, medie_rez, st_dev, medie_timp, total_timp);
                

            }
            else
            if(k == 3){
                

                double statistica_rez[40], statistica_timp[40];
                double st_dev = 0, maxim = -999999, minim = 999999;
                double medie_rez = 0, medie_timp = 0, total_timp = 0;


                ////
                int CURRENT_SELECTION_METHOD = 2;
                ////


                for(long long int i = 0; i < 40; i++){

                    //st_dev = 0, minim = -999999, maxim = 999999;
                    //medie_rez = 0, medie_timp = 0;

                    initializeRandom(i, generator);
                    auto start = chrono::high_resolution_clock::now();

                    statistica_rez[i] = Genetic_Algorithm(Rastrigin, nrParameters, a, b, popSize, nrGenerations, CURRENT_SELECTION_METHOD, tournamentSize_OR_nrElites, P_mut, sigma, P_cr);

                    auto end = std::chrono::high_resolution_clock::now();
                    chrono::duration<double, std::milli> duration = end - start;

                    statistica_timp[i] = (double)(duration.count()) * 0.001;
                    total_timp += statistica_timp[i];
                }

                
                for(long long int i = 0; i < 40; i++){
                    medie_rez += statistica_rez[i]; if(statistica_rez[i] < minim)minim = statistica_rez[i]; if(statistica_rez[i] > maxim)maxim = statistica_rez[i];
                    medie_timp += statistica_timp[i];
                }
                medie_rez /= 40; medie_timp /= 40;

                for(long long int i = 0; i < 40; i++){
                    st_dev += (statistica_rez[i]-medie_rez)*(statistica_rez[i]-medie_rez);
                }
                st_dev /= 39; st_dev = sqrt(st_dev);

                printf("\nREZULTATE GA: selectionMethod = %d:\n\nVALOARE MINIMA: %.5f\nVALOARE MAXIMA: %.5f\nMEDIE VALORI: %.5f\nDEVIATIE STANDARD: %.5f\nTIMP MEDIU: %.2f SECUNDE\nTOTAL TIMP EXECUTIE PROCES: %.2f\n*****************\n",CURRENT_SELECTION_METHOD ,minim, maxim, medie_rez, st_dev, medie_timp, total_timp);

                
            }
            else
            {
                

                double statistica_rez[40], statistica_timp[40];
                double st_dev = 0, maxim = -999999, minim = 999999;
                double medie_rez = 0, medie_timp = 0, total_timp = 0;


                ////
                int CURRENT_SELECTION_METHOD = 3;
                ////


                for(long long int i = 0; i < 40; i++){

                    //st_dev = 0, minim = -999999, maxim = 999999;
                    //medie_rez = 0, medie_timp = 0;

                    initializeRandom(i, generator);
                    auto start = chrono::high_resolution_clock::now();

                    statistica_rez[i] = Genetic_Algorithm(Rastrigin, nrParameters, a, b, popSize, nrGenerations, CURRENT_SELECTION_METHOD, tournamentSize_OR_nrElites, P_mut, sigma, P_cr);

                    auto end = std::chrono::high_resolution_clock::now();
                    chrono::duration<double, std::milli> duration = end - start;

                    statistica_timp[i] = (double)(duration.count()) * 0.001;
                    total_timp += statistica_timp[i];
                }

                
                for(long long int i = 0; i < 40; i++){
                    medie_rez += statistica_rez[i]; if(statistica_rez[i] < minim)minim = statistica_rez[i]; if(statistica_rez[i] > maxim)maxim = statistica_rez[i];
                    medie_timp += statistica_timp[i];
                }
                medie_rez /= 40; medie_timp /= 40;

                for(long long int i = 0; i < 40; i++){
                    st_dev += (statistica_rez[i]-medie_rez)*(statistica_rez[i]-medie_rez);
                }
                st_dev /= 39; st_dev = sqrt(st_dev);

                printf("\nREZULTATE GA: selectionMethod = %d:\n\nVALOARE MINIMA: %.5f\nVALOARE MAXIMA: %.5f\nMEDIE VALORI: %.5f\nDEVIATIE STANDARD: %.5f\nTIMP MEDIU: %.2f SECUNDE\nTOTAL TIMP EXECUTIE PROCES: %.2f\n*****************\n",CURRENT_SELECTION_METHOD ,minim, maxim, medie_rez, st_dev, medie_timp, total_timp);
            }

            exit(0);
        }
    }
    for (long long int p = 1; p <= 4; p++) wait(NULL);
    return 0;
}

/*
int main(){



    //FUNCTIE
    int nrParameters = 30;
    double a, b;
    //a = -100; b = 100; //Easom
    //a = -500, b = 500; //Schwefel
    //a=400;b=422;
    //a = 0, b = pi; //Michalewicz
    a = -5.12, b = 5.12; // DeJong/Rastrigin

    int popSize = 30;
    int nrGenerations = 20000;

    int selectionMethod = 3;
    int tournamentSize_OR_nrElites = 2;   double d = 5;

    double P_mut = 1.0/ceil(log2((b - a) * pow(10, d)));///(double)nrParameters;
    double sigma = 0.3;

    double P_cr = 0.01;

    auto start = std::chrono::high_resolution_clock::now();

    double sum = 0;
    for(int i = 0; i < 1; i++)sum+=Genetic_Algorithm(Rastrigin, nrParameters, a, b, popSize, nrGenerations, selectionMethod, tournamentSize_OR_nrElites, P_mut, sigma, P_cr);
    printf("%.5f",sum);


    auto end = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double, std::milli> duration = end - start;
    std::cout << endl << "Execution time: " << duration.count() * 0.001 << " seconds" << std::endl<< "**************"<< endl << endl;
    return 0;
}

*/
