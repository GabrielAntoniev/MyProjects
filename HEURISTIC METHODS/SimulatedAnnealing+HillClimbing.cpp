#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <cmath>
#include <time.h>
#include <cmath>
#include <random>
#include <chrono>
#include <iostream>
#include <vector>
#include <utility>
#include <thread>
#include <stdio.h>
#include <unistd.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <fcntl.h>
#include <string.h>
#include <stdbool.h>
#include <dirent.h>
#include <stdlib.h>
#include <signal.h>
#include <sys/types.h>
#include <errno.h>
#include <wait.h>
#include <sys/wait.h>

using namespace std;
const double pi = 3.14159265358979323846;
const double radicaldindoipi = 2.50662827463100050241;
//double temperatura, sigma;
double simpli1, simpli2, S;

/*
struct results{
    double minutes, output;
};
*/

double Rastrigin(const vector<double> x) {
    double s = 0;
    long long int length = x.size();
    for (long long int i = 0; i < length; i++) {
        s += x[i] * x[i] - 10 * cos(2 * pi * x[i]);
    }
    return 10 * length + s;
}

double Schwefel(const vector<double> x) {
    double s = 0;
    long long int length = x.size();
    for (long long int i = 0; i < length; i++) {
        s += (x[i] * sin(sqrt(abs(x[i]))));
    }
    return 418.9829 * length - s;
}

double Michalewicz(const vector<double> x) {
    double s = 0;
    double m = 10;
    long long int length = x.size();
    for (long long int i = 0; i < length; i++) {
        s += sin(x[i]) * pow(sin((i+1) * x[i] * x[i] / pi), 2 * m);
    }
    return -s;
}

double DeJong(const vector<double> x) {
    long long int length = x.size();
    double s = 0;
    for (long long int i = 0; i < length; i++) {
        s += (x[i] * x[i]);
    }
    return s;
}

double Easom(const vector<double> x){
    long long int length = x.size();
    double aux = 0, p = -1;
    for(long long int i = 0; i < length; i++){
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

uniform_real_distribution<> dis(0.0, 1.0);

mt19937_64 generator;
void initializeRandom(long long int rep, mt19937_64& generator) {
    mt19937_64 helper;
    helper.seed(time(NULL) + rep + clock() * 1000 + hash<thread::id>{}(this_thread::get_id()));
    helper.discard(23412 + rep);
    generator.seed(helper());
}

vector<bool> generateRandomBitString(mt19937_64& generator, const long long int& nrbits_adk_precizia) {
    long long int random = 0;
    vector<bool> bitString;
    for (long long int i = 0; i < nrbits_adk_precizia; i++) {
        if (random <= 3) {
            initializeRandom(i, generator);
            random = generator();
        }
        bitString.push_back(random & 1);
        random = random >> 1;
    }
    return bitString;
}

vector<double> Transform_into_Mathematical_Object(const vector<bool>& Candidate, const double& a, const double& b, const long long int& ParameterbitStringLength) {

    //se aplica pt fiecare parametru al Candidatului
    long long int length = Candidate.size();
    vector<double> rez;
    double b_minus_a_totul_supra_maximul_reprezentabil = (b - a) / ((1 << ParameterbitStringLength) - 1);
    long long int x, P, L = length - ParameterbitStringLength;
    for (long long int i = 0; i <= L; i += ParameterbitStringLength) {
        x = 0; P = ParameterbitStringLength + i - 1;
        for (long long int j = i; j < i + ParameterbitStringLength; j++) x = x + Candidate[j] * (1 << (P - j));//cout<<x<<" ";
        rez.push_back(a + x * b_minus_a_totul_supra_maximul_reprezentabil);
    }
    return rez;
}

long long int randomPosition(const long long int& length){
    return generator() % length;
}

long long int bestPosition(double (*FUNCTIE)(vector<double>), vector<bool>& Candidate, const double& a, const double& b, const long long int& ParameterbitStringLength){
    double best;// = 999999;
    Candidate[0] = 1 - Candidate[0];
    best = FUNCTIE(Transform_into_Mathematical_Object(Candidate, a, b, ParameterbitStringLength));
    Candidate[0] = 1 - Candidate[0];
    double eval = 999999;
    long long int pos = 0;
    long long int length = Candidate.size();
    for(long long int i = 1; i < length; i++){
        Candidate[i] = 1 - Candidate[i];
        eval = FUNCTIE(Transform_into_Mathematical_Object(Candidate, a, b, ParameterbitStringLength));
        Candidate[i] = 1 - Candidate[i];
        if(eval < best) best = eval, pos = i;
    }
    return pos;
}
/*
vector<bool> BestImprovement(double (*FUNCTIE)(vector<double>), vector<bool> Candidate, const double& a, const double& b, const long long int& ParameterbitStringLength) {

    long long int length = Candidate.size();
    double BestValue = FUNCTIE(Transform_into_Mathematical_Object(Candidate, a, b, ParameterbitStringLength));
    long long int pos = 0;
    double eval = 999999;
    for (long long int i = 0; i < length; i++) {
        Candidate[i] = 1 - Candidate[i];
        eval = FUNCTIE(Transform_into_Mathematical_Object(Candidate, a, b, ParameterbitStringLength));
        if (eval < BestValue) {BestValue = eval; pos = i;}
        Candidate[i] = 1 - Candidate[i];
    }
    Candidate[pos] = 1 - Candidate[pos];
    return Candidate;
}*/

double HillClimbing_BestImprovement(double (*FUNCTIE)(vector<double>), long long int& max_repeats, long long int& nrParameters, double& a, double& b) {

    long long int d = 5;// pt ca vrem 5 zecimale dupa virgula
    long long int ParameterbitStringLength = ceil(log2((b - a) * pow(10, d)));
    vector<bool> Candidate;
    double bestValue = 999999;
    double eval;
    long long int length = ParameterbitStringLength * nrParameters;
    long long int pos = 0;
    double CandidateVal, ImprovedCandidateVal;//=765;
    
    for (long long int t = 0; t < max_repeats; t++) {

        //initializeRandom(t, generator);
        Candidate = generateRandomBitString(generator, ParameterbitStringLength * nrParameters);
        CandidateVal = FUNCTIE(Transform_into_Mathematical_Object(Candidate, a, b, ParameterbitStringLength));
        pos = 0;
        ImprovedCandidateVal = CandidateVal;
        while (true) {
            
            //Candidate[0] = 1 - Candidate[0], ImprovedCandidateVal = FUNCTIE(Transform_into_Mathematical_Object(Candidate, a, b, ParameterbitStringLength)), Candidate[0] = 1 - Candidate[0];
            ImprovedCandidateVal = FUNCTIE(Transform_into_Mathematical_Object(Candidate, a, b, ParameterbitStringLength));
            pos = 0;
            for (long long int i = 0; i < length; i++) {
                Candidate[i] = 1 - Candidate[i];
                eval = FUNCTIE(Transform_into_Mathematical_Object(Candidate, a, b, ParameterbitStringLength));
                Candidate[i] = 1 - Candidate[i];
                if (eval < ImprovedCandidateVal) ImprovedCandidateVal = eval, pos = i;
            }
            if (ImprovedCandidateVal < CandidateVal) Candidate[pos] = 1 - Candidate[pos], CandidateVal = ImprovedCandidateVal;
            else break;
        }
        if (CandidateVal < bestValue) bestValue = CandidateVal;
    }
    return bestValue;
}

double HillClimbing_FirstImprovement(double (*FUNCTIE)(vector<double>), long long int& max_repeats, long long int& nrParameters, double& a, double& b) {

    long long int d = 5;// pt ca vrem 5 zecimale dupa virgula
    long long int ParameterbitStringLength = ceil(log2((b - a) * pow(10, d)));
    vector<bool> Candidate;
    double bestValue = 999999;
    double eval;
    long long int length = ParameterbitStringLength * nrParameters;
    double CandidateVal, ImprovedCandidateVal;
    double reserveCandidateVal;

    for (long long int t = 0; t < max_repeats; t++) {

        //initializeRandom(t, generator);
        Candidate = generateRandomBitString(generator, ParameterbitStringLength * nrParameters);
        CandidateVal = FUNCTIE(Transform_into_Mathematical_Object(Candidate, a, b, ParameterbitStringLength));
        ImprovedCandidateVal = CandidateVal;
        //reserveCandidateVal = CandidateVal;
        while (true) {
            
            //Candidate[0] = 1 - Candidate[0], ImprovedCandidateVal = FUNCTIE(Transform_into_Mathematical_Object(Candidate, a, b, ParameterbitStringLength)), Candidate[0] = 1 - Candidate[0];
            ImprovedCandidateVal = FUNCTIE(Transform_into_Mathematical_Object(Candidate, a, b, ParameterbitStringLength));
            for (long long int i = 0; i < length; i++) {
                Candidate[i] = 1 - Candidate[i];
                eval = FUNCTIE(Transform_into_Mathematical_Object(Candidate, a, b, ParameterbitStringLength));
                if (eval < ImprovedCandidateVal) {
                    ImprovedCandidateVal = eval;
                    break;
                }Candidate[i] = 1 - Candidate[i];
            }

            if (ImprovedCandidateVal < CandidateVal)CandidateVal = ImprovedCandidateVal;
            else break;
        }
        if (CandidateVal < bestValue) bestValue = CandidateVal;
    }
    return bestValue;
}

double HillClimbing_WorstImprovement(double (*FUNCTIE)(vector<double>), long long int& max_repeats, long long int& nrParameters, double& a, double& b) {

    long long int d = 5;// pt ca vrem 5 zecimale dupa virgula
    long long int ParameterbitStringLength = ceil(log2((b - a) * pow(10, d)));
    vector<bool> Candidate;
    double bestValue = 999999;
    double eval;
    long long int length = ParameterbitStringLength * nrParameters;
    long long int pos = 0;
    double CandidateVal, ImprovedCandidateVal;//=765;
    
    for (long long int t = 0; t < max_repeats; t++) {

        //initializeRandom(t, generator);
        Candidate = generateRandomBitString(generator, ParameterbitStringLength * nrParameters);
        CandidateVal = FUNCTIE(Transform_into_Mathematical_Object(Candidate, a, b, ParameterbitStringLength));
        pos = 0;
        ImprovedCandidateVal = CandidateVal;
        vector<double> HCW(length);
        while (true) {

            
            //Candidate[0] = 1 - Candidate[0], ImprovedCandidateVal = FUNCTIE(Transform_into_Mathematical_Object(Candidate, a, b, ParameterbitStringLength)), Candidate[0] = 1 - Candidate[0];
            //HCW[0] = ImprovedCandidateVal;
            for(long long int i = 0; i < HCW.size(); i++)HCW[i]=999999;
            ImprovedCandidateVal = FUNCTIE(Transform_into_Mathematical_Object(Candidate, a, b, ParameterbitStringLength));
            for (long long int i = 0; i < length; i++) {
                Candidate[i] = 1 - Candidate[i];
                eval = FUNCTIE(Transform_into_Mathematical_Object(Candidate, a, b, ParameterbitStringLength));
                Candidate[i] = 1 - Candidate[i];
                if (eval < ImprovedCandidateVal) HCW[i] = eval;
            }
            
            double maxi=-999999;
            for(long long int ii = 0; ii < length; ii++){if(HCW[ii]>maxi)maxi=HCW[ii], pos=ii;}
            if (maxi < CandidateVal) Candidate[pos] = 1 - Candidate[pos], CandidateVal = maxi;
            else break;
        }
        if (CandidateVal < bestValue) bestValue = CandidateVal;
    }
    return bestValue;
}

double SimulatedAnnealing(double (*FUNCTIE)(vector<double>), double& minTemp, long long int& nrParameters, double& a, double& b, double& startingTemperature){

    long long int d = 5;// pt ca vrem 5 zecimale dupa virgula
    long long int ParameterbitStringLength = ceil(log2((b - a) * pow(10, d)));
    double T = startingTemperature;
    long long int length = nrParameters * ParameterbitStringLength;

    
    vector<bool> Candidate = generateRandomBitString(generator, length);
    
    /*
    vector<bool> v;
    long long int mijloc_interval = ((a + b) / 2 - a) * ((1 << ParameterbitStringLength) - 1) / (b - a);
        while(mijloc_interval){
            v.push_back(mijloc_interval & 1);
            mijloc_interval = mijloc_interval >> 1;
        }

        if(v.size()<ParameterbitStringLength){
            int aux=ParameterbitStringLength-v.size();
            for(int i=1;i<=aux;i++)v.push_back(0);
        }

        for(int i=0;i<(v.size()/2);i++)swap(v[i],v[v.size()-i-1]);

    for(long long int i = 0; i < nrParameters; i++){        
        for(int j=0;j<v.size();j++)Candidate.push_back(v[j]);
    }
        */
    
    long long int randomIndex = 0;
    long long int stepswithoutimprovement = 0;
    double CandidateVal = 0, NeighbourVal = 0;
    double dif = pow(10, -d) / (b - a); //printf("\n%.20f\n", dif);
    double r = 0;
    
    do{
        T = Normal_Distribution(r);
        //T*=0.9;
        //T=logaritmic(r);
        while(true){

            CandidateVal = FUNCTIE(Transform_into_Mathematical_Object(Candidate, a, b, ParameterbitStringLength));

            randomIndex = randomPosition(length);
            Candidate[randomIndex] = 1 - Candidate[randomIndex];
            NeighbourVal = FUNCTIE(Transform_into_Mathematical_Object(Candidate, a, b, ParameterbitStringLength));

            if(!(dis(generator) < exp(-abs(CandidateVal - NeighbourVal) / T) || NeighbourVal < CandidateVal))
                Candidate[randomIndex] = 1 - Candidate[randomIndex];

            if(abs(CandidateVal - NeighbourVal) <= dif)break;
        }
        r++;
    }while(T > minTemp);

    
    double bestValue=999999, ImprovedCandidateVal,eval; int pos;
    CandidateVal = FUNCTIE(Transform_into_Mathematical_Object(Candidate, a, b, ParameterbitStringLength));
        pos = 0;
        ImprovedCandidateVal = CandidateVal;
        while (true) {
            
            //Candidate[0] = 1 - Candidate[0], ImprovedCandidateVal = FUNCTIE(Transform_into_Mathematical_Object(Candidate, a, b, ParameterbitStringLength)), Candidate[0] = 1 - Candidate[0];
            ImprovedCandidateVal = FUNCTIE(Transform_into_Mathematical_Object(Candidate, a, b, ParameterbitStringLength));
            pos = 0;
            for (long long int i = 0; i < length; i++) {
                Candidate[i] = 1 - Candidate[i];
                eval = FUNCTIE(Transform_into_Mathematical_Object(Candidate, a, b, ParameterbitStringLength));
                Candidate[i] = 1 - Candidate[i];
                if (eval < ImprovedCandidateVal) ImprovedCandidateVal = eval, pos = i;
            }
            if (ImprovedCandidateVal < CandidateVal) Candidate[pos] = 1 - Candidate[pos], CandidateVal = ImprovedCandidateVal;
            else break;
        }
        if (CandidateVal < bestValue) bestValue = CandidateVal;


    return bestValue;
}

void Simulated_Annealing_HC_BestImprovement(double (*FUNCTIE)(vector<double>), double& minTemp, long long int& nrParameters, double& a, double& b, double& startingTemperature){
    long long int d = 5;// pt ca vrem 5 zecimale dupa virgula
    long long int ParameterbitStringLength = ceil(log2((b - a) * pow(10, d)));
    double T = startingTemperature;
    vector<bool> Candidate = generateRandomBitString(generator, nrParameters * ParameterbitStringLength);
    long long int Index = 0;
    long long int stepswithoutimprovement = 0;
    double CandidateVal = 0, NeighbourVal = 0;
    double dif = pow(10, -d) / (b - a); //printf("\n%.20f\n", dif);
    double r = 0;
    
    do{
        T = Normal_Distribution(r);
        //T*=0.9;
        //T=logaritmic(r);
        while(true){

            CandidateVal = FUNCTIE(Transform_into_Mathematical_Object(Candidate, a, b, ParameterbitStringLength));

            Index = bestPosition(FUNCTIE, Candidate, a, b, ParameterbitStringLength);
            Candidate[Index] = 1 - Candidate[Index];
            NeighbourVal = FUNCTIE(Transform_into_Mathematical_Object(Candidate, a, b, ParameterbitStringLength));

            //if(!(dis(generator) < exp(-abs(CandidateVal - NeighbourVal) / T) || NeighbourVal < CandidateVal))
            if(NeighbourVal >= CandidateVal && dis(generator) >= exp(-abs(CandidateVal - NeighbourVal) / T))
                Candidate[Index] = 1 - Candidate[Index];

            if(abs(CandidateVal - NeighbourVal) < dif)break;
        }
        r++;
    }while(T > minTemp);

    printf("\n%.5f", CandidateVal);
}


int main(){

    double a, b;
    //a = -100; b = 100; //Easom
    a = -500, b = 500; //Schwefel
    //a=400;b=422;
    //a = 0, b = pi; //Michalewicz
    //a = -5.12, b = 5.12; // DeJong/Rastrigin

    //a=-0.01;b=500;


    long long int nrparam = 30, rep = 1000;
    double temperature = 640, halting_condition_temperature = pow(10, -4), sigma = 400;

    S = temperature;
    simpli1 = temperature * temperature / (sigma * radicaldindoipi);
    simpli2 = sigma * sigma * 2; 

    cout << "Schwefel - Simulated Annealing" << endl;
    //printf("%.2f %.1f\n", simpli1, simpli2);
    cout << "DOMENIUL DE DEFINITIE: [" << a << ", " << b << ']' << endl;
    cout << "NR. DE PARAMETRI FUNCTIE: " << nrparam << endl;
    printf("TEMPERATURA DE START: %.0f -> TEMPERATURA FINALA: %.4f ; SIGMA: %.0f\n", temperature, halting_condition_temperature, sigma);
    //cout << "NR. DE INCERCARI PT. GASIRE MINIM GLOBAL: " << rep << endl << endl;

    pid_t pid;

    for(int k=1;k<=4;k++){
    pid=fork();
    if(pid<0){ perror("fork failed");
            exit(1);}
    if(pid==0){

    auto start = std::chrono::high_resolution_clock::now();

    initializeRandom(0, generator);
    printf("%.5f",SimulatedAnnealing(Schwefel, halting_condition_temperature, nrparam, a, b, temperature));
    //printf("%.10f\n",HillClimbing_BestImprovement(Rastrigin, rep, nrparam, a, b));

    //Simulated_Annealing_HC_BestImprovement(Schwefel, halting_condition_temperature, nrparam, a, b, temperature);
    //HillClimbing_BestImprovement(Michalewicz, rep, nrparam, a, b);

    auto end = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double, std::milli> duration = end - start;
    std::cout << endl << "Execution time: " << duration.count() * 0.001 << " seconds" << std::endl<< "**************"<< endl << endl;
    exit(0);
    }

    }
    for (int p = 1; p <= 4; p++) wait(NULL);
    return 0;

}



/*
int main(int argc, char* argv[]) {
    
    double a, b;
    //a = -100; b = 100; //Easom
    //a = -500, b = 500; //Schwefel
    //a = 0, b = pi; //Michalewicz
    a = -5.12, b = 5.12; // DeJong/Rastrigin
    long long int nrparam = 5, rep = 1000;
    //double temperature = 868, halting_condition_temperature = pow(10, -10), sigma = 362;
    double temperature = 640, halting_condition_temperature = pow(10, -10), sigma = 400;

    S = temperature;
    simpli1 = temperature * temperature / (sigma * radicaldindoipi);
    simpli2 = sigma * sigma * 2; 

    //cout << "Rastrigin" << endl;
    //printf("%.2f %.1f\n", simpli1, simpli2);
    //cout << "DOMENIUL DE DEFINITIE: [" << a << ", " << b << ']' << endl;
    //cout << "NR. DE PARAMETRI FUNCTIE: " << nrparam << endl;
    //printf("TEMPERATURA DE START: %.0f -> TEMPERATURA FINALA: %.10f\n", temperature, halting_condition_temperature);
    //cout << "NR. DE INCERCARI PT. GASIRE MINIM GLOBAL: " << rep << endl << endl;


    pid_t pid; int k;

    for(k = 1; k <= 4; k++){
        pid=fork();
        if(pid<0){ perror("fork failed");exit(1);}

        if(pid == 0){

            if(k == 1){
                cout << "DeJong" << endl;
                //printf("%.2f %.1f\n", simpli1, simpli2);
                cout << "DOMENIUL DE DEFINITIE: [" << a << ", " << b << ']' << endl;
                cout << "NR. DE PARAMETRI FUNCTIE: " << nrparam << endl;
                //printf("TEMPERATURA DE START: %.0f -> TEMPERATURA FINALA: %.10f\n", temperature, halting_condition_temperature);
                cout << "NR. DE INCERCARI PT. GASIRE MINIM GLOBAL: " << rep << endl << endl;


                double statistica_rez[40], statistica_timp[40];
                double st_dev = 0, maxim = -999999, minim = 999999;
                double medie_rez = 0, medie_timp = 0, total_timp = 0;

                for(long long int i = 0; i < 40; i++){

                    //st_dev = 0, minim = -999999, maxim = 999999;
                    //medie_rez = 0, medie_timp = 0;

                    initializeRandom(i, generator);
                    auto start = chrono::high_resolution_clock::now();

                    statistica_rez[i] = HillClimbing_FirstImprovement(DeJong, rep, nrparam, a, b);

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

                printf("\nREZULTAT MEDIU - HC_FI :\n\nVALOARE MINIMA: %.5f\nVALOARE MAXIMA: %.5f\nMEDIE VALORI: %.5f\nDEVIATIE STANDARD: %.5f\nTIMP MEDIU: %.2f SECUNDE\nTOTAL TIMP EXECUTIE PROCES: %.2f\n*****************\n", minim, maxim, medie_rez, st_dev, medie_timp, total_timp);
            }
            else 
            if(k == 2){
                cout << "DeJong" << endl;
                //printf("%.2f %.1f\n", simpli1, simpli2);
                cout << "DOMENIUL DE DEFINITIE: [" << a << ", " << b << ']' << endl;
                cout << "NR. DE PARAMETRI FUNCTIE: " << nrparam << endl;
                //printf("TEMPERATURA DE START: %.0f -> TEMPERATURA FINALA: %.10f\n", temperature, halting_condition_temperature);
                cout << "NR. DE INCERCARI PT. GASIRE MINIM GLOBAL: " << rep << endl << endl;


                double statistica_rez[40], statistica_timp[40];
                double st_dev = 0, maxim = -999999, minim = 999999;
                double medie_rez = 0, medie_timp = 0, total_timp = 0;

                for(long long int i = 0; i < 40; i++){

                    //st_dev = 0, minim = -999999, maxim = 999999;
                    //medie_rez = 0, medie_timp = 0;

                    initializeRandom(i, generator);
                    auto start = chrono::high_resolution_clock::now();

                    statistica_rez[i] = HillClimbing_BestImprovement(DeJong, rep, nrparam, a, b);

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

                printf("\nREZULTAT MEDIU - HC_BI :\n\nVALOARE MINIMA: %.5f\nVALOARE MAXIMA: %.5f\nMEDIE VALORI: %.5f\nDEVIATIE STANDARD: %.5f\nTIMP MEDIU: %.2f SECUNDE\nTOTAL TIMP EXECUTIE PROCES: %.2f\n*****************\n", minim, maxim, medie_rez, st_dev, medie_timp, total_timp);
            }
            else
            if(k == 3){
                cout << "DeJong" << endl;
                //printf("%.2f %.1f\n", simpli1, simpli2);
                cout << "DOMENIUL DE DEFINITIE: [" << a << ", " << b << ']' << endl;
                cout << "NR. DE PARAMETRI FUNCTIE: " << nrparam << endl;
                //printf("TEMPERATURA DE START: %.0f -> TEMPERATURA FINALA: %.10f\n", temperature, halting_condition_temperature);
                cout << "NR. DE INCERCARI PT. GASIRE MINIM GLOBAL: " << rep << endl << endl;


                double statistica_rez[40], statistica_timp[40];
                double st_dev = 0, maxim = -999999, minim = 999999;
                double medie_rez = 0, medie_timp = 0, total_timp = 0;

                for(long long int i = 0; i < 40; i++){

                    //st_dev = 0, minim = -999999, maxim = 999999;
                    //medie_rez = 0, medie_timp = 0;

                    initializeRandom(i, generator);
                    auto start = chrono::high_resolution_clock::now();
//nrparam*=3;
                    //double aa=-5.12,bb=5.12;long long int nrpa=30;
                    statistica_rez[i] = HillClimbing_WorstImprovement(DeJong, rep, nrparam, a, b);

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

                printf("\nREZULTAT MEDIU - HC_WI :\n\nVALOARE MINIMA: %.5f\nVALOARE MAXIMA: %.5f\nMEDIE VALORI: %.5f\nDEVIATIE STANDARD: %.5f\nTIMP MEDIU: %.2f SECUNDE\nTOTAL TIMP EXECUTIE PROCES: %.2f\n*****************\n", minim, maxim, medie_rez, st_dev, medie_timp, total_timp);
            }
            else
            {
                cout << "DeJong" << endl;
                printf("%.2f %.1f\n", simpli1, simpli2);
                cout << "DOMENIUL DE DEFINITIE: [" << a << ", " << b << ']' << endl;
                cout << "NR. DE PARAMETRI FUNCTIE: " << nrparam << endl;
                printf("TEMPERATURA DE START: %.0f -> TEMPERATURA FINALA: %.10f ; SIGMA: %.0f\n", temperature, halting_condition_temperature, sigma);
                //cout << "NR. DE INCERCARI PT. GASIRE MINIM GLOBAL: " << rep << endl << endl;


                //initializeRandom(0,generator);printf("%.5f",SimulatedAnnealing(Schwefel, halting_condition_temperature, nrparam, a, b, temperature));exit(0);



                double statistica_rez[40], statistica_timp[40];
                double st_dev = 0, maxim = -999999, minim = 999999;
                double medie_rez = 0, medie_timp = 0, total_timp = 0;

                for(long long int i = 0; i < 40; i++){

                    //st_dev = 0, minim = -999999, maxim = 999999;
                    //medie_rez = 0, medie_timp = 0;

                    initializeRandom(i, generator);
                    auto start = chrono::high_resolution_clock::now();

                    statistica_rez[i] = SimulatedAnnealing(DeJong, halting_condition_temperature, nrparam, a, b, temperature);

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

                printf("\nREZULTAT MEDIU - SA :\n\nVALOARE MINIMA: %.5f\nVALOARE MAXIMA: %.5f\nMEDIE VALORI: %.5f\nDEVIATIE STANDARD: %.5f\nTIMP MEDIU: %.2f SECUNDE\nTOTAL TIMP EXECUTIE PROCES: %.2f\n*****************\n", minim, maxim, medie_rez, st_dev, medie_timp, total_timp);
            }

            exit(0);
        }
    }
    for (long long int p = 1; p <= 4; p++) wait(NULL);
    return 0;
}
*/


/*
int main(){
    double a,b, nrparam,d;

    //a=-5.12;b=5.12;
    a=-3;b=7;
    double r =(a+b)/2;
    nrparam=30;
    d=5;
    long long Parameterbitstringlength;
    //long long L=
    
    Parameterbitstringlength=ceil(log2((b - a) * pow(10, d)));

    int x = ((r-a)* ( (1<<Parameterbitstringlength)-1 ) )/(b-a);
    cout<<x<<endl;int cop=x;
    vector<bool> v;
    while(cop){
        v.push_back(cop%2);
        cop/=2;
    }



    if(v.size()<Parameterbitstringlength){
        int aux=Parameterbitstringlength-v.size();
        for(int i=1;i<=aux;i++)v.push_back(0);
    }
    for(int i=0;i<(v.size()/2);i++)swap(v[i],v[v.size()-i-1]);

    for(int i =0;i<v.size();i++)cout<<v[i];cout<<" "<<Parameterbitstringlength<< " "<<v.size() <<endl;
    vector<double>rez=Transform_into_Mathematical_Object(v,a,b,v.size());//cout<<endl;
    for(int i=0;i<rez.size();i++)printf("%.20f",rez[0]);

    cout<<endl;
    int x1=2;double x2=1.9999999999999999;
    if(x1==x2)cout<<"DA";else cout<<"NU";


    return 0;
}
*/



/*
int main(){

    double a, b;
    //a = -100; b = 100; //Easom
    //a = -500, b = 500; //Schwefel
    //a = 0, b = pi; //Michalewicz
    a = -5.12, b = 5.12; // DeJong/Rastrigin
    long long int nrparam = 30, rep = 1000;
    double temperature = 868, halting_condition_temperature = pow(10, -10), sigma = 362;

    S = temperature;
    simpli1 = temperature * temperature / (sigma * radicaldindoipi);
    simpli2 = sigma * sigma * 2; 

    cout << "Rastrigin" << endl;
    //printf("%.2f %.1f\n", simpli1, simpli2);
    cout << "DOMENIUL DE DEFINITIE: [" << a << ", " << b << ']' << endl;
    cout << "NR. DE PARAMETRI FUNCTIE: " << nrparam << endl;
    printf("TEMPERATURA DE START: %.0f -> TEMPERATURA FINALA: %.10f ; SIGMA: %.0f\n", temperature, halting_condition_temperature, sigma);
    //cout << "NR. DE INCERCARI PT. GASIRE MINIM GLOBAL: " << rep << endl << endl;

    pid_t pid;

    for(int k=1;k<=4;k++){
    pid=fork();
    if(pid<0){ perror("fork failed");
            exit(1);}
    if(pid==0){

    auto start = std::chrono::high_resolution_clock::now();

    initializeRandom(0, generator);
    printf("%.5f",SimulatedAnnealing(Rastrigin, halting_condition_temperature, nrparam, a, b, temperature));
    //printf("%.5f\n",HillClimbing_BestImprovement(Rastrigin, rep, nrparam, a, b));

    //Simulated_Annealing_HC_BestImprovement(Schwefel, halting_condition_temperature, nrparam, a, b, temperature);
    //HillClimbing_BestImprovement(Michalewicz, rep, nrparam, a, b);

    auto end = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double, std::milli> duration = end - start;
    std::cout << endl << "Execution time: " << duration.count() * 0.001 << " seconds" << std::endl<< "**************"<< endl << endl;
    exit(0);
    }

    }
    for (int p = 1; p <= 4; p++) wait(NULL);
    return 0;

}
*/


/*

int main() {
   
/*
    



    vector<bool>test(35);
    long long int f[35];for(int i=0;i<35;i++)f[i]=0;

auto start = std::chrono::high_resolution_clock::now();

initializeRandom(2,generator);
    for(int i=1;i<99999;i++){
        
    f[randomPosition(test)]++;
    }
   for(int i=0;i<35;i++)cout<<f[i]<<" ";

//long long int rez = 0;
//long long int cop = test.size() - 1; while(cop) cop = cop >> 1, rez++; //ceil(log2)
//cout<<rez;

   auto end = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double, std::milli> duration = end - start;
    std::cout << endl << "Execution time: " << duration.count() * 0.001 << " seconds" << std::endl<< "**************"<< endl << endl;

    return 0;
}














int main(){

    double a, b;
    //a = -100; b = 100; //Easom
    //a = -500, b = 500; //Schwefel
    //a = 0, b = pi; //Michalewicz
    a = -5.12, b = 5.12; // DeJong/Rastrigin
    long long int nrparam = 30, rep = 1000;
    double temperature = 571, halting_condition_temperature = pow(10, -10), sigma = 364;

    S = temperature;
    simpli1 = temperature * temperature / (sigma * radicaldindoipi);
    simpli2 = sigma * sigma * 2; 

    cout << "Rastrigin" << endl;
    //printf("%.2f %.1f\n", simpli1, simpli2);
    cout << "DOMENIUL DE DEFINITIE: [" << a << ", " << b << ']' << endl;
    cout << "NR. DE PARAMETRI FUNCTIE: " << nrparam << endl;
    printf("TEMPERATURA DE START: %.0f -> TEMPERATURA FINALA: %.10f\n", temperature, halting_condition_temperature);
    //cout << "NR. DE INCERCARI PT. GASIRE MINIM GLOBAL: " << rep << endl << endl;

    pid_t pid;

    for(int k=1;k<=4;k++){
    pid=fork();
    if(pid<0){ perror("fork failed");
            exit(1);}
    if(pid==0){

    auto start = std::chrono::high_resolution_clock::now();

    initializeRandom(0, generator);
    printf("%.5f",SimulatedAnnealing(Rastrigin, halting_condition_temperature, nrparam, a, b, temperature));
    //printf("%.5f\n",HillClimbing_BestImprovement(Rastrigin, rep, nrparam, a, b));

    //Simulated_Annealing_HC_BestImprovement(Schwefel, halting_condition_temperature, nrparam, a, b, temperature);
    //HillClimbing_BestImprovement(Michalewicz, rep, nrparam, a, b);

    auto end = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double, std::milli> duration = end - start;
    std::cout << endl << "Execution time: " << duration.count() * 0.001 << " seconds" << std::endl<< "**************"<< endl << endl;
    exit(0);
    }

    }
    for (int p = 1; p <= 4; p++) wait(NULL);
    return 0;

}





























int main(){

    double a, b;
    //a = -100; b = 100; //Easom
    //a = -500, b = 500; //Schwefel
    //a = 0, b = pi; //Michalewicz
    a = -5.12, b = 5.12; // DeJong/Rastrigin
    long long int nrparam = 30, rep = 1000;
    double temperature = 640, halting_condition_temperature = pow(10, -10), sigma = 400;

    S = temperature;
    simpli1 = temperature * temperature / (sigma * radicaldindoipi);
    simpli2 = sigma * sigma * 2; 

    cout << "Rastrigin" << endl;
    //printf("%.2f %.1f\n", simpli1, simpli2);
    cout << "DOMENIUL DE DEFINITIE: [" << a << ", " << b << ']' << endl;
    cout << "NR. DE PARAMETRI FUNCTIE: " << nrparam << endl;
    printf("TEMPERATURA DE START: %.0f -> TEMPERATURA FINALA: %.10f\n", temperature, halting_condition_temperature);
    //cout << "NR. DE INCERCARI PT. GASIRE MINIM GLOBAL: " << rep << endl << endl;

    pid_t pid;

    for(int k=1;k<=4;k++){
    pid=fork();
    if(pid<0){ perror("fork failed");
            exit(1);}
    if(pid==0){

    auto start = std::chrono::high_resolution_clock::now();

    double best_temperature, best_sigma, global_minimum=999999;



        for(long long int time=0;time<301;time++){
                initializeRandom(time + (long long int)((clock()*1000)), generator);

                temperature= randomPosition(1000);
                sigma= randomPosition(1000);                   


                simpli1 = temperature * temperature / (sigma * radicaldindoipi);
                simpli2 = sigma * sigma * 2; 

                    //halting_condition_temperature=pow(10,-h);
                    double rez=SimulatedAnnealing(Rastrigin, halting_condition_temperature, nrparam, a, b, temperature);
                    //printf("%.5f\n", rez);
                    if(rez<global_minimum){
                        global_minimum=rez;
                        best_temperature=temperature;
                        best_sigma=sigma;
                    }
                    printf("intermediar:temperatura: %.5f   sigma: %.5f    minim global: %.5f\n", temperature, sigma, rez);
                    //printf("intermediar:temperatura: %.5f   sigma: %.5f    minim global: %.5f\n", best_temperature, best_sigma, global_minimum);

        }

        printf("final:temperatura: %.5f   sigma: %.5f    minim global: %.5f\n", best_temperature, best_sigma, global_minimum);


    //initializeRandom(0, generator);
    //printf("%.5f",SimulatedAnnealing(Rastrigin, halting_condition_temperature, nrparam, a, b, temperature));





    //printf("%.5f\n",HillClimbing_BestImprovement(Rastrigin, rep, nrparam, a, b));

    //Simulated_Annealing_HC_BestImprovement(Schwefel, halting_condition_temperature, nrparam, a, b, temperature);
    //HillClimbing_BestImprovement(Michalewicz, rep, nrparam, a, b);

    auto end = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double, std::milli> duration = end - start;
    std::cout << endl << "Execution time: " << duration.count() * 0.001 << " seconds" << std::endl<< "**************"<< endl << endl;
    exit(0);
    }

    }
    for (int p = 1; p <= 4; p++) wait(NULL);
    return 0;

}











*/