#include <iostream>
#include <vector>
#include <cmath>
#include <fstream>
#define NMAX 30
using namespace std;

vector<double> BME;
vector <double> PME;
vector<double> PMT;
vector<double> ANALYTICAL;

double fun(double t){
    return (10*t*t +20.)/(t*t+1);
}

void funBME(double t, double dt){
    for(int i = 0; i<NMAX; i++) {
        BME.push_back(BME[BME.size() - 1] - fun(t) * (BME[BME.size() - 1] - 1) * dt);
        t+=dt;
    }
}
void funPME(double t, double dt){
    for(int i = 0; i<NMAX; i++) {
        PME.push_back((PME[PME.size() - 1] * ((t + dt) * (t + dt) + 1.) + (10. * (t + dt) * (t + dt) + 20.) * dt) /
                      ((t + dt) * (t + dt) + 1. + (10. * (t + dt) * (t + dt) + 20.) * dt));
        t+=dt;
    }
}
void funPMT(double t, double dt){
    for(int i = 0; i<NMAX; i++) {
        PMT.push_back((2 * PMT[PMT.size() - 1] + dt * (fun(t) * (1 - PMT[PMT.size() - 1]) + fun(t + dt))) /
                      (2 + fun(t + dt) * dt));
        t+=dt;
    }
}
void funANALYTICAL(double t, double dt){
    for(int i = 0; i<NMAX; i++) {
        ANALYTICAL.push_back(1. - exp(-10. * (t + atan(t))));
        t+=dt;
    }
}



void vectorsClear(){
    BME.clear();
    PME.clear();
    PMT.clear();
    ANALYTICAL.clear();
    BME.push_back(0);
    PME.push_back(0);
    PMT.push_back(0);
}
int main() {
    vectorsClear();
    double dt1 = 0.1;
    double dt2 = 0.25;
    ofstream plik("../analityczne1.txt");
    ofstream plik0("../analityczne2.txt");
    ofstream plik1("../stabilneBME.txt");
    ofstream plik2("../niestabilneBME.txt");
    ofstream plik3("../PME.txt");
    ofstream plik4("../PMT.txt");
    ofstream plik5("../bladStabilneBME.txt");
    ofstream plik6("../bladNiestabilneBME.txt");
    ofstream plik7("../bladPME.txt");
    ofstream plik8("../bladPMT.txt");
    double t = 0;
    funANALYTICAL(t, dt1);
    funBME(t, dt1);
    funPME(t, dt1);
    funPMT(t,  dt1);
    for(int i = 0; i<NMAX; i++){
        plik<<t<<' '<<ANALYTICAL[i]<<endl;
        plik1<<t<<' '<<BME[i]<<endl;
        plik3<<t<<' '<<PME[i]<<endl;
        plik4<<t<<' '<<PMT[i]<<endl;
        plik5<<log10(t)<<' '<<BME[i]-ANALYTICAL[i]<<endl;
        plik7<<log10(t)<<' '<<PME[i]-ANALYTICAL[i]<<endl;
        plik8<<log10(t)<<' '<<PMT[i]-ANALYTICAL[i]<<endl;
        t +=dt1;
    }
    vectorsClear();
    t=0;
    funANALYTICAL(t, dt2);
    funBME(t, dt2);
    for(int i = 0; i<NMAX; i++){
        plik0<<t<<' '<<ANALYTICAL[i]<<endl;
        plik2<<t<<' '<<BME[i]<<endl;
        plik6<<log10(t)<<' '<<BME[i]-ANALYTICAL[i]<<endl;
        t +=dt2;
    }
}
