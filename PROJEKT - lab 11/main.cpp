#include <iostream>
#include <cmath>
#include <vector>
#include <cfloat>
#include <fstream>
#include <ctime>

using namespace std;

const double TMIN = 0, TMAX = 0.5, D = 1, XMIN = 0, XMAX = 1;
const double ALFA = 0.0, BETA = 1.0, GAMMA = 0;
const double FI = 0.0, PSI = 1.0, TETA = 0.0;
const double eps = 10e-12;
vector<vector<double>> kmbU;
vector<vector<double>> AnalyticalU;
vector<vector<double>> laasonenThomasU;
vector<vector<double>> laasonenLuU;

double Ufun(double t, double x) {
    if (t == 0) {
        return sin(M_PI * x);
    } else if (x == 0 || x == 1) {
        return 0;
    } else {
        return exp(-M_PI * M_PI * D * t) * sin(M_PI * x);
    }
}

void init(int N, int M, double h, vector<vector<double>> &tab) {
    for (int k = 0; k < M; k++) {
        vector<double> pusty(N, 0);
        for (int i = 0; i < N; i++) {
            tab.push_back(pusty);
        }
    }
    double x = h;
    for (int i = 1; i < N - 1; i++) {
        tab[0][i] = Ufun(0, x);
        x += h;
    }
}

bool Thomas(int N, const double L[], const double U[], double DW[], double BR[], int k) {
    for (int i = 2; i <= N; i++) {
        if (fabs(DW[i]) < eps)
            return false;
        DW[i] = DW[i] - L[i] / DW[i - 1] * U[i - 1];
    }
    for (int i = 2; i <= N; i++) {
        if (fabs(DW[i - 1]) < eps)
            return false;
        BR[i] = BR[i] - L[i] / DW[i - 1] * BR[i - 1];
    }

    laasonenThomasU[k][N] = 1.0 / DW[N] * BR[N];
    for (int i = N - 1; i > 0; i--) {
        laasonenThomasU[k][i] = 1.0 / DW[i] * (BR[i] - U[i] * laasonenThomasU[k][i + 1]);
    }
    return true;
}

bool Aeliminacja(int N, vector<vector<double> > &U, vector<vector<double> > &L, int swapping[]) {
    for (int i = 0; i < N - 1; i++) {
        for (int j = i + 1; j < N; j++) {
            if (fabs(U[i][i]) < eps) {
                double maxi = -DBL_MAX;
                int i_maxi = i;
                for (int k = i + 1; k < N; k++) {
                    if (U[k][i] > maxi && fabs(U[k][i]) > eps) {
                        maxi = U[k][i];
                        i_maxi = k;
                    }
                }
                if (fabs(maxi) < eps)
                    return false;
                for (int k = i; k < N; k++)
                    swap(U[i][k], U[i_maxi][k]);
                swapping[i] = i_maxi;
                if (i > 0)
                    swap(L[i][i - 1], L[i_maxi][i - 1]);
            }
            double iloraz = -U[j][i] / U[i][i];
            if (fabs(iloraz) < eps) {
                iloraz = 0;
            }
            L[j][i] -= iloraz;
            U[j][i] = 0;
            for (int k = i + 1; k < N; k++) {
                U[j][k] += iloraz * U[i][k];
            }
        }
    }
    return true;
}

void Beliminacja(int N, double B[], vector<vector<double> > &L, const int swapping[]) {
    for (int i = 0; i < N - 1; i++) {
        swap(B[i], B[swapping[i]]);

        for (int j = i + 1; j < N; j++) {
            swap(L[j][i], L[swapping[j]][i]);
            B[j] -= L[j][i] * B[i];
        }
        for (int j = i + 1; j < N; j++) {
            swap(L[j][i], L[swapping[j]][i]);
        }
    }
}

void obliczYX(int N, double X[], double Y[], const double B[], vector<vector<double> > &L, vector<vector<double> > &U) {
    for (int i = 0; i < N; i++) {
        Y[i] = B[i];
        for (int j = 0; j < i; j++) {
            Y[i] -= L[j][i] * Y[j];
        }
    }
    for (int i = N - 1; i >= 0; i--) {
        X[i] = Y[i] / U[i][i];
        for (int j = N - 1; j > i; j--) {
            X[i] -= U[i][j] * X[j] / U[i][i];
        }
    }
}

void LaasonenThomas(int N, int M, double h, double dt) {

    double lambda = D * dt / (h * h);
    double tU[N], tDW[N + 1], tL[N + 1], tBR[N + 1];

    init(N, M, h, laasonenThomasU);

    for (int k = 0; k < M - 1; k++) {
        tU[1] = ALFA / h;
        tDW[1] = BETA - ALFA / h;
        tBR[1] = -GAMMA;

        for (int i = 2; i < N - 1; i++) {
            tL[i] = lambda;
            tDW[i] = -(1 + 2 * lambda);
            tU[i] = lambda;
            tBR[i] = -laasonenThomasU[k][i];
        }

        tL[N] = -FI / h;
        tDW[N] = FI / h + PSI;
        tBR[N] = -TETA;
        Thomas(N, tL, tU, tDW, tBR, k);
    }
}

void LaasonenLU(int N, int M, double h, double dt){

    double lambda = D * dt / (h * h);

    vector<vector<double>> luA, luL, luU;
    double luX[N], luY[N], luB[N];
    int swapping[N];

    init(N, M, h, laasonenLuU);

    for (int k = 0; k < N; k++) {
        vector<double> pusty(N, 0);
        for (int i = 0; i < N; i++) {
            luA.push_back(pusty);
            luL.push_back(pusty);
            luU.push_back(pusty);
        }
    }
    for(int i =0; i<N; i++){
        for(int j =0; j<N; j++){
            if(j == i-1 || j == i+1){
                luA[i][j]= lambda;
            }else if(j == i){
                luA[i][j] = -(1+2*lambda);
            }else {
                luA[i][j] = 0;
            }
        }
    }
    luA[0][0] = BETA - ALFA/h;
    luA[0][1] = ALFA/h;
    luA[N-1][N-2] = -FI/h;
    luA[N-1][N-1] = FI/h+PSI;

    for (int i = 0; i < N; i++) {
        for (int j = 0; j < N; j++) {
            luU[i][j] = luA[i][j];
        }
    }

    for (int i = 0; i < N; i++) {
        luY[i] = 0;
        luX[i] = 0;
        swapping[i] = i;
        for (int j = 0; j < N; j++) {
            if (i == j)
                luL[i][j] = 1;
            else
                luL[i][j] = 0;
        }
    }
    Aeliminacja(N,luU, luL, swapping);

    for (int k = 0; k < M - 1; k++) {

        luB[0] = -GAMMA;
        for(int i = 1; i<N-1; i++) {
            luB[i] = -laasonenThomasU[k][i];
        }
        luB[N-1]= -TETA;

        Beliminacja(N, luB, luL, swapping);
        obliczYX(N, luX, luY, luB, luL, luU);
    }
}


void KMB(int N, int M, double h, double dt) {
    double lambda = D * dt / (h * h);
    init(N, M, h, kmbU);
    for (int k = 0; k < M - 1; k++) {
        for (int i = 1; i < N - 1; i++) {
            kmbU[k + 1][i] = lambda * kmbU[k][i - 1] + (1 - 2 * lambda) * kmbU[k][i] + lambda * kmbU[k][i + 1];
        }
    }
}
void Analytical(int N, int M, double h, double dt) {
    init(N, M, h, AnalyticalU);
    double t = dt;
    for (int k = 0; k < M - 1; k++) {
        double x = h;
        for (int i = 1; i < N - 1; i++) {
            kmbU[k + 1][i] = Ufun(t,x);
            x+=h;
        }
        t+=dt;
    }
}
double Liczenie1(vector<vector<double> > & tab, int N, int M){
    double maksBlad = 0;
    for(int i = 0; i< N; i++){
        double tmpBlad = fabs(tab[M][i] - AnalyticalU[M][i]);
        if(tmpBlad>maksBlad){
            maksBlad = tmpBlad;
        }
    }
    return maksBlad;
}

int main() {

    ofstream plik("../KMB_vs_Analytical.txt");
    ofstream plik1("../LaasonenThomas_vs_Analytical.txt");
    ofstream plik2("../LaasonenLU_vs_Analytical.txt");
    //for(double h = 0.001; h<(XMAX-XMIN)/2; h+=0.001) {
    double h = 0.0025;
        double dt = 0.4 * h * h;
        int N = (int) ((XMAX - XMIN) / h + 1.0);
        int M = (int) ((TMAX - TMIN) / dt + 1.0);
        cout << h << ' ' << dt << ' ' << N << ' ' << M << endl;
        time_t timer = time(NULL);
        KMB(N, M, h, dt);
        cout<<"time: "<<time(NULL)-timer<<endl;

        timer = time(NULL);
        Analytical(N, M, h, dt);
        cout<<"time: "<<time(NULL)-timer<<endl;

        timer = time(NULL);
        LaasonenThomas(N, M, h, dt);
        cout<<"time: "<<time(NULL)-timer<<endl;

        timer = time(NULL);
        LaasonenLU(N, M, h, dt);
        cout<<"time: "<<time(NULL)-timer<<endl;

        plik << log10(h) <<' '<< log10(Liczenie1(kmbU, N, M)) << endl;
        plik1 << log10(h) <<' '<< log10(Liczenie1(laasonenThomasU, N, M)) << endl;
        plik2 << log10(h) <<' '<< log10(Liczenie1(laasonenLuU, N, M)) << endl;
    //}
    return 0;
}
