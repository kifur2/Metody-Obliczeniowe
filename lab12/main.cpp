// lab12.cpp : This file contains the 'main' function. Program execution begins and ends there.
//

#include <iostream>
#include <vector>
#include <cmath>
#include <fstream>

#define N 60
#define N1 60
using namespace std;

constexpr double x_pocz = -1;
constexpr double x_kon = 1;

vector<double> argumenty[2];
vector<double> ilorazy[2][N1];


double function(double x) {
    return 1.0 / (1 + 10 * x * x * x * x * x * x);
}

//obliczanie wielomianu interpolacyjnego
double horner0(double x) {
    double wynik = ilorazy[0][N - 1][0] * (x - argumenty[0][N - 2]) + ilorazy[0][N - 2][0];
    for (int i = N - 3; i >= 0; i--) {
        wynik *= (x - argumenty[0][i]);
        wynik += ilorazy[0][i][0];
    }
    return wynik;
}

double horner1(double x) {
    double wynik = ilorazy[1][N1 - 1][0] * (x - argumenty[1][N1 - 2]) + ilorazy[1][N1 - 2][0];
    for (int i = N1 - 3; i >= 0; i--) {
        wynik *= (x - argumenty[1][i]);
        wynik += ilorazy[1][i][0];
    }
    return wynik;
}

int main() {
    double h = (x_kon - x_pocz) / (N - 1);
    double x0 = x_pocz;
    double x1;
    //obliczanie argumentów i wartości funkcji
    for (int i = 0; i < N; i++) {
        argumenty[0].push_back(x0);
        ilorazy[0][0].push_back(function(x0));
        x0 += h;
    }
    for (int i = 0; i < N1; i++) {

        x1 = (x_kon + x_pocz) / 2 + (x_kon - x_pocz) / 2 * cos((2. * i + 1) / (2 * N1 + 2) * M_PI);
        argumenty[1].push_back(x1);
        ilorazy[1][0].push_back(function(x1));

    }

    // obliczanie reszty tabelki
    int j = N - 1;
    while (j > 0) {
        for (int i = 0; i < j; i++) {
            ilorazy[0][N - j].push_back((ilorazy[0][N - j - 1][i + 1] - ilorazy[0][N - j - 1][i]) /
                                        (argumenty[0][i + N - j] - argumenty[0][i]));
            //cout<<"ilorazy[0]["<<N-j<<"]["<<i<<"] = "<<ilorazy[0][N-j][i]<<endl;
        }
        j--;
    }

    j = N1 - 1;
    while (j > 0) {
        for (int i = 0; i < j; i++) {
            ilorazy[1][N1 - j].push_back((ilorazy[1][N1 - j - 1][i + 1] - ilorazy[1][N1 - j - 1][i]) /
                                         (argumenty[1][i + N1 - j] - argumenty[1][i]));
        }
        j--;
    }

    ofstream wielomianR("../wielomianRownoodlegle.txt");
    ofstream wielomianC("../wielomianCzybyszew.txt");
    ofstream funkcjaR("../funkcjaRownoodlegle.txt");
    ofstream funkcjaC("../funkcjaCzybyszew.txt");


    //double h1 = (x_kon - x_pocz) / (N1 - 1);
    x0 = x_pocz;
    for (int i = 0; i < N; i++) {

        cout << x0 << ' ' << horner0(x0) << "\t";
        cout << function(x0) << endl;

        wielomianR << x0 << ' ' << horner0(x0) << endl;
        funkcjaR << x0 << ' ' << function(x0) << endl;
        x0 += h;
    }
    for (int i = 0; i < N1; i++) {

        x1 = (x_kon + x_pocz) / 2 + (x_kon - x_pocz) / 2 * cos((2. * i + 1) / (2 * N1 + 2) * M_PI);

        wielomianC << x1 << ' ' << horner1(x1) << endl;
        funkcjaC << x1 << ' ' << function(x1) << endl;
    }

}

