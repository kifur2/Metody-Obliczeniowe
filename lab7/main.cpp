#include <iostream>
#include <cmath>
#include <vector>
#include <iomanip>

#define N 4
#define W (0.5)
#define NMAX 50
#define BLAD 1e-12
using namespace std;
struct dane {
    long double x[N], e[N], r[N];
};
vector<dane> jacobXn;
vector<dane> gauss1Xn;
vector<dane> gauss2Xn;
vector<dane> sorXn;
dane xn;
long double A[N][N] = {
        100, -1, 2, -3,
        1, 200, -4, 5,
        -2, 4, 300, -6,
        3, -5, 6, 400
};
long double jacobM[N][N];
long double gauss1M[N][N];
long double gauss2M[N][N];
long double sorM[N][N];
long double b[N] = {116, -226, 912, -1174};

void printMatrix(long double tab[][N], const string &message) {

    cout << message << endl;

    for (int i = 0; i < N; i++) {
        for (int j = 0; j < N; j++) {
            cout << setw(15) << tab[i][j] << ' ';
        }
        cout << endl;
    }
}

bool printAll(long double estymator, long double reziduum){
    bool out = false;
    if (estymator < BLAD) {
        for (int i = 0; i < 3; i++)
            cout << setw(15) << "TOLX";
        cout << endl << endl;
        out = true;
    } else if (reziduum < BLAD) {
        for (int i = 0; i < 3; i++)
            cout << setw(15) << "TOLF";
        cout << endl << endl;
        out = true;

    } else {

        for (int i = 0; i < N; i++) {
            cout << "[" << setw(15) << xn.x[i] << "]";
            cout << "[" << setw(15) << xn.e[i] << "]";
            cout << "[" << setw(15) << xn.r[i] << "]"<<endl;
        }
        cout << endl;
    }
    return out;
}

void jacob() {
    cout << "jacob: " << endl;
    cout << setw(15) << "x" << setw(15) << "estym" << setw(15) << "reziduum" << endl;

    for (int k = 1; k <= NMAX; k++) {

        cout << k << ": " << endl;
        long double estymator = 0, reziduum = 0;

        for (int i = 0; i < N; i++) {
            xn.x[i] = 0;
            for (int j = 0; j < N; j++) {
                xn.x[i] += jacobM[i][j] * jacobXn[jacobXn.size() - 1].x[j];
            }
            xn.x[i] += b[i] / A[i][i] ;
        }

        for (int i = 0; i < N; i++) {
            xn.e[i] = fabs(jacobXn[jacobXn.size() - 1].x[i] - xn.x[i]);
            estymator = max(xn.e[i], estymator);
        }
        for (int i = 0; i < N; i++) {
            xn.r[i] = 0;
            for (int j = 0; j < N; j++) {
                xn.r[i] += A[i][j] * xn.x[j];
            }
            xn.r[i]-=b[i];
            xn.r[i] = fabs(xn.r[i]);
            reziduum = max(xn.r[i], reziduum);
        }
        jacobXn.push_back(xn);

        if(printAll(estymator, reziduum))
            break;
    }
}

void gauss1() {

    cout << "gauss1: " << endl;
    cout << setw(15) << "x" << setw(15) << "estym" << setw(15) << "reziduum" << endl;

    for (int k = 1; k <= NMAX; k++) {

        cout << k << ": " << endl;

        long double estymator = 0, reziduum = 0;
        for (int i = 0; i < N; i++) {
            xn.x[i] = 0;
            for (int j = i + 1; j < N; j++) {
                xn.x[i] += gauss1M[i][j] * gauss1Xn[gauss1Xn.size() - 1].x[j];
            }
            xn.x[i] += b[i];
        }

        for (int i = 0; i < N; i++) {
            for (int j = 0; j < i; j++) {
                xn.x[i] -= xn.x[j] * A[i][j];
            }
            xn.x[i] /= A[i][i];
        }
        for (int i = 0; i < N; i++) {
            xn.e[i] = fabs(gauss1Xn[gauss1Xn.size() - 1].x[i] - xn.x[i]);
            estymator = max(xn.e[i], estymator);
        }

        for (int i = 0; i < N; i++) {
            xn.r[i] = 0;
            for (int j = 0; j < N; j++) {
                xn.r[i] += A[i][j] * xn.x[j];
            }
            xn.r[i]-=b[i];
            xn.r[i] = fabs(xn.r[i]);
            reziduum = max(xn.r[i], reziduum);
        }

        gauss1Xn.push_back(xn);

        if(printAll(estymator, reziduum))
            break;
    }
}

void gauss2() {

    cout << "gauss2: " << endl;
    cout << setw(15) << "x" << setw(15) << "estym" << setw(15) << "reziduum" << endl;

    for (int k = 1; k <= NMAX; k++) {

        cout << k << ": " << endl;

        long double estymator = 0, reziduum = 0;
        for (int i = 0; i < N; i++) {
            xn.x[i] = 0;
            for (int j = 0; j < i; j++) {
                xn.x[i] += gauss2M[i][j] * gauss2Xn[gauss2Xn.size() - 1].x[j];
            }
            xn.x[i] += b[i];
        }
        for (int i = N - 1; i >= 0; i--) {
            for (int j = N - 1; j > i; j--) {
                xn.x[i] -= xn.x[j] * A[i][j];
            }
            xn.x[i] /= A[i][i];
        }
        for (int i = 0; i < N; i++) {
            xn.e[i] = fabs(gauss2Xn[gauss2Xn.size() - 1].x[i] - xn.x[i]);
            estymator = max(xn.e[i], estymator);
        }

        for (int i = 0; i < N; i++) {
            xn.r[i] = 0;
            for (int j = 0; j < N; j++) {
                xn.r[i] += A[i][j] * xn.x[j];
            }
            xn.r[i]-=b[i];
            xn.r[i] = fabs(xn.r[i]);
            reziduum = max(xn.r[i], reziduum);
        }
        gauss2Xn.push_back(xn);

        if(printAll(estymator, reziduum))
            break;
    }
}

void sor() {

    cout << "sor: " << endl;
    cout << setw(15) << "x" << setw(15) << "estym" << setw(15) << "reziduum" << endl;

    for (int k = 1; k <= NMAX; k++) {

        cout << k << ": " << endl;

        long double estymator = 0, reziduum = 0;
        for (int i = 0; i < N; i++) {
            xn.x[i] = 0;
            for (int j = i; j < N; j++) {
                    xn.x[i] += sorM[i][j] * sorXn[sorXn.size() - 1].x[j];
            }
            xn.x[i] += b[i];
        }

        for (int i = 0; i < N; i++) {
            for (int j = 0; j < i; j++) {
                xn.x[i] -= xn.x[j] * A[i][j];
            }
            xn.x[i] /= (A[i][i]/W);
        }
        for (int i = 0; i < N; i++) {
            xn.e[i] = fabs(sorXn[sorXn.size() - 1].x[i] - xn.x[i]);
            estymator = max(xn.e[i], estymator);
        }
        for (int i = 0; i < N; i++) {
            xn.r[i] = 0;
            for (int j = 0; j < N; j++) {
                xn.r[i] += A[i][j] * xn.x[j];
            }
            xn.r[i]-=b[i];
            xn.r[i] = fabs(xn.r[i]);
            reziduum = max(xn.r[i], reziduum);
        }

        sorXn.push_back(xn);

        if(printAll(estymator, reziduum))
            break;
    }
}

void allCountM() {

    for (int i = 0; i < N; i++)
        for (int j = 0; j < N; j++) {

            /* Metoda Jacobiego */
            if (i != j)
                jacobM[i][j] = -1. / A[i][i] * A[i][j];

            /* Metoda Gaussa - Seidela 1 */
            if (i < j)
                gauss1M[i][j] = -A[i][j];

            /* Metoda Gaussa - Seidela 2 */
            if (i > j)
                gauss2M[i][j] = -A[i][j];

            /* Metoda sukcesywnej nadrelaksacji */
            if (i == j)
                sorM[i][j] = -((1. - 1.0 / W) * A[i][j]);
            if (i < j)
                sorM[i][j] = -A[i][j];

        }
    printMatrix(jacobM, "jacobM: ");
    printMatrix(gauss1M, "gauss1M: ");
    printMatrix(gauss2M, "gauss2M: ");
    printMatrix(sorM, "sorM: ");
}

void setAll() {
    for (long double &i: xn.x)
        i = 2;
    jacobXn.push_back(xn);
    gauss1Xn.push_back(xn);
    gauss2Xn.push_back(xn);
    sorXn.push_back(xn);
}

int main() {
    setAll();
    allCountM();
    jacob();
    gauss1();
    gauss2();
    sor();
}

