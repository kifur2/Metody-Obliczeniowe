#include <iostream>
#include <cmath>
#include <iomanip>
#include <float.h>

using namespace std;

const int n = 4;
double B[n]={35,
             104,
             -366,
             -354};
double L[n][n];
double U[n][n];
int swapping[n];

bool Aeliminacja(){
   double eps = 10e-12;
   for(int i = 0; i < n-1; i++){

       if(fabs(U[i][i])<eps){

           int i_max = i;
           double maxi = 0;
           for(int j = i+1; j<n; j++){
               if(fabs(U[j][i])>maxi){
                   i_max = j;
                   maxi = fabs(U[j][i]);
               }
           }

           for(int j = 0; j<n; j++){
               swap(U[i][j], U[i_max][j]);
           }
           swapping[i] = swapping[i_max];
           for(int j = 0; j<i; j++){
               swap(L[i][j], L[i_max][j]);
           }
       }

       for(int j = i+1; j<n; j++){
           L[j][i] += U[j][i]/U[i][i];
           U[j][i] = 0;
           for(int k = i+1; k<n; k++){
               U[j][k] -= U[i][k]*L[j][i];
           }
       }
   }
   return true;
}
/*
 * Aby roziwazac rownanie:
 * A * x = b
 * Korzystamy z tego ze:
 * A = L * U
 * Oraz tworzymy uklad rowanan:
 * L * y = b
 * y = U * x
 */

/**
 * Wyznaczenie y potrzebnego do rozwiazania powyzszego ukladu rownan
 * @param macierzL maczierz L (trojkatna dolna)
 * @param wektorB wektor b
 * @param index tablica indexow
 * @param n rozmiar macierzy
 */
void wyznaczY() {
  double suma = 0.0;
  for(int i = 0; i< n; i++){
      swap(B[i], B[swapping[i]]);
  }

  // Poruszanie sie po kolumnie, zaczyna od lewego gornego
  for (int i = 0; i < n; i++) {
    // Poruszanie sie wierszu w prawo
    for (int j = 0; j < i; j++)
      suma += L[i][j] * B[j];

    B[i] = (B[i] - suma) / 1.0;     //Na glownej przekatnej sa 1, bo L jest polaczone z U

    suma = 0.0;
  }
}

/**
 * Wyznaczenie x potrzebnego do rozwiazania powyzszego ukladu rownan
 * @param macierzU maczierz U (trojkatna gorna)
 * @param wektorB wektor b
 * @param index tablica indexow
 * @param n rozmiar macierzy
 */
void wyznaczX() {
  double suma = 0.0;
  // Zaczyna od prawego dolnego rogu
  for (int i = n; i >= 0; i--) {
    for (int j = i + 1; j <= n; j++)
      suma += U[i][j] * B[j];

//#ifdef opcja2
    if(U[i][i] == 0)
      std::cout << "DZIELENIE PRZEZ 0\n";
//#endif

    B[i] = (B[i] - suma) / (U[i][i]);     //UWAGA MOŻLIWE DZIELENIE PRZEZ 0

    suma = 0.0;
  }
}

int main()
{
    double A[n][n] ={{1,-20,30,-4},
                     {2,-40,-6,50},
                     {9,-180,11,-12},
                     {-16,15,-140,13}};

    for(int i=0; i<n; i++){
        for(int j=0; j<n; j++){
            U[i][j]=A[i][j];
        }
    }
    for(int i=0; i<n; i++){
        swapping[i]=i;
        for(int j=0; j<n; j++){
            L[i][j]=0;
            if(i==j){
                L[i][j]=1;
            }
        }
    }
    Aeliminacja();
    cout<<"A: "<<endl;
    for(int i=0; i<n; i++){
        for(int j=0; j<n; j++){
            cout<<setw(10)<<A[i][j]<<" ";
        }
        cout<<endl;
    }
    cout<<"L: "<<endl;
    for(int i=0; i<n; i++){
        for(int j=0; j<n; j++){
            cout<<setw(10)<<L[i][j]<<" ";
        }
        cout<<endl;
    }
    cout<<"U: "<<endl;
    for(int i=0; i<n; i++){
        for(int j=0; j<n; j++){
            cout<<setw(10)<<U[i][j]<<" ";
        }
        cout<<endl;
    }
    wyznaczY();
    cout<<"B: "<<endl;
    for(int j=0; j<n; j++){
        cout<<setw(10)<<B[j]<<endl;
    }
    wyznaczX();
    cout<<"B: "<<endl;
    for(int j=0; j<n; j++){
        cout<<setw(10)<<B[j]<<endl;
    }
    cout<<"swap: "<<endl;
    for(int j=0; j<n; j++){
        cout<<setw(10)<<swapping[j]<<endl;
    }
    return 0;
}
