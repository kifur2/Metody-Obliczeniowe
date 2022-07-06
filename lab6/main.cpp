#include <iostream>
#include <cmath>
#include <iomanip>

using namespace std;

const int n = 6;
double U[n]={0,0.5,0.25,1./6.,0.125,0.1};
double DW[n+1]={0, 10,20,30,30,20,10};
double L[n+1]={0,0, 1./3.,0.2,1./7.,1./9.,1./11.};
double BR[n+1]={0,31,165./4.,917./30.,851./28.,3637./90.,332./11.};
double X[n+1];


bool eliminacja1(){
    double eps = 10e-12;
    for (int i=2; i<=n; i++)
    {
        if(fabs(DW[i])<eps)
            return false;
        DW[i]=DW[i]-L[i]/DW[i-1]*U[i-1];
    }
    return true;
}
bool eliminacja2()
{
    double eps = 10e-12;
    for (int i=2; i<=n; i++)
    {
        if(fabs(DW[i-1])<eps)
            return false;
        BR[i]=BR[i]-L[i]/DW[i-1]*BR[i-1];
    }

    X[n]=1.0/DW[n]*BR[n];
    for(int i = n-1; i>0; i--)
    {
        X[i]=1.0/DW[i]*(BR[i]-U[i]*X[i+1]);
    }

    return true;
}

void wypiszA(){
    cout<<"Macierz A: "<<endl;
    cout<<setw(10)<<DW[1]<<setw(10)<<U[1];
    for(int i=0; i<n-2; i++){
        cout<<setw(10)<<0;
    }
    cout<<endl;
    for(int i=2; i<n; i++){
        for(int j=2;j<i;j++){
            cout<<setw(10)<<0;
        }
        cout<<setw(10)<<L[i]<<setw(10)<<DW[i]<<setw(10)<<U[i];
        for(int j=i+2; j<n; j++){
            cout<<setw(10)<<0;
        }
        cout<<endl;
    }
    for(int i=0; i<n-2; i++){
        cout<<setw(10)<<0;
    }
    cout<<setw(10)<<L[n]<<setw(10)<<DW[n]<<endl;
}

void wypiszW(){
    cout<<"Macierz W: "<<endl;
    cout<<setw(10)<<DW[1]<<setw(10)<<U[1];
    for(int i=0; i<n-2; i++){
        cout<<setw(10)<<0;
    }
    cout<<endl;
    for(int i=2; i<n; i++){
        for(int j=1;j<i;j++){
            cout<<setw(10)<<0;
        }
        cout<<setw(10)<<DW[i]<<setw(10)<<U[i];
        for(int j=i+2; j<n; j++){
            cout<<setw(10)<<0;
        }
        cout<<endl;
    }
    for(int i=0; i<n-1; i++){
        cout<<setw(10)<<0;
    }
    cout<<setw(10)<<DW[n]<<endl;
}

void wypiszB(const string& wiadomosc){
    cout<<wiadomosc<<endl;
    for(int i=1; i<=n; i++){
        cout<<setw(10)<<BR[i]<<endl;
    }
}
void wypiszX(){
    cout<<"Macierz X: "<<endl;
    for(int i=1; i<=n; i++){
        cout<<setw(10)<<X[i]<<endl;
    }
}

int main()
{
    wypiszA();
    wypiszB("Macierz B: ");
    if(!eliminacja1()){
        cout<<"wystapil blad w 1 ----> exit";
        exit(1);
    }
    wypiszW();
    if(!eliminacja2()){
        cout<<"wystapil blad w 2 ----> exit";
        exit(1);
    }
    wypiszB("Macierz R: ");
    wypiszX();
    return 0;
}