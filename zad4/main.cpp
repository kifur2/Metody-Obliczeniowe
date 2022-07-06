#include <iostream>
#include <math.h>
#include <vector>
#include <iomanip>
#include <fstream>
#define NMAX 25
#define BLAD 1e-60
using namespace std;
struct dane
{
    long double x,y,z,f_1,f_2,f_3, e_1, e_2, e_3;
};
long double f1(long double x, long double y, long double z){
    return x*x+y*y+z*z-2;
}
long double f1dx(long double x, long double y, long double z){
    return 2*x;
}
long double f1dy(long double x, long double y, long double z){
    return 2*y;
}
long double f1dz(long double x, long double y, long double z){
    return 2*z;
}
long double f2(long double x, long double y, long double z){
    return x*x+y*y-1;
}
long double f2dx(long double x, long double y, long double z){
    return 2*x;
}
long double f2dy(long double x, long double y, long double z){
    return 2*y;
}
long double f2dz(long double x, long double y, long double z){
    return 0;
}
long double f3(long double x, long double y, long double z){
    return x*x-y;
}
long double f3dx(long double x, long double y, long double z){
    return 2*x;
}
long double f3dy(long double x, long double y, long double z){
    return -1;
}
long double f3dz(long double x, long double y, long double z){
    return 0;
}
long double Jacob[3][3];
long double macierz_funkcji[3];
void eliminacja(){
    long double iloraz;
    if(Jacob[0][0]!=0)
    {
        iloraz = Jacob[1][0]/Jacob[0][0];
        Jacob[1][0]=0;
        Jacob[1][1]-=Jacob[0][1]*iloraz;
        Jacob[1][2]-=Jacob[0][2]*iloraz;
        macierz_funkcji[1]-=macierz_funkcji[0]*iloraz;
        iloraz = Jacob[2][0]/Jacob[0][0];
        Jacob[2][0]=0;
        Jacob[2][1]-=Jacob[0][1]*iloraz;
        Jacob[2][2]-=Jacob[0][2]*iloraz;
        macierz_funkcji[2]-=macierz_funkcji[0]*iloraz;
    }
    else
    {
        if(Jacob[1][0]!=0)
        {
            swap(Jacob[0][0],Jacob[1][0]);
            swap(Jacob[0][1],Jacob[1][1]);
            swap(Jacob[0][2],Jacob[1][2]);
            swap(macierz_funkcji[0],macierz_funkcji[1]);
            iloraz = Jacob[2][0]/Jacob[0][0];
            Jacob[2][0]=0;
            Jacob[2][1]-=Jacob[0][1]*iloraz;
            Jacob[2][2]-=Jacob[0][2]*iloraz;
            macierz_funkcji[2]-=macierz_funkcji[0]*iloraz;
        }
        else
        {
            swap(Jacob[0][0],Jacob[2][0]);
            swap(Jacob[0][1],Jacob[2][1]);
            swap(Jacob[0][2],Jacob[2][2]);
            swap(macierz_funkcji[0],macierz_funkcji[2]);
        }
    }
    if(Jacob[1][1]!=0)
    {
        iloraz = Jacob[2][1]/Jacob[1][1];
        Jacob[2][1]=0;
        Jacob[2][2]-=Jacob[1][2]*iloraz;
        macierz_funkcji[2]-=macierz_funkcji[1]*iloraz;
    }
    else
    {
        swap(Jacob[1][0],Jacob[2][0]);
        swap(Jacob[1][1],Jacob[2][1]);
        swap(Jacob[1][2],Jacob[2][2]);
        swap(macierz_funkcji[1],macierz_funkcji[2]);
    }
}
vector <dane > newton_vector;
string e_newton;
dane xn;
void newton_fun()
{
    for(int i=0; i<NMAX; i++)
    {
        long double x = newton_vector[i].x;
        long double y = newton_vector[i].y;
        long double z = newton_vector[i].z;
        macierz_funkcji[0]=newton_vector[i].f_1;
        macierz_funkcji[1]=newton_vector[i].f_2;
        macierz_funkcji[2]=newton_vector[i].f_3;
        Jacob[0][0]= f1dx(x,y,z);
        Jacob[0][1]= f1dy(x,y,z);
        Jacob[0][2]= f1dz(x,y,z);
        Jacob[1][0]= f2dx(x,y,z);
        Jacob[1][1]= f2dy(x,y,z);
        Jacob[1][2]= f2dz(x,y,z);
        Jacob[2][0]= f3dx(x,y,z);
        Jacob[2][1]= f3dy(x,y,z);
        Jacob[2][2]= f3dz(x,y,z);
        eliminacja();
        long double x_3 = macierz_funkcji[2]/Jacob[2][2];
        long double x_2 =(macierz_funkcji[1]-Jacob[1][2]*x_3)/Jacob[1][1];
        long double x_1 =(macierz_funkcji[1]-Jacob[0][2]*x_3-Jacob[0][1]*x_2)/Jacob[0][0];
        xn.x=x-x_1;
        xn.y=y-x_2;
        xn.z=z-x_3;
        xn.f_1=f1(xn.x, xn.y,xn.z);
        xn.f_2=f2(xn.x, xn.y,xn.z);
        xn.f_3=f3(xn.x, xn.y,xn.z);
        xn.e_1=x_1;
        xn.e_2=x_2;
        xn.e_3=x_3;
        if(fabs(xn.e_1)<BLAD && fabs(xn.e_2)<BLAD && fabs(xn.e_3)<BLAD)
        {
            e_newton="TOLX";
            break;
        }
        if(fabs(xn.f_1)<BLAD && fabs(xn.f_2)<BLAD && fabs(xn.f_3)<BLAD )
        {
            e_newton="TOLF";
            break;
        }
        newton_vector.push_back(xn);
    }
}

void set_all()
{
    newton_vector.clear();
    xn.x=-0.8;
    xn.y=0.6;
    xn.z=1.0;
    xn.f_1=f1(xn.x, xn.y,xn.z);
    xn.f_2=f2(xn.x, xn.y,xn.z);
    xn.f_3=f3(xn.x, xn.y,xn.z);
    newton_vector.push_back(xn);
}

int main()
{
    set_all();
    newton_fun();
    ofstream plik("../newton");
    plik<<setw(15)<<"x"<<setw(15)<<"estym"<<setw(15)<<"reziduum";
    plik<<setw(15)<<"y"<<setw(15)<<"estym"<<setw(15)<<"reziduum";
    plik<<setw(15)<<"z"<<setw(15)<<"estym"<<setw(15)<<"reziduum"<<endl;
    for(int i=0; i<NMAX; i++)
    {
        if(i<newton_vector.size())
        {
            plik<<setw(15)<<newton_vector[i].x<<setw(15)<<newton_vector[i].e_1<<setw(15)<<newton_vector[i].f_1;
            plik<<setw(15)<<newton_vector[i].y<<setw(15)<<newton_vector[i].e_2<<setw(15)<<newton_vector[i].f_2;
            plik<<setw(15)<<newton_vector[i].z<<setw(15)<<newton_vector[i].e_3<<setw(15)<<newton_vector[i].f_3;
        }
        else
            plik<<setw(15)<<e_newton<<setw(15)<<e_newton<<setw(15)<<e_newton<<setw(15)<<e_newton<<setw(15)<<e_newton<<setw(15)<<e_newton<<setw(15)<<e_newton<<setw(15)<<e_newton<<setw(15)<<e_newton;
        plik<<endl;
    }
}
