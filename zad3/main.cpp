#include <iostream>
#include <math.h>
#include <vector>
#include <iomanip>
#include <fstream>
#define NMAX 50
#define BLAD 1e-16
using namespace std;
vector <pair<long double,pair<long double, long double> > > sin_picard;
vector <pair<long double,pair<long double, long double> > > tan_picard;
vector <pair<long double,pair<long double, long double> > > sin_bisekcja;
vector <pair<long double,pair<long double, long double> > > tan_bisekcja;
vector <pair<long double,pair<long double, long double> > > sin_newton;
vector <pair<long double,pair<long double, long double> > > tan_newton;
vector <pair<long double,pair<long double, long double> > > sin_sieczne;
vector <pair<long double,pair<long double, long double> > > tan_sieczne;

pair<long double,pair<long double, long double> > xn;
string e_sin_picard;
string e_tan_picard;
string e_sin_bisekcja;
string e_tan_bisekcja;
string e_sin_newton;
string e_tan_newton;
string e_sin_sieczne;
string e_tan_sieczne;

int i1=0, i2=0;
void clear_all()
{
    sin_picard.clear();
    tan_picard.clear();
    sin_bisekcja.clear();
    tan_bisekcja.clear();
    sin_newton.clear();
    tan_newton.clear();
    sin_sieczne.clear();
    tan_sieczne.clear();
    pair<long double,pair<long double, long double> > x;
    x.first=0.5;
    sin_picard.push_back(x);
    tan_picard.push_back(x);
    sin_newton.push_back(x);
    tan_newton.push_back(x);
    sin_sieczne.push_back(x);
    tan_sieczne.push_back(x);
    x.first=1.5;
    x.second.first=1;
    sin_sieczne.push_back(x);
    tan_sieczne.push_back(x);
}
void picard_sin()
{
    for(int i=0; i<NMAX; i++)
    {
        long double x = sin_picard[i].first;
        xn.first = sin(x / 4) * sin(x / 4);
        xn.second.first = xn.first - x;
        xn.second.second = sin(xn.first / 4) * sin(xn.first / 4);
        sin_picard.push_back(xn);
        if(fabs(xn.second.first)<BLAD)
        {
            e_sin_picard="TOLX";
            break;
        }
        if(fabs(xn.second.second)<BLAD)
        {
            e_sin_picard="TOLF";
            break;
        }
    }
}
void picard_tan()
{
    for(int i=0; i<NMAX; i++)
    {
        long double x = tan_picard[i].first;
        xn.first = tan(2 * x) - 1;
        xn.second.first = xn.first - x;
        xn.second.second = tan(2 * xn.first) - 1;
        tan_picard.push_back(xn);
        if(fabs(xn.second.first)<BLAD)
        {
            e_tan_picard="TOLX";
            break;
        }
        if(fabs(xn.second.second)<BLAD)
        {
            e_tan_picard="TOLF";
            break;
        }
    }
}

void bisekcja_sin(long double a, long double b)
{
    i1++;
    if(i1<=NMAX)
    {
        //cout<<i1<<endl;
        xn.first = (a + b) / 2;
        xn.second.first=(b-a)/2;
        xn.second.second = sin(xn.first / 4) * sin(xn.first / 4)- xn.first;
        sin_bisekcja.push_back(xn);
        long double wynika = sin(a / 4) * sin(a / 4)-a;
        long double wynikxn = sin(xn.first / 4) * sin(xn.first / 4)-xn.first;
        long double wynikb = sin(b / 4) * sin(b / 4)-b;
        if(fabs(xn.second.first)<BLAD)
        {
            e_sin_bisekcja="TOLX";
            i1=NMAX;
            cout<<"x";
        }
        if(fabs(xn.second.second)<BLAD)
        {
            e_sin_bisekcja="TOLF";
            i1=NMAX;
            cout<<"d";
        }
        /*cout<<wynika<<endl;
        cout<<wynikb<<endl;
        cout<<wynikxn<<endl;
        cout<<wynika*wynikxn<<endl;
        cout<<(wynika*wynikxn<0) <<endl;*/
        if (wynika * wynikxn < 0)
            bisekcja_sin(a, xn.first);
        else if (wynikb * wynikxn < 0)
            bisekcja_sin(xn.first, b);
    }
}
void bisekcja_tan(long double a, long double b)
{
    i2++;
    if(i2<=NMAX)
    {
        xn.first = (a + b) / 2;
        xn.second.first=(b-a)/2;
        xn.second.second = tan(2 * xn.first)- xn.first - 1;
        tan_bisekcja.push_back(xn);
        long double wynika = tan(2*a)-a-1;
        long double wynikxn = tan(2 * xn.first)-xn.first - 1;
        long double wynikb = tan(2 * b)-b - 1;

        if(fabs(xn.second.first)<BLAD)
        {
            e_tan_bisekcja="TOLX";
            i2=NMAX;
        }
        if(fabs(xn.second.second)<BLAD)
        {
            e_tan_bisekcja="TOLF";
            i2=NMAX;
        }

        if (wynika * wynikxn < 0)
            bisekcja_tan(a, xn.first);
        else if (wynikb * wynikxn < 0)
            bisekcja_tan(xn.first, b);
    }
}
void newton_sin()
{
    for(int i=0; i<NMAX; i++)
    {
        long double x = sin_newton[i].first;
        xn.first=x - (sin(x / 4) * sin(x / 4) - x) / (sin(x / 2) / 4 - 1);
        xn.second.first = xn.first - x;
        xn.second.second = sin(xn.first / 4) * sin(xn.first / 4);
        sin_newton.push_back(xn);
        if(fabs(xn.second.first)<BLAD)
        {
            e_sin_newton="TOLX";
            break;
        }
        if(fabs(xn.second.second)<BLAD)
        {
            e_sin_newton="TOLF";
            break;
        }
    }
}
void newton_tan()
{
    for(int i=0; i<NMAX; i++)
    {
        long double x = tan_newton[i].first;
        xn.first = x - (tan(2 * x) - x - 1) / (2 / (cos(2 * x) * cos(2 * x)) - 1);
        xn.second.first = xn.first - x;
        xn.second.second = tan(2 * xn.first) - 1;
        tan_newton.push_back(xn);
        if(fabs(xn.second.first)<BLAD)
        {
            e_tan_newton="TOLX";
            break;
        }
        if(fabs(xn.second.second)<BLAD)
        {
            e_tan_newton="TOLF";
            break;
        }
    }
}
void sieczne_sin()
{
    for(int i=1; i<NMAX; i++)
    {
        long double x = sin_sieczne[i-1].first;
        long double x1 = sin_sieczne[i].first;
        xn.first = x1 - (sin(x1 / 4) * sin(x1 / 4) - x1) /
                        (((sin(x1 / 4) * sin(x1 / 4) - x1) - (sin(x / 4) * sin(x / 4) - x)) / (x1 - x));
        xn.second.first = xn.first - x1;
        xn.second.second = sin(xn.first / 4) * sin(xn.first / 4);
        sin_sieczne.push_back(xn);
        if(fabs(xn.second.first)<BLAD)
        {
            e_sin_sieczne="TOLX";
            break;
        }
        if(fabs(xn.second.second)<BLAD)
        {
            e_sin_sieczne="TOLF";
            break;
        }
    }
}
void sieczne_tan()
{
    for(int i=1; i<NMAX; i++)
    {
        long double x = tan_sieczne[i-1].first;
        long double x1 = tan_sieczne[i].first;
        xn.first = x1-(tan(2*x1)-x1-1)/
                      (((tan(2*x1)-x1-1)-((tan(2*x)-x-1)))/(x1-x));
        xn.second.first = xn.first - x1;
        xn.second.second = tan(2 * xn.first) - 1;
        tan_sieczne.push_back(xn);
        if(fabs(xn.second.first)<BLAD)
        {
            e_tan_sieczne="TOLX";
            break;
        }
        if(fabs(xn.second.second)<BLAD)
        {
            e_tan_sieczne="TOLF";
            break;
        }
    }
}
int main()
{
    clear_all();
    picard_sin();
    picard_tan();
    bisekcja_sin(-0.5,10);
    bisekcja_tan(-0.5,5);
    newton_sin();
    newton_tan();
    sieczne_sin();
    sieczne_tan();

    ofstream plik("../sinus");
    if(!plik.is_open())
        exit(1);
    plik<<setw(15)<<"picard"<<setw(15)<<"estym"<<setw(15)<<"reziduum";
    plik<<setw(15)<<"bisekcja"<<setw(15)<<"estym"<<setw(15)<<"reziduum";
    plik<<setw(15)<<"newton"<<setw(15)<<"estym"<<setw(15)<<"reziduum";
    plik<<setw(15)<<"sieczne"<<setw(15)<<"estym"<<setw(15)<<"reziduum"<<endl;
    for(int i=0; i<NMAX; i++)
    {
        if(i<sin_picard.size())
            plik<<setw(15)<<sin_picard[i].first<<setw(15)<<sin_picard[i].second.first<<setw(15)<<sin_picard[i].second.second;
        else
            plik<<setw(15)<<e_sin_picard<<setw(15)<<e_sin_picard<<setw(15)<<e_sin_picard;
        if(i<sin_bisekcja.size())
            plik<<setw(15)<<sin_bisekcja[i].first<<setw(15)<<sin_bisekcja[i].second.first<<setw(15)<<sin_bisekcja[i].second.second;
        else
            plik<<setw(15)<<e_sin_bisekcja<<setw(15)<<e_sin_bisekcja<<setw(15)<<e_sin_bisekcja;
        if(i<sin_newton.size())
            plik<<setw(15)<<sin_newton[i].first<<setw(15)<<sin_newton[i].second.first<<setw(15)<<sin_newton[i].second.second;
        else
            plik<<setw(15)<<e_sin_newton<<setw(15)<<e_sin_newton<<setw(15)<<e_sin_newton;
        if(i<sin_sieczne.size())
            plik<<setw(15)<<sin_sieczne[i].first<<setw(15)<<sin_sieczne[i].second.first<<setw(15)<<sin_sieczne[i].second.second;
        else
            plik<<setw(15)<<e_sin_sieczne<<setw(15)<<e_sin_sieczne<<setw(15)<<e_sin_sieczne;
        plik<<endl;
    }
    plik.close();
    ofstream plik1("../tangens");
    if(!plik1.is_open())
        exit(1);
    plik1<<setw(15)<<"picard"<<setw(15)<<"estym"<<setw(15)<<"reziduum";
    plik1<<setw(15)<<"bisekcja"<<setw(15)<<"estym"<<setw(15)<<"reziduum";
    plik1<<setw(15)<<"newton"<<setw(15)<<"estym"<<setw(15)<<"reziduum";
    plik1<<setw(15)<<"sieczne"<<setw(15)<<"estym"<<setw(15)<<"reziduum"<<endl;
    for(int i=0; i<NMAX; i++)
    {
        if(i<tan_picard.size())
            plik1<<setw(15)<<tan_picard[i].first<<setw(15)<<tan_picard[i].second.first<<setw(15)<<tan_picard[i].second.second;
        else
            plik1<<setw(15)<<e_tan_picard<<setw(15)<<e_tan_picard<<setw(15)<<e_tan_picard;
        if(i<tan_bisekcja.size())
            plik1<<setw(15)<<tan_bisekcja[i].first<<setw(15)<<tan_bisekcja[i].second.first<<setw(15)<<tan_bisekcja[i].second.second;
        else
            plik1<<setw(15)<<e_tan_bisekcja<<setw(15)<<e_tan_bisekcja<<setw(15)<<e_tan_bisekcja;
        if(i<tan_newton.size())
            plik1<<setw(15)<<tan_newton[i].first<<setw(15)<<tan_newton[i].second.first<<setw(15)<<tan_newton[i].second.second;
        else
            plik1<<setw(15)<<e_tan_newton<<setw(15)<<e_tan_newton<<setw(15)<<e_tan_newton;
        if(i<tan_sieczne.size())
            plik1<<setw(15)<<tan_sieczne[i].first<<setw(15)<<tan_sieczne[i].second.first<<setw(15)<<tan_sieczne[i].second.second;
        else
            plik1<<setw(15)<<e_tan_sieczne<<setw(15)<<e_tan_sieczne<<setw(15)<<e_tan_sieczne;
        plik1<<endl;
    }
    plik1.close();
}
