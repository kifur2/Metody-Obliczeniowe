#include <iostream>
#include <vector>
#include <fstream>
#include <cmath>
using namespace std;

vector<long double> x_dane;
vector<long double> fx_dane;
vector<long double> log_dane;
vector<long double> log_fx_wynik;
vector<long double> log_fx_taylora_wynik;
unsigned long long liczba_danych=0;

void save(vector<long double> *wynik,const string &nazwa_pliku);

long double taylor(long double x, int n);

int main()
{
    ifstream plik("../dane.txt");
    string str;
    if (plik.is_open())
    {
        for(int i=0; i<3900;i++)
        {
            plik>>str;
            (log_dane).push_back(stod(str));
            plik>>str;
            (x_dane).push_back(stod(str));
            plik>>str;
            (fx_dane).push_back(stod(str));
        }
        plik.close();
    }

    liczba_danych= x_dane.size();
    for (unsigned long long i = 0; i < liczba_danych; ++i)
    {
        long double fun = (1.0 - exp(-x_dane[i])) / x_dane[i];
        log_fx_wynik.push_back(log10(fabs((fx_dane[i] - ((1.0 - exp(-x_dane[i])) / x_dane[i])) / fx_dane[i])));
    }

    save(&log_fx_wynik,"../wykres1.txt");

    for (unsigned long long i = 0; i < liczba_danych; ++i)
    {
        log_fx_taylora_wynik.push_back(log10(fabs((fx_dane[i] - taylor(fx_dane.at(i), 99999)) / fx_dane[i])));
    }

    save(&log_fx_taylora_wynik, "../wykres2.txt");
}

void save(vector<long double> *wynik,const string &nazwa_pliku)
{
    ofstream plik(nazwa_pliku);
    if (plik.is_open())
    {
        for (unsigned long long i = 0; i < liczba_danych; ++i)
        {
            plik << wynik->at(i) << ' '<< log_dane[i] << endl;
        }
        plik.close();
    }
}

long double taylor(long double x, int n)
{
    long double suma = 0, iloczyn = 1, silnia = 1;
    for (int i = 1; i <= n; ++i)
    {
        silnia *= i;
        if (i % 2 == 0)
        {
            suma += iloczyn / silnia;
        }
        else
        {
            suma -= iloczyn / silnia;
        }
        iloczyn *= x;
    }
    return suma;

}