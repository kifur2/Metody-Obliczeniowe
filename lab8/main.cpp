#include <iostream>
#include <cmath>
#include <fstream>

using namespace std;
template <class A> A wewnetrznaProgresywna(A x, A h){
    return fabs(cos(x) - (sin(x+h)-sin(x))/h);
}
template <class A> A wewnetrznaWsteczna(A x, A h){
    return fabs(cos(x) - (sin(x)-sin(x-h))/h);
}
template <class A> A wewnetrznaCentralna(A x, A h){
    return fabs(cos(x) - (sin(x+h)-sin(x-h))/2/h);
}
template <class A> A dwupunktowaPoczatek(A x, A h){
    return fabs(cos(x) - (sin(x+h)-sin(x))/h);
}
template <class A> A trzypunktowaPoczatek(A x, A h){
    return fabs(cos(x) - (-(A)1.5 * sin(x) + 2 * sin(x+h)- (A)0.5*sin(x+h+h))/h);
}
template <class A> A dwupunktowaKoniec(A x, A h){
    return fabs(cos(x) - (sin(x)-sin(x-h))/h);
}
template <class A> A trzypunktowaKoniec(A x, A h){
    return fabs(cos(x) - ((A)0.5 * sin(x-2*h) - 2 * sin(x-h) + (A)1.5*sin(x))/h);
}

int main() {
    ofstream plik("../wewnetrznaProgresywnaDouble.txt");
    ofstream plik1("../wewnetrznaWstecznaDouble.txt");
    ofstream plik2("../wewnetrznaCentralnaDouble.txt");
    ofstream plik3("../dwupunktowaPoczatekDouble.txt");
    ofstream plik4("../trzypunktowaPoczatekDouble.txt");
    ofstream plik5("../dwupunktowaKoniecDouble.txt");
    ofstream plik6("../trzypunktowaKoniecDouble.txt");
    for (double h = (double) (M_PI / 400000); h < M_PI / 4; h += (double) (M_PI / 400000)) {
        plik << log10(h) <<' '<< log10(wewnetrznaProgresywna<double>(M_PI / 4, h)) << endl;
        plik1 << log10(h) <<' '<< log10(wewnetrznaWsteczna<double>(M_PI / 4, h)) << endl;
        plik2 << log10(h) <<' '<< log10(wewnetrznaCentralna<double>(M_PI / 4, h)) << endl;
        plik3 << log10(h) <<' '<< log10(dwupunktowaPoczatek<double>(0, h)) << endl;
        plik4 << log10(h) <<' '<< log10(trzypunktowaPoczatek<double>(0, h)) << endl;
        plik5 << log10(h) <<' '<< log10(dwupunktowaKoniec<double>(M_PI / 2, h)) << endl;
        plik6 << log10(h) <<' '<< log10(trzypunktowaKoniec<double>(M_PI / 2, h)) << endl;
        cout<<"h: "<<h<<endl;
    }
    plik.close();
    plik1.close();
    plik2.close();
    plik3.close();
    plik4.close();
    plik5.close();
    plik6.close();
    ofstream plik7("../wewnetrznaProgresywnaFloat.txt");
    ofstream plik8("../wewnetrznaWstecznaFloat.txt");
    ofstream plik9("../wewnetrznaCentralnaFloat.txt");
    ofstream plik10("../dwupunktowaPoczatekFloat.txt");
    ofstream plik11("../trzypunktowaPoczatekFloat.txt");
    ofstream plik12("../dwupunktowaKoniecFloat.txt");
    ofstream plik13("../trzypunktowaKoniecFloat.txt");
    for (float h = (float) (M_PI / 400000); h < M_PI / 4; h+= (float) (M_PI / 400000)) {
        plik7 << log10(h) <<' '<< log10(wewnetrznaProgresywna<float>(M_PI / 4, h)) << endl;
        plik8 << log10(h) <<' '<< log10(wewnetrznaWsteczna<float>(M_PI / 4, h)) << endl;
        plik9 << log10(h) <<' '<< log10(wewnetrznaCentralna<float>(M_PI / 4, h)) << endl;
        plik10 << log10(h) <<' '<< log10(dwupunktowaPoczatek<float>(0, h)) << endl;
        plik11 << log10(h) <<' '<< log10(trzypunktowaPoczatek<float>(0, h)) << endl;
        plik12 << log10(h) <<' '<< log10(dwupunktowaKoniec<float>(M_PI / 2, h)) << endl;
        plik13 << log10(h) <<' '<< log10(trzypunktowaKoniec<float>(M_PI / 2, h)) << endl;
        cout<<"h: "<<h<<endl;
    }
    plik7.close();
    plik8.close();
    plik9.close();
    plik10.close();
    plik11.close();
    plik12.close();
    plik13.close();
}
