#include <iostream>
#include <fstream>
#include <cmath>
#include <vector>
#include <string>
#include <valarray>

using namespace std;
int mySign(double v) {
    if (v <= 0) return -1;
    if (v > 0) return 1;
    return 0;
}

double fun1( double y)
{
    return y;
}

double fun2(double x, double y,double L ,double Cd2m, double k2m,double g2m){
    if (x<L){
        return 9.81 - mySign(y)* (Cd2m) * y*y ;
    }else return 9.81 - mySign(y) * (Cd2m) * y*y - ((k2m) * (x - L) + (g2m) * y);
}

int main(int argc, char** argv) {
    double L =30 ;
    double Lcut=50;
    // Read design variables
    std::ifstream in("task.dat", std::ios::in);
    int n;
    in >> n;
    std::valarray<double> b;
    b.resize(n);
    b = 0;
    for (int i = 0; i < n; i++) in >> b[i];
    in.close();
    //
    double k2m = b[0];
    double g2m =b[1];
    double Cd2m = b[2];
    double t0 =  b[3];


    vector<double> y(2), aux(2);
    vector<vector<double> > ak(2, vector<double>(4));

    double time = 0;
    double deltat = 0.05;
    int neqs = 2;
    y[0] = 0;
    y[1] = 0;


    int tl=0 ;
    for (int ktimestep = 1; ktimestep < 2000; ktimestep++) {
        //Step 1:
        ak[0][0] = deltat * fun1(y[1]);
        ak[1][0] = deltat * fun2(y[0], y[1] ,L , Cd2m,  k2m, g2m);
        aux[0] = y[0] + 0.5e0 * ak[0][0];
        aux[1] = y[1] + 0.5e0 * ak[1][0];
        //Step 2:
        ak[0][1] = deltat * fun1(aux[1]);
        ak[1][1] = deltat * fun2(aux[0], aux[1],L , Cd2m,  k2m, g2m);
        aux[0] = y[0] + 0.5e0 * ak[0][1];
        aux[1] = y[1] + 0.5e0 * ak[1][1];
        //Step 3:
        ak[0][2] = deltat * fun1(aux[1]);
        ak[1][2] = deltat * fun2(aux[0], aux[1],L , Cd2m,  k2m, g2m);
        aux[0] = y[0] + ak[0][2];
        aux[1] = y[1] + ak[1][2];
        //Step 4:
        ak[0][3] = deltat * fun1(aux[1]);
        ak[1][3] = deltat * fun2(aux[0], aux[1],L , Cd2m,  k2m, g2m);
        //Synthesis:
        y[0] = y[0] + (ak[0][0] + 2.e0 * ak[0][1] + 2.e0 * ak[0][2] + ak[0][3]) / 6.e0;
        y[1] = y[1] + (ak[1][0] + 2.e0 * ak[1][1] + 2.e0 * ak[1][2] + ak[1][3]) / 6.e0;


        if (sin(0.6283 * ((time+t0))) <= 0 && y[0] >Lcut) {
            tl++;
        } else if (sin(0.628 * (time+t0)) > 0 && y[0] >= Lcut)

            break;


        cout << time+t0 << "\t";
        for (int i = 0; i < neqs; i++)
            cout << y[i]<<"    " << "\t" ;
           // cout << y[i]<<"    " << "\t" <<"    "<<50*abs(sin(0.6283 * (time+t0)))/(sin(0.6283 * (time+t0)));
        cout << endl;
        time = time + deltat;

    }
    double f1=-tl*deltat;
    double f2=-Cd2m;
    cout << "time under Lcut is " << f1 << endl;

    // Return
    std::ofstream out("task.res", std::ios::out);
    out.precision(14);
    out << std::scientific << f1 << std::endl;
    out<<std::scientific<<f2<<std::endl;
    out.close();
    //

    return 0;
}


//// με comment out τις γραμμές 96 και 103 είναι ο κώδικας του soo