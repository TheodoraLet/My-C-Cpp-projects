#include <iostream>
#include <vector>
#include <cmath>
#include <stdlib.h>
#include <valarray>
#include <cstdlib>
#include<fstream>

using namespace std;

void UELEF0 (double EPSF,double gamma,vector<double> &XA,vector<double> &XB,vector<double> &XO,vector<double> &U ) {
    const double PI4 = 4 * 3.14159365;
    vector<double> RA(3),RB(3),AB(3), AL(3);
    double TEST, HE, FUN, AR, DAL, DRA, DRB, DAB, H, C;;
    for (int K = 0; K < 3; K++) {
        RA[K] = XO[K] - XA[K];
        RB[K] = XO[K] - XB[K];
        AB[K] = XB[K] - XA[K];
    }
    AL[0] = RA[1] * RB[2] - RA[2] * RB[1]; // πρωτη συνιστωστα εξ γινομενου
    AL[1] = RA[2] * RB[0] - RA[0] * RB[2]; /// δευτερη συνιστωσα
    AL[2] = RA[0] * RB[1] - RA[1] * RB[0]; /// τριτη συνιστωσα
    DRA = sqrt(RA[0] * RA[0] + RA[1] * RA[1] + RA[2] * RA[2]); /// μετρο r1
    DRB = sqrt(RB[0] * RB[0] + RB[1] * RB[1] + RB[2] * RB[2]); /// μετρο r2
    DAB = sqrt(AB[0] * AB[0] + AB[1] * AB[1] + AB[2] * AB[2]); /// μετρο r0
    DAL = AL[0] * AL[0] + AL[1] * AL[1] + AL[2] * AL[2];
    TEST = DAL * DRA * DRB;
    AR = sqrt(DAL); /// μετρο εξωτερικου γινομενου |r1xr2|

   C = ((RA[0] * AB[0] + RA[1] * AB[1] + RA[2] * AB[2]) / DRA -
         (RB[0] * AB[0] + RB[1] * AB[1] + RB[2] * AB[2]) / DRB) / (PI4 * AR * DAB);
    H = AR / DAB; /// |r1xr2|/|r0|
        if (H < 3 * EPSF) {
        U[0] = 0.0;
        U[1] = 0.0;
        U[2] = 0.0;
    } else {
        C = C / H;
        U[0] = AL[0] * C * gamma + U[0];
        U[1] = AL[1] * C * gamma + U[1];
        U[2] = AL[2] * C * gamma + U[2];
    }
}

int main() {
    ofstream outdata;
    int n=100; // αριθμός δινοπετάλων και κέντρων
    int m=3*n; // αριθμός δινοσωλήνων
    double b=10; // εκπέτασμα
    double ang=0.0873; // γωνία 5 μοιρών σε rad
    double c=2.0; // χορδή πτέρυγας
    double h=20.0; // το διπλάσιο ύψος που απέχει απο το έφαδος
    double dx=b/n;
    double Uinf=10.0; //m/s
    double Uind=0.0;
    int l=0;
    double EPSF =0.01;
    double a0=-0.0873;
    double a=0.0 ; // σε γωνία πτήσης 5 μοιρών
    double temp;
    double rho=1.2;
    vector <double> XA(3,0), XB(3,0),XO(3,0),U(3),B(n),GAMMA(n), RA(3,0),RB(3,0),AB(3,0),AL(3,0),W(n),L(n),D(n), y(n) ;
    double A1[n][n], A2[n][n], A[n][n], Aug[n][n+1], Aw1[n][n], Aw2[n][n], Aw[n][n];

    for(int i=0;i<n;i++){
        for(int j=0;j<n;j++){
            A[i][j]=0.0;
        }
    }
    for(int i=0;i<n;i++) {
        for (int j = 0; j < n; j++) {
            A1[i][j] = 0.0;
            A2[i][j] = 0.0;
        }
    }
          for( int i=0;i<n;i++){
        for(int j=0;j<=n;j++){
            Aug[i][j]=0.0;
        }
    }
          for (int i=0;i<3;i++){
              U[i]=0.0;
          }

    outdata.open("Y.txt");

    for (int i=0;i<n;i++) {
        if (i < n / 2) {
            XO[1] = dx * (i + 0.5);
            XO[0] = -tan(ang) * (b / 2 - dx * (i + 0.5));
            XO[2] = 0;
        }if (i >= n / 2){
            XO[1] = -dx * (i-n/2 + 0.5);
        XO[0] = -tan(ang) * (b/2 - dx * (i -n/2 +0.5));
        XO[2] = 0;
        }
        double g=1.0;
        for (int k=1;k<=m;k++) {
            if (k <= m / 2) {
                if (k % 3 == 1) {
                    XA[1] = dx * (g - 1);
                    XA[0] = 10 * c;
                    XA[2] = 0;
                    XB[1] = dx * (g - 1);
                    XB[0] = -tan(ang) * (b / 2 - dx * (g - 1));
                    XB[2] = 0;
                }


               
                 if (k % 3 == 2) {
                    XA[1] = dx * (g - 1);
                    XA[0] = -tan(ang) * (b / 2 - dx * (g - 1));
                    XA[2] = 0;
                    XB[1] = dx * g;
                    XB[0] = -tan(ang) * (b / 2 - dx * g);
                    XB[2] = 0;
                }
                if (k % 3 == 0) {
                    XA[1] = dx * g;
                    XA[0] = -tan(ang) * (b / 2 - dx * g);
                    XA[2] = 0;
                    XB[1] = dx * g;
                    XB[0] = 10 * c;
                    XB[2] = 0;
                    g++;
                }
            }

            if (k > m / 2) {
                if (k % 3 == 1) {
                    XB[1] = -dx * (g - n / 2 - 1);
                    XB[0] = 10 * c;
                    XB[2] = 0;
                    XA[1] = -dx * (g - n / 2 - 1);
                    XA[0] = -tan(ang) * (b / 2 - dx * (g - n / 2 - 1));
                    XA[2] = 0;
                }
                if (k % 3 == 2) {
                    XB[1] = -dx * (g - n / 2 - 1);
                    XB[0] =- tan(ang) * (b / 2 - dx * (g - n / 2 - 1));
                    XB[2] = 0;
                    XA[1] = -dx * (g - n / 2);
                    XA[0] = -tan(ang) * (b / 2 - dx * (g - n / 2));
                    XA[2] = 0;
                }
                if (k % 3 == 0) {
                    XB[1] = -dx * (g - n / 2);
                    XB[0] = -tan(ang) * (b / 2 - dx * (g - n / 2));
                    XB[2] = 0;
                    XA[1] = -dx * (g - n / 2);
                    XA[0] = 10 * c;
                    XA[2] = 0;
                    g++;
                }

            }
            UELEF0(EPSF, 1.0, XA, XB, XO, U);
            if (k % 3 != 0) {
                Uind += U[2];
            }
               



            if (k % 3 == 0) {
                    Uind += U[2];
                    l = (k / 3)-1 ;
                    A1[i][l] = Uind / Uinf;
                    Aw1[i][l]=Uind;
                    Uind = 0;
                }
                U[0]=0.0; U[1]=0.0; U[2]=0.0;
            }
           outdata << XO[1] << "\n";
           y[i]=XO[1];

        }

outdata.close();

    ////// το κατοπτρικο
    Uind=0.0;
    for (int i=0;i<n;i++) {
        if (i < n / 2) {
            XO[1] = dx * (i + 0.5);
            XO[0] =- tan(ang) * (b/2 -dx * (i + 0.5));
            XO[2] = 0;
        }
        if (i >= n / 2) {
            XO[1] = -dx * (i - n / 2 + 0.5);
            XO[0] =- tan(ang) * (b / 2 - dx * (i - n / 2 + 0.5));
            XO[2] = 0;
        }
        double g = 1.0;
        for (int k = 1; k <= m; k++) {
            if (k <= m / 2) {
                if (k % 3 == 1) {
                    XB[1] = dx * (g - 1);
                    XB[0] = 10 * c;
                    XB[2] = -h;
                    XA[1] = dx * (g - 1);
                    XA[0] = -tan(ang) * (b/2 -dx * (g - 1));
                    XA[2] = -h;
                }
                if (k % 3 == 2) {
                    XB[1] = dx * (g - 1);
                    XB[0] = -tan(ang) * (b/2 -dx * (g - 1));
                    XB[2] = -h;
                    XA[1] = dx * g;
                    XA[0] = -tan(ang) * (b/2 -dx * g);
                    XA[2] = -h;
                }
                if (k % 3 == 0) {
                    XB[1] = dx * g;
                    XB[0] = -tan(ang) * (b/2 -dx * g);
                    XB[2] = -h;
                    XA[1] = dx * g;
                    XA[0] = 10 * c;
                    XA[2] = -h;
                    g++;
                }
            }

            
              if (k > m / 2) {
                if (k % 3 == 1) {
                    XA[1] = -dx * (g - n / 2 - 1);
                    XA[0] = 10 * c;
                    XA[2] = -h;
                    XB[1] = -dx * (g - n / 2 - 1);
                    XB[0] = -tan(ang) * (b / 2 - dx * (g - n / 2 - 1));
                    XB[2] = -h;
                }
                if (k % 3 == 2) {
                    XA[1] = -dx * (g - n / 2 - 1);
                    XA[0] = -tan(ang) * (b / 2 - dx * (g - n / 2 - 1));
                    XA[2] = -h;
                    XB[1] = -dx * (g - n / 2);
                    XB[0] = -tan(ang) * (b / 2 - dx * (g - n / 2));
                    XB[2] = -h;
                }
                if (k % 3 == 0) {
                    XA[1] = -dx * (g - n / 2);
                    XA[0] = -tan(ang) * (b / 2 - dx * (g - n / 2));
                    XA[2] = -h;
                    XB[1] = -dx * (g - n / 2);
                    XB[0] = 10 * c;
                    XB[2] = -h;
                    g++;
                }
            }
            UELEF0(EPSF, 1.0, XA, XB, XO, U);
            if (k % 3 != 0) {
                Uind += U[2];
            }
            if (k % 3 == 0) {
                Uind += U[2];
                l = (k / 3) -1;
                A2[i][l] = Uind / Uinf;
                Aw2[i][l]=Uind;
                Uind = 0;

            }
            U[0] = 0.0;
            U[1] = 0.0;
            U[2] = 0.0;

        }
    }

    for(int i=0;i<n;i++) {
        for (int j = 0; j < n; j++) {
            A[i][j] = A1[i][j] + A2[i][j];
            Aw[i][j]= Aw1[i][j] +Aw2[i][j];
        }
    }



        
    for (int i=0;i<n;i++){
            A[i][i]=A[i][i] -1/(Uinf*3.1415*c);
            B[i]=a0-a;
        }

        for(int i=0;i<n;i++) {
            for (int j = 0; j < n ; j++) {
                Aug[i][j] = A[i][j];
            }

            Aug[i][n]=B[i];
        }

/// Gauss Elimination
    for(int j=0; j<n-1; j++) {
        for (int i = j + 1; i < n; i++) {
            temp = Aug[i][j] / Aug[j][j];

            for (int k = 0; k < n + 1; k++) {
                Aug[i][k] -= Aug[j][k] * temp;
            }
        }
    }
        for(int i=n-1; i>=0; i--)
        {
           double  s=0;
            for(int j=i+1; j<n; j++) {
                s += Aug[i][j] * GAMMA[j];
              //  GAMMA[i] = (Aug[i][n] - s) / Aug[i][i];
            }
            GAMMA[i] = (Aug[i][n] - s) / Aug[i][i];
        }



    outdata.open("Gamma.txt");
        for(int i=0;i<n;i++){
            cout<<GAMMA[i]<<endl;
            outdata<<GAMMA[i]<<"\n";
        }
    outdata.close();

        for(int i=0;i<n;i++){
            for(int j=0;j<n;j++){
                W[i] += Aw[i][j]*GAMMA[j];
            }
        }

        outdata.open("W.txt");
    for(int i=0;i<n;i++){
        cout<<W[i]<<endl;
        outdata<<W[i]<<"\n";
    }
    outdata.close();

        
for (int i=0;i<n;i++){
    L[i] = rho * Uinf * GAMMA[i] * dx;
    D[i] = rho * W[i] * GAMMA[i] * dx;
    }

        

        outdata.open("Lift.txt");
    for(int i=0;i<n;i++){
        cout<<L[i]<<endl;
        outdata<<L[i]<<"\n";
    }
    outdata.close();

    outdata.open("Drag.txt");
    for(int i=0;i<n;i++){
        cout<<D[i]<<endl;
        outdata<<D[i]<<"\n";
    }
    outdata.close();

 return 0;
}
