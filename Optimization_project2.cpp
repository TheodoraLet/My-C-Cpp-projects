#include <iostream>
#include <fstream>
#include <cmath>
#include <vector>
#include <string>
#include <valarray>
#include <iomanip>
#include<fstream>


using namespace std;
using std:: ofstream;



double M1=0;
double M2 =6;

double tfun(double t,double x0,double x1,double x2,double x3,double x4,double x5, double xt){
   return (t*t*t*t*t*(x5-5*x4+10*x3-10*x2+5*x1-x0)+t*t*t*t*(5*x4-20*x3+30*x2-20*x1+5*x0)+t*t*t*(10*x3-30*x2+30*x1-10*x0)+t*t*(10*x2-20*x1+10*x0)+t*(5*x1-5*x0)+x0)-xt;
};
double dtfun(double t, double x0,double x1,double x2,double x3,double x4,double x5){
  return  5*t*t*t*t*(x5-5*x4+10*x3-10*x2+5*x1-x0)+4*t*t*t*(5*x4-20*x3+30*x2-20*x1+5*x0)+3*t*t*(10*x3-30*x2+30*x1-10*x0)+2*t*(10*x2-20*x1+10*x0)+5*x1-5*x0;
};

double newton(double t0, double er, int nmax, double tr, double ea, double fr,double x0,double x1,double x2,double x3,double x4,double x5,double xt)
{
//--  initial values
    int niter=0;
    double told = t0;
    fr=tfun(told, x0, x1, x2,x3,x4,x5,xt);
//
//___start of iterations __________
//
   // cout<<" niter tr ea fr\n";
    for(niter;niter<nmax;niter++)
    {

        double fgr = dtfun(told,x0,x1,x2,x3,x4,x5);
        if(fgr!=0)
        {
            tr = told - fr/fgr;
        }
        else
        {
            cout<<"error!  zero derivative"<<endl;
            return 0;
        }
//-- xr is a roote?
        fr=tfun(tr,x0,x1,x2,x3,x4,x5,xt);
        if(fr==0)
        {
            ea=0;
          return tr;
        }
//-- relatice (or absolute) error
        if(tr!=0)
        {
            ea = fabs((tr-told)/tr);
        }
        else
        {
            ea = fabs(tr-told);
        }
//-- monitoring (optional)
     //   cout<<niter<<" "<<tr<<" "<<ea<<" "<<fr<<"\n";
//-- exit checks
        if(ea<=er)	return tr;
        if(niter>=nmax)
        {
            cout<<"warning! iterations limit"<<endl;
            return 0;
        }
        told=tr;
    }
    return tr;
}


double S (double t,double b0,double b1,double b2,double b3,double b4,double b5){

  return t*t*t*t*t*(b5-5*b4+10*b3-10*b2+5*b1-b0)+t*t*t*t*(5*b4-20*b3+30*b2-20*b1+5*b0)+t*t*t*(10*b3-30*b2+30*b1-10*b0)+t*t*(10*b2-20*b1+10*b0)+t*(5*b1-5*b0)+b0;
};

double fun(double T, double x,double S){
    return 0.01*T*T*S -32*x ;

}


double F(double T, double T0){
    double obj=T-T0;
    return pow(obj,2);
}

double trap(int n, vector <double> xg, vector <double> Tg, double &sum)
{
    sum=0.e0;
    for(int i=1;i<n+1;i++)
    {
        if(xg[i]<xg[i-1])
        {
            cout<<"x NOT in INCREASING order"<<endl;
            break;
        }
        sum=sum+.5e0*(xg[i]-xg[i-1])*(Tg[i]+Tg[i-1]);
    }
    return sum ;
}

double Fb1fun(double t){
    return 5*t*t*t*t*t-20*t*t*t*t+30*t*t*t-20*t*t+5*t;
}

double Fb2fun(double t){
    return -10*t*t*t*t*t+30*t*t*t*t-30*t*t*t+10*t*t;
}

double Fb3fun(double t){
    return 10*t*t*t*t*t-20*t*t*t*t+10*t*t*t;
}

double Fb4fun(double t){
    return -5*t*t*t*t*t+5*t*t*t*t;
}

double fae(double Tf,double T0,double S,double Psi){
    return 2*(Tf-T0)-0.02*Tf*S*Psi;
}

void tdma(vector<double> &a, vector<double> &b,vector<double> &c,vector<double> &d,int n){
    n--;
    c[0] /= b[0];
    d[0] /= b[0];

    for (int i = 1; i < n; i++) {
        c[i] /= b[i] - a[i]*c[i-1];
        d[i] = (d[i] - a[i]*d[i-1]) / (b[i] - a[i]*c[i-1]);
    }

    d[n] = (d[n] - a[n]*d[n-1]) / (b[n] - a[n]*c[n-1]);

    for (int i = n; i-- > 0;) {
        d[i] -= c[i]*d[i+1];
    }
}


double DirectDif(double T,double DerivTb,double y,double DerivS,double x){
    return 0.02*T*DerivTb*y+0.01*T*T*DerivS-32*x;
}


int main() {
     ofstream fdL;
     ofstream fd1;
     ofstream fd2;
     ofstream fd3;
     ofstream fd4;
    double dx;
    double L =4;
    int n=200;
   vector<double> T(1), aux(3), x(n+1),taf(n+1),Ttot(n+1),y(n+1),Tg(n+1),Psi(n+1),tempPsi(1),termLfun(n+1),deltaFL(2),termxfun(n+1),deltaFb1(2);
   vector<double> termS1fun(n+1),termS2fun(n+1),termS3fun(n+1),termS4fun(n+1),derivFb1(n+1),derivFb2(n+1),derivFb3(n+1),derivFb4(n+1),deltaFb2(2),deltaFb3(2),deltaFb4(2);
   vector<double> TL(1),Tc1(1),Tc2(1),Tc3(1),Tc4(1), TtL(n+1),Ttc1(n+1),Ttc2(n+1),Ttc3(n+1),Ttc4(n+1),TgL(n+1),Tgc1(n+1),Tgc2(n+1),Tgc3(n+1),Tgc4(n+1),tafn(n+1),xn(n+1);
    vector<double>TLm(1),Tc1m(1),Tc2m(1),Tc3m(1),Tc4m(1), TtLm(n+1),Ttc1m(n+1),Ttc2m(n+1),Ttc3m(n+1),Ttc4m(n+1),TgLm(n+1),Tgc1m(n+1),Tgc2m(n+1),Tgc3m(n+1),Tgc4m(n+1),tafnm(n+1),xnm(n+1);
   double dRdb[n+1][4], dFdb[1][4], dTdb[n+1][4],dRdT[n+1][n+1],dFdbnew[1][4];
   vector<double> DirTb1(n+1,0),DirTb2(n+1,0),DirTb3(n+1,0),DirTb4(n+1,0),tempDirTb1(1),tempDirTb2(1),tempDirTb3(1),tempDirTb4(1),DDfb1(n+1),DDfb2(n+1),DDfb3(n+1),DDfb4(n+1);
   vector<double> b(n+1),a(n+1),c(n+1),d(n+1),bF(n+1),aF(n+1),cF(n+1),dF(n+1),dF2(n+1),dF3(n+1),dF4(n+1);
    vector<double> bF2(n+1),aF2(n+1),cF2(n+1),bF3(n+1),aF3(n+1),cF3(n+1),bF4(n+1),aF4(n+1),cF4(n+1);
   vector<double> dTdb1(n+1),dTdb2(n+1),dTdb3(n+1),dTdb4(n+1), dFdb1(n+1),dFdb2(n+1),dFdb3(n+1),dFdb4(n+1),dFdT(n+1),dRb1(n+1),dRb2(n+1),dRb3(n+1),dRb4(n+1),DDFb1(1),DDFb2(1),DDFb3(1),DDFb4(1);;
   //T[0]=3;
    double x0=0;
    double x1=0.5;
    double x2=1.2;
    double x3=1.5;
    double x4=1.9;
    double x5=L;
    double b1=1.5;
    double b2=7;
    double b3=2;
    double b4=5;
    double b5=3;
    double b0=1;
    double t0=0.0;
    double er=pow(10,-4);
     x[0]=0.0;
    double k1, k2 ,k3 , k4;
    Ttot[0]=T[0];
    int nmax=20;
    double ea;
    double tr;
    double h=L/n;
    double fr;
    double fgr;
    double sum, sumL,sum1,sum2,sum3,sum4,sumn;
    double sumLm,sum1m,sum2m,sum3m,sum4m;
    double sumb1,sumb2,sumb3,sumb4;
    double result;
    double resultx,results1,results2,results3,results4;
    double trn; double trnm;
    double epsilon=pow(10,-3);
    double Ln; double Lnm;
    double dxn;double dxnm;
    double x5n; double x5nm;

    deltaFL[0]=0;deltaFL[1]=0.1; deltaFb1[0]=0;deltaFb1[1]=0.1;deltaFb2[0]=0;deltaFb2[1]=0.1;deltaFb3[0]=0;deltaFb3[1]=0.1;deltaFb4[0]=0;deltaFb4[1]=0.1;

   while(abs(deltaFL[1]-deltaFL[0])>0.001 &&abs(deltaFb1[1]-deltaFb1[0])>0.001 && abs(deltaFb2[1]-deltaFb2[0])>0.001 &&abs(deltaFb3[1]-deltaFb3[0])>0.001 &&abs(deltaFb4[1]-deltaFb4[0])>0.001 && (2<=L) && (4>=L)){
       Ln=L+epsilon;
       Lnm=L-epsilon;
        x5 = L;
        x5n=Ln;
        x5nm=Lnm;
        dx = L / n;
        dxn=Ln/n;
        dxnm=Lnm/n;
        t0 = 0;
        x[0] = 0.0;
        xn[0]=0.0;
        xnm[0]=0.0;
        taf[0] = newton(t0, er, nmax, tr, ea, fr, x0, x1, x2, x3, x4, x5, x[0]);
        tafn[0]=newton(t0,er,nmax,trn,ea,fr,x0,x1,x2,x3,x4,x5n,xn[0]);
        tafnm[0]=newton(t0,er,nmax,trnm,ea,fr,x0,x1,x2,x3,x4,x5nm,xnm[0]);
        T[0] = pow(L, -1) - 25 * pow(L, -2);
        Tc1[0]=T[0];
        Tc2[0]=T[0];
        Tc3[0]=T[0];
        Tc4[0]=T[0];
        Tc1m[0]=T[0];
        Tc2m[0]=T[0];
        Tc3m[0]=T[0];
        Tc4m[0]=T[0];
        TL[0] = pow(Ln + epsilon, -1) - 25 * pow(Ln+epsilon, -2);
        TLm[0] = pow(Ln - epsilon, -1) - 25 * pow(Ln-epsilon, -2);
        Ttot[0] = T[0];
        TtL[0]=TL[0];
        TtLm[0]=TLm[0];
        Ttc1[0] = Tc1[0];
        Ttc2[0] = Tc2[0];
        Ttc3[0] = Tc3[0];
        Ttc4[0] = Tc4[0];
        Ttc1m[0] = Tc1m[0];
        Ttc2m[0] = Tc2m[0];
        Ttc3m[0] = Tc3m[0];
        Ttc4m[0] = Tc4m[0];
        y[0] = S(taf[0], b0, b1, b2, b3, b4, b5);
        for (int i = 1; i < n + 1; i++) {
            x[i] = x[i - 1] + dx;
            xn[i]=xn[i-1]+dxn;
            xnm[i]=xnm[i-1]+dxnm;
            taf[i] = newton(t0, er, nmax, tr, ea, fr, x0, x1, x2, x3, x4, x5, x[i]);
            tafn[i]=newton(t0,er,nmax,trn,ea,fr,x0,x1,x2,x3,x4,x5n,xn[i]);
            tafnm[i]=newton(t0,er,nmax,trn,ea,fr,x0,x1,x2,x3,x4,x5n,xn[i]);

            double slope = fun(T[0], x[i], S(taf[i], b0, b1, b2, b3, b4, b5));      /// Λύνω το primal για το CA και τα 10 primal για τις FD
            T[0] = T[0] + dx * slope;
            double slopeL = fun(TL[0], xn[i], S(tafn[i], b0, b1, b2, b3, b4, b5));
            TL[0] = TL[0] + dxn * slopeL;
            double slopeLm = fun(TLm[0], xnm[i], S(tafnm[i], b0, b1, b2, b3, b4, b5));
            TLm[0] = TLm[0] + dxnm * slopeLm;
            double slopec1 = fun(Tc1[0], x[i], S(taf[i], b0, b1+epsilon, b2, b3, b4, b5));
            Tc1[0] = Tc1[0] + dx * slopec1;
            double slopec1m = fun(Tc1m[0], x[i], S(taf[i], b0, b1-epsilon, b2, b3, b4, b5));
            Tc1m[0] = Tc1m[0] + dx * slopec1m;
            double slopec2 = fun(Tc2[0], x[i], S(taf[i], b0, b1, b2+epsilon, b3, b4, b5));
            Tc2[0] = Tc2[0] + dx * slopec2;
            double slopec2m = fun(Tc2m[0], x[i], S(taf[i], b0, b1, b2-epsilon, b3, b4, b5));
            Tc2m[0] = Tc2m[0] + dx * slopec2m;
            double slopec3 = fun(Tc3[0], x[i], S(taf[i], b0, b1, b2, b3+epsilon, b4, b5));
            Tc3[0] = Tc3[0] + dx * slopec3;
            double slopec3m = fun(Tc3m[0], x[i], S(taf[i], b0, b1, b2, b3-epsilon, b4, b5));
            Tc3m[0] = Tc3m[0] + dx * slopec3m;
            double slopec4 = fun(Tc4[0], x[i], S(taf[i], b0, b1, b2, b3, b4+epsilon, b5));
            Tc4[0] = Tc4[0] + dx * slopec4;
            double slopec4m = fun(Tc4m[0], x[i], S(taf[i], b0, b1, b2, b3, b4-epsilon, b5));
            Tc4m[0] = Tc4m[0] + dx * slopec4m;
            Ttot[i] = T[0];
            TtL[i] = TL[0];
            TtLm[i] = TLm[0];
            Ttc1[i] = Tc1[0];
            Ttc1m[i] = Tc1m[0];
            Ttc2[i] = Tc2[0];
            Ttc2m[i] = Tc2m[0];
            Ttc3[i] = Tc3[0];
            Ttc3m[i] = Tc3m[0];
            Ttc4[i] = Tc4[0];
            Ttc4m[i] = Tc4m[0];
            y[i] = S(taf[i], b0, b1, b2, b3, b4, b5);
        }
       ofstream outdata;
       outdata.open("TemperaturePrimal.txt");
       if (!outdata){
           cout<<"Error file could not open"<<"\t";
           exit(1);
       }

        for (int i = 0; i < n + 1; i++) {
            cout << Ttot[i] << "\t";
            outdata<<Ttot[i]<<endl;
            cout << taf[i] << "\t";
            cout << x[i] << "\t";
            cout << y[i] << "\t";
        }

        outdata.close();
       outdata.open("XPrimal.txt");
       if (!outdata){
           cout<<"Error file could not open"<<"\t";
           exit(1);
       }
       for(int i=0;i<n+1;i++){
           outdata<<x[i]<<endl;
       }
        outdata.close();
       outdata.open("SPrimal.txt");
       if (!outdata){
           cout<<"Error file could not open"<<"\t";
           exit(1);
       }
       for(int i=0;i<n+1;i++){
           outdata<<y[i]<<endl;
       }
       outdata.close();

        for (int i = 0; i < n + 1; i++) {  /// Υπολογίζω όρους για την objective function
            Tg[i] = F(Ttot[i], Ttot[0]);
            TgL[i] = F(TtL[i], TtL[0]);
            TgLm[i] = F(TtLm[i], TtLm[0]);
            Tgc1[i] = F(Ttc1[i], Ttc1[0]);
            Tgc1m[i] = F(Ttc1m[i], Ttc1m[0]);
            Tgc2[i] = F(Ttc2[i], Ttc2[0]);
            Tgc2m[i] = F(Ttc2m[i], Ttc2m[0]);
            Tgc3[i] = F(Ttc3[i], Ttc3[0]);
            Tgc3m[i] = F(Ttc3m[i], Ttc3m[0]);
            Tgc4[i] = F(Ttc4[i], Ttc4[0]);
            Tgc4m[i] = F(Ttc4m[i], Ttc4m[0]);
            //cout<<Tg[i]<<"\t";
        }

        double ObjFun = trap(n, x, Tg, sum);  /// Υπολογίζω την objective function
       double ObjFL = trap(n, xn, TgL, sumL);
       double ObjFLm = trap(n, xnm, TgLm, sumLm);
       double ObjFc1 = trap(n, x, Tgc1, sum1);
       double ObjFc1m = trap(n, x, Tgc1m, sum1m);
       double ObjFc2 = trap(n, x, Tgc2, sum2);
       double ObjFc2m = trap(n, x, Tgc2m, sum2m);
       double ObjFc3 = trap(n, x, Tgc3, sum3);
       double ObjFc3m = trap(n, x, Tgc3m, sum3m);
       double ObjFc4 = trap(n, x, Tgc4, sum4);
       double ObjFc4m = trap(n, x, Tgc4m, sum4m);
        // cout<<ObjFun<<"\n";

        Psi[n] = 0;  /// Λύνω με Euler με αρνητικό βήμα για να βρω το Ψ
        tempPsi[0] = Psi[n];
        for (int i = n - 1; i > -1; i--) {
            double slopePsi = fae(Ttot[i], Ttot[0], y[i], tempPsi[0]);
            tempPsi[0] = tempPsi[0] - dx * slopePsi;
            Psi[i] = tempPsi[0];
              //cout<<Psi[i]<<"\t";
        }

//        for (int i = 0; i < n + 1; i++) {
//            cout << Psi[i] << "\t";
//        }
        for (int i = 0; i < n + 1; i++) {   /// Υπολογίζω όρους για να βρω τις παραγώγους
            termLfun[i] = 2 * (Ttot[i] - Ttot[0]) * (50 * pow(L, -3) - pow(L, -2));
            termxfun[i] = 32* Psi[i] / n;
            derivFb1[i] = Fb1fun(taf[i]);
            derivFb2[i] = Fb2fun(taf[i]);
            derivFb3[i] = Fb3fun(taf[i]);
            derivFb4[i] = Fb4fun(taf[i]);
            termS1fun[i] = 0.01 * Psi[i] * Ttot[i] * Ttot[i] * derivFb1[i];
            termS2fun[i] = 0.01 * Psi[i] * Ttot[i] * Ttot[i] * derivFb2[i];
            termS3fun[i] = 0.01 * Psi[i] * Ttot[i] * Ttot[i] * derivFb3[i];
            termS4fun[i] = 0.01 * Psi[i] * Ttot[i] * Ttot[i] * derivFb4[i];
        }

        double termL = trap(n, x, termLfun, result); /// Υπολογίζω όρους για να βρω τις παραγώγους
        double termx = trap(n, x, termxfun, resultx);
        double termS1 = trap(n, x, termS1fun, results1);
        double termS2 = trap(n, x, termS2fun, results2);
        double termS3 = trap(n, x, termS3fun, results3);
        double termS4 = trap(n, x, termS4fun, results4);
//        cout<<termL<<"\t";
//        cout<<termx<<"\t";
//        cout<<termS1<<"\t";
//        cout<<termS2<<"\t";
//        cout<<termS3<<"\t";
//        cout<<termS4<<"\t";

       double fdL = (ObjFL - ObjFLm) / (2*epsilon);  /// Υπολογίζω τις παραγώγους για τις πεπερασμένες διαφορές
       double fdc1 = (ObjFc1 - ObjFc1m) / (2*epsilon);
       double fdc2 = (ObjFc2 - ObjFc2m) / (2*epsilon);
       double fdc3 = (ObjFc3 - ObjFc3m) / (2*epsilon);
       double fdc4 = (ObjFc4 - ObjFc4m) / (2*epsilon);

       cout<<"Finite Differences"<<"\t";
       cout << fdL << "\n";
       cout << fdc1 << "\n";
       cout << fdc2 << "\n";
       cout << fdc3 << "\n";
       cout << fdc4 << "\n";

       outdata.open("FdL.text",ios_base::app);
       outdata<<fdL<<endl;
       outdata.close();
       outdata.open("Fd1.text",ios_base::app);
       outdata<<fdc1<<endl;
       outdata.close();
       outdata.open("Fd2.text",ios_base::app);
       outdata<<fdc2<<endl;
       outdata.close();
       outdata.open("Fd3.text",ios_base::app);
       outdata<<fdc3<<endl;
       outdata.close();
       outdata.open("Fd4.text",ios_base::app);
       outdata<<fdc4<<endl;
       outdata.close();



        deltaFL[0] = deltaFL[1]; /// Όροι για να βρω παραγώγους για το CA
        deltaFL[1] = pow(Ttot[n] - Ttot[0], 2) - termL - Psi[0] * (50 * pow(L, -3) - pow(L, -2)) + termx;
        deltaFb1[0] = deltaFb1[1];
        deltaFb1[1] =   - termS1;
        deltaFb2[0] = deltaFb2[1];
        deltaFb2[1] =   - termS2;
        deltaFb3[0] = deltaFb3[1];
        deltaFb3[1] =   - termS3;
        deltaFb4[0] = deltaFb4[1];
        deltaFb4[1] =   - termS4;
        cout << deltaFL[1] << "\n";  /// παράγωγοι CA
        cout << deltaFb1[1] << "\n";
        cout << deltaFb2[1] << "\n";
        cout << deltaFb3[1] << "\n";
        cout << deltaFb4[1] << "\n";
       outdata.open("DfL.text",ios_base::app);
        outdata << deltaFL[1] << "\n";
       outdata.close();
       outdata.open("Df1.text",ios_base::app);
        outdata << deltaFb1[1] << "\n";
       outdata.close();
       outdata.open("Df2.text",ios_base::app);
        outdata<< deltaFb2[1] << "\n";
       outdata.close();
       outdata.open("Df3.text",ios_base::app);
        outdata << deltaFb3[1] << "\n";
       outdata.close();
       outdata.open("Df4.text",ios_base::app);
        outdata<< deltaFb4[1] << "\n";
        outdata.close();

        double hta = pow(10,-5);

        double Lold=L;
        L = L - hta * deltaFL[1];
        if(L<2 || L>4){
            L=Lold;
            break;
        }else
           /// Βελτιστοποίηση με Steepest Descent
            b1 = b1 - hta * deltaFb1[1];
            b2 = b2 - hta * deltaFb2[1];
            b3 = b3 - hta * deltaFb3[1];
            b4 = b4 - hta * deltaFb4[1];

            cout << L << "\n";
            cout << b1 << "\n";
            cout << b2 << "\n";
            cout << b3 << "\n";
            cout << b4 << "\n";
       outdata.open("L.text",ios_base::app);
       outdata<<L<<endl;
       outdata.close();
       outdata.open("b1.text",ios_base::app);
       outdata<<b1<<endl;
       outdata.close();
       outdata.open("b2.text",ios_base::app);
       outdata<<b2<<endl;
       outdata.close();
       outdata.open("b3.text",ios_base::app);
       outdata<<b3<<endl;
       outdata.close();
       outdata.open("b4.text",ios_base::app);
       outdata<<b4<<endl;
       outdata.close();


   }

///////////////// Discrete
    b[0]=1.e0;    ///// Λύνω με τριδιαγώνιο για να βρω το Ψ
    d[0]=-dx*(n-2)-(Ttot[n]-Ttot[0])*dx;
    a[0]=0.e0;
    a[1]=0.e0;

    for (int i = 1; i < n; i++) {
        b[i] = -0.02* Ttot[i] *y[i];
        d[i] = 2*(Ttot[i]-Ttot[0])*dx;
    }

    for (int i=2;i<n+1;i++){
        a[i] = pow(2 * dx, -1);
    }
    for(int i=0;i<n-1;i++){
        c[i] = -pow(2 * dx, -1);
    }

    int m=n+1;
    c[n-1]=-pow(dx,-1);
    c[n]=0.e0;
    b[n] =-0.02* Ttot[n]* y[n]+1/dx;
    d[n] = (Ttot[n]-Ttot[0])*dx;
    tdma( a,b,c,d,m);
    for (int i = 0; i < n; i++) {
        cout << d[i] << endl;
    }

        dRdb[0][0]=0.0;
        dRdb[0][1]=0.0;
        dRdb[0][2]=0.0;
        dRdb[0][3]=0.0;
    for (int i=1;i<n+1;i++){    /// Υπολογίζω το θR/θb
        dRdb[i][0]=-0.01*Ttot[i]*Ttot[i]*Fb1fun(taf[i]);
        dRdb[i][1]=-0.01*Ttot[i]*Ttot[i]*Fb2fun(taf[i]);
        dRdb[i][2]=-0.01*Ttot[i]*Ttot[i]*Fb3fun(taf[i]);
        dRdb[i][3]=-0.01*Ttot[i]*Ttot[i]*Fb4fun(taf[i]);
    }

    for(int i=0;i<n+1;i++){
        cout<<dRdb[i][0]<<"\t";
        cout<<dRdb[i][1]<<"\t";
        cout<<dRdb[i][2]<<"\t";
        cout<<dRdb[i][3]<<"\t";
    }

    for(int i=0;i<1;i++)
        for (int j = 0; j < 4; j++)
        {
            dFdb[i][j] = 0.0;
        }

    for(int i=0;i<1;i++)
        for (int j = 0; j < 4; j++)
            for (int k = 0; k < n + 1; k++)
            {
                dFdb[i][j] += -d[k] * dRdb[k][j];    /// Με πολλαπλασιασμό πινάκων βρίσκω το θF/θb
            }


    for(int i=0;i<1;i++)
        for (int j = 0; j < 4; j++)
        {
            cout << "Discrete Derivatives" << "\t";
            cout << dFdb[i][j] << "\t";
        }



//////////////////////////////////// Direct

     for(int i=0;i<n+1;i++){
         dRb1[i]=dRdb[i][0];
         dRb2[i]=dRdb[i][1];
         dRb3[i]=dRdb[i][2];
         dRb4[i]=dRdb[i][3];
     }

     dFdT[0]=0;                    /// Υπολογίζω με 4 τριδιαγώνιους τα θT/θb
     for (int i=1;i<n;i++){
         dFdT[i]=2*(Ttot[i]-Ttot[0])*dx;
     }
     dFdT[n]=(Ttot[n]-Ttot[0])*dx;
/// 1os
    bF[0] = 1.e0;
    dF[0] = -dRb1[0];
    cF[0] = 0;
    aF[0] = 0;

    for (int i = 1; i < n; i++) {
        bF[i] = -0.02 * Ttot[i] * y[i];
        dF[i] = -dRb1[i];
    }
    for (int i = 1; i < n; i++) {
        aF[i] = -pow(2 * dx, -1);
    }
    for (int i = 1; i < n; i++) {
        cF[i] = pow(2 * dx, -1);
    }

    aF[n] = -pow(dx, -1);
    bF[n] = -0.02 * Ttot[n] * y[n] + 1 / dx;
    dF[n] = -dRb1[n];
    cF[n]=0;
    tdma(aF, bF, cF, dF, m);
    for (int i = 0; i < n + 1; i++) {
        dTdb1[i] = dF[i];
        //  cout<<dTdb1[i]<<"\t";
    }
/// 2os
    bF2[0] = 1.e0;
    dF2[0] = -dRb2[0];
    cF2[0] = 0;
    aF2[0] = 0;

    for (int i = 1; i < n; i++) {
        bF2[i] = -0.02 * Ttot[i] * y[i];
        dF2[i] = -dRb2[i];
    }
    for (int i = 1; i < n; i++) {
        aF2[i] = -pow(2 * dx, -1);
    }
    for (int i = 1; i < n; i++) {
        cF2[i] = pow(2 * dx, -1);
    }

    aF2[n] = -pow(dx, -1);
    bF2[n] = -0.02 * Ttot[n] * y[n] + 1 / dx;
    dF2[n] = -dRb2[n];
    cF2[n]=0;
    tdma(aF2, bF2, cF2, dF2, m);
    for (int i = 0; i < n + 1; i++) {
        dTdb2[i] = dF2[i];
        // cout<<dTdb2[i]<<"\t";
    }
/// 3os
    bF3[0] = 1.e0;
    dF3[0] = -dRb3[0];
    cF3[0] = 0;
    aF3[0] = 0;

    for (int i = 1; i < n; i++) {
        bF3[i] = -0.02 * Ttot[i] * y[i];
        dF3[i] = -dRb3[i];
    }
    for (int i = 1; i < n; i++) {
        aF3[i] = -pow(2 * dx, -1);
    }
    for (int i = 1; i < n; i++) {
        cF3[i] = pow(2 * dx, -1);
    }

    aF3[n] = -pow(dx, -1);
    bF3[n] = -0.02 * Ttot[n] * y[n] + 1 / dx;
    dF3[n] = -dRb3[n];
    cF3[n]=0;
    tdma(aF3, bF3, cF3, dF3, m);
    for (int i = 0; i < n + 1; i++) {
        dTdb3[i] = dF3[i];
        //cout<<dTdb3[i]<<"\t";
    }
/// 4os
    bF4[0] = 1.e0;
    dF4[0] = -dRb4[0];
    cF4[0] = 0;
    aF4[0] = 0;

    for (int i = 1; i < n; i++) {
        bF4[i] = -0.02 * Ttot[i] * y[i];
        dF4[i] = -dRb4[i];
    }
    for (int i = 1; i < n; i++) {
        aF4[i] = -pow(2 * dx, -1);
    }
    for (int i = 1; i < n; i++) {
        cF4[i] = pow(2 * dx, -1);
    }

    aF4[n] = -pow(dx, -1);
    bF4[n] = -0.02 * Ttot[n] * y[n] + 1 / dx;
    dF4[n] = -dRb4[n];
    cF4[n]=0;
    tdma(aF4, bF4, cF4, dF4, m);
    for (int i = 0; i < n + 1; i++) {
        dTdb4[i] = dF4[i];
        // cout<<dTdb4[i]<<"\t";
    }

    //for (int j=0;j<1;j++)
    // for(int i=0;i<4;i++)
    for (int k = 0; k < n + 1; k++)  /// Με πολλαπλασιασμό πινάκων θF/θΤ *θΤ/θb βρίσκω τις παραγώγους
    {
        DDFb1[0] += dFdT[k] * dTdb1[k];
        DDFb2[0] += dFdT[k] * dTdb2[k];
        DDFb3[0] += dFdT[k] * dTdb3[k];
        DDFb4[0] += dFdT[k] * dTdb4[k];
    }

    //for(int i=0;i<1;i++)
    //   for(int j=0;j<4;j++)
    cout<<DDFb1[0]<<"\t";
    cout<<DDFb2[0]<<"\t";
    cout<<DDFb3[0]<<"\t";
    cout<<DDFb4[0]<<"\t";



    return 0;
}

