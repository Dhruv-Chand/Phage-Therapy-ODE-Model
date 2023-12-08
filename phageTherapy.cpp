#include <iostream>
#include <cmath>
#include "cpgplot.h"

double dxdt(double x, double v, double a, double b) {
    return (a * x) - (b * v * x);
}

double dydt(double y, double v, double x, double a, double b, double k) {
    return (a * y) + (b * v * x) - (k * y);
}

double dvdt(double v, double y, double x, double k, double L, double m, double b) {
    return (k * L * y) - (b * v * x) - (m * v);
}

int rk4(double a, double b, double k, double L, double m) {
    //intilize variables
    int n;
    double tstart;
    double tend; 
    double t;
    double kx1 = 0, kx2 = 0, kx3 = 0, kx4 = 0;
    double ky1 = 0, ky2 = 0, ky3 = 0, ky4 = 0;
    double kv1 = 0, kv2 = 0, kv3 = 0, kv4 = 0;
    tstart = 0; //fixed for now
    tend = 70; //fixed for now
    n = 1000000; //set n to 1000000 for now
    double deltaT = (tend - tstart)/n;
    t = tstart; 
    double x = 3 * pow(10, 4); // uninfected
    double y = 1 * pow(10, 4); // infected
    double v = 100;
    int i = 0;
    float tp[n+1];
    float yp[n+1];
    float xp[n+1];
    float vp[n+1];
    tp[0] = t;
    yp[0] = y;
    xp[0] = x;
    vp[0] = v;
    //rk4 estimation Loop
    while (i < n) {
        //std::cout << "t = " << tp[i] << " x = " << xp[i]<< " y = " << yp[i] << " v = " << vp[i] << "\n";
        kx1 = deltaT * (dxdt(x, v, a, b));
        ky1 = deltaT * (dydt(y, v, x, a, b, k));
        kv1 = deltaT * (dvdt(v, y, x, k, L, m, b));
        //k2
        kx2 = deltaT * (dxdt(x + 0.5 * kx1, v + 0.5 * kv1, a, b));
        ky2 = deltaT * (dydt(y + 0.5 * ky1, v + 0.5 * kv1, x + 0.5 * kx1, a, b, k));
        kv2 = deltaT * (dvdt(v + 0.5 * kv1, y + 0.5 * ky1, x + 0.5 * kx1, k, L, m, b));
        //k3
        kx3 = deltaT * (dxdt(x + 0.5 * kx2, v + 0.5 * kv2, a, b));
        ky3 = deltaT * (dydt(y + 0.5 * ky2, v + 0.5 * kv2, x + 0.5 * kx2, a, b, k));
        kv3 = deltaT * (dvdt(v + 0.5 * kv2, y + 0.5 * ky2, x + 0.5 * kx2, k, L, m, b));
        //k4
        kx4 = deltaT * (dxdt(x + kx3, v + kv3, a, b));
        ky4 = deltaT * (dydt(y + ky3, v + kv3, x + kx3, a, b, k));
        kv4 = deltaT * (dvdt(v + kv3, y + ky3, x + kx3, k, L, m, b));
        //iterate values
        t = t + deltaT;
        x = x + (1.0/6.0) * (kx1 + 2*kx2 + 2*kx3 + kx4);
        y = y + (1.0/6.0) * (ky1 + 2*ky2 + 2*ky3 + ky4);
        v = v + (1.0/6.0) * (kv1 + 2*kv2 + 2*kv3 + kv4);
        i++;
        tp[i] = t;
        yp[i] = y;
        xp[i] = x;
        vp[i] = v;
    }
    for (int pop = 0; pop < n; pop= pop + 1000) {
        std::cout << "t = " << tp[pop] << " x = " << xp[pop]<< " y = " << yp[pop] << " v = " << vp[pop] << "\n";
        //yp[pop] = yp[pop] + xp[pop];
    }
    if (!cpgopen("/XWINDOW")) return 1;
    cpgenv(tstart, tend, 0, 10000, 0, 0);
    cpglab("Inoculation Time (hours)", "Severity (Bacteria Amount)", "Bacteriophage Graph (Green = Infected Bac., Pink = Uninfected Bac., Blue = Phage Bac.)");
    cpgsci(3);//green y = infected
    cpgline(n + 1, tp, yp);
    cpgsci(6);//pink x = uninfected
    cpgline(n + 1, tp, xp);
    cpgsci(12);//purple v = free phage
    cpgline(n + 1, tp, vp);
    cpgclos();
    return 0;
    
}


int main() {
    std::cout << "Please input values for a, b, k, L, m (or input default to apply default values)";
    double a = 0.3;
    double b = pow(10, (-6));
    double k = 1.2;
    double L = 100;
    double m = 1.8;
    double x = 3 * pow(10, 4); // uninfected
    double y = 1 * pow(10, 4); // infected
    double v = 100;// phage
    rk4(a, b, k, L, m);
    return 0;
}
