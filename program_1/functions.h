//
//  functions.h
//  program_1
//
//  Created by Kozlov Mikhail on 11.03.16.
//  Copyright © 2016 Kozlov Mikhail. All rights reserved.
//

#ifndef functions_h
#define functions_h

#include <stdlib.h>
#include <math.h>

double norm1(double *y2, double *y1){
    double max = fabs(y2[0]-y1[0]);
    
    for (int i = 1; i <= 10; i++){
        if (fabs(y2[i]-y1[i]) > max){
            max = fabs(y2[i]-y1[i]);
        }
    }
    
    return max;
}

double norm2(double *y2, double *y1){
    double sum = 0;
    
    for (int i = 1; i <= 10; i++){
        sum += fabs(y2[i]-y1[i]);
    }
    
    return sum;
}

double norm3(double *y2, double *y1){
    double sum = 0;
    
    for (int i = 1; i <= 10; i++){
        sum += pow(y2[i]-y1[i], 2);
    }
    
    return pow(sum, 0.5);
}

inline double k1_x(double x){ //x<x0
    return pow(x, 2)+0.5;
}

inline double k2_x(double x){ //x>x0
    return pow(x, 2)+0.5;
}

inline double q1_x(double x){ //x<x0
    return exp(-pow(x, 2));
}

inline double q2_x(double x){ //x>x0
    return 1;
}

inline double f1_x(double x){ //x<x0
    return cos(x);
}

inline double f2_x(double x){ //x>x0
    return 1;
}

inline double a1_x(double x, double h){
    return k1_x(x+h/2);
}

inline double b1_x(double x, double h){
    return -(k1_x(x+h/2)+k1_x(x-h/2)+q1_x(x)*pow(h, 2));
}

inline double c1_x(double x, double h){
    return k1_x(x-h/2);
}

inline double d1_x(double x, double h){
    return -f1_x(x)*pow(h, 2);
}

inline double a2_x(double x, double h){
    return k2_x(x+h/2);
}

inline double b2_x(double x, double h){
    return -(k2_x(x+h/2)+k2_x(x-h/2)+q2_x(x)*pow(h, 2));
}

inline double c2_x(double x, double h){
    return k2_x(x-h/2);
}

inline double d2_x(double x, double h){
    return -f2_x(x)*pow(h, 2);
}

inline double a1(double x0){
    return k1_x(x0);
}

inline double b1(double x0, double h){
    return -2*k1_x(x0)-q1_x(x0)*pow(h, 2);
}

inline double c1(double x0){
    return k1_x(x0);
}

inline double d1(double x0, double h){
    return -f1_x(x0)*pow(h, 2);
}

inline double a2(double x0){
    return k2_x(x0);
}

inline double b2(double x0, double h){
    return -2*k2_x(x0)-q2_x(x0)*pow(h, 2);
}

inline double c2(double x0){
    return k2_x(x0);
}

inline double d2(double x0, double h){
    return -f2_x(x0)*pow(h, 2);
}

inline double lambda1(double x0){
    return sqrt(q1_x(x0)/k1_x(x0));
}

inline double mu1(double x0){
    return f1_x(x0)/q1_x(x0);
}

inline double lambda2(double x0){
    return sqrt(q2_x(x0)/k2_x(x0));
}

inline double mu2(double x0){
    return f2_x(x0)/q2_x(x0);
}

void tma_mod(double *u, double h, double x0, int n, double xInit){
    double alpha[10*n+1];
    double beta[10*n+1];
    double uBuf;
    int la = 0, lb = 0;
    
    alpha[1] = -a1(x0)/b1(x0, h);
    beta[1] = (d1(x0, h)-c1(x0)*u[0])/b1(x0, h);
    alpha[10*n-1] = -c2(x0)/b2(x0, h);
    beta[10*n-1] = (d2(x0, h)-c2(x0)*u[10])/b2(x0, h);
    
    for (int l = 1; xInit+l*h < x0; l++){
        la += 1;
    }
    
    lb = la+1;
    
    for (int l = 2; l < la; l++){
        alpha[l] = -a1(x0)/(b1(x0, h)+c1(x0)*alpha[l-1]);
        beta[l] = (d1(x0, h)-c1(x0)*beta[l-1])/(b1(x0, h)+c1(x0)*alpha[l-1]);
    }
    
    for (int l = 10*n-2; l > lb; l--){
        alpha[l] = -c2(x0)/(b2(x0, h)+a2(x0)*alpha[l+1]);
        beta[l] = (d2(x0, h)-a2(x0)*beta[l+1])/(b2(x0, h)+c2(x0)*alpha[l+1]);
    }
    
    uBuf = (k1_x(x0)*beta[la-1]+k2_x(x0)*beta[lb+1])/(k1_x(x0)*(1-alpha[la-1])+k2_x(x0)*(1-alpha[lb+1])); //u(la)
    
    for (int l = la; l >= 2; l--){
        uBuf = alpha[l-1]*uBuf+beta[l-1];
        
        if ((l-1)%n == 0){
            u[(l-1)/n] = uBuf;
        }
    }
    
    uBuf = (k1_x(x0)*beta[la-1]+k2_x(x0)*beta[lb+1])/(k1_x(x0)*(1-alpha[la-1])+k2_x(x0)*(1-alpha[lb+1])); //u(lb)
    
    for (int l = lb; l <= 10*n-2; l++){
        uBuf = alpha[l+1]*uBuf+beta[l+1];
        
        if ((l+1)%n == 0){
            u[(l+1)/n] = uBuf;
        }
    }
}

void tma(double *u, double h, double x0, int n, double xInit){
    double alpha[10*n+1];
    double beta[10*n+1];
    double uBuf;
    int la = 0, lb = 0;
    
    alpha[1] = -a1_x(xInit+h, h)/b1_x(xInit+h, h);
    beta[1] = (d1_x(xInit+h, h)-c1_x(xInit+h, h)*u[0])/b1_x(xInit+h, h);
    alpha[10*n-1] = -c2_x(xInit+h*(10*n-1), h)/b2_x(xInit+h*(10*n-1), h);
    beta[10*n-1] = (d2_x(xInit+h*(10*n-1), h)-c2_x(xInit+h*(10*n-1), h)*u[10])/b2_x(xInit+h*(10*n-1), h);
    
    for (int l = 1; xInit+l*h < x0; l++){
        la += 1;
    }
    
    lb = la+1;
    
    for (int l = 2; l < la; l++){
        alpha[l] = -a1_x(xInit+h*l, h)/(b1_x(xInit+h*l, h)+c1_x(xInit+h*l, h)*alpha[l-1]);
        beta[l] = (d1_x(xInit+h*l, h)-c1_x(xInit+h*l, h)*beta[l-1])/(b1_x(xInit+h*l, h)+c1_x(xInit+h*l, h)*alpha[l-1]);
    }
    
    for (int l = 10*n-2; l > lb; l--){
        alpha[l] = -c2_x(xInit+h*l, h)/(b2_x(xInit+h*l, h)+a2_x(xInit+h*l, h)*alpha[l+1]);
        beta[l] = (d2_x(xInit+h*l, h)-a2_x(xInit+h*l, h)*beta[l+1])/(b2_x(xInit+h*l, h)+c2_x(xInit+h*l, h)*alpha[l+1]);
    }
    
    uBuf = (k1_x(xInit+h*la)*beta[la-1]+k2_x(xInit+h*lb)*beta[lb+1])/(k1_x(xInit+h*la)*(1-alpha[la-1])+k2_x(xInit+h*lb)*(1-alpha[lb+1])); //u(la)
    
    for (int l = la; l >= 2; l--){
        uBuf = alpha[l-1]*uBuf+beta[l-1];
        
        if ((l-1)%n == 0){
            u[(l-1)/n] = uBuf;
        }
    }
    
    uBuf = (k1_x(xInit+h*la)*beta[la-1]+k2_x(xInit+h*lb)*beta[lb+1])/(k1_x(xInit+h*la)*(1-alpha[la-1])+k2_x(xInit+h*lb)*(1-alpha[lb+1])); //u(lb)
    
    for (int l = lb; l <= 10*n-2; l++){
        uBuf = alpha[l+1]*uBuf+beta[l+1];
        
        if ((l+1)%n == 0){
            u[(l+1)/n] = uBuf;
        }
    }
}


#endif /* functions_h */
