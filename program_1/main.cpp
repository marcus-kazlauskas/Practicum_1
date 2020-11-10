//
//  main.cpp
//  program_1 Работа №4, Вариант №3, Задание №10
//
//  Created by Kozlov Mikhail on 26.02.16.
//  Copyright © 2016 Kozlov Mikhail. All rights reserved.
//

#include <iostream>
#include <iomanip>

#include "functions.h"

using namespace std;

int main(int argc, const char * argv[]) {
    // insert code here...
    
    double xInit = 0;                                   // начало отрезка
    double xEnd = 1;                                    // конец отрезка
    double x0 = 1/sqrt(2);                              // точка разрыва
    double e = 0.0001;                                  // относительная погрешность
    double u_mod[11], u1_mod[11], u2_mod[11];
    double u1[11], u2[11];
    double h = (xEnd-xInit)/10;                         // начальное значение шага вычислений (отрезок делится на 10 частей, чтобы посчитать значения ф-ции в 11 точках)
    double hConst = h;                                  // шаг для модельной задачи с постоянными коэффициентами
    u_mod[0] = 0; u1_mod[0] = 0; u2_mod[0] = 0;
    u_mod[10] = 1; u1_mod[10] = 1; u2_mod[10] = 1;
    u1[0] = 0; u2[0] = 0;
    u1[10] = 1; u2[10] = 1;
    int n = 1;                                          // отрезок при вычислениях делится на 10*n частей
//    cout.setf(ios::scientific);
    
    for (int l = 1; l < 10; l++){
        u1_mod[l] = 0;
        u2_mod[l] = 1;
        u1[l] = 0;
        u2[l] = 1;
    }
    
    while ((norm1(u2_mod, u1_mod)/1) > e){       // численное решение модельной задачи; магическое число 1 = xEnd-xInit
        tma_mod(u2_mod, h, x0, n, xInit);
        tma_mod(u1_mod, h/2, x0, n*2, xInit);
        
        h /= 2;
        n *= 2;
    }
    
    double A11, A12, A21, A22;                          // аналитическое решение модельной задачи
    double B1, B2;
    double C1, C2, C3, C4;
    
    A11 = exp(-lambda1(x0)*x0)-exp(lambda1(x0)*x0);
    A12 = exp(lambda2(x0)*(2-x0))-exp(lambda2(x0)*x0);
    A21 = k1_x(x0)*lambda1(x0)*(exp(lambda1(x0)*x0)+exp(-lambda1(x0)*x0));
    A22 = k2_x(x0)*lambda2(x0)*(exp(lambda2(x0)*(2-x0))+exp(lambda2(x0)*x0));
    B1 = mu2(x0)-mu1(x0)+(mu1(x0)-u_mod[0])*exp(lambda1(x0)*x0)-(mu2(x0)-u_mod[10])*exp(lambda2(x0)*(1-x0));
    B2 = k1_x(x0)*lambda1(x0)*(u_mod[0]-mu1(x0))*exp(lambda1(x0)*x0)+k2_x(x0)*lambda2(x0)*(u_mod[10]-mu2(x0))*exp(lambda2(x0)*(1-x0));
    C1 = (((u_mod[0]-mu1(x0))*A11-B1)*A22-((u_mod[0]-mu1(x0))*A21-B2)*A12)/(A11*A22-A12*A21);               // C1, C2 - коэф. общего решения
    C2 = (B1*A22-B2*A12)/(A11*A22-A12*A21);                                                                 // C4, C3 - коэф. общего решения
    C3 = (B2*A11-B1*A21)/(A11*A22-A12*A21);                                                                 // C1..C4 вычисляются по ф-ле (7) со стр. 68
    C4 = (u_mod[10]-mu2(x0))*exp(lambda2(x0))-C3*exp(2*lambda2(x0));
    
    for (int l = 1; l <= 7; l++){   // 1/sqrt(2) ~= 0.7071 > 7
        u_mod[l] = C1*exp(lambda1(x0)*(xInit+hConst*l))+C2*exp(-lambda1(x0)*(xInit+hConst*l))+mu1(x0);      // аналитическое реш. модельной зад. при x < x0
    }                                                                                                       // u = u_общее + u_частное
    
    for (int l = 9; l >= 8; l--){   // 1/sqrt(2) ~= 0.7071 < 8
        u_mod[l] = C3*exp(lambda2(x0)*(xInit+hConst*l))+C4*exp(-lambda2(x0)*(xInit+hConst*l))+mu2(x0);      // аналитическое реш. модельной зад. при x > x0
    }                                                                                                       // u = u_общее + u_частное
    
    cout << "Model task solution:" << endl;
    
    for (int i = 0; i <= 10; i++){
        cout << setprecision(6) << u2_mod[i] << " " << u1_mod[i] << " " << u_mod[i] << endl;
    }
    
    cout << h << endl;
    h = hConst;
    n = 1;
    
    while ((norm1(u2, u1)/1) > e){               // численное решение основной задачи; магическое число 1 = xEnd-xInit
        tma(u2, h, x0, n, xInit);
        tma(u1, h/2, x0, n*2, xInit);
        
        h /= 2;
        n *= 2;
    }
    
    cout << "Solution of the task with variable coefficients:" << endl;
    
    for (int i = 0; i <= 10; i++){
        cout << u2[i] << " " << u1[i] << endl;
    }
    
    cout << h << endl;
    
    return 0;
}
