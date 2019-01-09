#include "gauss.h"
void ulamek(){
    int i, t, x, y, h=0, m;  // liczniki
    std::ifstream fp("macierze.txt");
    ofstream abcXCzas;
    abcXCzas.open ("FactorCzasyABCX.txt");
    abcXCzas.precision(PRECD);

    ofstream aXCzas;
    aXCzas.open ("FactorCzasyAX.txt");
    aXCzas.precision(PRECD);

    ofstream abcCzas;
    abcCzas.open ("FactorCzasyABC.txt");
    abcCzas.precision(PRECD);

    ofstream abcXWyn;
    abcXWyn.open ("FactorWynikiABCX.txt");
    abcXWyn.precision(PRECD);

    ofstream aXWyn;
    aXWyn.open ("FactorWynAX.txt");
    aXWyn.precision(PRECD);

    ofstream abcWyn;
    abcWyn.open ("FactorWynABC.txt");
    abcWyn.precision(PRECD);


    for(i=10;i<MAXU;i = i+ 10){
        Fraction f2dArrayA[i][i];
        Fraction f2dArrayB[i][i];
        Fraction f2dArrayC[i][i];
        Fraction fVector[i][1];
        Fraction fWyn[i][i];
        Fraction fWyn2[i][i];
        Fraction fWyn3[i][i];
        Fraction fWynVector[i][1];
        Fraction fWynVector2[i][1];
        Fraction fWynAbc[i][i];

        for ( x=0; x<i;x++){
            fWynVector[x][0].set_num(0);
            fWynVector[x][0].set_denum(1);
            fWynVector2[x][0].set_num(0);
            fWynVector2[x][0].set_denum(1);
            for(y=0; y<i;y++){
                fWyn3[x][y].set_num(0);
                fWyn3[x][y].set_denum(1);
                fWyn2[x][y].set_num(0);
                fWyn2[x][y].set_denum(1);
                fWyn[x][y].set_num(0);
                fWyn[x][y].set_denum(1);
                fWynAbc[x][y].set_num(0);
                fWynAbc[x][y].set_denum(1);
                unsigned  long long r1;
                fp >> r1;
                unsigned  long long r2;
                fp >> r2;
                f2dArrayA[x][y].set_num(r1);
                f2dArrayA[x][y].set_denum(r2);
                fp >> r1;
                fp >> r2;
                f2dArrayB[x][y].set_num(r1);
                f2dArrayB[x][y].set_denum(r2);
                fp >> r1;
                fp >> r2;
                f2dArrayC[x][y].set_num(r1);
                f2dArrayC[x][y].set_denum(r2);
            }
        }

        for (t=0; t<i; t++) { // wektor
            unsigned long long r1;
            fp >> r1;
            unsigned long long r2;
            fp >> r2;
            fVector[t][0].set_num(r1);
            fVector[t][0].set_denum(r2);
        }

///////////////////////////(A + B + C) * X//////////////////////
        abcXCzas << i << endl;
        auto start1 = std::chrono::high_resolution_clock::now();
        for ( x=0; x<i;x++){
            for(y=0; y<i;y++){
                fWyn[x][y] = (f2dArrayA[x][y]+f2dArrayB[x][y]+f2dArrayC[x][y])*fVector[h][0];
                h++;
                fWynVector[x][0]=fWynVector[x][0] + fWyn[x][y];
            }
            h=0;
        }
        auto end1 = std::chrono::high_resolution_clock::now();
        std::chrono::duration<double> time1;

        time1 = std::chrono::duration_cast<std::chrono::duration<double>>(end1 - start1);
        abcXCzas << time1.count() << endl;

        abcXWyn << "i: "<< i << endl;
        for ( x=0; x<i;x++){
            double out = (double) fWynVector[x][0].get_num() / (double) fWynVector[x][0].get_denum();
            abcXWyn << out << endl;
        }


///////////////////////////A * X//////////////////////
        aXCzas << i << endl;
        auto start2 = std::chrono::high_resolution_clock::now();
        for ( x=0; x<i;x++){
            for(y=0; y<i;y++){
                fWyn2[x][y] = f2dArrayA[x][y]*fVector[h][0];
                h++;
                fWynVector2[x][0]=fWynVector2[x][0]+fWyn2[x][y];
            }
            h=0;
        }
        auto end2 = std::chrono::high_resolution_clock::now();
        double res1Time = std::chrono::duration_cast<std::chrono::nanoseconds>(end2 - start2).count();
        //res1Time -= 1000000;
        //res1Time /= 1000000;

        aXCzas << res1Time<< endl;

        aXWyn << "i: "<< i << endl;
        for ( x=0; x<i;x++){
            double out = (double) fWynVector2[x][0].get_num() / (double) fWynVector2[x][0].get_denum();
            aXWyn << out << endl;
        }
///////////////////////////A *(B*C)//////////////////////
        abcCzas << i << endl;
        auto start3 = std::chrono::high_resolution_clock::now();

        for ( x=0; x<i;x++)
            for(y=0; y<i;y++)
                for(m=0;m<i;m++)
                    fWyn3[x][y] = fWyn3[x][y] + (f2dArrayB[x][m] * f2dArrayC[m][y]);

         for ( x=0; x<i;x++)
            for(y=0; y<i;y++)
                for(m=0;m<i;m++){
                    fWynAbc[x][y] = fWynAbc[x][y] + (f2dArrayA[x][m] * fWyn3[m][y]);
                }
        auto end3 = std::chrono::high_resolution_clock::now();
        std::chrono::duration<double> time3;

        time3 = std::chrono::duration_cast<std::chrono::duration<double>>(end3 - start3);
        abcCzas << time3.count() << endl;

        abcWyn << "i: "<< i << endl;
        for ( x=0; x<i;x++){
            for(y=0; y<i;y++){
                    double out = (double) fWynAbc[x][y].get_num() / (double) fWynAbc[x][y].get_denum();
                    abcWyn << out;
                    abcWyn << " ";
                }
                abcWyn << endl;
            }
///////////////////////////////////////////////////////////

    }
    abcXCzas.close();
    aXCzas.close();
    abcCzas.close();

    abcXWyn.close();
    aXWyn.close();
    abcWyn.close();
}