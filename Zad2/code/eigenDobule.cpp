#include "gauss.h"
void eigenDouble(){
    int i, t, x, y;
    double num, denum;
    std::ifstream fp("macierze.txt");
    ofstream abcXCzas;
    abcXCzas.open ("DEigenCzasyABCX.txt");
    abcXCzas.precision(PRECD);

    ofstream aXCzas;
    aXCzas.open ("DEigenCzasyAX.txt");
    aXCzas.precision(PRECD);

    ofstream abcCzas;
    abcCzas.open ("DEigenCzasyABC.txt");
    abcCzas.precision(PRECD);

    ofstream abcXWyn;
    abcXWyn.open ("DEigenWynikiABCX.txt");
    abcXWyn.precision(PRECD);

    ofstream aXWyn;
    aXWyn.open ("DEigenWynAX.txt");
    aXWyn.precision(PRECD);

    ofstream abcWyn;
    abcWyn.open ("DEigenWynABC.txt");
    abcWyn.precision(PRECD);

    ofstream fullGaussWyn;
    fullGaussWyn.open ("DEigenfullGaussWyn.txt");
    fullGaussWyn.precision(PRECD);

    ofstream partGaussWyn;
    partGaussWyn.open ("DEigenpartGaussWyn.txt");
    partGaussWyn.precision(PRECD);

    ofstream fullGaussCzas;
    fullGaussCzas.open ("DEigenfullGaussCzas.txt");
    fullGaussCzas.precision(PRECD);

    ofstream partGaussCzas;
    partGaussCzas.open ("DEigenpartGaussCzas.txt");
    partGaussCzas.precision(PRECD);

        for(i=10;i<MAX;i = i+ 10){
        MatrixXd emA(i,i);
        MatrixXd emB(i,i);
        MatrixXd emC(i,i);
        VectorXd evV(i);
        VectorXd evW(i);
        MatrixXd emW(i,i);
        MatrixXd emW2(i,i);

        for ( x=0; x<i;x++){
            evW(x) = 0.0;
            for(y=0; y<i;y++){
                emW(x,y) = 0.0;
                emW2(x,y) = 0.0;
                fp >> num;
                fp >> denum;
                double r1 = num / denum;
                emA(x,y) = r1;
                fp >> num;
                fp >> denum;
                r1 = num / denum;
                emB(x,y) = r1;
                fp >> num;
                fp >> denum;
                r1 = num / denum;
                emC(x,y) = r1;
            }
        }

        for (t=0; t<i; t++) { // wektor
            double num;
            fp >> num;
            double denum;
            fp >> denum;
            double wyn = num / denum;
            evV(t) = wyn;
        }

///////////////// A * X /////////////////
        aXCzas << i << endl;
        auto start1 = std::chrono::high_resolution_clock::now();

        evW = emA * evV;

        auto end1 = std::chrono::high_resolution_clock::now();
        std::chrono::duration<double> time1;

        time1 = std::chrono::duration_cast<std::chrono::duration<double>>(end1 - start1);
        aXCzas << time1.count() << endl;

        aXWyn << "i: "<< i << endl;
        aXWyn << evW << endl;
///////////////////////////(A + B + C) * X//////////////////////
        abcXCzas << i << endl;
        auto start2 = std::chrono::high_resolution_clock::now();

        evW = (emA + emB + emC) * evV;

        auto end2 = std::chrono::high_resolution_clock::now();
        std::chrono::duration<double> time2;

        time2 = std::chrono::duration_cast<std::chrono::duration<double>>(end2 - start2);
        abcXCzas << time2.count() << endl;

        abcXWyn << "i: "<< i << endl;
        abcXWyn << evW << endl;
///////////////////////////A *(B*C)//////////////////////
        abcCzas << i << endl;
        auto start3 = std::chrono::high_resolution_clock::now();

        emW = emB * emC;
        emW2 = emA * emW;

        auto end3 = std::chrono::high_resolution_clock::now();
        std::chrono::duration<double> time3;

        time3 = std::chrono::duration_cast<std::chrono::duration<double>>(end3 - start3);
        abcCzas << time3.count() << endl;

        abcWyn << "i: "<< i << endl;
        abcWyn << emW2 << endl;
/////////////////////////////////////////////////////////

                    // fullPivLu
        fullGaussCzas << i << endl;
        auto start4 = std::chrono::high_resolution_clock::now();

        VectorXd fullGauss = emA.fullPivLu().solve(evV);

        auto end4 = std::chrono::high_resolution_clock::now();
        std::chrono::duration<double> time4;
        time4 = std::chrono::duration_cast<std::chrono::duration<double>>(end4 - start4);
        fullGaussCzas << time4.count() << endl;

        fullGaussWyn << i << endl;
        fullGaussWyn << fullGauss << endl;

            // PartialPivLu
        auto start5 = std::chrono::high_resolution_clock::now();

        VectorXd partialGauss = emA.partialPivLu().solve(evV);

        auto end5 = std::chrono::high_resolution_clock::now();
        std::chrono::duration<double> time5;
        time5 = std::chrono::duration_cast<std::chrono::duration<double>>(end5 - start5);
        partGaussCzas << time5.count() << endl;

        partGaussWyn << i << endl;
        partGaussWyn << partialGauss << endl;

    }
    abcXCzas.close();
    aXCzas.close();
    abcCzas.close();

    abcXWyn.close();
    aXWyn.close();
    abcWyn.close();

    fullGaussWyn.close();
    partGaussWyn.close();

    fullGaussCzas.close();
    partGaussCzas.close();
}