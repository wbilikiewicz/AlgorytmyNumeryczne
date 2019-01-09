#include "gauss.h"
template<typename T>
void matrixFD(){
    int i, t, x, y, m, v;  // liczniki
    T num, denum;
    std::ifstream fp("macierze.txt");

    ofstream abcXCzas;
    abcXCzas.open ("FMatrixCzasyABCX.txt");//D
    abcXCzas.precision(PRECD);

    ofstream aXCzas;
    aXCzas.open ("FMatrixCzasyAX.txt");//D
    aXCzas.precision(PRECD);

    ofstream abcCzas;
    abcCzas.open ("FMatrixCzasyABC.txt");//D
    abcCzas.precision(PRECD);

    ofstream abcXWyn;
    abcXWyn.open ("FMatrixWynikiABCX.txt");//D
    abcXWyn.precision(PRECF);//

    ofstream aXWyn;
    aXWyn.open ("FMatrixWynAX.txt");//D
    aXWyn.precision(PRECF);//

    ofstream abcWyn;
    abcWyn.open ("FMatrixWynABC.txt");//D
    abcWyn.precision(PRECF);//

    ofstream fullGaussWyn;
    fullGaussWyn.open ("FMatrixfullGaussWyn.txt");//D
    fullGaussWyn.precision(PRECF);//

    ofstream partGaussWyn;
    partGaussWyn.open ("FMatrixpartGaussWyn.txt");//D
    partGaussWyn.precision(PRECF);//

    ofstream noneGaussWyn;
    noneGaussWyn.open ("FMatrixnoneGaussWyn.txt");//D
    noneGaussWyn.precision(PRECF);//

    ofstream fullGaussCzas;
    fullGaussCzas.open ("FMatrixfullGaussCzas.txt");//D
    fullGaussCzas.precision(PRECD);

    ofstream partGaussCzas;
    partGaussCzas.open ("FMatrixpartGaussCzas.txt");//D
    partGaussCzas.precision(PRECD);

    ofstream noneGaussCzas;
    noneGaussCzas.open ("FMatrixnoneGaussCzas.txt");//D
    noneGaussCzas.precision(PRECD);

    for(i=10;i<MAX;i = i+ 10){
        MyMatrix<T> matrixA(i, i, 0.0);
        MyMatrix<T> matrixB(i, i, 0.0);
        MyMatrix<T> matrixC(i, i, 0.0);
        vector<T> mVector(i);
        vector<T> mVectorWyn(i);
        vector<T> mVectorWyn2(i);
        MyMatrix<T> mWyn(i, i, 0.0);
        MyMatrix<T> mWyn2(i, i, 0.0);
        MyMatrix<T> mWyn3(i, i, 0.0);
        MyMatrix<T> matrixAV(i, i+1, 1.0);
        MyMatrix<T> matrixAV2(i, i+1, 1.0);
        MyMatrix<T> matrixAV3(i, i+1, 1.0);

        for ( x=0; x<i;x++){

            for(y=0; y<i;y++){
                fp >> num;
                fp >> denum;
                T r1 = num / denum;
                matrixA.matrix[x][y] = r1;
                matrixAV.matrix[x][y] = r1;
                matrixAV2.matrix[x][y] = r1;
                matrixAV3.matrix[x][y] = r1;
                fp >> num;
                fp >> denum;
                r1 = num / denum;
                matrixB.matrix[x][y] = r1;
                fp >> num;
                fp >> denum;
                r1 = num / denum;
                matrixC.matrix[x][y] = r1;
            }
            mVector[x] = 0;
            mVectorWyn[x] = 0;
            mVectorWyn2[x] = 0;
        }
        for (t=0; t<i; t++) { // wektor
            T num;
            fp >> num;
            T denum;
            fp >> denum;
            T wyn = num / denum;
            mVector[t] = wyn;
        }
        for(v =0; v<matrixAV.getRowCount(); v++){
                matrixAV(v,matrixAV.getColCount()-1)=mVector[v];
                matrixAV2(v,matrixAV2.getColCount()-1)=mVector[v];
                matrixAV3(v,matrixAV3.getColCount()-1)=mVector[v];
        }

//////////// A * X /////////////
        aXCzas << i << endl;
        auto start1 = std::chrono::high_resolution_clock::now();

        mVectorWyn = matrixA *  mVector;

        auto end1 = std::chrono::high_resolution_clock::now();
        std::chrono::duration<double> time1;
        time1 = std::chrono::duration_cast<std::chrono::duration<double>>(end1 - start1);
        aXCzas << time1.count() << endl;

        aXWyn << "i: "<< i << endl;
        for(x=0; x<i;x++){
            aXWyn << mVectorWyn[x] << endl;
        }
///////////////////////////(A + B + C) * X//////////////////////
        abcXCzas << i << endl;
        auto start2 = std::chrono::high_resolution_clock::now();

        for ( x=0; x<i;x++)
            for(y=0; y<i;y++){
                 mWyn.matrix[x][y] = matrixA.matrix[x][y] + matrixB.matrix[x][y] + matrixC.matrix[x][y];
                }
            mVectorWyn = mWyn * mVector;

        auto end2 = std::chrono::high_resolution_clock::now();
        std::chrono::duration<double> time2;
        time2 = std::chrono::duration_cast<std::chrono::duration<double>>(end2 - start2);
        abcXCzas << time2.count() << endl;

        abcXWyn << "i: "<< i << endl;
        for(x=0; x<i;x++){
            abcXWyn << mVectorWyn[x] << endl;
        }
///////////////////////////A *(B*C)//////////////////////
        abcCzas << i << endl;
        auto start3 = std::chrono::high_resolution_clock::now();

              for ( x=0; x<i;x++)
                for(y=0; y<i;y++)
                  for(m=0;m<i;m++)
                    mWyn3.matrix[x][y] = mWyn3.matrix[x][y] + (matrixB.matrix[x][m] * matrixC.matrix[m][y]);

              for ( x=0; x<i;x++)
                for(y=0; y<i;y++)
                  for(m=0;m<i;m++)
                    mWyn2.matrix[x][y] = mWyn2.matrix[x][y] + (matrixA.matrix[x][m] * mWyn3.matrix[m][y]);

        auto end3 = std::chrono::high_resolution_clock::now();
        std::chrono::duration<double> time3;
        time3 = std::chrono::duration_cast<std::chrono::duration<double>>(end3 - start3);
        abcCzas << time3.count() << endl;

        abcWyn << "i: "<< i << endl;
        for ( x=0; x<i;x++){
            for(y=0; y<i;y++){
                abcWyn << mWyn2.matrix[x][y];
                abcWyn << " ";
            }
            abcWyn << endl;
        }
//////////////////////////////////////////////////////////
        vector<T> gaussFullOut(i);
        vector<T> gaussPartOut(i);
        vector<T> gaussNoneOut(i);
///////////////////////////GAUSS FULL///////////////////////////

        fullGaussCzas << i << endl;
        auto start4 = std::chrono::high_resolution_clock::now();

        gaussFullOut = matrixAV.Gauss1();

        auto end4 = std::chrono::high_resolution_clock::now();
        std::chrono::duration<double> time4;
        time4 = std::chrono::duration_cast<std::chrono::duration<double>>(end4 - start4);
        fullGaussCzas << time4.count() << endl;

        fullGaussWyn << i << endl;
        for ( x=0; x<i;x++)
            fullGaussWyn << gaussFullOut[x] << endl;

/////////////////////////GAUSS PARTIAL/////////////////////////
        partGaussCzas << i << endl;
        auto start5 = std::chrono::high_resolution_clock::now();

        gaussPartOut = matrixAV2.Gauss2();

        auto end5 = std::chrono::high_resolution_clock::now();
        std::chrono::duration<double> time5;
        time5 = std::chrono::duration_cast<std::chrono::duration<double>>(end5 - start5);
        partGaussCzas << time5.count() << endl;

        partGaussWyn << i << endl;
        for ( x=0; x<i;x++)
            partGaussWyn << gaussPartOut[x] << endl;


////////////////////////////GAUSS NONE/////////////////////////////
        noneGaussCzas << i << endl;
        auto start6 = std::chrono::high_resolution_clock::now();

        gaussNoneOut = matrixAV3.Gauss3();

        auto end6 = std::chrono::high_resolution_clock::now();
        std::chrono::duration<double> time6;
        time6 = std::chrono::duration_cast<std::chrono::duration<double>>(end6 - start6);
        noneGaussCzas << time6.count() << endl;

        noneGaussWyn << i << endl;
        for ( x=0; x<i;x++)
            noneGaussWyn << gaussNoneOut[x] << endl;

        }
    abcXCzas.close();
    aXCzas.close();
    abcCzas.close();

    abcXWyn.close();
    aXWyn.close();
    abcWyn.close();

    fullGaussWyn.close();
    noneGaussWyn.close();
    partGaussWyn.close();

    fullGaussCzas.close();
    noneGaussCzas.close();
    partGaussCzas.close();
};