#include "gauss.h"
void matrixNorm(){
    int x,t,y;
    std::ifstream fMf("FMatrixfullGaussWyn.txt");
    std::ifstream fMp("FMatrixpartGaussWyn.txt");
    std::ifstream fMn("FMatrixnoneGaussWyn.txt");
  /*std::ifstream dMf("DMatrixfullGaussWyn.txt");
    std::ifstream dMp("DMatrixpartGaussWyn.txt");
    std::ifstream dMn("DMatrixnoneGaussWyn.txt");*/

    std::ifstream fEf("FEigenfullGaussWyn.txt");
    std::ifstream fEp("FEigenpartGaussWyn.txt");
  //std::ifstream dEf("DEigenfullGaussWyn.txt");
  //std::ifstream dEp("DEigenpartGaussWyn.txt");


    ofstream FfullGaussNorm;
    FfullGaussNorm.open ("FfullGaussNorm.txt");
    FfullGaussNorm.precision(PRECF);
    ofstream FpartGaussNorm;
    FpartGaussNorm.open ("FpartGaussNorm.txt");
    FpartGaussNorm.precision(PRECF);
    ofstream FnoneGaussNorm;
    FnoneGaussNorm.open ("FnoneGaussNorm.txt");
    FnoneGaussNorm.precision(PRECF);
  /*ofstream DfullGaussNorm;
    DfullGaussNorm.open ("DfullGaussNorm.txt");
    DfullGaussNorm.precision(PRECF);
    ofstream DpartGaussNorm;
    DpartGaussNorm.open ("DpartGaussNorm.txt");
    DpartGaussNorm.precision(PRECF);
    ofstream DnoneGaussNorm;
    DnoneGaussNorm.open ("DnoneGaussNorm.txt");
    DnoneGaussNorm.precision(PRECF);*/


    for(x=10;x<MAX;x=x+10){
        float ef, efF, mf, xff=0.0, xfp=0.0, xfn=0.0;
        double ed, edF, md, xdf=0.0, xdp=0.0, xdn=0.0;
        for(y=0;y<x;y++){
            fMf >> t; fMp >> t; fMn >> t;
          //dMf >> t; dMp >> t; dMn >> t; 

            fEf >> t; fEp >> t;
          //dEf >> t; dEp >> t;

            fMf >> mf;
            fEf >> efF;
            xff += pow((mf-efF), 2);

            fMp >> mf;
            fEp >> ef;
            xfp += pow((mf-ef), 2);

            fMn >> mf;
            xfn += pow((mf-efF), 2);
          /*dMf >> md;
            dEf >> edF;
            xdf += pow((md-edF),2);

            dMp >> md;
            dEp >> ed;
            xdp += pow((md-ed),2);

            dMn >> md;
            xdn += pow((md-edF), 2);*/
        }
        FfullGaussNorm << sqrt(xff) << endl;
        FpartGaussNorm << sqrt(xfp) << endl;
        FnoneGaussNorm << sqrt(xfn) << endl;
      /*DfullGaussNorm << sqrt(xdf) << endl;
        DpartGaussNorm << sqrt(xdp) << endl;
        DnoneGaussNorm << sqrt(xdn) << endl;*/
    }
    FfullGaussNorm.close();
    FpartGaussNorm.close();
    FnoneGaussNorm.close();
  /*DfullGaussNorm.close();
    DpartGaussNorm.close();
    DnoneGaussNorm.close();*/
}