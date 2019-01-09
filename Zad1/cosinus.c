#include <stdio.h>
#include <math.h>

    double pot(double x, long n){
        int i;
        double out = x;
        if (n==0) return 1;
        for(i=0;i<n-1;i++)
        {
            out *= x;
        }
        return out;
    }

     long long silnia(long div){
         long long out =1;
         while(div>0){
             out *= div;
             div--;
         }
         return out;
     }

     double absV(double x, double y){
        if(x > y) return x-y;
        else return y-x;
    }

int main(){
        int i,n;
        long double cosMath;
        double x = 2;
        long double cosx = 0;
        cosMath = cos(x);
        double licznik = 0;
        double jedynka = 0;
        unsigned long long mianownik = 0;
        double div = 0;

        for(n = 0; absV(cosx,cosMath) >= 0.000000000000001; n++){
            mianownik = silnia(2*n);
            licznik = pot(x, 2*n);
            jedynka = pot(-1, n);
            div = licznik / mianownik;
            cosx += jedynka * div;
            //printf("%i . %.15Lf\n", n, cosx);
            printf("%.15Lf\n", cosx);
        }



        printf("\n\n%.15Lf\n", cosx);
        printf("%.15Lf", cosMath);

}
