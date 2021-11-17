#include "stdafx.h"
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "interpolation.h"

using namespace alglib;
void function_debt_func(const real_1d_array& c, const real_1d_array& x, double& func, void* ptr)
{
    //
    // this callback calculates f(c,x)=c[0]*(1+c[1]*(pow(x[0]-1999,c[2])-1))
    //
    
    func = 2 * (c[0] * (1 - exp(-x[0] / c[1])) + c[2]);// c[3] * x[0] + c[4] * pow(x[0], 2));

    //func = 2 * (c[0] * (1 - exp(-x[0] / c[1])) + c[2]+ c[3] * x[0] + c[4] * pow(x[0], 2));
    //func = 2 * (c[0] * (1 - exp(-x[0] / c[1])) + c[2] + c[3] * pow(x[0], 2));


}
int main(int argc, char** argv)
{
   
    int acquisition_frequency = 10;
    double Dt = 1. / acquisition_frequency;
    int nb_images = 2499;
    int nb_wave_vector = 128;
    double temp = nb_images / acquisition_frequency;
    int max_t_ad = nb_images / acquisition_frequency +1;
    int size = nb_images * nb_wave_vector;
    double* tab_dinamics = new double[size];
   // real_1d_array q1_dinamics="[1 , 5]";// = new double[nb_images];
    double* q1_dinamics  = new double[max_t_ad];
    double* xx_d = new double[max_t_ad];

    FILE* images_file=fopen("C:\\samples\\output\\test_ORCA_20082_dinamics.mat", "rb");
    if (images_file == NULL)
        std::cout << "ERROR" << std::endl;

    uint16_t sz;
    fread(&sz, sizeof(uint16_t), 1, images_file);
    uint16_t dimy;
    fread(&dimy, sizeof(uint16_t), 1, images_file);
    uint16_t dimx;
    fread(&dimx, sizeof(uint16_t), 1, images_file);

   int nb= fread(tab_dinamics, sizeof(double), dimx*dimy, images_file);

   fclose(images_file);

   int max_t_adbis = 150;

   double* q1_dinamicsbis = new double[max_t_adbis];
   double* xxbis = new double[max_t_adbis];
   //for (int j = 30; j < 31; j++) {  //nb_wave_vector

       for (int i = 0; i < 150; i++) { //max_t_ad

           xxbis[i] = i + 1;               //tab_dinamics[j + i * nb_wave_vector] / 4.3563e+03;
           q1_dinamicsbis[i] = 2 * (5.12345 * (1 - exp(-xxbis[i] / 12.352489787)) + 0.05);
           //std::cout << i<<"   "<<q1_dinamics[i] << std::endl;
       }
   //}
   real_1d_array yy;
   yy.setcontent(max_t_adbis, q1_dinamicsbis);

   
   real_2d_array xx;
   for (int i = 0; i < max_t_ad; i++) {
       xx_d[i] = Dt * (i+1);
      // std::cout << i+1<<"   "<<xx_d[i]<<"    "<<q1_dinamics[i] << std::endl;

   }

   xx.setcontent(max_t_adbis, 1, xxbis);

   real_1d_array cc = "[20, 1, 0]";
   //real_1d_array ss = "[2, 10, 204.14, 0.05, 0]";
   real_1d_array bndl_ = "[0, 0, 0]";
   real_1d_array bndu_ = "[100, 50, 5]";

    //real_2d_array x = "[[2000],[2001],[2002],[2003],[2004],[2005],[2006],[2007],[2008]]";
   

    //real_1d_array y = "[4323239600000.0, 4560913100000.0, 5564091500000.0, 6743189300000.0, 7284064600000.0, 7050129600000.0, 7092221500000.0, 8483907600000.0, 8625804400000.0]";
    //real_1d_array c = "[10, 1, 1]";
    double epsx = 1.0e-5;
    //real_1d_array bndl = "[-inf, -10, 0.1]";
    //real_1d_array bndu = "[+inf, +10, 2.0]";
    //real_1d_array s = "[1, 1, 1]";
    ae_int_t maxits = 0;
    ae_int_t info;
    lsfitstate state;
    lsfitreport rep;
    double diffstep = 1.0e-5;

    lsfitcreatef(xx, yy, cc, diffstep, state);
    //lsfitsetcond(state, epsx, maxits);
    lsfitsetbc(state, bndl_, bndu_);
    lsfitsetcond(state, epsx, maxits);
    //lsfitsetscale(state, ss);
    alglib::lsfitfit(state, function_debt_func);
    lsfitresults(state, info, cc, rep);
    //printf("%d\n", int(info)); // EXPECTED: 2
    printf("%s\n", cc.tostring(-2).c_str()); // EXPECTED: [4.142560E+12, 0.434240, 0.565376]
   // double  calc;
   // double yy_;
    //for (int i = 0; i < max_t_ad; i++) {
        //yy_ = 2 * (cc[0] * (1 - exp(-xx[i][0] / cc[1])) + cc[2] + cc[3] * xx[i][0] + cc[4] * pow(xx[i][0], 2));
          //  2 * (cc[0] * (1 - exp(-xx[i][0] / cc[1])) + cc[2] + cc[3] * pow(xx[i][0], 2));
        //calc = x[i][0] -1999;
        //std::cout << exp(10) << std::endl;
        //std::cout << i <<"   "<<xx[i][0]<<"   "<<yy[i]<<"   "<< 2 * (cc[0] * (1 - exp(-xx[i][0] / cc[1])) + cc[2] + cc[3] * xx[i][0] + cc[4] * pow(xx[i][0], 2))<< std::endl;
        //std::cout << i << "   " << yy[i] << "   " << yy_ <<"    "<< abs(yy_-yy[i])<<std::endl;

   // }
    return 0;
}
