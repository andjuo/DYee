/*

Minnesota BLUE Results Combiner
Version 1.0
May 10, 2011

J. Pastika and J. Mans

This tar contains version 1.0 of the UMN BLUE results combiner.

This will combine the results of 2 independant measurements each with
their own covariance matrix.  The method is summarized in [1].  The
method is implemented as a root script which accepts the input files
of the seperate measurements as text files. 

[1] Valassi, A. "Combining correlated measurements of several
different physical quantities."  Nuclear instruments & methods in
physics research. Section A, Accelerators, spectrometers, detectors
and associated equipment 500.1-3 (2003):391

 */

#include "TMatrixT.h"
#include "TVectorT.h"
#include <stdio.h>
#include <math.h>
#include <cmath>
#ifndef VECTOR
#define VECTOR
#include <vector>
#endif

#include <iostream>

void readDataVector(const char * name, TVectorD& dv, double binLimits[][2] = 0, int ftr = 1, double* ccCovar = NULL);
void readCovMatrix(const char * name, TMatrixT<double>& dv);
void simpleMatrixDivider(TMatrixT<double>& mat, double factor, bool multiply);

//const int theNBins = 2;
//Double_t mBins[theNBins+1] = {15, 20,25};
const int theNBins = 40;
Double_t mBins[theNBins+1] = {15, 20, 25, 30, 35, 40, 45, 50, 55, 60, 64, 68, 72, 76, 81, 86, 91, 96, 101, 106, 110, 115, 120, 126, 133, 141, 150, 160, 171, 185, 200, 220, 243, 273, 320, 380, 440, 510, 600, 1000, 1500};

void resultCombiner_2D(const char * name_dv1, const char * name_cm1, const char* slength1, const char * name_dv2, const char * name_cm2, const char* slength2, int ftr = 1, const char * outputNameStub = "output")
{
    const int length1 = 40; //atoi(slength1);
    const int length2 = 40; //atoi(slength2);

    TVectorD dv1(length1), dv2(length2);
    TMatrixT<double> cm1(length1, length1), cm2(length2, length2);
    double binLimits1[length1][2], binLimits2[length2][2], ccCov1[length1], ccCov2[length1];

    readDataVector(name_dv1, dv1, binLimits1, ftr, ccCov1);
    printf("Read data vector 1 (%d)\n",length1);
    //for(int i = 0;i<length1;++i){
      //std::cout<<" i = "<<i<<" dv1[i] = "<<dv1(i)<<std::endl;
    //}

    readDataVector(name_dv2, dv2, binLimits2, ftr, ccCov2);
    printf("Read data vector 2 (%d)\n",length2);

    readCovMatrix(name_cm1, cm1);
    printf("Read covariance matrix 1\n");
    simpleMatrixDivider(cm1, 1000000000000000., true);

    readCovMatrix(name_cm2, cm2);
    printf("Read covariance matrix 2\n");
    simpleMatrixDivider(cm2, 1000000000000000., true);

    std::vector<double*> binLimits;
    std::vector<std::vector<int > > preU;
    int i1 = 0, i2 = 0;
   std::cout<<"00000"<<std::endl;
    while(i1 < length1 || i2 < length2)
    {
        if(i1 < length1 && i2 < length2)
        {
            //std::cout<<" i1 = "<<i1<<" i2 = "<<i2<<" L1[i1][1] = "<<binLimits1[i1][1]<<" L1[i1][0] = "<<binLimits1[i1][0]<<" L2[i2][1] = "<<binLimits2[i2][1]<<" L2[i2][0] = "<<binLimits2[i2][0]<<std::endl;
            if((binLimits1[i1][1] + binLimits1[i1][0])/2 > binLimits2[i2][0] && (binLimits1[i1][1] + binLimits1[i1][0])/2 < binLimits2[i2][1])
            {
                binLimits.push_back(binLimits1[i1]);
                std::vector<int> tmp;
                tmp.push_back(i1);
                tmp.push_back(i2);
                preU.push_back(tmp);
                i1++;
                i2++;
            }
            else if((binLimits1[i1][1] + binLimits1[i1][0])/2 <= binLimits2[i2][0])
            {
                binLimits.push_back(binLimits1[i1]);
                std::vector<int> tmp;
                tmp.push_back(i1);
                tmp.push_back(-1);
                preU.push_back(tmp);
                i1++;
            }
            else
            {
                binLimits.push_back(binLimits2[i2]);
                std::vector<int> tmp;
                tmp.push_back(-1);
                tmp.push_back(i2);
                preU.push_back(tmp);
                i2++;
            }
        }
        else if(i1 < length1 && i2 >= length2)
        {
            binLimits.push_back(binLimits1[i1]);
            std::vector<int> tmp;
            tmp.push_back(i1);
            tmp.push_back(-1);
            preU.push_back(tmp);
            i1++;
        }
        else
        {
            binLimits.push_back(binLimits2[i2]);
            std::vector<int> tmp;
            tmp.push_back(-1);
            tmp.push_back(i2);
            preU.push_back(tmp);
            i2++;
        }
    }
//

   ofstream myfile;
   std::string myName = "test_preU.dat";
   myfile.open (myName.c_str());
    for (int i =0;i<preU.size();++i){
      for (int j =0;j<2;++j){
        myfile<<" i = "<<i<<" j = "<<j<<" preU[i][j] = "<<preU[i][j]<<"\n";
      }
    }

    myfile.close();
    std::cout<<" tests written to "<<myName<<std::endl;

//

    TVectorD dv(length1 + length2);
    for(int i = 0; i < length1 + length2; i++)
    {
        dv[i] = (i < length1) ? dv1[i] : dv2[i - length1];
    }

    TMatrixT<double> cm(length1 + length2, length1 + length2), U(length1 + length2, preU.size());
    for(int i = 0; i < length1; i++)
    {
        for(int j = 0; j < length1; j++)
        {
            cm[i][j] = cm1[i][j];
        }
    }
    for(int ii = length1; ii < length1 + length2; ii++)
    {
        for(int j = length1; j < length1 + length2; j++)
        {
            cm[ii][j] = cm2[ii - length1][j - length1];
        }
    }
    
    //for(unsigned int i = 0; i < preU.size(); i++)
    for(int i = 0; i < preU.size(); i++)
    {
        if(preU[i][0] >= 0) U[preU[i][0]][i] = 1;
        if(preU[i][1] >= 0) U[preU[i][1] + length1][i] = 1;
        if(ftr > 1 && preU[i][0] >= 0 && preU[i][1] >= 0)  cm[preU[i][0]][preU[i][1] + length1] = cm[preU[i][1] + length1][preU[i][0]] = ccCov1[preU[i][0]]*ccCov2[preU[i][1]];
    }
    
    //        cm.Print();
    cout << "Covariance matrix" << endl;
    cm.Print();

    std::cout<<" Transpose U"<<std::endl;
    cout << "U matrix" << endl;
    U.Print();
    TMatrixT<double> Ut(U);
    Ut.T();

    std::cout<<" Invert cm"<<std::endl;
    TMatrixT<double> cmInv(cm);
    cmInv.Invert();
    TMatrixT<double> step1 = Ut * cmInv * U;
    TMatrixT<double> step2 = Ut * cmInv;
    std::cout<<" Invert Ut * cmInv * U"<<std::endl;
    TMatrixT<double> lambda = step1.Invert() * step2;
    cout << "dv1 initial vector" << endl;
    dv1.Print();
    cout << "dv2 initial vector" << endl;
    dv2.Print();
    cout << "Lambda multipliers" << endl;
    lambda.Print();
    TVectorD bV = lambda*dv;
    cout << "bV final result vector" << endl;
    bV.Print();
    //simpleDataVectorDivider(bV,100000000.);
    std::cout<<" Invert Ut * cmInv * U"<<std::endl;
    TMatrixT<double> bcm = (Ut * cmInv * U).Invert();
    cout << "bcm final covariance matrix" << endl;
    simpleMatrixDivider(bcm, 1000000000000000., false);
    bcm.Print();

    bV.Draw();
    //dv.Draw();
    gStyle->SetMarkerStyle(26);
    dv1.Draw("samep");
    gStyle->SetMarkerStyle(28);
    dv2.Draw("samep");
    //lambda.Draw("colz");
    bcm.Draw("colz");
    //U.Draw("colz");
    printf("Done with combination.\n");


    //FIXME addition: look for combination
    cout << "Doing chi2" << endl;
    double largeNumber = 100000000.;
    TMatrixT<double> cov_invert(bcm);
    cov_invert*=largeNumber;
    cov_invert.Invert();
    cov_invert*=largeNumber;
    TVectorT<double> data_TV_tmp(dv1 - dv2);
    data_TV_tmp *= cov_invert;
    double chi2 = data_TV_tmp*(dv1 - dv2);
    cout << "XXX Chi2 " << chi2/40. << endl;
/*
    // write the provided input in readable format (like the output)
    simpleMatrixDivider(cm1, 1000000000000000., false);
    simpleMatrixDivider(cm2, 1000000000000000., false);
    FILE *file_input;
    char bVinName[128], CMinName[128];
    sprintf(bVinName, "input_data1.txt");

    file_input = fopen(bVinName, "w");
    if(file_input)
    {
        fprintf(file_input, "#\n#%9s %9s %9s %15s %15s %15s\n", "Bin", "Y_min", "Y_max", "Value", "Uncertainty", "Uncertainty (%)");
        for(int i = 0; i < dv1.GetNoElements(); i++)
        {
            fprintf(file_input, " %9i %9.2f %9.2f %15e %15e     %4.2f \n", i + 1, binLimits[i][0], binLimits[i][1], dv1[i], sqrt(cm1[i][i]), 100.*sqrt(cm1[i][i])/dv1[i]);
        }
        fclose(file_input);
    }

    FILE *file_input;
    //char bVinName[128], CMinName[128];
    sprintf(bVinName, "input_data2.txt");

    file_input = fopen(bVinName, "w");
    if(file_input)
    {
        fprintf(file_input, "#\n#%9s %9s %9s %15s %15s %15s\n", "Bin", "Y_min", "Y_max", "Value", "Uncertainty", "Uncertainty (%)");
        for(int i = 0; i < dv1.GetNoElements(); i++)
        {
            fprintf(file_input, " %9i %9.2f %9.2f %15e %15e     %4.2f \n", i + 1, binLimits[i][0], binLimits[i][1], dv2[i], sqrt(cm2[i][i]), 100.*sqrt(cm2[i][i])/dv2[i]);
        }
        fclose(file_input);
    }

    //write output
    FILE *file;
    char bVoutName[128], CMoutName[128];
    sprintf(bVoutName, "%s_data.txt", outputNameStub);

    file = fopen(bVoutName, "w");
    if(file)
    {
        fprintf(file, "#\n#%9s %9s %9s %15s %15s\n", "Bin", "Y_min", "Y_max", "Value", "Uncertainty");
        for(int i = 0; i < bV.GetNoElements(); i++)
        {
            fprintf(file, " %9i %9.2f %9.2f %15e %15e\n", i + 1, binLimits[i][0], binLimits[i][1], bV[i],sqrt(bcm[i][i]));
        }
        fclose(file);
    }

    sprintf(CMoutName, "%s_covMat.txt", outputNameStub);

    file = fopen(CMoutName, "w");
    if(file)
    {
        fprintf(file, "#\n#%9s %9s %15s\n", "Bin i", "Bin j", "Value");
        for(int i = 0; i < bcm.GetNrows(); i++)
        {
            for(int j = 0; j < bcm.GetNcols(); j++)
            {
                fprintf(file, " %9i %9i %15e\n", i + 1, j + 1, bcm[i][j]);
            }
        }
        fclose(file);
    }
    printf("Output complete.\n");
*/
}

void readDataVector(const char * name, TVectorD& dv, double binLimits[][2], int ftr, double* ccCovar)
{
    //read data vectors
    FILE *file;
    char *c, buff[4096], *k;
    int i;
    double d1, d2, d3, f;

    file = fopen(name, "r");
    if(file)
    {
        while(!feof(file))
        {
            c = fgets(buff, 4096, file);
            std::cout<<" buff = "<<buff<<std::endl;
            for(k = strchr(buff, ','); k != 0; k = strchr(buff, ',')) *k = ' ';

            if(c != NULL && buff[0] != '#')
            {
                if(sscanf(buff, "%d %lf %lf %lf", &i, &d1, &d2, &d3) == 4)
                {
		  
                    if(binLimits != 0)
                    {
                        binLimits[i - 1][0] = d1;
                        binLimits[i - 1][1] = d2;
                    }
                    dv[i - 1] = d3;
                }
                if(ccCovar != NULL && ftr > 1)
                {
                	ccCovar[i] = 0.0;
                	char delim[128];
                	sprintf(delim, "%%*d");
                	for(int j = 0; j <= ftr; j++)
                	{
                		sprintf(delim, "%s %%*f", delim);
                	}
                	sprintf(delim, "%s %%f", delim);
                	if(sscanf(buff, delim, &f) == 1) 
                	{
                		ccCovar[i - 1] = f;
                	}
                }
            }
        }
        fclose(file);
    }

    return;
}

void simpleMatrixDivider(TMatrixT<double>& mat, double factor, bool multiply) { 

     for(int i = 0; i < mat.GetNrows(); i++)
      {
         for(int j = 0; j < mat.GetNcols(); j++)
            {
                if (multiply) {
                     mat[i][j] *= factor;
                } else {
                     mat[i][j] /= factor;
                }
           }
       }
}

void readCovMatrix(const char * name, TMatrixT<double>& cm)
{
    //read covariance matrix
    FILE *file;
    char *k, *c, buff[4096];
    int i = 0, j = 0;
    double d;

    file = fopen(name, "r");
    if(file)
    {
        while(!feof(file))
        {
            c = fgets(buff, 4096, file);
            for(k = strchr(buff, ','); k != 0; k = strchr(buff, ',')) *k = ' ';
            if(c != NULL && buff[0] != '#')
            {
                if(sscanf(buff, "%d %d %lf", &i, &j, &d) == 3)
                {
		  if (i<=cm.GetNrows() && j<=cm.GetNcols())
                    cm[i - 1][j - 1] = d;
                }
            }
        }
        fclose(file);
    }

    return;
}
