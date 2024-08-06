//This file generates all terms of the spin correlation that is expressed via pseudo-fermion vertices.
//The terms are obtained by expressing spin correlations as in Eq. (43) of Phys. Rev. B 109, 174414.
//Different symmetry scenarios can be considered, namely models with or without time-reversal symmetry or global U(1) spin rotation symmetry.
//Some aspects of this code may become more clearer when considered together with the corresponding PFFRG code.

#include <stdio.h>
#include <math.h>
#include <iostream>
#include <fstream>
//#include <gsl/gsl_errno.h>
//#include <gsl/gsl_matrix.h>
//#include <gsl/gsl_odeiv2.h>
#include <complex>
#include <array>
#include <stdlib.h>
#include <iomanip>
#include <algorithm>
#include <list>
#include<vector>

using namespace std;

complex<int> iu(0,1);

//The $\bm{\beta}_{a}$ matrices are defined by the equation
// \bm{\sigma}^{\mu} \bm{\sigma}^{\nu} = \sum_{a=0,x,y,z} \beta^{\mu\nu}_{a} \bm{\sigma}^{a},
//with $\bm{\sigma}^{\mu}$ being Pauli matrices for $\mu=x,y,z$ and the identity matrix for $\mu=0$.
complex<int> beta[4][4][4] = {{{1,0,0,0},{0,1,0,0},{0,0,1,0},{0,0,0,1}},
                              {{0,1,0,0},{1,0,0,0},{0,0,0,iu},{0,0,-iu,0}},
                              {{0,0,1,0},{0,0,0,-iu},{1,0,0,0},{0,iu,0,0}},
                              {{0,0,0,1},{0,0,iu,0},{0,-iu,0,0},{1,0,0,0}}};

int main(int argc, char *argv[])
{
    bool TermsSimplified = true;
    //Specify symmetry properties of the considered model
    bool U1Symmetric = true;
    bool TRSymmetric = true;

    string termsA = "(";
    string termsB = "(";
    string termsC = "";

    //Specify for which matrix element of the spin correlation terms are given out
    int mu = 2;
    int nu = 1;

    string tempTermA1;
    string tempTermB1;
    for(int a=0;a<4;a++){
        //Propagators only have 0 and z components in case of a global U(1) spin rotation symmetry
        if(U1Symmetric){a=3*a;}
        if(a>0 && TRSymmetric){continue;}

        if(TermsSimplified){
            tempTermA1 += "+ga"+to_string(a)+"*(";
            tempTermB1 += "+ga"+to_string(a)+"*(";
        }else{
            tempTermA1 += "+iGLam("+to_string(a/3)+",w,G_vec)*(";
            tempTermB1 += "+iGLam("+to_string(a/3)+",w,G_vec)*(";
        }
        string tempTermA2;
        string tempTermB2;
        for(int b=0;b<4;b++){
            if(U1Symmetric){b=3*b;}
            if(b>0 && TRSymmetric){continue;}

            if(TermsSimplified){
                tempTermA2 += "+gb2_"+to_string(b)+"*(";
                tempTermB2 += "+gb2_"+to_string(b)+"*(";
            }else{
                tempTermA2 += "+iGLam("+to_string(b/3)+",w2,G_vec)*(";
                tempTermB2 += "+iGLam("+to_string(b/3)+",w2,G_vec)*(";
            }
            string tempTermA3;
            string tempTermB3;
            for(int c=0;c<4;c++){
                if(U1Symmetric){c=3*c;}
                if(c>0 && TRSymmetric){continue;}

                if(TermsSimplified){
                    tempTermA3 += "+ga2_"+to_string(c)+"*(";
                    tempTermB3 += "+ga2_"+to_string(c)+"*(";
                }else{
                    tempTermA3 += "+iGLam("+to_string(c/3)+",w,G_vec)*(";
                    tempTermB3 += "+iGLam("+to_string(c/3)+",w,G_vec)*(";
                }

                string tempTermA4;
                string tempTermB4;
                for(int d=0;d<4;d++){
                    if(U1Symmetric){d=3*d;}
                    if(d>0 && TRSymmetric){continue;}

                    complex<int> propagatorPrefactor = 1;
                    if(a==0){propagatorPrefactor *= -iu;}
                    if(b==0){propagatorPrefactor *= -iu;}
                    if(c==0){propagatorPrefactor *= -iu;}
                    if(d==0){propagatorPrefactor *= -iu;}

                    if(TermsSimplified){
                        tempTermA4 += "+gb"+to_string(d)+"*(";
                        tempTermB4 += "+gb"+to_string(d)+"*(";
                    }else{
                        tempTermA4 += "+iGLam("+to_string(d/3)+",w2,G_vec)*(";
                        tempTermB4 += "+iGLam("+to_string(d/3)+",w2,G_vec)*(";
                    }

                    string tempTermASub = "";
                    string tempTermBSub = "";

                    for(int g=0;g<4;g++){
                        for(int h=0;h<4;h++){
                            bool vertexImaginary = false;
                            if((g==0 && h!=0) || (g!=0 && h==0))
                            {
                                vertexImaginary = true;
                            }
                            for(int k=0;k<4;k++){
                                for(int l=0;l<4;l++){
                                    for(int m=0;m<4;m++){
                                        for(int p=0;p<4;p++){
                                            //First term (without Kronecker delta)
                                            //Evaluation of the traces
                                            complex<int> f1 = beta[k][mu][c]*beta[l][k][g]*beta[0][l][a];
                                            complex<int> f2 = beta[m][nu][d]*beta[p][m][h]*beta[0][p][b];

                                            bool s1 = (f1*f2 != 0);

                                            if(s1){
                                                string prefactorSymbol="";
                                                complex<int> f12=propagatorPrefactor*f1*f2;
                                                if(vertexImaginary)
                                                {
                                                    f12 *=iu;
                                                }
                                                if(f12 == 1){
                                                    prefactorSymbol="+";
                                                }else if(f12 == -1){
                                                    prefactorSymbol="-";
                                                }else if(f12 == iu){
                                                    prefactorSymbol="+i";
                                                }else if(f12 == -iu){
                                                    prefactorSymbol="-i";
                                                }else if(f12 == 0){
                                                    continue;
                                                }
                                                if(TermsSimplified){
                                                    tempTermASub += prefactorSymbol+"G"+to_string(g)+to_string(h);
                                                }else{
                                                    tempTermASub += prefactorSymbol+"getIntpolG(G_vec,1,"+to_string(g)+","+to_string(h)+", R, pw_wpw2,1,pw_wmw2)";
                                                }
                                            }

                                            for(int q=0;q<4;q++){
                                                for(int r=0;r<4;r++){
                                                    //Second term (with Kronecker delta)
                                                    //Evaluation of the trace
                                                    complex<int> f3 = beta[k][mu][c]*beta[l][k][g]*beta[m][l][b]*beta[p][m][nu]*beta[q][p][d]*beta[r][q][h]*beta[0][r][a];

                                                    bool s3 = (f3 != 0);

                                                    if(s3){
                                                        string prefactorSymbol="";
                                                        f3 *= propagatorPrefactor;
                                                        if(vertexImaginary)
                                                        {
                                                            f3 *=iu;
                                                        }
                                                        if(f3 == 1){
                                                            prefactorSymbol="+";
                                                        }else if(f3 == -1){
                                                            prefactorSymbol="-";
                                                        }else if(f3 == iu){
                                                            prefactorSymbol="+i";
                                                        }else if(f3 == -iu){
                                                            prefactorSymbol="-i";
                                                        }else if(f3 == 0){
                                                            continue;
                                                        }

                                                        if(TermsSimplified){
                                                            tempTermBSub += prefactorSymbol+"G"+to_string(g)+to_string(h);
                                                        }else{
                                                            tempTermBSub += prefactorSymbol+"getIntpolG(G_vec,1,"+to_string(g)+","+to_string(h)+", R0, pw_wpw2,pw_wmw2,1)";
                                                        }
                                                    }
                                                }}}}

                                }}}}
                    tempTermA4 += tempTermASub;
                    tempTermB4 += tempTermBSub;
                    tempTermA4 += ")";
                    tempTermB4 += ")";
                    //termsA += tempTermA4;
                    //termsB += tempTermB4;
                }
                tempTermA3 += tempTermA4;
                tempTermB3 += tempTermB4;
                tempTermA3 += ")";
                tempTermB3 += ")";
                //termsA += tempTermA4;
                //termsB += tempTermB4;
            }
            tempTermA2 += tempTermA3;
            tempTermB2 += tempTermB3;
            tempTermA2 += ")";
            tempTermB2 += ")";
        }
        tempTermA1 += tempTermA2;
        tempTermB1 += tempTermB2;
        tempTermA1 += ")";
        tempTermB1 += ")";
    }

    termsA += tempTermA1+")";
    termsB += tempTermB1+")";

    //Give out all terms
    cout << "First terms(includes summation over lattice sites): " << endl << termsA << endl;
    cout << "Second terms: " << endl << termsB << endl;


    //Generate terms for the spin correlation that are quadratic in propagators

    for(int a=0;a<4;a++){
        //Propagators only have 0 and z components in case of a global U(1) spin rotation symmetry
        if(U1Symmetric){a=3*a;}
        if(a>0 && TRSymmetric){continue;}
        for(int b=0;b<4;b++){
            //Propagators only have 0 and z components in case of a global U(1) spin rotation symmetry
            if(U1Symmetric){b=3*b;}
            if(b>0 && TRSymmetric){continue;}
            complex<int> propagatorPrefactor = 1;
            if(a==0){propagatorPrefactor *= -iu;}
            if(b==0){propagatorPrefactor *= -iu;}

            for(int l=0;l<4;l++){
                for(int m=0;m<4;m++){

                    //Evaluation of the trace
                    complex<int> f = beta[l][mu][a]*beta[m][l][nu]*beta[0][m][b];

                    if(f != 0){
                        string prefactorSymbol="";
                        complex<int> prefactor=-propagatorPrefactor*f;
                        if(prefactor == 1){
                            prefactorSymbol="+";
                        }else if(prefactor == -1){
                            prefactorSymbol="-";
                        }else if(prefactor == iu){
                            prefactorSymbol="+i";
                        }else if(prefactor == -iu){
                            prefactorSymbol="-i";
                        }else if(prefactor == 0){
                            continue;
                        }
                        termsC += prefactorSymbol+"ga" + to_string(a) + "*ga" + to_string(b);
                    }
                }}

        }}

    cout << "Terms quadratic in propagators: " << endl;
    cout << termsC << endl;



    return 0;
}
