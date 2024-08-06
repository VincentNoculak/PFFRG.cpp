//This file generates all terms of the two-particle vertex flow equation for a general model without any continuous spin rotation symmetries and with finite magnetic fields.
//The terms are obtained by expressing the flow equation as done in Eq. (35) of Phys. Rev. B 109, 174414.
//The generated terms are inserted in the PFFRG code for a Heisenberg model on a triangular lattice in a magnetic field.
//Some aspects of this code may become more clearer when considered together with the corresponding PFFRG code.

#include <stdio.h>
#include <math.h>
#include <iostream>
#include <fstream>
#include <complex>
#include <array>
#include <stdlib.h>
#include <iomanip>
#include <algorithm>
#include <list>
#include<vector>

using namespace std;

//Decide whether the flow equations should be generated for a model with a global U(1) spin rotation symmetry (the "for" loops for the propagator variables will also have to be adjusted in this case)
bool U1Symmetric = false;


int TermCounter[5] = {0,0,0,0,0};

complex<int> iu(0,1);

//The $\bm{\beta}_{a}$ matrices are defined by the equation
// \bm{\sigma}^{\mu} \bm{\sigma}^{\nu} = \sum_{a=0,x,y,z} \beta^{\mu\nu}_{a} \bm{\sigma}^{a},
//with $\bm{\sigma}^{\mu}$ being Pauli matrices for $\mu=x,y,z$ and the identity matrix for $\mu=0$.
complex<int> beta[4][4][4] = {{{1,0,0,0},{0,1,0,0},{0,0,1,0},{0,0,0,1}},
                              {{0,1,0,0},{1,0,0,0},{0,0,0,iu},{0,0,-iu,0}},
                              {{0,0,1,0},{0,0,0,-iu},{1,0,0,0},{0,iu,0,0}},
                              {{0,0,0,1},{0,0,iu,0},{0,-iu,0,0},{1,0,0,0}}};

int subl = 1; //Sublattice-dependent values subl=0,1,2 that are relevant to the RPA-channel summation over lattice sites. This variable is applied in the PFFRG implementation of the triangular lattice Heisenberg model in a magnetic field, where multiple sublattices inequivalent by symmetries exist.

int main(int argc, char *argv[])
{
    string Channel[21][4][4]; //Gathers terms in the form of strings
    string SDETerms[2];

    int maxComponent = 3;
    int prFactor = 1;
    if(U1Symmetric){
        maxComponent = 1;
        prFactor = 3;
    }

    for(int mu=0;mu<4;mu++){
        for(int nu=0;nu<4;nu++){
            bool vertexImaginary = false;
            if((mu==0 && nu!=0) || (mu!=0 && nu==0))
            {
                vertexImaginary = true;
            }

            //Channel A (s channel)
            for(int prA=0;prA<=maxComponent;prA++){
                for(int prB=0;prB<=maxComponent;prB++){
                    /*if(prB==0 && prA==1){
                    continue;
                }*/
                    Channel[0][mu][nu] += "+(";
                    for(int a=0;a<4;a++){
                        for(int b=0;b<4;b++){
                            for(int d=0;d<4;d++){
                                for(int c=0;c<4;c++){
                                    for(int g=0;g<4;g++){
                                        for(int h=0;h<4;h++){
                                            if(U1Symmetric)
                                            {
                                                if(!((a==b) || (a==1 && b==2) || (a==2 && b==1) || (a==0 && b==3) || (a==3 && b==0))){
                                                    continue;
                                                }
                                                if(!((c==d) || (c==1 && d==2) || (c==2 && d==1) || (c==0 && d==3) || (c==3 && d==0))){
                                                    continue;
                                                }
                                            }
                                            //Additional factor, because vertices can be imaginary
                                            complex<int> extraFactor=1;
                                            if((a==0 && b!=0) || (a!=0 && b==0))
                                            {
                                                extraFactor*=iu;
                                            }
                                            if((c==0 && d!=0) || (c!=0 && d==0))
                                            {
                                                extraFactor*=iu;
                                            }
                                            if(prA==0){
                                                extraFactor*=-iu;
                                            }
                                            if(prB==0){
                                                extraFactor*=-iu;
                                            }

                                            complex<int> prefactor = extraFactor*beta[g][a][prFactor*prA]*beta[mu][g][c]*beta[h][b][prFactor*prB]*beta[nu][h][d];
                                            //There should actually be an additional minus sign, which can be implemented as in the following lines that are commented out.
                                            //However, the lines are commented out, beause the minus sign is implemented directly in the PFFRG code and not here in the code generator.

                                            //"Why one extra (-1) Factor?": It comes from a frequency transformation omega_4->-omega_4
                                            //prB  has to belong to the propagator with argument 4 for for the following implementation to be correct (this is important)
                                            /*if(prB==0){
                                                prefactor *= -1;
                                                }*/

                                            if(vertexImaginary)
                                            {
                                                prefactor *=-iu;
                                            }
                                            string prefactorSymbol="";
                                            if(prefactor == 1){
                                                prefactorSymbol=" +";
                                            }else if(prefactor == -1){
                                                prefactorSymbol=" -";
                                            }else if(prefactor == iu){
                                                prefactorSymbol="+i";
                                            }else if(prefactor == -iu){
                                                prefactorSymbol="-i";
                                            }else if(prefactor == 0){
                                                continue;
                                            }
                                            //if(a==b && c==d){
                                            string tempTerm = prefactorSymbol+"G"+to_string(a)+to_string(b)+"*"+"G"+to_string(c)+to_string(d)+" ";
                                            string tempTermCodeFormatShort = prefactorSymbol+"Ch1A1_"+to_string(a)+to_string(b)+"*Ch1A2_"+to_string(c)+to_string(d);
                                            string tempTermCodeFormat = prefactorSymbol+"getIntpolG(G_vec,1,"+to_string(a)+","+to_string(b)+", R, ns,pw_1a,pw_1b)*getIntpolG(G_vec,1,"+to_string(c)+","+to_string(d)+", R, ns,pw_2a,pw_2b)";
                                            string tempTermCodeFormatDiag = prefactorSymbol+"getIntpolG(G_vec,1,"+to_string(a)+", R, ns,pw_1a,pw_1b)*getIntpolG(G_vec,1,"+to_string(c)+", R, ns,pw_2a,pw_2b)";
                                            string tempTermCodeFormatTwoLoop = prefactorSymbol+"(getIntpolG(DG_TwoLoop_T,1,"+to_string(a)+","+to_string(b)+",0, R, ns,pw_1a,pw_1b)+getIntpolG(DG_TwoLoop_U,1,"+to_string(a)+","+to_string(b)+",0, R, ns,pw_1a,pw_1b))"+"*getIntpolG(G_vec,1,"+to_string(c)+","+to_string(d)+",0, R, ns,pw_2a,pw_2b)"
                                                    +prefactorSymbol+"getIntpolG(G_vec,1,"+to_string(a)+","+to_string(b)+",0, R, ns,pw_1a,pw_1b)*(getIntpolG(DG_TwoLoop_T,1,"+to_string(c)+","+to_string(d)+",0, R, ns,pw_2a,pw_2b)+getIntpolG(DG_TwoLoop_U,1,"+to_string(c)+","+to_string(d)+",0, R, ns,pw_2a,pw_2b))";
                                            string tempTermCodeFormatTwoLoopShort = prefactorSymbol+"Ch1A1_"+to_string(a)+to_string(b)+"TL*Ch1A2_"+to_string(c)+to_string(d)
                                                    +prefactorSymbol+"Ch1A1_"+to_string(a)+to_string(b)+"*Ch1A2_"+to_string(c)+to_string(d)+"TL";

                                            Channel[0][mu][nu] += tempTermCodeFormatShort;
                                            TermCounter[0]++;
                                            //Channel[0][mu][nu] += tempTermCodeFormatTwoLoop;
                                            //}
                                        }
                                    }
                                }
                            }}}
                    Channel[0][mu][nu] += ")*Pt"+to_string(prA)+to_string(prB);
                }}

            //Channel B (t channel, RPA channel)
            for(int prA=0;prA<=maxComponent;prA++){
                for(int prB=0;prB<=maxComponent;prB++){
                    for(int b=0;b<4;b++){
                        for(int c=0;c<4;c++){
                            for(int j=0;j<4;j++){
                                if(U1Symmetric)
                                {
                                    if(!((mu==b) || (mu==1 && b==2) || (mu==2 && b==1) || (mu==0 && b==3) || (mu==3 && b==0))){
                                        continue;
                                    }
                                    if(!((c==nu) || (c==1 && nu==2) || (c==2 && nu==1) || (c==0 && nu==3) || (c==3 && nu==0))){
                                        continue;
                                    }
                                }
                                complex<int> extraFactor=1;
                                if((mu==0 && b!=0) || (mu!=0 && b==0))
                                {
                                    extraFactor*=iu;
                                }
                                if((c==0 && nu!=0) || (c!=0 && nu==0))
                                {
                                    extraFactor*=iu;
                                }
                                if(prA==0){
                                    extraFactor*=-iu;
                                }
                                if(prB==0){
                                    extraFactor*=-iu;
                                }
                                complex<int> prefactor = extraFactor*beta[j][b][prFactor*prA]*beta[prFactor*prB][j][c];
                                if(vertexImaginary)
                                {
                                    prefactor *=-iu;
                                }

                                //Neglect factor 2 here, because we include the factor directly in the setG method of the PFFRG code
                                string prefactorSymbol="";
                                if(prefactor == 1){
                                    prefactorSymbol="-";
                                }else if(prefactor == -1){
                                    prefactorSymbol="+";
                                }else if(prefactor == iu){
                                    prefactorSymbol="-i*";
                                }else if(prefactor == -iu){
                                    prefactorSymbol="+i*";
                                }else if(prefactor == 0){
                                    continue;
                                }
                                //if(mu==b &&c==nu){
                                string tempTerm = prefactorSymbol+"G"+to_string(mu)+to_string(b)+"*"+"G"+to_string(c)+to_string(nu)+" ";
                                string tempTermCodeFormatShort = prefactorSymbol+"Ch2A1_"+to_string(mu)+to_string(b)+"*Ch2A2_"+to_string(c)+to_string(nu);
                                //string optimizedTerms = prefactorSymbol+"vertexProductBuffer["+to_string(64*mu+16*b+4*c+nu)+"]";
                                string optimizedTerms = "vertexProductBuffer["+to_string(64*mu+16*b+4*c+nu)+"]";
                                string tempTermCodeFormat = prefactorSymbol+"getIntpolG(G_vec,1,"+to_string(mu)+","+to_string(b)+",R1j,pw_1a,nt,pw_1b)*getIntpolG(G_vec,1,"+to_string(c)+","+to_string(nu)+", Rj2, pw_2a,nt,pw_2b)";
                                string tempTermCodeFormatDiag = prefactorSymbol+"getIntpolG(G_vec,1,"+to_string(mu)+",R1j,pw_1a,nt,pw_1b)*getIntpolG(G_vec,1,"+to_string(c)+", {-Rj2.a1,-Rj2.a2,Rj2.b}, pw_2a,nt,pw_2b)";
                                string tempTermCodeFormatTwoLoop = prefactorSymbol+"(getIntpolG(DG_TwoLoop_S,1,"+to_string(mu)+","+to_string(b)+",0,R1j,pw_1a,nt,pw_1b)+getIntpolG(DG_TwoLoop_U,1,"+to_string(mu)+","+to_string(b)+",0,R1j,pw_1a,nt,pw_1b))*getIntpolG(G_vec,1,"+to_string(c)+","+to_string(nu)+",0, {-Rj2.a1,-Rj2.a2,Rj2.b}, pw_2a,nt,pw_2b)"
                                        +prefactorSymbol+"getIntpolG(G_vec,1,"+to_string(mu)+","+to_string(b)+",0,R1j,pw_1a,nt,pw_1b)*(getIntpolG(DG_TwoLoop_S,1,"+to_string(c)+","+to_string(nu)+",0, {-Rj2.a1,-Rj2.a2,Rj2.b}, pw_2a,nt,pw_2b)+getIntpolG(DG_TwoLoop_U,1,"+to_string(c)+","+to_string(nu)+",0, {-Rj2.a1,-Rj2.a2,Rj2.b}, pw_2a,nt,pw_2b))";
                                string tempTermCodeFormatTwoLoopShort = prefactorSymbol+"Ch2A1_"+to_string(mu)+to_string(b)+"TL*Ch2A2_"+to_string(c)+to_string(nu)
                                        +prefactorSymbol+"Ch2A1_"+to_string(mu)+to_string(b)+"*Ch2A2_"+to_string(c)+to_string(nu)+"TL";
                                //                        if(((m==3 || m==2)&&!((nu==3 || nu==2)))||(!(m==3 || m==2)&&((nu==3 || nu==2)))){
                                //                             tempTermCodeFormatShort = prefactorSymbol+"Ch2A1_"+to_string(mu)+to_string(m)+"*(-1)*Ch2A2_"+to_string(m)+to_string(nu);
                                //            tempTermCodeFormat = prefactorSymbol+"getIntpolG(G_vec,1,"+to_string(mu)+","+to_string(m)+",R1j,pw_1a,nt,pw_1b)*(-1)*getIntpolG(G_vec,1,"+to_string(m)+","+to_string(nu)+", Rj2, pw_2a,nt,pw_2b)";
                                //        }


                                string PropagatorTerms = prefactorSymbol+"Pt["+ to_string(subl) +"][" + to_string(4*prA+prB) +"]";
                                TermCounter[1]++;
                                //Channel[5+prA*4+prB][mu][nu] += optimizedTerms;//tempTermCodeFormatShort;
                                Channel[5+b*4+c][mu][nu] +=  PropagatorTerms;//tempTermCodeFormatShort;
                                //Channel[1][mu][nu] += tempTermCodeFormatTwoLoop;
                                //}
                            }
                        }}
                }}
            //Channel C (t channel)
            for(int prA=0;prA<=maxComponent;prA++){
                for(int prB=0;prB<=maxComponent;prB++){
                    Channel[2][mu][nu] += "+(";
                    for(int b=0;b<4;b++){
                        for(int c=0;c<4;c++){
                            for(int d=0;d<4;d++){
                                for(int j=0;j<4;j++){
                                    for(int k=0;k<4;k++){
                                        for(int l=0;l<4;l++){
                                            if(U1Symmetric)
                                            {
                                                if(!((mu==b) || (mu==1 && b==2) || (mu==2 && b==1) || (mu==0 && b==3) || (mu==3 && b==0))){
                                                    continue;
                                                }
                                                if(!((c==d) || (c==1 && d==2) || (c==2 && d==1) || (c==0 && d==3) || (c==3 && d==0))){
                                                    continue;
                                                }
                                            }
                                            complex<int> extraFactor=1;
                                            if((mu==0 && b!=0) || (mu!=0 && b==0))
                                            {
                                                extraFactor*=iu;
                                            }
                                            if((d==0 && c!=0) || (d!=0 && c==0))
                                            {
                                                extraFactor*=iu;
                                            }
                                            if(prA==0){
                                                extraFactor*=-iu;
                                            }
                                            if(prB==0){
                                                extraFactor*=-iu;
                                            }


                                            complex<int> prefactor = extraFactor*beta[j][d][prFactor*prB]*beta[k][j][b]*beta[l][k][prFactor*prA]*beta[nu][l][c];
                                            if(vertexImaginary)
                                            {
                                                prefactor *=-iu;
                                            }

                                            string prefactorSymbol="";
                                            if(prefactor == 1){
                                                prefactorSymbol=" +";
                                            }else if(prefactor == -1){
                                                prefactorSymbol=" -";
                                            }else if(prefactor == iu){
                                                prefactorSymbol="+i*";
                                            }else if(prefactor == -iu){
                                                prefactorSymbol="-i*";
                                            }else if(prefactor == 0){
                                                continue;
                                            }
                                            //if(mu==b && c==d){
                                            string tempTerm = prefactorSymbol+"G"+to_string(mu)+to_string(b)+"*"+"G"+to_string(c)+to_string(d)+" ";
                                            //                            string tempTermCodeFormat = prefactorSymbol+"getIntpolG(G_vec,1,"+to_string(mu)+","+to_string(ml)+",0, R , pw_1a,nt,pw_1b)*invertFactor(R.b,"+to_string(mm)+","+to_string(mr)+")*getIntpolG(G_vec,1,"+to_string(mm)+","+to_string(mr)+",0, R0, pw_2a,pw_2b,nt)";
                                            //string tempTermCodeFormatShort = prefactorSymbol+"Ch3A1_"+to_string(mu)+to_string(ml)+"*Ch3A2_"+to_string(mm)+to_string(mr);
                                            //string tempTermCodeFormatShort = prefactorSymbol+"Ch3A1_"+to_string(mu)+to_string(b)+"*invertFactor(R.b,"+to_string(c)+","+to_string(d)+")*Ch3A2_"+to_string(c)+to_string(d);
                                            string tempTermCodeFormatShort = prefactorSymbol+"Ch3A1_"+to_string(mu)+to_string(b)+"*Ch3A2_"+to_string(c)+to_string(d);
                                            string tempTermCodeFormat = prefactorSymbol+"getIntpolG(G_vec,1,"+to_string(mu)+","+to_string(b)+", R , pw_1a,nt,pw_1b)*getIntpolG(G_vec,1,"+to_string(c)+","+to_string(d)+", R0, pw_2a,pw_2b,nt)";
                                            string tempTermCodeFormatDiag = prefactorSymbol+"getIntpolG(G_vec,1,"+to_string(mu)+", R , pw_1a,nt,pw_1b)*getIntpolG(G_vec,1,"+to_string(c)+", R0, pw_2a,pw_2b,nt)";
                                            string tempTermCodeFormatTwoLoop = prefactorSymbol+"(getIntpolG(DG_TwoLoop_S,1,"+to_string(mu)+","+to_string(b)+",0, R , pw_1a,nt,pw_1b)+getIntpolG(DG_TwoLoop_U,1,"+to_string(mu)+","+to_string(b)+",0, R , pw_1a,nt,pw_1b))*getIntpolG(G_vec,1,"+to_string(c)+","+to_string(d)+",0, R0, pw_2a,pw_2b,nt)"
                                                    +prefactorSymbol+"getIntpolG(G_vec,1,"+to_string(mu)+","+to_string(b)+",0, R , pw_1a,nt,pw_1b)*(getIntpolG(DG_TwoLoop_S,1,"+to_string(c)+","+to_string(d)+",0, R0, pw_2a,pw_2b,nt)+getIntpolG(DG_TwoLoop_T,1,"+to_string(c)+","+to_string(d)+",0, R0, pw_2a,pw_2b,nt))";
                                            string tempTermCodeFormatTwoLoopShort = prefactorSymbol+"Ch3A1_"+to_string(mu)+to_string(b)+"TL*Ch3A2_"+to_string(c)+to_string(d)
                                                    +prefactorSymbol+"Ch3A1_"+to_string(mu)+to_string(b)+"*Ch3A2_"+to_string(c)+to_string(d)+"TL";
                                            Channel[2][mu][nu] += tempTermCodeFormatShort;
                                            //Channel[2][mu][nu] += tempTermCodeFormatTwoLoop;
                                            //}
                                            TermCounter[2]++;
                                        }
                                    }
                                }
                            }}}
                    Channel[2][mu][nu] += ")*Pt3_"+to_string(prA)+to_string(prB);
                }}

            //Channel D (t channel)
            for(int prA=0;prA<=maxComponent;prA++){
                for(int prB=0;prB<=maxComponent;prB++){
                    Channel[3][mu][nu] += "+(";
                    for(int c=0;c<4;c++){
                        for(int b=0;b<4;b++){
                            for(int a=0;a<4;a++){
                                for(int j=0;j<4;j++){
                                    for(int k=0;k<4;k++){
                                        for(int l=0;l<4;l++){
                                            if(U1Symmetric)
                                            {
                                                if(!((a==b) || (a==1 && b==2) || (a==2 && b==1) || (a==0 && b==3) || (a==3 && b==0))){
                                                    continue;
                                                }
                                                if(!((c==nu) || (c==1 && nu==2) || (c==2 && nu==1) || (c==0 && nu==3) || (c==3 && nu==0))){
                                                    continue;
                                                }
                                            }
                                            complex<int> extraFactor=1;
                                            if((nu==0 && c!=0) || (nu!=0 && c==0))
                                            {
                                                extraFactor*=iu;
                                            }
                                            if((a==0 && b!=0) || (a!=0 && b==0))
                                            {
                                                extraFactor*=iu;
                                            }
                                            if(prA==0){
                                                extraFactor*=-iu;
                                            }
                                            if(prB==0){
                                                extraFactor*=-iu;
                                            }

                                            complex<int> prefactor = extraFactor*beta[j][a][prFactor*prA]*beta[k][j][c]*beta[l][k][prFactor*prB]*beta[mu][l][b];
                                            if(vertexImaginary)
                                            {
                                                prefactor *=-iu;
                                            }

                                            string prefactorSymbol="";
                                            if(prefactor == 1){
                                                prefactorSymbol=" +";
                                            }else if(prefactor == -1){
                                                prefactorSymbol=" -";
                                            }else if(prefactor == iu){
                                                prefactorSymbol="+i*";
                                            }else if(prefactor == -iu){
                                                prefactorSymbol="-i*";
                                            }else if(prefactor == 0){
                                                continue;
                                            }
                                            //if(a==b && c==nu){
                                            string tempTerm = prefactorSymbol+"G"+to_string(a)+to_string(b)+"*"+"G"+to_string(c)+to_string(nu)+" ";
                                            string tempTermCodeFormatShort = prefactorSymbol+"Ch4A1_"+to_string(a)+to_string(b)+"*Ch4A2_"+to_string(c)+to_string(nu);
                                            string tempTermCodeFormat = prefactorSymbol+"getIntpolG(G_vec,1,"+to_string(a)+","+to_string(b)+", R0, pw_1a,pw_1b,nt)*getIntpolG(G_vec,1,"+to_string(c)+","+to_string(nu)+", R , pw_2a,nt,pw_2b)";
                                            string tempTermCodeFormatDiag = prefactorSymbol+"getIntpolG(G_vec,1,"+to_string(a)+", R0, pw_1a,pw_1b,nt)*getIntpolG(G_vec,1,"+to_string(c)+", R , pw_2a,nt,pw_2b)";
                                            string tempTermCodeFormatTwoLoop = prefactorSymbol+"(getIntpolG(DG_TwoLoop_S,1,"+to_string(a)+","+to_string(b)+",0, R0, pw_1a,pw_1b,nt)+getIntpolG(DG_TwoLoop_T,1,"+to_string(a)+","+to_string(b)+",0, R0, pw_1a,pw_1b,nt))*getIntpolG(G_vec,1,"+to_string(c)+","+to_string(nu)+",0, R , pw_2a,nt,pw_2b)"
                                                    +prefactorSymbol+"getIntpolG(G_vec,1,"+to_string(a)+","+to_string(b)+",0, R0, pw_1a,pw_1b,nt)*(getIntpolG(DG_TwoLoop_S,1,"+to_string(c)+","+to_string(nu)+",0, R , pw_2a,nt,pw_2b)+getIntpolG(DG_TwoLoop_U,1,"+to_string(c)+","+to_string(nu)+",0, R , pw_2a,nt,pw_2b))";
                                            string tempTermCodeFormatTwoLoopShort = prefactorSymbol+"Ch4A1_"+to_string(a)+to_string(b)+"TL*Ch4A2_"+to_string(c)+to_string(nu)
                                                    +prefactorSymbol+"Ch4A1_"+to_string(a)+to_string(b)+"*Ch4A2_"+to_string(c)+to_string(nu)+"TL";
                                            Channel[3][mu][nu] += tempTermCodeFormatShort;
                                            //Channel[3][mu][nu] += tempTermCodeFormatTwoLoop;
                                            //}
                                            TermCounter[3]++;
                                        }
                                    }
                                }
                            }}}
                    Channel[3][mu][nu] += ")*Pt4_"+to_string(prA)+to_string(prB);
                }}


            //Channel E (u channel)
            for(int prA=0;prA<=maxComponent;prA++){
                for(int prB=0;prB<=maxComponent;prB++){
                    Channel[4][mu][nu] += "+(";
                    for(int b=0;b<4;b++){
                        for(int a=0;a<4;a++){
                            for(int c=0;c<4;c++){
                                for(int d=0;d<4;d++){
                                    for(int g=0;g<4;g++){
                                        for(int h=0;h<4;h++){
                                            if(U1Symmetric)
                                            {
                                                if(!((a==b) || (a==1 && b==2) || (a==2 && b==1) || (a==0 && b==3) || (a==3 && b==0))){
                                                    continue;
                                                }
                                                if(!((c==d) || (c==1 && d==2) || (c==2 && d==1) || (c==0 && d==3) || (c==3 && d==0))){
                                                    continue;
                                                }
                                            }
                                            //Additional factor, because vertices can be imaginary
                                            complex<int> extraFactor=1;
                                            if((b==0 && a!=0) || (b!=0 && a==0))
                                            {
                                                extraFactor*=iu;
                                            }
                                            if((d==0 && c!=0) || (d!=0 && c==0))
                                            {
                                                extraFactor*=iu;
                                            }
                                            if(prA==0){
                                                extraFactor*=-iu;
                                            }
                                            if(prB==0){
                                                extraFactor*=-iu;
                                            }

                                            complex<int> prefactor = extraFactor*beta[g][d][prFactor*prB]*beta[mu][g][b]*beta[h][a][prFactor*prA]*beta[nu][h][c];
                                            if(vertexImaginary)
                                            {
                                                prefactor *=-iu;
                                            }
                                            //This is a consequence of the \omega_4->-\omega_4 transformation
                                            if(prA != prB){
                                                prefactor *=-1;
                                            }

                                            string prefactorSymbol="";
                                            if(prefactor == 1){
                                                prefactorSymbol=" +";
                                            }else if(prefactor == -1){
                                                prefactorSymbol=" -";
                                            }else if(prefactor == iu){
                                                prefactorSymbol="+i*";
                                            }else if(prefactor == -iu){
                                                prefactorSymbol="-i*";
                                            }else if(prefactor == 0){
                                                continue;
                                            }
                                            //if(a == b && c == d){
                                            string tempTerm = prefactorSymbol+"G"+to_string(b)+to_string(a)+"*"+"G"+to_string(d)+to_string(c)+" ";
                                            string tempTermCodeFormatShort = prefactorSymbol+"Ch5A1_"+to_string(b)+to_string(a)+"*Ch5A2_"+to_string(d)+to_string(c);
                                            string tempTermCodeFormat = prefactorSymbol+"getIntpolG(G_vec,1,"+to_string(b)+","+to_string(a)+", R, pw_1a,pw_1b,nu)*getIntpolG(G_vec,1,"+to_string(d)+","+to_string(c)+", R, pw_2a,pw_2b,nu)";
                                            string tempTermCodeFormatDiag = prefactorSymbol+"getIntpolG(G_vec,1,"+to_string(a)+", R, pw_1a,pw_1b,nu)*getIntpolG(G_vec,1,"+to_string(c)+", R, pw_2a,pw_2b,nu)";
                                            string tempTermCodeFormatTwoLoop = prefactorSymbol+"(getIntpolG(DG_TwoLoop_S,1,"+to_string(b)+","+to_string(a)+",0, R, pw_1a,pw_1b,nu)*getIntpolG(DG_TwoLoop_T,1,"+to_string(b)+","+to_string(a)+",0, R, pw_1a,pw_1b,nu))*getIntpolG(G_vec,1,"+to_string(d)+","+to_string(c)+",0, R, pw_2a,pw_2b,nu)"
                                                    +prefactorSymbol+"getIntpolG(G_vec,1,"+to_string(b)+","+to_string(a)+",0, R, pw_1a,pw_1b,nu)*(getIntpolG(DG_TwoLoop_S,1,"+to_string(d)+","+to_string(c)+",0, R, pw_2a,pw_2b,nu)+getIntpolG(DG_TwoLoop_T,1,"+to_string(d)+","+to_string(c)+",0, R, pw_2a,pw_2b,nu))";
                                            string tempTermCodeFormatTwoLoopShort = prefactorSymbol+"Ch5A1_"+to_string(b)+to_string(a)+"TL*Ch5A2_"+to_string(d)+to_string(c)
                                                    +prefactorSymbol+"Ch5A1_"+to_string(b)+to_string(a)+"*Ch5A2_"+to_string(d)+to_string(c)+"TL";
                                            Channel[4][mu][nu] += tempTermCodeFormatShort;
                                            //Channel[4][mu][nu] += tempTermCodeFormatTwoLoop;
                                            //}
                                            TermCounter[4]++;
                                        }
                                    }
                                }
                            }}}
                    Channel[4][mu][nu] += ")*Pt"+to_string(prA)+to_string(prB);
                }}


        }
    }


    //Important: Verify that the generation of SDE terms is correct in case you plan to apply it.
    for(int mu=0;mu<2;mu++){
        SDETerms[mu] += "BV*(";
        //Generate Schwinger-Dyson equation terms (only works for Heisenberg interactions at the moment)
        for(int prA=0;prA<=maxComponent;prA++){
            SDETerms[mu] += "+(";
            string tempTermsA = "";
            for(int prB=0;prB<=maxComponent;prB++){
                SDETerms[mu] += "+(";
                string tempTermsB = "";
                for(int prC=0;prC<=maxComponent;prC++){
                    SDETerms[mu] += "+(";
                    //SDETerms[0][mu] += "+(";
                    //SDETerms[1][mu] += "+(";
                    string tempTerms = "";
                    for(int a=1;a<4;a++){
                        for(int b=0;b<4;b++){
                            for(int c=0;c<4;c++){
                                for(int k=0;k<4;k++){
                                    for(int l=0;l<4;l++){
                                        for(int m=0;m<4;m++){
                                            if(U1Symmetric)
                                            {
                                                if(!((c==b) || (c==1 && b==2) || (c==2 && b==1) || (c==0 && b==3) || (c==3 && b==0))){
                                                    continue;
                                                }
                                            }
                                            //Additional factor, because vertices can be imaginary
                                            complex<int> extraFactor=1;
                                            if((c==0 && b!=0) || (c!=0 && b==0))
                                            {
                                                extraFactor*=iu;
                                            }
                                            if(prA==0){
                                                extraFactor*=-iu;
                                            }
                                            if(prB==0){
                                                extraFactor*=-iu;
                                            }
                                            if(prC==0){
                                                extraFactor*=-iu;
                                            }
                                            if(mu==0){
                                                extraFactor*=iu;
                                            }

                                            complex<int> prefactor1 = extraFactor*beta[k][a][prFactor*prA]*beta[prFactor*mu][k][b]*beta[l][prFactor*prC][a]*beta[m][l][prFactor*prB];
                                            complex<int> prefactor2 = extraFactor*beta[k][a][prFactor*prB]*beta[prFactor*mu][k][b]*beta[l][prFactor*prC][a]*beta[m][l][prFactor*prA];

                                            if(m != c){
                                                prefactor1 = 0;
                                                prefactor2 = 0;
                                            }
                                            bool noTerms1 = false;
                                            bool noTerms2 = false;
                                            string prefactorSymbol1="";
                                            if(prefactor1 == 1){
                                                prefactorSymbol1="+";
                                            }else if(prefactor1 == -1){
                                                prefactorSymbol1="-";
                                            }else if(prefactor1 == iu){
                                                prefactorSymbol1="+i";
                                            }else if(prefactor1 == -iu){
                                                prefactorSymbol1="-i";
                                            }else if(prefactor1 == 0){
                                                noTerms1 = true;
                                            }
                                            string prefactorSymbol2="";
                                            if(prefactor2 == 1){
                                                prefactorSymbol2="+";
                                            }else if(prefactor2 == -1){
                                                prefactorSymbol2="-";
                                            }else if(prefactor2 == iu){
                                                prefactorSymbol2="+i";
                                            }else if(prefactor2 == -iu){
                                                prefactorSymbol2="-i";
                                            }else if(prefactor2 == 0){
                                                noTerms2 = true;
                                            }

                                            //BV: bare vertex, FV: full vertex
                                            string TermCodeFormatShort = prefactorSymbol1+"FV1_"+to_string(b)+to_string(c);
                                            if(!noTerms1){
                                                tempTerms += TermCodeFormatShort;
                                            }
                                            TermCodeFormatShort = prefactorSymbol2+"FV2_"+to_string(b)+to_string(c);
                                            if(!noTerms2){
                                                tempTerms += TermCodeFormatShort;
                                            }
                                        }
                                    }
                                }
                            }}}
                    if(tempTerms != ""){
                        tempTermsB += tempTerms +")*PrC_"+to_string(prC);
                        SDETerms[mu] += tempTerms +")*PrC_"+to_string(prC);
                    }
                }
                if(tempTermsB != ""){
                    tempTermsA += ")*PrB_"+to_string(prB);
                    SDETerms[mu] += ")*PrB_"+to_string(prB);
                }
            }
            if(tempTermsA != ""){
                SDETerms[mu] += ")*PrA_"+to_string(prA);
            }
        }
        SDETerms[mu] += ")*weight2*weight3";
    }

    cout << "SDE terms : " << endl
         << "0-component: " << SDETerms[0] << endl << "z-component: " << SDETerms[1] << endl << endl;


    //Output

    /*   for(int channelIndex=0;channelIndex<5;channelIndex++){
        cout << "Channel no "+to_string(channelIndex)+": " << endl;
    for(int mu=0;mu<4;mu++){
        for(int nu=0;nu<4;nu++){
            cout << "["+to_string(mu)+"]"+"["+to_string(nu)+"]: " << Channel[channelIndex][mu][nu] << ";" << endl;
        }
        }
        cout << endl;
        }
return 0;*/



    //Print out flow equation terms

    //s-channel terms
    int channelIndex=0;
    cout << "Channel no "+to_string(channelIndex)+": " << endl;
    for(int mu=0;mu<4;mu++){
        for(int nu=0;nu<4;nu++){
            if(U1Symmetric){
                if(!((mu==nu) || (mu==1 && nu==2) || (mu==2 && nu==1) || (mu==0 && nu==3) || (mu==3 && nu==0))){
                    continue;
                }
            }
            //if(mu==nu){
            cout << "addG(( "
                 << Channel[channelIndex][mu][nu]
                    << ")/(2*pi) , DG_vec,1,"+to_string(mu)+","+to_string(nu)+", R, ns,nt,nu);"
                    << endl;
            //}
        }
    }
    cout << endl;

    //RPA summations
    channelIndex=1;
    //Belongs to the RPA section
    //Print out expressions for t(mu)(nu)
    for(int mu=0;mu<4;mu++){
        for(int nu=0;nu<4;nu++){
            cout << "double t"+to_string(mu) + to_string(nu) + "i" + to_string(subl) + " = ";
            for(int b=0;b<4;b++){
                for(int c=0;c<4;c++){
                    //Important: Note the "shift" summand needs to dissapear for subl=0 and becomes "shift2" for subl=2
                    cout << " +("+ Channel[5+b*4+c][mu][nu] + ")*vertexProduct[shift+"+ to_string(64*mu+16*b+4*c+nu) + "]";
                }}
            cout << ";"<<endl;
        }}
    cout << endl;

    //Set all t(mu)(nu) to zero
    /*for(int mu=0;mu<4;mu++){
        for(int nu=0;nu<4;nu++){
            for(int a=0;a<4;a++){
                for(int b=0;b<4;b++){
          cout << "t"+to_string(mu) + to_string(nu) +"_" + to_string(a) + to_string(b) + "=0,";
                }}
        }}*/


    cout << endl;


    //t-channel terms
    channelIndex=2;
    cout << "Channel no "+to_string(channelIndex)+": " << endl;
    for(int mu=0;mu<4;mu++){
        for(int nu=0;nu<4;nu++){
            if(U1Symmetric){
                if(!((mu==nu) || (mu==1 && nu==2) || (mu==2 && nu==1) || (mu==0 && nu==3) || (mu==3 && nu==0))){
                    continue;
                }
            }
            //if(mu==nu){
            cout << "addG(( 2*t"+to_string(mu)+to_string(nu)+" "
                 << Channel[2][mu][nu] << Channel[3][mu][nu]
                    << ")/(2*pi) , DG_vec,1,"+to_string(mu)+","+to_string(nu)+", R, ns,nt,nu);"
                    << endl;
            //}
        }
    }
    cout << endl;


    //u-channel terms
    channelIndex=4;
    cout << "Channel no "+to_string(channelIndex)+": " << endl;
    for(int mu=0;mu<4;mu++){
        for(int nu=0;nu<4;nu++){
            if(U1Symmetric){
                if(!((mu==nu) || (mu==1 && nu==2) || (mu==2 && nu==1) || (mu==0 && nu==3) || (mu==3 && nu==0))){
                    continue;
                }
            }
            //if(mu==nu){
            cout << "addG(( " << Channel[channelIndex][mu][nu]
                    << ")/(2*pi) , DG_vec,1,"+to_string(mu)+","+to_string(nu)+", R, ns,nt,nu);"
                    << endl;
            //}
        }
    }
    cout << endl;



    for(int i=0; i<5;i++){
        cout << "Number of terms in channel " << i << ": " << TermCounter[i] << endl;
    }
    return 0;
}
