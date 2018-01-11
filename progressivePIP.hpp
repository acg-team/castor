/*******************************************************************************
 * Licensed Materials - Property of Lorenzo Gatti & Massimo Maiolo
 *
 *
 * Copyright (C) 2015-2017 by Lorenzo Gatti & Massimo Maiolo
 *******************************************************************************
 *
 * This file is part of miniJATI
 *
 * miniJATI is a free software: you can redistribute it and/or modify it
 * under the terms of the GNU Lesser General Public License as published
 * by the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * miniJATI is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public
 * License along with miniJATI. If not, see <http://www.gnu.org/licenses/>.
 *******************************************************************************/

#ifndef MINIJATI_PROGRESSIVEPIP_HPP
#define MINIJATI_PROGRESSIVEPIP_HPP

#include <string>
#include <map>
#include <vector>
#include <iostream>
#include <sstream>
#include <cstdlib>
#include <algorithm>
#include <Eigen/Core>
#include <Eigen/Dense>
#include <unsupported/Eigen/MatrixFunctions>
#include <random>
#include <set>
#include <chrono>

#include <Alignment.hpp>
#include <TreeRearrangment.hpp>

using namespace tshlib;

namespace progressivePIP{

#define ERR_STATE 0
#define MATCH_STATE 1
#define GAP_X_STATE 2
#define GAP_Y_STATE 3
#define STOP_STATE 4

#define MATCH_CHAR '1'
#define GAP_X_CHAR '2'
#define GAP_Y_CHAR '3'

//DP-PIP
    struct CompareFirst
    {
        CompareFirst(std::string val) : val_(val) {}
        bool operator()(const std::pair<std::string,std::string>& elem) const {
            return val_ == elem.first;
        }
    private:
        std::string val_;
    };

    const char mytable[256] = { -1, -1, -1, -1, -1, -1, -1, -1, -1,
                                -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,
                                -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,
                                4, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,
                                -1, -1, 2, -1, 1, -1, -1, -1, 3, -1, -1, -1, -1, -1, -1, -1, -1, -1,
                                -1, -1, -1, 0, 0, -1, -1, 5, -1, -1, -1, -1, -1, -1, -1, -1, 2, -1, 1,
                                -1, -1, -1, 3, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, 0, 0,
                                -1, -1, 5, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,
                                -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,
                                -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,
                                -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,
                                -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,
                                -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,
                                -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,
                                -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1 };

    const char mytableAA[256] = { -1, -1, -1, -1, -1, -1, -1, -1, -1,
                                  -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,
                                  -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,
                                  -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,
                                  -1, -1, 0, 20, 1, 2, 3, 4, 5, 6, 7, 20, 8, 9, 10, 11, 20, 12, 13, 14,
                                  15, 16, 20, 17, 18, 20, 19, 20, -1, -1, -1, -1, -1, -1, 0, 20, 1, 2, 3,
                                  4, 5, 6, 7, 20, 8, 9, 10, 11, 20, 12, 13, 14, 15, 16, 20, 17, 18, 20,
                                  19, 20, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,
                                  -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,
                                  -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,
                                  -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,
                                  -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,
                                  -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,
                                  -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,
                                  -1, -1, -1, -1, -1, -1, -1, -1, -1 };


    struct ProgressivePIPResult;

    struct ProgressivePIPResult{
        std::vector<std::pair<std::string,std::string > > MSA;
        std::vector< std::pair<std::string,std::string> > MSAs;
        std::string traceback_path;
        double score;
        //Eigen::VectorXd Pc;
        std::map<std::string,Eigen::VectorXd> fv_map;
        double lk_gap;
        double pc0;
        const bpp::Alphabet *alphabet;
        int alphabetSize;
    };
//DP-PIP
/*    void dealloc_3D_matrix(long double ***mat,int depth,int height){

        for(int k=(depth-1); k>=0; k--){
            for(int j=(height-1); j>=0; j--){
                delete[] mat[k][j];
            }
            delete[] mat[k];
        }
        delete[] mat;

        mat=NULL;
    }*/

    Eigen::VectorXd fv_observed(std::string &s,int &idx,int alphabetSize,const bpp::Alphabet *alphabet){
        int ii;

        char ch=s[idx];

        ii=alphabet->charToInt(&ch);


        Eigen::VectorXd fv = Eigen::VectorXd::Zero(alphabetSize+1);
        //TODO: Integrate mytable selection in Alignment object
        //ii=(int)mytable[(int)s[idx]];
        ii=ii<0?alphabetSize:ii;

        fv[ii]=1.0;
        idx++;

        return fv;
    }
//DP-PIP
    Eigen::VectorXd go_down(VirtualNode *tree,std::string &s,int &idx,int alphabetSize,const bpp::Alphabet *alphabet){
        Eigen::VectorXd fv;
        Eigen::VectorXd fvL;
        Eigen::VectorXd fvR;

        if(tree->isTerminalNode()){

            fv=fv_observed(s,idx,alphabetSize,alphabet);

        }else{

            fvL=go_down(tree->getNodeLeft(),s,idx,alphabetSize,alphabet);
            fvR=go_down(tree->getNodeRight(),s,idx,alphabetSize,alphabet);

            fv=(tree->getNodeLeft()->getPr()*fvL).cwiseProduct(tree->getNodeRight()->getPr()*fvR);

        }

        return fv;
    }
//DP-PIP
    void reset_corner(int &up_corner_i,int &up_corner_j,int &bot_corner_i, int &bot_corner_j,int h,int w){
        int delta;

        if(up_corner_j>=w){
            delta=up_corner_j-w+1;
            up_corner_j-=delta;
            up_corner_i+=delta;
        }
        if(bot_corner_i>=h){
            delta=bot_corner_i-h+1;
            bot_corner_i-=delta;
            bot_corner_j+=delta;
        }

    }
//DP-PIP
    std::string create_col_MSA_gap(int len){

        std::string colMSA (len,'-');

        return colMSA;
    }
//DP-PIP
    int index_of_max(double m, double x, double y,double epsilon,
                     std::default_random_engine &generator,
                     std::uniform_real_distribution<double> &distribution){
        double random_number;

        int ERR=-1;

        if(not(std::isinf(m)) & not(std::isinf(x)) & (fabs(m-x)<epsilon)){
            x=m;
        }

        if(not(std::isinf(m)) & not(std::isinf(y)) & (fabs(m-y)<epsilon)){
            y=m;
        }

        if(not(std::isinf(x)) & not(std::isinf(y)) & (fabs(x-y)<epsilon)){
            y=x;
        }

        if(m>x){
            if(m>y){
                return int(MATCH_STATE);
            }else if (y>m){
                return int(GAP_Y_STATE);
            }else{
                if(abs(m-y)<epsilon){
                    //m or y
                    random_number  = distribution(generator);
                    if(random_number < (1.0/2.0) ){
                        return int(MATCH_STATE);
                    }else{
                        return int(GAP_Y_STATE);
                    }
                }else{
                    perror("ERROR in index_of_max\n");
                    exit(EXIT_FAILURE);
                }
            }
        }else if (x>m){
            if(x>y){
                return int(GAP_X_STATE);
            }else if (y>x){
                return int(GAP_Y_STATE);
            }else{
                if(abs(x-y)<epsilon){
                    //x or y
                    random_number  = distribution(generator);
                    if(random_number < (1.0/2.0) ){
                        return int(GAP_X_STATE);
                    }else{
                        return int(GAP_Y_STATE);
                    }
                }else{
                    perror("ERROR in index_of_max_3\n");
                    exit(EXIT_FAILURE);
                }
            }
        }else{

            double mx=x;
            if(mx>y){
                //m or x
                random_number  = distribution(generator);
                if(random_number < (1.0/2.0) ){
                    return int(MATCH_STATE);
                }else{
                    return int(GAP_X_STATE);
                }
            }else if (y>mx){
                return int(GAP_Y_STATE);
            }else{
                if(abs(mx-y)<epsilon){
                    //m or x or y
                    random_number  = distribution(generator);
                    if(random_number < (1.0/3.0)){
                        return int(MATCH_STATE);
                    }else if(random_number < (2.0/3.0)){
                        return int(GAP_X_STATE);
                    }else{
                        return int(GAP_Y_STATE);
                    }
                }else{
                    perror("ERROR in index_of_max_3\n");
                    exit(EXIT_FAILURE);
                }
            }
        }

        return ERR;
    }
//DP-PIP
    double max_of_three(double a, double b, double c,double epsilon){

        //------------------------------------
        if(fabs(a)<epsilon){
            a=-INFINITY;
        }
        if(fabs(b)<epsilon){
            b=-INFINITY;
        }
        if(fabs(c)<epsilon){
            c=-INFINITY;
        }
        //------------------------------------

        if(std::isinf(a) && std::isinf(b) && std::isinf(c)){
            perror("max_of_three: all inf\n");
            exit(EXIT_FAILURE);
        }

        if(a>b){
            if(a>c){
                return a;
            }
            return c;
        }else{
            if(b>c){
                return b;
            }
            return c;
        }

    }
//DP-PIP
    double compute_nu(double tau,double lambda,double mu){

        if(fabs(mu)<1e-8){
            perror("ERROR in compute_nu: mu too small");
        }

        return lambda*(tau+1/mu);
    }
//DP-PIP
    bool is_inside(int x0,int y0,int xf,int yf,int xt,int yt){

        if((xt<x0) || (yt>y0) || (xt>xf) || (yt<yf)){

            return false;
        }

        if( (y0-yt)>(xt-x0) ){
            return false;
        }

        return true;
    }
//DP-PIP
    void set_indeces_M(int &up_corner_i,int &up_corner_j,int &bot_corner_i, int &bot_corner_j,int level,int h,int w){

        if(level==0){
            up_corner_i=0;
            up_corner_j=0;
            bot_corner_i=0;
            bot_corner_j=0;
        }else{
            up_corner_i=1+level-std::min(w-1,level);
            up_corner_j=std::min(w-1,level);
            bot_corner_i=std::min(h-1,level);
            bot_corner_j=1+level-std::min(h-1,level);
        }

    }
//DP-PIP
    void set_indeces_X(int &up_corner_i,int &up_corner_j,int &bot_corner_i, int &bot_corner_j,int level,int h,int w){

        if(level==0){
            up_corner_i=0;
            up_corner_j=0;
            bot_corner_i=0;
            bot_corner_j=0;
        }else{
            up_corner_i=1+level-1-std::min(w-1,level-1);
            up_corner_j=std::min(w-1,level-1);
            bot_corner_i=std::min(h-1,level);
            bot_corner_j=level-std::min(h-1,level);
        }

    }
//DP-PIP
    void set_indeces_Y(int &up_corner_i,int &up_corner_j,int &bot_corner_i, int &bot_corner_j,int level,int h,int w){

        if(level==0){
            up_corner_i=0;
            up_corner_j=0;
            bot_corner_i=0;
            bot_corner_j=0;
        }else{
            up_corner_i=level-std::min(w-1,level);
            up_corner_j=std::min(w-1,level);
            bot_corner_i=std::min(h-1,level-1);
            bot_corner_j=1+level-1-std::min(h-1,level-1);
        }

    }
//DP-PIP
    void set_indeces_T(int &up_corner_i,int &up_corner_j,int &bot_corner_i, int &bot_corner_j,int level,int h,int w){

        int up_corner_i_x;
        int up_corner_i_y;

        int up_corner_j_x;
        int up_corner_j_y;

        int bot_corner_i_x;
        int bot_corner_i_y;

        int bot_corner_j_x;
        int bot_corner_j_y;

        set_indeces_X(up_corner_i_x,up_corner_j_x,bot_corner_i_x,bot_corner_j_x,level,h,w);

        set_indeces_Y(up_corner_i_y,up_corner_j_y,bot_corner_i_y,bot_corner_j_y,level,h,w);

        int delta_i,delta_j;

        delta_i=bot_corner_i_x-up_corner_i_y;
        delta_j=up_corner_j_y-bot_corner_j_x;

        if(delta_i>delta_j){
            up_corner_i=up_corner_i_y;
            up_corner_j=up_corner_j_y;
            bot_corner_i=up_corner_i_y+delta_i;
            bot_corner_j=up_corner_j_y-delta_i;
        }else{
            up_corner_i=bot_corner_i_x-delta_j;
            up_corner_j=bot_corner_j_x+delta_j;
            bot_corner_i=bot_corner_i_x;
            bot_corner_j=bot_corner_j_x;
        }

    }
//DP-PIP
    int get_indices_M(int nx,int ny,int up_corner_i,int up_corner_j,int bot_corner_i,int bot_corner_j,int m,int h,int w){

        int idx;

        set_indeces_M(up_corner_i,up_corner_j,bot_corner_i,bot_corner_j,m,h,w);

        if(is_inside(up_corner_i,up_corner_j,bot_corner_i,bot_corner_j,nx,ny)){

            int dx,sx;

            dx=nx-up_corner_i+1;

            sx=((dx+1)*dx/2)-1;

            idx=sx+(ny-up_corner_j);
        }else{
            idx=-999;
        }

        return idx;

    }
//DP-PIP
    int get_indices_X(int nx,int ny,int up_corner_i,int up_corner_j,int bot_corner_i,int bot_corner_j,int m,int h,int w){

        int idx;

        set_indeces_X(up_corner_i,up_corner_j,bot_corner_i,bot_corner_j,m,h,w);

        if(is_inside(up_corner_i,up_corner_j,bot_corner_i,bot_corner_j,nx,ny)){

            int dx,sx;

            dx=nx-up_corner_i+1;

            sx=((dx+1)*dx/2)-1;

            idx=sx+(ny-up_corner_j);
        }else{
            idx=-999;
        }

        return idx;

    }
//DP-PIP
    int get_indices_Y(int nx,int ny,int up_corner_i,int up_corner_j,int bot_corner_i,int bot_corner_j,int m,int h,int w){

        int idx;

        set_indeces_Y(up_corner_i,up_corner_j,bot_corner_i,bot_corner_j,m,h,w);

        if(is_inside(up_corner_i,up_corner_j,bot_corner_i,bot_corner_j,nx,ny)){

            int dx,sx;

            dx=nx-up_corner_i+1;

            sx=((dx+1)*dx/2)-1;

            idx=sx+(ny-up_corner_j);
        }else{
            idx=-999;
        }

        return idx;

    }
//DP-PIP
    int get_indices_T(int nx,int ny,int up_corner_i,int up_corner_j,int bot_corner_i,int bot_corner_j,int m,int h,int w){

        set_indeces_T(up_corner_i,up_corner_j,bot_corner_i,bot_corner_j,m,h,w);

        reset_corner(up_corner_i,up_corner_j,bot_corner_i,bot_corner_j,h,w);

        int idx;
        int dx,sx;

        dx=nx-up_corner_i+1;

        sx=((dx+1)*dx/2)-1;

        idx=sx+(ny-up_corner_j);

        return idx;

    }
//DP-PIP
    void allgaps(VirtualNode *tree,std::string &s,int &idx,bool &flag){

        if(tree->isTerminalNode()){
            char ch=s[idx];

            idx++;

            if(ch!='-'){
                flag=false;
            }

        }else{
            allgaps(tree->getNodeLeft(),s,idx,flag);
            allgaps(tree->getNodeRight(),s,idx,flag);
        }

    }
//DP-PIP
    double compute_lk_gap_down(VirtualNode *tree,std::string &s,Eigen::VectorXd const &pi,int alphabetSize,const bpp::Alphabet *alphabet){

        double pr=0;
        double pL=0;
        double pR=0;
        int idx;
        Eigen::VectorXd fvL;
        Eigen::VectorXd fvR;
        Eigen::VectorXd fv;
        double fv0;

        if(tree->isTerminalNode()){
            idx=0;
            fv=go_down(tree,s,idx,alphabetSize,alphabet);
            fv0=fv.dot(pi);
            pr=tree->getIota()-tree->getIota()*tree->getBeta()+tree->getIota()*tree->getBeta()*fv0;

            return pr;
        }else{
            idx=0;
            fv=go_down(tree,s,idx,alphabetSize,alphabet);
            fv0=fv.dot(pi);
            pr=tree->getIota()-tree->getIota()*tree->getBeta()+tree->getIota()*tree->getBeta()*fv0;

            bool flagL=true;
            bool flagR=true;
            idx=0;
            allgaps(tree->getNodeLeft(),s,idx,flagL);
            int ixx=idx;
            allgaps(tree->getNodeRight(),s,idx,flagR);
            int len;

            std::string sL;//=stringFromSequence(s);
            len=ixx;
            sL=s.substr(0,len);
            pL=compute_lk_gap_down(tree->getNodeLeft(),sL,pi,alphabetSize,alphabet);

            std::string sR;//=stringFromSequence(s);
            sR=s.substr(ixx);
            pR=compute_lk_gap_down(tree->getNodeRight(),sR,pi,alphabetSize,alphabet);
        }

        return pr+pL+pR;
    }
//DP-PIP
    bool checkboundary(int up_corner_i,int up_corner_j,int bot_corner_i,int bot_corner_j,int h,int w){

        if( (up_corner_i  >=0) & (up_corner_i  <h) &\
	   (up_corner_j  >=0) & (up_corner_j  <w) &\
	   (bot_corner_i >=0) & (bot_corner_i <h) &\
	   (bot_corner_j >=0) & (bot_corner_j <w)){
            return true;
        }

        return false;
    }
//DP-PIP
    std::string create_col_MSA(std::vector<std::pair<std::string,std::string>> &result,int index){
        std::string colMSA;

        for(unsigned int i=0;i<result.size();i++){
            colMSA.append(result.at(i).second,index,1);
        }

        return colMSA;
    }
//DP-PIP
    std::vector<std::pair<std::string,std::string>> align_seq_left(	std::vector<std::pair<std::string,std::string>> &MSA_in,
                                                                       std::string &traceback_path){

        unsigned int idx;

        std::vector<std::pair<std::string,std::string>> MSA_out;

        typedef typename std::vector<std::pair<std::string,std::string>>::iterator vect_iterator;
        for (vect_iterator iter = MSA_in.begin(); iter != MSA_in.end(); iter++){

            std::pair<std::string,std::string> seq = (*iter);

            std::string seq_name=seq.first;
            std::string seq_not_aligned=seq.second;
            std::string seq_aligned(traceback_path.size(),'-');

            idx=0;
            for(unsigned int j=0;j<traceback_path.size();j++){

                if(traceback_path.at(j)==MATCH_CHAR){

                    seq_aligned[j]=seq_not_aligned.at(idx);
                    idx++;

                }else if(traceback_path.at(j)==GAP_X_CHAR){

                    seq_aligned[j]=seq_not_aligned.at(idx);
                    idx++;

                }else if(traceback_path.at(j)==GAP_Y_CHAR){

                    seq_aligned[j]='-';

                }else{
                    perror("ERROR in align_seq_left");
                }
            }

            MSA_out.push_back(std::make_pair(seq_name,seq_aligned));
        }

        return MSA_out;
    }
//DP-PIP
    std::vector<std::pair<std::string,std::string>> align_seq_right(std::vector<std::pair<std::string,std::string>> &result,
                                                                    std::string &traceback_path){
        unsigned int idx;

        std::vector<std::pair<std::string,std::string>> MSA_out;

        typedef typename std::vector<std::pair<std::string,std::string>>::iterator vect_iterator;
        for (vect_iterator iter = result.begin(); iter != result.end(); iter++){

            std::pair<std::string,std::string> seq = (*iter);

            std::string seq_name=seq.first;
            std::string seq_not_aligned=seq.second;
            std::string seq_aligned(traceback_path.size(),'-');

            idx=0;
            for(unsigned int j=0;j<traceback_path.size();j++){

                if(traceback_path.at(j)==MATCH_CHAR){

                    seq_aligned[j]=seq_not_aligned.at(idx);
                    idx++;

                }else if(traceback_path.at(j)==GAP_X_CHAR){

                    seq_aligned[j]='-';

                }else if(traceback_path.at(j)==GAP_Y_CHAR){

                    seq_aligned[j]=seq_not_aligned.at(idx);
                    idx++;

                }else{
                    perror("ERROR in align_seq");
                }
            }

            MSA_out.push_back(std::make_pair(seq_name,seq_aligned));
        }

        return MSA_out;
    }
//DP-PIP
    static std::vector<std::pair<std::string,std::string>> build_MSA(std::string traceback_path,std::vector<std::pair<std::string,std::string>> &MSA_L,std::vector<std::pair<std::string,std::string>> &MSA_R){
        std::vector<std::pair<std::string,std::string>> MSA;
        std::vector<std::pair<std::string,std::string>> MSA_L_out;
        std::vector<std::pair<std::string,std::string>> MSA_R_out;

        MSA_L_out=align_seq_left(MSA_L,traceback_path);
        MSA_R_out=align_seq_right(MSA_R,traceback_path);

        for (unsigned int i=0;i<MSA_L_out.size();i++){
            MSA.push_back(MSA_L_out.at(i));
        }

        for (unsigned int i=0;i<MSA_R_out.size();i++){
            MSA.push_back(MSA_R_out.at(i));
        }

        return MSA;
    }
//DP-PIP
    double compute_lk_down(VirtualNode *tree,std::string &s,Eigen::VectorXd &pi,int alphabetSize,const bpp::Alphabet *alphabet){

        double pr;
        int idx;
        Eigen::VectorXd fvL;
        Eigen::VectorXd fvR;
        Eigen::VectorXd fv;
        double fv0;

        if(tree->isTerminalNode()){

            idx=0;
            fv=go_down(tree,s,idx,alphabetSize,alphabet);
            fv0=fv.dot(pi);
            pr=tree->getIota()*tree->getBeta()*fv0;

            return pr;

        }else{

            idx=0;
            fv=go_down(tree,s,idx,alphabetSize,alphabet);
            fv0=fv.dot(pi);
            pr=tree->getIota()*tree->getBeta()*fv0;

            bool flagL=true;
            bool flagR=true;
            idx=0;
            allgaps(tree->getNodeLeft(),s,idx,flagL);
            int ixx=idx;
            allgaps(tree->getNodeRight(),s,idx,flagR);

            int len;
            if(flagR){
                std::string sL;//=stringFromSequence(s);
                len=ixx;
                sL=s.substr(0,len);
                return pr + compute_lk_down(tree->getNodeLeft(),sL,pi,alphabetSize,alphabet);
            }

            if(flagL){
                std::string sR;//=stringFromSequence(s);
                sR=s.substr(ixx);
                return pr + compute_lk_down(tree->getNodeRight(),sR,pi,alphabetSize,alphabet);
            }

        }

        return pr;
    }
//DP-PIP
    double compute_pr_gap_all_edges_s(	  VirtualNode *tree,
                                            std::string &sL,
                                            std::string &sR,
                                            Eigen::VectorXd &pi,
                                            int alphabetSize,
                                            const bpp::Alphabet *alphabet){

        double fv0;
        double pr;
        Eigen::VectorXd fvL;
        Eigen::VectorXd fvR;
        Eigen::VectorXd fv;
        int idx;

        idx=0;
        fvL=go_down(tree->getNodeLeft(),sL,idx,alphabetSize,alphabet);

        idx=0;
        fvR=go_down(tree->getNodeRight(),sR,idx,alphabetSize,alphabet);

        fv=(tree->getNodeLeft()->getPr()*fvL).cwiseProduct(tree->getNodeRight()->getPr()*fvR);

        fv0=fv.dot(pi);

        if(tree->isRootNode()){
            pr=(tree->getIota()*fv0);
        }else{
            pr=(tree->getIota() - tree->getIota()*tree->getBeta() + tree->getIota()*tree->getBeta()*fv0);
        }


        double pL,pR;
        pL=compute_lk_gap_down(tree->getNodeLeft(),sL,pi,alphabetSize,alphabet);
        pR=compute_lk_gap_down(tree->getNodeRight(),sR,pi,alphabetSize,alphabet);

        pr=pr+pL+pR;

        //*************************************************************************************************
        VirtualNode *p_tree=tree;

        //TODO:controllare se torna null
        while(p_tree->getNodeUp()){

            fv=(p_tree->getPr()*fv);

            fv0=fv.dot(pi);

            if(p_tree->getNodeUp()->getNodeUp()==NULL){
                pr+=(p_tree->getNodeUp()->getIota()*fv0);
            }else{
                //TODO:controllare parentesi
                pr+=(p_tree->getNodeUp()->getIota() -
                     p_tree->getNodeUp()->getIota()*p_tree->getNodeUp()->getBeta() +
                     p_tree->getNodeUp()->getIota()*p_tree->getNodeUp()->getBeta()*fv0);
            }

            p_tree=p_tree->getNodeUp();
        }
        //*************************************************************************************************

        return pr;
    }
//DP-PIP
    double compute_pr_gap_local_tree_s(VirtualNode *tree,
                                       std::string &sL,
                                       std::string &sR,
                                       Eigen::VectorXd &pi,
                                       int alphabetSize,
                                       const bpp::Alphabet *alphabet){

        double fv0;
        double pr;
        Eigen::VectorXd fvL;
        Eigen::VectorXd fvR;
        Eigen::VectorXd fv;
        int idx;

        idx=0;
        fvL=go_down(tree->getNodeLeft(),sL,idx,alphabetSize,alphabet);

        idx=0;
        fvR=go_down(tree->getNodeRight(),sR,idx,alphabetSize,alphabet);

        fv=(tree->getNodeLeft()->getPr()*fvL).cwiseProduct(tree->getNodeRight()->getPr()*fvR);

        fv0=fv.dot(pi);

        pr=tree->getIota()*fv0;

        double pL,pR;
        pL=compute_lk_gap_down(tree->getNodeLeft(),sL,pi,alphabetSize,alphabet);
        pR=compute_lk_gap_down(tree->getNodeRight(),sR,pi,alphabetSize,alphabet);

        pr=pr+pL+pR;

        return pr;
    }
//DP-PIP
    double computeLK_M_all_edges_s_opt(double valM,
                                       double valX,
                                       double valY,
                                       double nu,
                                       VirtualNode *tree,
                                       std::string &sL,
                                       std::string &sR,
                                       Eigen::VectorXd &pi,
                                       int m,
                                       std::map<std::string,double> &lkM,
                                       int alphabetSize,
                                       const bpp::Alphabet *alphabet){


        double pr;
        double val;
        //------------------------------------------------------------------------------------------
        typedef typename std::map<std::string,double>::iterator MapIterator;
        std::string s;

        s.append(sL);
        s.append(sR);

        MapIterator it=lkM.find(s);
        if(it == lkM.end()){

            //------------------------------------------------------------------------------------------
            //DOWN
            Eigen::VectorXd fvL;
            Eigen::VectorXd fvR;
            Eigen::VectorXd fv;
            double fv0;
            int idx;

            idx=0;
            fvL=go_down(tree->getNodeLeft(),sL,idx,alphabetSize,alphabet);
            idx=0;
            fvR=go_down(tree->getNodeRight(),sR,idx,alphabetSize,alphabet);

            fv=(tree->getNodeLeft()->getPr()*fvL).cwiseProduct(tree->getNodeRight()->getPr()*fvR);

            fv0=fv.dot(pi);

            pr=tree->getIota()*tree->getBeta()*fv0;
            //------------------------------------------------------------------------------------------
            //UP
            VirtualNode *p_tree=tree;

            //TODO:controllare se torna null
            while(p_tree->getNodeUp()!=NULL){

                fv=(p_tree->getPr()*fv);

                fv0=fv.dot(pi);

                pr+=(p_tree->getNodeUp()->getIota()*p_tree->getNodeUp()->getBeta()*fv0);

                p_tree=p_tree->getNodeUp();
            }
            //------------------------------------------------------------------------------------------

            pr=log(pr);

            lkM[s]=pr;

        }else{
            pr=it->second;
        }

        val=-log(double(m))+log(nu)+pr+max_of_three(valM,valX,valY,(double)DBL_EPSILON);

        return val;
    }
//DP-PIP
    double computeLK_M_local_tree_s_opt(double valM,
                                        double valX,
                                        double valY,
                                        double nu,
                                        VirtualNode *tree,
                                        std::string &sL,
                                        std::string &sR,
                                        Eigen::VectorXd &pi,
                                        int m,
                                        std::map<std::string,double> &lkM,
                                        int alphabetSize,
                                        const bpp::Alphabet *alphabet){


        double pr;
        double val;
        //------------------------------------------------------------------------------------------
        typedef typename std::map<std::string,double>::iterator MapIterator;
        std::string s;
        Eigen::VectorXd fvL;
        Eigen::VectorXd fvR;
        Eigen::VectorXd fv;
        double fv0;
        int idx;

        s.append(sL);
        s.append(sR);

        MapIterator it=lkM.find(s);
        if(it == lkM.end()){

            //------------------------------------------------------------------------------------------


            idx=0;
            fvL=go_down(tree->getNodeLeft(),sL,idx,alphabetSize,alphabet);
            idx=0;
            fvR=go_down(tree->getNodeRight(),sR,idx,alphabetSize,alphabet);

            fv=(tree->getNodeLeft()->getPr()*fvL).cwiseProduct(tree->getNodeRight()->getPr()*fvR);

            fv0=fv.dot(pi);

            pr=tree->getIota()*tree->getBeta()*fv0;

            //------------------------------------------------------------------------------------------

            pr=log(pr);

            lkM[s]=pr;

        }else{
            pr=it->second;
        }
        //------------------------------------------------------------------------------------------

#ifdef VERBOSE
        std::cout<<"prM("<<sL<<":"<<sR<<")"; printf(" %18.16lf\n",pr);
#endif

        val=-log(double(m))+log(nu)+pr+max_of_three(valM,valX,valY,(double)DBL_EPSILON);

        return val;
    }
//DP-PIP
    double computeLK_X_all_edges_s_opt(double valM,
                                       double valX,
                                       double valY,
                                       double nu,
                                       VirtualNode *tree,
                                       std::string &sL,
                                       std::string &col_gap_R,
                                       Eigen::VectorXd &pi,
                                       int m,
                                       std::map<std::string,double> &lkX,
                                       int alphabetSize,
                                       const bpp::Alphabet *alphabet){


        double pr;
        double val;
        //------------------------------------------------------------------------------------------
        typedef typename std::map<std::string,double>::iterator MapIterator;
        std::string s;

        s.append(sL);
        s.append(col_gap_R);

        MapIterator it=lkX.find(s);
        if(it == lkX.end()){

            //------------------------------------------------------------------------------------------
            Eigen::VectorXd fvL;
            Eigen::VectorXd fvR;
            Eigen::VectorXd fv;
            int idx;
            double fv0;

            idx=0;
            fvL=go_down(tree->getNodeLeft(),sL,idx,alphabetSize,alphabet);

            idx=0;
            fvR=go_down(tree->getNodeRight(),col_gap_R,idx,alphabetSize,alphabet);

            fv=(tree->getNodeLeft()->getPr()*fvL).cwiseProduct(tree->getNodeRight()->getPr()*fvR);
            fv0=fv.dot(pi);

            pr=tree->getIota()*tree->getBeta()*fv0;

            double pL;

            //------------------------------------------------------------------------------------------
            pL=compute_lk_down(tree->getNodeLeft(),sL,pi,alphabetSize,alphabet);
            //------------------------------------------------------------------------------------------

            pr+=pL;
            //------------------------------------------------------------------------------------------


            //******************************************************************************************
            VirtualNode *p_tree=tree;

            //TODO: controllare se torna null
            while(p_tree->getNodeUp()!=NULL){

                fv=(p_tree->getPr()*fv);

                fv0=fv.dot(pi);

                pr+=(p_tree->getNodeUp()->getIota()*p_tree->getNodeUp()->getBeta()*fv0);

                p_tree=p_tree->getNodeUp();
            }
            //******************************************************************************************

            pr=log(pr);

            lkX[s]=pr;

        }else{
            pr=it->second;
        }
        //------------------------------------------------------------------------------------------

        val=-log(double(m))+log(nu)+pr+max_of_three(valM,valX,valY,(double)DBL_EPSILON);

        return val;
    }
//DP-PIP
    double computeLK_X_local_tree_s_opt(double valM,
                                        double valX,
                                        double valY,
                                        double nu,
                                        VirtualNode *tree,
                                        std::string &sL,
                                        std::string &col_gap_R,
                                        Eigen::VectorXd &pi,
                                        int m,
                                        std::map<std::string,double> &lkX,
                                        int alphabetSize,
                                        const bpp::Alphabet *alphabet){


        double pr;
        double val;
        //------------------------------------------------------------------------------------------
        typedef typename std::map<std::string,double>::iterator MapIterator;
        std::string s;
        Eigen::VectorXd fvL;
        Eigen::VectorXd fvR;
        Eigen::VectorXd fv;
        int idx;
        double fv0;

        s.append(sL);
        s.append(col_gap_R);

        MapIterator it=lkX.find(s);
        if(it == lkX.end()){

            //------------------------------------------------------------------------------------------


            idx=0;
            fvL=go_down(tree->getNodeLeft(),sL,idx,alphabetSize,alphabet);

            idx=0;
            fvR=go_down(tree->getNodeRight(),col_gap_R,idx,alphabetSize,alphabet);

            fv=(tree->getNodeLeft()->getPr()*fvL).cwiseProduct(tree->getNodeRight()->getPr()*fvR);
            fv0=fv.dot(pi);

            pr=tree->getIota()*tree->getBeta()*fv0;

            double pL;

            //------------------------------------------------------------------------------------------
            pL=compute_lk_down(tree->getNodeLeft(),sL,pi,alphabetSize,alphabet);
            //------------------------------------------------------------------------------------------

            pr+=pL;
            //------------------------------------------------------------------------------------------

            pr=log(pr);

            lkX[s]=pr;

        }else{
            pr=it->second;
        }
        //------------------------------------------------------------------------------------------

        val=-log(double(m))+log(nu)+pr+max_of_three(valM,valX,valY,(double)DBL_EPSILON);

        return val;
    }
//DP-PIP
    double computeLK_Y_all_edges_s_opt(double valM,
                                       double valX,
                                       double valY,
                                       double nu,
                                       VirtualNode *tree,
                                       std::string &col_gap_L,
                                       std::string &sR,
                                       Eigen::VectorXd &pi,
                                       int m,
                                       std::map<std::string,double> &lkY,
                                       int alphabetSize,
                                       const bpp::Alphabet *alphabet){


        double pr;
        double val;
        //------------------------------------------------------------------------------------------
        typedef typename std::map<std::string,double>::iterator MapIterator;
        std::string s;

        s.append(col_gap_L);
        s.append(sR);

        MapIterator it=lkY.find(s);
        if(it == lkY.end()){

            //------------------------------------------------------------------------------------------
            Eigen::VectorXd fvL;
            Eigen::VectorXd fvR;
            Eigen::VectorXd fv;
            double fv0;
            int idx=0;
            fvL=go_down(tree->getNodeLeft(),col_gap_L,idx,alphabetSize,alphabet);

            idx=0;
            fvR=go_down(tree->getNodeRight(),sR,idx,alphabetSize,alphabet);

            fv=(tree->getNodeLeft()->getPr()*fvL).cwiseProduct(tree->getNodeRight()->getPr()*fvR);
            //------------------------------------------------------------------------------------------

            fv0=fv.dot(pi);

            pr=tree->getIota()*tree->getBeta()*fv0;

            double pR;

            //------------------------------------------------------------------------------------------
            pR=compute_lk_down(tree->getNodeRight(),sR,pi,alphabetSize,alphabet);
            //------------------------------------------------------------------------------------------

            pr+=pR;

            //*******************************************************************************************
            VirtualNode *p_tree=tree;

            //TODO:controllare che torni null
            while(p_tree->getNodeUp()!=NULL){

                fv=(p_tree->getPr()*fv);

                fv0=fv.dot(pi);

                pr+=(p_tree->getNodeUp()->getIota()*p_tree->getNodeUp()->getBeta()*fv0);

                p_tree=p_tree->getNodeUp();
            }
            //*******************************************************************************************

            pr=log(pr);

            lkY[s]=pr;

        }else{
            pr=it->second;
        }
        //------------------------------------------------------------------------------------------

        val=-log(double(m))+log(nu)+pr+max_of_three(valM,valX,valY,(double)DBL_EPSILON);

        return val;
    }
//DP-PIP
    double computeLK_Y_local_tree_s_opt(double valM,
                                        double valX,
                                        double valY,
                                        double nu,
                                        VirtualNode *tree,
                                        std::string &col_gap_L,
                                        std::string &sR,
                                        Eigen::VectorXd &pi,
                                        int m,
                                        std::map<std::string,double> &lkY,
                                        int alphabetSize,
                                        const bpp::Alphabet *alphabet){


        double pr;
        double val;
        //------------------------------------------------------------------------------------------
        typedef typename std::map<std::string,double>::iterator MapIterator;
        std::string s;
        Eigen::VectorXd fvL;
        Eigen::VectorXd fvR;
        Eigen::VectorXd fv;
        double fv0;
        int idx=0;

        s.append(col_gap_L);
        s.append(sR);

        MapIterator it=lkY.find(s);
        if(it == lkY.end()){

            //------------------------------------------------------------------------------------------

            fvL=go_down(tree->getNodeLeft(),col_gap_L,idx,alphabetSize,alphabet);

            idx=0;
            fvR=go_down(tree->getNodeRight(),sR,idx,alphabetSize,alphabet);

            fv=(tree->getNodeLeft()->getPr()*fvL).cwiseProduct(tree->getNodeRight()->getPr()*fvR);
            //------------------------------------------------------------------------------------------

            fv0=fv.dot(pi);

            pr=tree->getIota()*tree->getBeta()*fv0;

            double pR;

            //------------------------------------------------------------------------------------------
            pR=compute_lk_down(tree->getNodeRight(),sR,pi,alphabetSize,alphabet);
            //------------------------------------------------------------------------------------------

            pr+=pR;

            pr=log(pr);

            lkY[s]=pr;

        }else{
            pr=it->second;
        }
        //------------------------------------------------------------------------------------------


        val=-log(double(m))+log(nu)+pr+max_of_three(valM,valX,valY,(double)DBL_EPSILON);


        return val;
    }
//DP-PIP
/*    void fill_scores(double lk,double &max_lk,double &prev_max_lk,int &level_max_lk,
                     int &last_d,int m,bool &flag_exit,double *scores,int &counter,
                     bool CENTER,int num_subopt,int &index0){

        if (lk>max_lk){
            prev_max_lk=max_lk;
            max_lk=lk;
            level_max_lk=m;

            if (CENTER){
                if(std::isinf(prev_max_lk)){
                    scores[0]=max_lk;
                    counter=1;
                    index0=m;
                }else{
                    if(num_subopt>2){
                        scores[0]=prev_max_lk;
                        scores[1]=max_lk;
                        counter=2;
                        index0=m-1;
                    }else{
                        scores[0]=max_lk;
                        counter=1;
                        index0=m;
                    }
                }
            }else{
                scores[0]=max_lk;
                counter=1;
                index0=m;
            }

        }else{

            if (counter>=num_subopt){
                flag_exit=true;
                last_d=m;
                return;
            }

            scores[counter]=lk;
            counter=counter+1;

        }

    }*/
//DP-PIP
    bool check_uniform_len_s(std::vector<std::pair<std::string,std::string>> &result){
        unsigned int len;

        if(result.size()==0){
            return false;
        }

        len=result.at(0).second.length();

        for(unsigned int i=1;i<result.size();i++){
            if(len!=result.at(i).second.length()){
                return false;
            }
        }

        return true;
    }
//DP-PIP
    int get_length_seq_s(std::vector<std::pair<std::string,std::string>> &result){

        if(result.size()==0){
            return 0;
        }

        if(!check_uniform_len_s(result)){
            perror("ERROR in get_length_seq: non aligned");
        }

        return result.at(0).second.length();
    }
//DP-PIP
    static ProgressivePIPResult compute_DP3D_PIP(ProgressivePIPResult &result_L,
                                                 ProgressivePIPResult &result_R,
                                                 VirtualNode *tree,
                                                 Likelihood *likelihood,
                                                 Alignment *alignment,
                                                 const bpp::Alphabet *alphabet,
                                                 bpp::SiteContainer *sites,
                                                 double gamma_rate,
                                                 bool local){


        //sites->getSite(0).getAlphabet()->charToInt("A");



        ProgressivePIPResult result;


        double tau;
        double nu;
        //int alphabetSize;
        //char *mapping_table; //alignment.table;


        result.alphabet=alphabet;
        result.alphabetSize=alphabet->getSize();

        //alphabetSize=alignment->align_alphabetsize;

        //@gamma_distribution
        double lambda_gamma=likelihood->lambda*gamma_rate;
        double mu_gamma=likelihood->mu*gamma_rate;

        if(local){
            //TODO: calcola tau dal nodo attuale
            tau=tree->computeTotalTreeLength();
            nu=compute_nu(tau,lambda_gamma,mu_gamma);
            //TODO: tau and nu are already set?
            //tree->set_tau(tau);
            //tree_>set_nu(nu);
            //@gamma_distribution
            //tree.set_iota_local(tau,mu_gamma);
            //tree.set_beta_local(tau,mu_gamma);
        }else{
            tau=likelihood->tau;
            nu=likelihood->nu;
        }


        int up_corner_i;
        int up_corner_j;
        int bot_corner_i;
        int bot_corner_j;
        int lw;
        int h,w;

        h=get_length_seq_s(result_L.MSAs)+1;
        w=get_length_seq_s(result_R.MSAs)+1;
        int d=(h-1)+(w-1)+1;

        //TODO: vectorXd to matrixXd
        Eigen::VectorXd pi = likelihood->pi;

        double pc0;
        //.-----.------.------.-------.------.-------.-------.--------.-------.--------.//
        std::string sLs;
        std::string sRs;
        std::string col_gap_Ls;
        std::string col_gap_Rs;
        //.-----.------.------.-------.------.-------.-------.--------.-------.--------.//
        col_gap_Ls=create_col_MSA_gap(result_L.MSAs.size());
        col_gap_Rs=create_col_MSA_gap(result_R.MSAs.size());
        //.-----.------.------.-------.------.-------.-------.--------.-------.--------.//
//	std::cout<<"random generator ON\n";
        unsigned int seed = std::chrono::system_clock::now().time_since_epoch().count();
        std::default_random_engine generator(seed);
        std::uniform_real_distribution<double> distribution(0.0,1.0);
        //.-----.------.------.//
//	std::cout<<"random generator OFF\n";
//	unsigned seed = 0;
//	std::default_random_engine generator(seed);
//	std::uniform_real_distribution<double> distribution(0.0,1.0);
        //.-----.------.------.-------.------.-------.-------.--------.-------.--------.//

        double epsilon=DBL_EPSILON;

        //***************************************************************************************
        //***************************************************************************************
        if(local){
            pc0=compute_pr_gap_local_tree_s(tree,
                                            col_gap_Ls,
                                            col_gap_Rs,
                                            pi,
                                            result.alphabetSize,alphabet);
        }else{
            pc0=compute_pr_gap_all_edges_s(tree,
                                           col_gap_Ls,
                                           col_gap_Rs,
                                           pi,
                                           result.alphabetSize,alphabet);
        }

        //***************************************************************************************
        //***************************************************************************************

        double** LogM = new double*[2];
        double** LogX = new double*[2];
        double** LogY = new double*[2];

        int** TR = new int*[d];

        LogM[0] = new double[int((w*(h+1))/2)];
        LogX[0] = new double[int((w*(h+1))/2)];
        LogY[0] = new double[int((w*(h+1))/2)];
        LogM[1] = new double[int((w*(h+1))/2)];
        LogX[1] = new double[int((w*(h+1))/2)];
        LogY[1] = new double[int((w*(h+1))/2)];

        LogM[0][0]=nu*(pc0-1.0);
        LogX[0][0]=nu*(pc0-1.0);
        LogY[0][0]=nu*(pc0-1.0);

        TR[0] = new int[1]();
        TR[0][0]=STOP_STATE;

        double max_of_3;
        double max_lk=-INFINITY;
        double prev_max_lk=-INFINITY;
        int level_max_lk=INT_MIN;
        double val;
        int m_binary_this;
        int m_binary_prev;

        double valM;
        double valX;
        double valY;

        int idx;

        int coordSeq_1;
        int coordSeq_2;
        int coordTriangle_this_i;
        int coordTriangle_this_j;
        int coordTriangle_prev_i;
        int coordTriangle_prev_j;

        int counter;

        //double *scores = new double[num_subopt];
        double scores;// = new double[num_subopt];
        //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        //bool CENTER = true;
        //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        int start_depth;
        int depth;


        bool flag_exit=false;
        int last_d=d-1;
        int size_tr,tr_up_i,tr_up_j,tr_down_i,tr_down_j;
        std::map<std::string,double> lkM;
        std::map<std::string,double> lkX;
        std::map<std::string,double> lkY;

        counter=0;
        for(int m=1;m<d;m++){

            if(flag_exit){
                break;
            }

            m_binary_this=m%2;
            m_binary_prev=(m+1)%2;
            //***************************************************************************************
            //***************************************************************************************
            set_indeces_M(up_corner_i,up_corner_j,bot_corner_i,bot_corner_j,m,h,w);

            if(checkboundary(up_corner_i,up_corner_j,bot_corner_i,bot_corner_j,h,w)){

                lw=0;
                for(int i=up_corner_i;i<=bot_corner_i;i++){

                    coordTriangle_this_i=i;
                    coordSeq_1=coordTriangle_this_i-1;
                    coordTriangle_prev_i=coordTriangle_this_i-1;
                    sLs=create_col_MSA(result_L.MSAs,coordSeq_1);


                    for(int j=0;j<=lw;j++){

                        coordTriangle_this_j=up_corner_j-j;
                        coordSeq_2=coordTriangle_this_j-1;
                        coordTriangle_prev_j=coordTriangle_this_j-1;
                        sRs=create_col_MSA(result_R.MSAs,coordSeq_2);

                        idx=get_indices_M(coordTriangle_prev_i,coordTriangle_prev_j,up_corner_i,up_corner_j,bot_corner_i,bot_corner_j,m-1,h,w);
                        if(idx>=0){
                            valM=LogM[m_binary_prev][idx];
                        }else{
                            valM=-INFINITY;
                        }

                        idx=get_indices_X(coordTriangle_prev_i,coordTriangle_prev_j,up_corner_i,up_corner_j,bot_corner_i,bot_corner_j,m-1,h,w);
                        if(idx>=0){
                            valX=LogX[m_binary_prev][idx];
                        }else{
                            valX=-INFINITY;
                        }

                        idx=get_indices_Y(coordTriangle_prev_i,coordTriangle_prev_j,up_corner_i,up_corner_j,bot_corner_i,bot_corner_j,m-1,h,w);
                        if(idx>=0){
                            valY=LogY[m_binary_prev][idx];
                        }else{
                            valY=-INFINITY;
                        }

                        if(std::isinf(valM) && std::isinf(valX) && std::isinf(valY)){
                            exit(EXIT_FAILURE);
                        }

                        if(local){
                            val=computeLK_M_local_tree_s_opt(	valM,
                                                                 valX,
                                                                 valY,
                                                                 nu,
                                                                 tree,
                                                                 sLs,sRs,
                                                                 pi,
                                                                 m,
                                                                 lkM,
                                                                 result.alphabetSize,alphabet);
                        }else{
                            val=computeLK_M_all_edges_s_opt(	valM,
                                                                valX,
                                                                valY,
                                                                nu,
                                                                tree,
                                                                sLs,sRs,
                                                                pi,
                                                                m,
                                                                lkM,
                                                                result.alphabetSize,alphabet);
                        }

                        if(std::isinf(val)){
                            exit(EXIT_FAILURE);
                        }

                        if(std::isnan(val)){
                            exit(EXIT_FAILURE);
                        }

                        idx=get_indices_M(coordTriangle_this_i,coordTriangle_this_j,up_corner_i,
                                          up_corner_j,bot_corner_i,bot_corner_j,m,h,w);

                        LogM[m_binary_this][idx]=val;
                    }
                    lw++;
                }

            }
            //***************************************************************************************
            //***************************************************************************************
            set_indeces_X(up_corner_i,up_corner_j,bot_corner_i,bot_corner_j,m,h,w);
            tr_down_i=bot_corner_i;
            tr_down_j=bot_corner_j;
            if(checkboundary(up_corner_i,up_corner_j,bot_corner_i,bot_corner_j,h,w)){

                lw=0;
                for(int i=up_corner_i;i<=bot_corner_i;i++){

                    coordTriangle_this_i=i;
                    coordTriangle_prev_i=coordTriangle_this_i-1;
                    coordSeq_1=coordTriangle_this_i-1;
                    sLs=create_col_MSA(result_L.MSAs,coordSeq_1);

                    for(int j=0;j<=lw;j++){

                        coordTriangle_this_j=up_corner_j-j;
                        coordTriangle_prev_j=coordTriangle_this_j;

                        idx=get_indices_M(coordTriangle_prev_i,coordTriangle_prev_j,up_corner_i,
                                          up_corner_j,bot_corner_i,bot_corner_j,m-1,h,w);
                        if(idx>=0){
                            valM=LogM[m_binary_prev][idx];
                        }else{
                            valM=-INFINITY;
                        }

                        idx=get_indices_X(coordTriangle_prev_i,coordTriangle_prev_j,up_corner_i,
                                          up_corner_j,bot_corner_i,bot_corner_j,m-1,h,w);
                        if(idx>=0){
                            valX=LogX[m_binary_prev][idx];
                        }else{
                            valX=-INFINITY;
                        }

                        idx=get_indices_Y(coordTriangle_prev_i,coordTriangle_prev_j,up_corner_i,
                                          up_corner_j,bot_corner_i,bot_corner_j,m-1,h,w);
                        if(idx>=0){
                            valY=LogY[m_binary_prev][idx];
                        }else{
                            valY=-INFINITY;
                        }

                        if(std::isinf(valM) && std::isinf(valX) && std::isinf(valY)){
                            exit(EXIT_FAILURE);
                        }

                        if(local){
                            val=computeLK_X_local_tree_s_opt(valM,
                                                             valX,
                                                             valY,
                                                             nu,
                                                             tree,
                                                             sLs,col_gap_Rs,
                                                             pi,
                                                             m,
                                                             lkX,
                                                             result.alphabetSize,alphabet);
                        }else{
                            val=computeLK_X_all_edges_s_opt(valM,
                                                            valX,
                                                            valY,
                                                            nu,
                                                            tree,
                                                            sLs,col_gap_Rs,
                                                            pi,
                                                            m,
                                                            lkX,
                                                            result.alphabetSize,alphabet);
                        }

                        if(std::isinf(val)){
                            exit(EXIT_FAILURE);
                        }

                        if(std::isnan(val)){
                            exit(EXIT_FAILURE);
                        }

                        idx=get_indices_X(coordTriangle_this_i,coordTriangle_this_j,up_corner_i,
                                          up_corner_j,bot_corner_i,bot_corner_j,m,h,w);

                        LogX[m_binary_this][idx]=val;
                    }
                    lw++;
                }

            }
            //***************************************************************************************
            //***************************************************************************************
            set_indeces_Y(up_corner_i,up_corner_j,bot_corner_i,bot_corner_j,m,h,w);
            tr_up_i=up_corner_i;
            tr_up_j=up_corner_j;
            if(checkboundary(up_corner_i,up_corner_j,bot_corner_i,bot_corner_j,h,w)){

                lw=0;
                for(int i=up_corner_i;i<=bot_corner_i;i++){
                    coordTriangle_this_i=i;
                    coordTriangle_prev_i=coordTriangle_this_i;
                    for(int j=0;j<=lw;j++){

                        coordTriangle_this_j=up_corner_j-j;
                        coordTriangle_prev_j=coordTriangle_this_j-1;
                        coordSeq_2=coordTriangle_this_j-1;
                        sRs=create_col_MSA(result_R.MSAs,coordSeq_2);

                        idx=get_indices_M(coordTriangle_prev_i,coordTriangle_prev_j,up_corner_i,
                                          up_corner_j,bot_corner_i,bot_corner_j,m-1,h,w);
                        if(idx>=0){
                            valM=LogM[m_binary_prev][idx];
                        }else{
                            valM=-INFINITY;
                        }

                        idx=get_indices_X(coordTriangle_prev_i,coordTriangle_prev_j,up_corner_i,
                                          up_corner_j,bot_corner_i,bot_corner_j,m-1,h,w);
                        if(idx>=0){
                            valX=LogX[m_binary_prev][idx];
                        }else{
                            valX=-INFINITY;
                        }

                        idx=get_indices_Y(coordTriangle_prev_i,coordTriangle_prev_j,up_corner_i,
                                          up_corner_j,bot_corner_i,bot_corner_j,m-1,h,w);
                        if(idx>=0){
                            valY=LogY[m_binary_prev][idx];
                        }else{
                            valY=-INFINITY;
                        }

                        if(std::isinf(valM) && std::isinf(valX) && std::isinf(valY)){
                            exit(EXIT_FAILURE);
                        }

                        if(local){
                            val=computeLK_Y_local_tree_s_opt(valM,
                                                             valX,
                                                             valY,
                                                             nu,
                                                             tree,
                                                             col_gap_Ls,sRs,
                                                             pi,
                                                             m,
                                                             lkY,
                                                             result.alphabetSize,alphabet);
                        }else{
                            val=computeLK_Y_all_edges_s_opt(valM,
                                                            valX,
                                                            valY,
                                                            nu,
                                                            tree,
                                                            col_gap_Ls,sRs,
                                                            pi,
                                                            m,
                                                            lkY,
                                                            result.alphabetSize,alphabet);
                        }

                        if(std::isinf(val)){
                            exit(EXIT_FAILURE);
                        }


                        if(std::isnan(val)){
                            exit(EXIT_FAILURE);
                        }

                        idx=get_indices_Y(coordTriangle_this_i,coordTriangle_this_j,up_corner_i,
                                          up_corner_j,bot_corner_i,bot_corner_j,m,h,w);

                        LogY[m_binary_this][idx]=val;
                    }
                    lw++;
                }

            }
            //***************************************************************************************
            //***************************************************************************************
            size_tr=int(ceil((tr_down_i-tr_up_i+1)*(tr_up_j-tr_down_j+1+1)/2));
            TR[m] = new int[size_tr](); /*TODO: optimize size TR*/
            memset(TR[m],0,size_tr*sizeof(TR[m][0]));
            set_indeces_T(up_corner_i,up_corner_j,bot_corner_i,bot_corner_j,m,h,w);

            if(checkboundary(up_corner_i,up_corner_j,bot_corner_i,bot_corner_j,h,w)){

                lw=0;
                for(int i=up_corner_i;i<=bot_corner_i;i++){
                    coordTriangle_this_i=i;
                    for(int j=0;j<=lw;j++){
                        coordTriangle_this_j=up_corner_j-j;

                        double mval;
                        double xval;
                        double yval;

                        idx=get_indices_M(coordTriangle_this_i,coordTriangle_this_j,up_corner_i,
                                          up_corner_j,bot_corner_i,bot_corner_j,m,h,w);
                        if(idx>=0){
                            mval=LogM[m_binary_this][idx];
                        }else{
                            mval=-INFINITY;
                        }

                        idx=get_indices_X(coordTriangle_this_i,coordTriangle_this_j,up_corner_i,
                                          up_corner_j,bot_corner_i,bot_corner_j,m,h,w);
                        if(idx>=0){
                            xval=LogX[m_binary_this][idx];
                        }else{
                            xval=-INFINITY;
                        }

                        idx=get_indices_Y(coordTriangle_this_i,coordTriangle_this_j,up_corner_i,
                                          up_corner_j,bot_corner_i,bot_corner_j,m,h,w);
                        if(idx>=0){
                            yval=LogY[m_binary_this][idx];
                        }else{
                            yval=-INFINITY;
                        }

                        mval=fabs(mval)<epsilon?-INFINITY:mval;
                        xval=fabs(xval)<epsilon?-INFINITY:xval;
                        yval=fabs(yval)<epsilon?-INFINITY:yval;

                        //.-----.------.------.-------.------.-------.-------.--------.-------.--------.//
                        int ttrr;

                        ttrr=index_of_max(mval,xval,yval,epsilon,generator,distribution);

                        idx=get_indices_T(coordTriangle_this_i,coordTriangle_this_j,up_corner_i,
                                          up_corner_j,bot_corner_i,bot_corner_j,m,h,w);

                        if(TR[m][idx]!=0){
                            exit(EXIT_FAILURE);
                        }

                        TR[m][idx]=ttrr;

                        //.-----.------.------.-------.------.-------.-------.--------.-------.--------.//

                        //.-----.------.------.-------.------.-------.-------.--------.-------.--------.//
                        if( (coordTriangle_this_i==(h-1)) & (coordTriangle_this_j==(w-1)) ){

                            max_of_3=max_of_three(mval,xval,yval,epsilon);

                            //TODO: check this part
                            //int num_subopt=1;
                            //fill_scores(max_of_3,max_lk,prev_max_lk,level_max_lk,last_d,m,flag_exit,
                            //            scores,counter,CENTER,num_subopt,start_depth);

                            prev_max_lk=max_lk;
                            max_lk=max_of_3;
                            level_max_lk=m;


                        }
                        //.-----.------.------.-------.------.-------.-------.--------.-------.--------.//



                    }
                    lw++;
                }
            }
        }

        //std::vector<ProgressivePIPResult> result_v;
        //for(int k=0;k<num_subopt;k++){
        //depth=start_depth+k;
        depth=start_depth;

        //if(depth>=d){
        //    break;
        //}


        //result.score=scores[k];
        result.score=scores;

        //.................................
//        out_score->precision(15);
//        *out_score<<result.score<<"\n";
        //.................................


        //std::string traceback_path(depth,ALPHABET::unknow);
        std::string traceback_path (depth, ' ');
        int id1=h-1;
        int id2=w-1;
        for(int lev=depth;lev>0;lev--){
            set_indeces_T(up_corner_i,up_corner_j,bot_corner_i,bot_corner_j,lev,h,w);
            idx=get_indices_T(id1,id2,up_corner_i,up_corner_j,bot_corner_i,bot_corner_j,lev,h,w);
            switch(TR[lev][idx]){
                case MATCH_STATE:
                    id1=id1-1;
                    id2=id2-1;
                    traceback_path[lev-1]=MATCH_CHAR;
                    break;
                case GAP_X_STATE:
                    id1=id1-1;
                    traceback_path[lev-1]=GAP_X_CHAR;
                    break;
                case GAP_Y_STATE:
                    id2=id2-1;
                    traceback_path[lev-1]=GAP_Y_CHAR;
                    break;
                default:
                    perror("ERROR in alignment_reconstruction !!!");
                    exit(EXIT_FAILURE);
            }
        }
        result.traceback_path=traceback_path;
        result.MSAs=build_MSA(traceback_path,result_L.MSAs,result_R.MSAs);
        //result_v.push_back(result);

        //}
        //delete scores;
        //.-----.------.------.-------.------.-------.-------.--------.-------.--------.//

        free(LogM[1]);
        free(LogM[0]);
        free(LogM);

        free(LogX[1]);
        free(LogX[0]);
        free(LogX);

        free(LogY[1]);
        free(LogY[0]);
        free(LogY);

        for(int i=last_d;i>=0;i--){
            free(TR[i]);
        }
        free(TR);
        //.-----.------.------.-------.------.-------.-------.--------.-------.--------.//

        //return result_v;
        return result;
    }
//DP-PIP
//    ProgressivePIPResult compute_DP3D_PIP_cross(ProgressivePIPResult &result_L,
//                                                ProgressivePIPResult &result_R,
//                                                VirtualNode *tree,
//                                                Likelihood &likelihood,
//                                                double gamma_rate,
//                                                int num_subopt,
//                                                bool local_tree,
//                                                bool stoch_backtracking_flag,
//                                                int ugglyflag,
//                                                std::ostream* &out_score,int alphabetSize){
//
//        ProgressivePIPResult result;
//
//        result=compute_DP3D_PIP(result_L,result_R,
//                                tree,likelihood,tau,
//                                nu,gamma_rate,num_subopt,local_tree,ugglyflag,out_score,alphabetSize);
//
//        return result;
//    }
//DP-PIP
/*    void add_sequence_to_vector_string(std::vector< std::pair<std::string,std::string> > &result,
                                       std::string name,
                                       std::vector< std::pair<std::string,std::string> > &sequences){

        typedef typename std::vector<std::pair<std::string,std::string>> ::iterator vect_iter;
        vect_iter iter = std::find_if(sequences.begin(),sequences.end(),CompareFirst(name));
        //std::string s=stringFromSequence(iter->second);
        std::string s=iter->second;
        result.push_back(std::make_pair(iter->first,s));

    }*/

    void add_sequence_to_alignment(ProgressivePIPResult &result,
                                   VirtualNode *tree,
                                   Alignment *alignment){

        int index;
        bool is_found=false;
        for(int i=0;i<alignment->align_dataset.size();i++){
            if(tree->vnode_name.compare(alignment->align_dataset.at(i)->seq_name)==0){
                index=i;
                is_found=true;
                break;
            }
        }

        if(!is_found){
            perror("ERROR: sequence not found\n");
            exit(EXIT_FAILURE);
        }

        result.MSAs.push_back(std::make_pair(alignment->align_dataset.at(index)->seq_name,alignment->align_dataset.at(index)->seq_data));

    }

//DP-PIP
    ProgressivePIPResult compute_DP3D_PIP_tree_cross(VirtualNode *tree,
                                                     Likelihood *likelihood,
                                                     Alignment *alignment,
                                                     const bpp::Alphabet *alphabet,
                                                     bpp::SiteContainer *sites,
                                                     double gamma_rate,
                                                     bool local_tree) {


        ProgressivePIPResult result;

        if (tree->isTerminalNode()) {

            add_sequence_to_alignment(result, tree, alignment);

        } else{

            ProgressivePIPResult result_L = compute_DP3D_PIP_tree_cross(tree->getNodeLeft(), likelihood, alignment, alphabet,sites,
                                                                        gamma_rate, local_tree);
            ProgressivePIPResult result_R = compute_DP3D_PIP_tree_cross(tree->getNodeRight(), likelihood, alignment, alphabet, sites,
                                                                        gamma_rate, local_tree);

            result = compute_DP3D_PIP(result_L, result_R, tree, likelihood, alignment, alphabet,sites, gamma_rate, local_tree);

        }

        return result;
    }

}
#endif //MINIJATI_PROGRESSIVEPIP_HPP
