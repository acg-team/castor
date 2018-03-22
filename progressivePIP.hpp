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
#include <cmath>
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


#include <limits>

#include "TSHTopologySearch.hpp"

double inf = std::numeric_limits<double>::infinity();


using namespace tshlib;

namespace progressivePIP{

#define ERR_STATE (-999)

#define MATCH_STATE 1
#define GAP_X_STATE 2
#define GAP_Y_STATE 3
#define STOP_STATE 4

#define MATCH_CHAR '1'
#define GAP_X_CHAR '2'
#define GAP_Y_CHAR '3'


    struct tree_msas{
        Utree *utree;
        //std::vector<bpp::SiteContainer *> msas;
        //std::vector<double> lks;
        double marginal_lk;
    };
    struct ProgressivePIPResult;

    struct ProgressivePIPResult{
        //std::vector<std::pair<std::string,std::string > > MSA_t;
        std::vector< std::pair<std::string,std::string> > MSAs;
        std::string traceback_path;
        double score;
        //std::map<std::string,Eigen::VectorXd> fv_map;
        //double lk_gap;
        //double pc0;
    };
    //DP-PIP
    Eigen::VectorXd fv_observed(std::string &s, unsigned long &idx,int alphabetSize,const bpp::Alphabet *alphabet){
        int ii;
        char ch=s[idx];
        Eigen::VectorXd fv = Eigen::VectorXd::Zero(alphabetSize+1);

        ii=alphabet->charToInt(&ch);
        ii=ii<0?alphabetSize:ii;

        fv[ii]=1.0;
        idx++;

        return fv;
    }
    //DP-PIP
    Eigen::VectorXd go_down(VirtualNode *tree,std::string &s, unsigned long &idx,int alphabetSize,const bpp::Alphabet *alphabet){
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
    void reset_corner(unsigned long &up_corner_i,unsigned long &up_corner_j,unsigned long &bot_corner_i,
                      unsigned long &bot_corner_j,unsigned long h,unsigned long w){
        unsigned long delta;

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
    std::string create_col_MSA_gap(unsigned long len){

        std::string colMSA (len,'-');

        return colMSA;
    }
    //DP-PIP
    int index_of_max(double m, double x, double y,double epsilon,
                     std::default_random_engine &generator,
                     std::uniform_real_distribution<double> &distribution){

        double random_number;

        if(not(std::isinf(m)) & not(std::isinf(x)) & (fabs((long double)(m-x))<epsilon)){
            x=m;
        }

        if(not(std::isinf(m)) & not(std::isinf(y)) & (fabs((long double)(m-y))<epsilon)){
            y=m;
        }

        if(not(std::isinf(x)) & not(std::isinf(y)) & (fabs((long double)(x-y))<epsilon)){
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
                    perror("ERROR in getIndexOfMax\n");
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

    }
    //DP-PIP
    double max_of_three(double a, double b, double c,double epsilon){

        //------------------------------------
        if(fabs((long double)a)<epsilon){
            //a=-INFINITY;
            a=-inf;
        }
        if(fabs((long double)b)<epsilon){
            //b=-INFINITY;
            b=-inf;
        }
        if(fabs((long double)c)<epsilon){
            //c=-INFINITY;
            c=-inf;
        }
        //------------------------------------

        if(std::isinf(a) && std::isinf(b) && std::isinf(c)){
            perror("getMaxOfThree: all inf\n");
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

        if(fabs((long double)mu)<1e-8){
            //perror("ERROR in compute_nu: mu too small");
            LOG(WARNING) << "[Parameter value] The parameter mu is too small! mu = " << mu;
        }

        return lambda*(tau+1/mu);
    }
    //DP-PIP
    bool is_inside(unsigned long x0,unsigned long y0,unsigned long xf,unsigned long yf,unsigned long xt,unsigned long yt){

        if((xt<x0) || (yt>y0) || (xt>xf) || (yt<yf)){
            return false;
        }

        if( (y0-yt)>(xt-x0) ){
            return false;
        }

        return true;
    }
    //DP-PIP
    void set_indeces_M(unsigned long &up_corner_i,unsigned long &up_corner_j,unsigned long &bot_corner_i,
                       unsigned long &bot_corner_j,unsigned long level,unsigned long h,unsigned long w){

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
    void set_indeces_X(unsigned long &up_corner_i,unsigned long &up_corner_j,unsigned long &bot_corner_i,
                       unsigned long &bot_corner_j,unsigned long level,unsigned long h,unsigned long w){

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
    void set_indeces_Y(unsigned long &up_corner_i,unsigned long &up_corner_j,unsigned long &bot_corner_i,
                       unsigned long &bot_corner_j,unsigned long level,unsigned long h,unsigned long w){

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
    void set_indeces_T(unsigned long &up_corner_i,unsigned long &up_corner_j,unsigned long &bot_corner_i,
                       unsigned long &bot_corner_j,unsigned long level,unsigned long h,unsigned long w){

        unsigned long up_corner_i_x;
        unsigned long up_corner_i_y;

        unsigned long up_corner_j_x;
        unsigned long up_corner_j_y;

        unsigned long bot_corner_i_x;
        unsigned long bot_corner_i_y;

        unsigned long bot_corner_j_x;
        unsigned long bot_corner_j_y;

        set_indeces_X(up_corner_i_x,up_corner_j_x,bot_corner_i_x,bot_corner_j_x,level,h,w);

        set_indeces_Y(up_corner_i_y,up_corner_j_y,bot_corner_i_y,bot_corner_j_y,level,h,w);

        unsigned long delta_i,delta_j;

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
    unsigned long get_indices_M(unsigned long nx,unsigned long ny,unsigned long up_corner_i,unsigned long up_corner_j,
                      unsigned long bot_corner_i,unsigned long bot_corner_j,unsigned long m,unsigned long h,unsigned long w){

        unsigned long idx;

        set_indeces_M(up_corner_i,up_corner_j,bot_corner_i,bot_corner_j,m,h,w);

        if(is_inside(up_corner_i,up_corner_j,bot_corner_i,bot_corner_j,nx,ny)){

            unsigned long dx,sx;

            dx=nx-up_corner_i+1;

            sx=((dx+1)*dx/2)-1;

            idx=sx+(ny-up_corner_j);
        }else{
            idx=ERR_STATE;
        }

        return idx;

    }
    //DP-PIP
    unsigned long get_indices_X(unsigned long nx,unsigned long ny,unsigned long up_corner_i,unsigned long up_corner_j,
                      unsigned long bot_corner_i,unsigned long bot_corner_j,unsigned long m,unsigned long h,unsigned long w){

        unsigned long idx;

        set_indeces_X(up_corner_i,up_corner_j,bot_corner_i,bot_corner_j,m,h,w);

        if(is_inside(up_corner_i,up_corner_j,bot_corner_i,bot_corner_j,nx,ny)){

            unsigned long dx,sx;

            dx=nx-up_corner_i+1;

            sx=((dx+1)*dx/2)-1;

            idx=sx+(ny-up_corner_j);
        }else{
            idx=ERR_STATE;
        }

        return idx;

    }
    //DP-PIP
    unsigned long get_indices_Y(unsigned long nx,unsigned long ny,unsigned long up_corner_i,unsigned long up_corner_j,
                      unsigned long bot_corner_i,unsigned long bot_corner_j,unsigned long m,unsigned long h,unsigned long w){

        unsigned long idx;

        set_indeces_Y(up_corner_i,up_corner_j,bot_corner_i,bot_corner_j,m,h,w);

        if(is_inside(up_corner_i,up_corner_j,bot_corner_i,bot_corner_j,nx,ny)){

            unsigned long dx,sx;

            dx=nx-up_corner_i+1;

            sx=((dx+1)*dx/2)-1;

            idx=sx+(ny-up_corner_j);
        }else{
            idx=ERR_STATE;
        }

        return idx;

    }
    //DP-PIP
    unsigned long get_indices_T(unsigned long nx,unsigned long ny,unsigned long up_corner_i,unsigned long up_corner_j,
                      unsigned long bot_corner_i,unsigned long bot_corner_j,unsigned long m,unsigned long h,unsigned long w){

        set_indeces_T(up_corner_i,up_corner_j,bot_corner_i,bot_corner_j,m,h,w);

        reset_corner(up_corner_i,up_corner_j,bot_corner_i,bot_corner_j,h,w);

        unsigned long idx;
        unsigned long dx,sx;

        dx=nx-up_corner_i+1;

        sx=((dx+1)*dx/2)-1;

        idx=sx+(ny-up_corner_j);

        return idx;

    }
    //DP-PIP
    void allgaps(VirtualNode *tree,std::string &s, unsigned long &idx,bool &flag){

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
        unsigned long idx;
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
            unsigned long ixx=idx;
            allgaps(tree->getNodeRight(),s,idx,flagR);
            unsigned long len;

            std::string sL;
            len=ixx;
            sL=s.substr(0,len);
            pL=compute_lk_gap_down(tree->getNodeLeft(),sL,pi,alphabetSize,alphabet);

            std::string sR;
            sR=s.substr(ixx);
            pR=compute_lk_gap_down(tree->getNodeRight(),sR,pi,alphabetSize,alphabet);
        }

        return pr+pL+pR;
    }
    //DP-PIP
    bool checkboundary(unsigned long up_corner_i,unsigned long up_corner_j,unsigned long bot_corner_i,
                       unsigned long bot_corner_j,unsigned long h,unsigned long w){

        if( (up_corner_i  >=0) & (up_corner_i  <h) &\
	   (up_corner_j  >=0) & (up_corner_j  <w) &\
	   (bot_corner_i >=0) & (bot_corner_i <h) &\
	   (bot_corner_j >=0) & (bot_corner_j <w)){
            return true;
        }

        return false;
    }
    //DP-PIP
    std::string create_col_MSA(std::vector<std::pair<std::string,std::string>> &result, unsigned long index){
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

            MSA_out.emplace_back(std::make_pair(seq_name,seq_aligned));
        }

        return MSA_out;
    }
    //DP-PIP
    std::vector<std::pair<std::string,std::string>> align_seq_right(std::vector<std::pair<std::string,std::string>> &result,
                                                                    std::string &traceback_path){
        unsigned int idx;

        std::vector<std::pair<std::string,std::string>> MSA_out;

        //typedef typename std::vector<std::pair<std::string,std::string>>::iterator vect_iterator;
        for (auto iter = result.begin(); iter != result.end(); iter++){

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

            MSA_out.emplace_back(std::make_pair(seq_name,seq_aligned));
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
        unsigned long idx;
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
            unsigned long ixx=idx;
            allgaps(tree->getNodeRight(),s,idx,flagR);

            unsigned long len;
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
        unsigned long idx;

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

            if(p_tree->getNodeUp()->getNodeUp()== nullptr){
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
        unsigned long idx;

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
                                       unsigned long m,
                                       std::map<std::string,double> &lkM,
                                       int alphabetSize,
                                       const bpp::Alphabet *alphabet){


        double pr;
        double val;
        //------------------------------------------------------------------------------------------
        //typedef typename std::map<std::string,double>::iterator MapIterator;
        std::string s;

        s.append(sL);
        s.append(sR);

        //MapIterator it=lkM.find(s);
        auto it=lkM.find(s);
        if(it == lkM.end()){

            //------------------------------------------------------------------------------------------
            //DOWN
            Eigen::VectorXd fvL;
            Eigen::VectorXd fvR;
            Eigen::VectorXd fv;
            double fv0;
            unsigned long idx;

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
            while(p_tree->getNodeUp()!= nullptr){

                fv=(p_tree->getPr()*fv);

                fv0=fv.dot(pi);

                pr+=(p_tree->getNodeUp()->getIota()*p_tree->getNodeUp()->getBeta()*fv0);

                p_tree=p_tree->getNodeUp();
            }
            //------------------------------------------------------------------------------------------

            pr=log((long double)pr);

            lkM[s]=pr;

        }else{
            pr=it->second;
        }

        val=-log((long double)m)+log((long double)nu)+pr+max_of_three(valM,valX,valY,DBL_EPSILON);

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
                                        unsigned long m,
                                        std::map<std::string,double> &lkM,
                                        int alphabetSize,
                                        const bpp::Alphabet *alphabet){


        double pr;
        double val;
        //------------------------------------------------------------------------------------------
        //typedef typename std::map<std::string,double>::iterator MapIterator;
        std::string s;
        Eigen::VectorXd fvL;
        Eigen::VectorXd fvR;
        Eigen::VectorXd fv;
        double fv0;
        unsigned long idx;

        s.append(sL);
        s.append(sR);

        auto it=lkM.find(s);
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

            pr=log((long double)pr);

            lkM[s]=pr;

        }else{
            pr=it->second;
        }
        //------------------------------------------------------------------------------------------

        val=-log((long double)m)+log((long double)nu)+pr+max_of_three(valM,valX,valY,DBL_EPSILON);

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
                                       unsigned long m,
                                       std::map<std::string,double> &lkX,
                                       int alphabetSize,
                                       const bpp::Alphabet *alphabet){


        double pr;
        double val;
        //------------------------------------------------------------------------------------------
        //typedef typename std::map<std::string,double>::iterator MapIterator;
        std::string s;

        s.append(sL);
        s.append(col_gap_R);

        auto it=lkX.find(s);
        if(it == lkX.end()){

            //------------------------------------------------------------------------------------------
            Eigen::VectorXd fvL;
            Eigen::VectorXd fvR;
            Eigen::VectorXd fv;
            unsigned long idx;
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
            while(p_tree->getNodeUp()!= nullptr){

                fv=(p_tree->getPr()*fv);

                fv0=fv.dot(pi);

                pr+=(p_tree->getNodeUp()->getIota()*p_tree->getNodeUp()->getBeta()*fv0);

                p_tree=p_tree->getNodeUp();
            }
            //******************************************************************************************

            pr=log((long double)pr);

            lkX[s]=pr;

        }else{
            pr=it->second;
        }
        //------------------------------------------------------------------------------------------

        val=-log((long double)m)+log((long double)nu)+pr+max_of_three(valM,valX,valY,DBL_EPSILON);

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
                                        unsigned long m,
                                        std::map<std::string,double> &lkX,
                                        int alphabetSize,
                                        const bpp::Alphabet *alphabet){


        double pr;
        double val;
        //------------------------------------------------------------------------------------------
        //typedef typename std::map<std::string,double>::iterator MapIterator;
        std::string s;
        Eigen::VectorXd fvL;
        Eigen::VectorXd fvR;
        Eigen::VectorXd fv;
        unsigned long idx;
        double fv0;

        s.append(sL);
        s.append(col_gap_R);

        auto it=lkX.find(s);
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

            pr=log((long double)pr);

            lkX[s]=pr;

        }else{
            pr=it->second;
        }
        //------------------------------------------------------------------------------------------

        val=-log((long double)m)+log((long double)nu)+pr+max_of_three(valM,valX,valY,DBL_EPSILON);

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
                                       unsigned long m,
                                       std::map<std::string,double> &lkY,
                                       int alphabetSize,
                                       const bpp::Alphabet *alphabet){


        double pr;
        double val;
        //------------------------------------------------------------------------------------------
        //typedef typename std::map<std::string,double>::iterator MapIterator;
        std::string s;

        s.append(col_gap_L);
        s.append(sR);

        auto it=lkY.find(s);
        if(it == lkY.end()){

            //------------------------------------------------------------------------------------------
            Eigen::VectorXd fvL;
            Eigen::VectorXd fvR;
            Eigen::VectorXd fv;
            double fv0;
            unsigned long idx=0;
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
            while(p_tree->getNodeUp()!= nullptr){

                fv=(p_tree->getPr()*fv);

                fv0=fv.dot(pi);

                pr+=(p_tree->getNodeUp()->getIota()*p_tree->getNodeUp()->getBeta()*fv0);

                p_tree=p_tree->getNodeUp();
            }
            //*******************************************************************************************

            pr=log((long double)pr);

            lkY[s]=pr;

        }else{
            pr=it->second;
        }
        //------------------------------------------------------------------------------------------

        val=-log((long double)m)+log((long double)nu)+pr+max_of_three(valM,valX,valY,DBL_EPSILON);

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
                                        unsigned long m,
                                        std::map<std::string,double> &lkY,
                                        int alphabetSize,
                                        const bpp::Alphabet *alphabet){


        double pr;
        double val;
        //------------------------------------------------------------------------------------------
        //typedef typename std::map<std::string,double>::iterator MapIterator;
        std::string s;
        Eigen::VectorXd fvL;
        Eigen::VectorXd fvR;
        Eigen::VectorXd fv;
        double fv0;
        unsigned long idx=0;

        s.append(col_gap_L);
        s.append(sR);

        auto it=lkY.find(s);
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

            pr=log((long double)pr);

            lkY[s]=pr;

        }else{
            pr=it->second;
        }
        //------------------------------------------------------------------------------------------


        val=-log((long double)m)+log((long double)nu)+pr+max_of_three(valM,valX,valY,DBL_EPSILON);


        return val;
    }
//DP-PIP
    bool check_uniform_len_s(std::vector<std::pair<std::string,std::string>> &result){
        unsigned long len;

        if(result.empty()){
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
    unsigned long get_length_seq_s(std::vector<std::pair<std::string,std::string>> &result){

        if(result.empty()){
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
                                                 VirtualNode *node,
                                                 bpp::Tree *tree,
                                                 UtreeBppUtils::treemap *tm,
                                                 Eigen::VectorXd pi,
                                                 double lambda,
                                                 double mu,
                                                 bpp::SequenceContainer *sequences,
                                                 const bpp::Alphabet *alphabet,
                                                 double gamma_rate,
                                                 bool local){

        ProgressivePIPResult result;


        int originalAlphabetSize=alphabet->getSize()-1;

        //@gamma_distribution
        double nu;
        double tau;

        double lambda_gamma = lambda * gamma_rate;
        double mu_gamma = mu * gamma_rate;

        if(local){
            //TODO: calcola tau dal nodo attuale
            bpp::Tree *subtree = tree->cloneSubtree(tm->right.at(node));
            tau = subtree->getTotalLength();
            //tau=tree->computeTotalTreeLength();
            nu=compute_nu(tau,lambda_gamma,mu_gamma);
            //TODO: tau and nu are already set?
            //tree->set_tau(tau);
            //tree_>set_nu(nu);
            //@gamma_distribution
            //tree.set_iota_local(tau,mu_gamma);
            //tree.set_beta_local(tau,mu_gamma);
        }else{

            //const Tree tree = tree->getTree();
            tau = tree->getTotalLength();
            nu = compute_nu(tau, lambda, mu);
        }


        unsigned long up_corner_i;
        unsigned long up_corner_j;
        unsigned long bot_corner_i;
        unsigned long bot_corner_j;
        unsigned long lw;
        unsigned long h,w;

        h=get_length_seq_s(result_L.MSAs)+1;
        w=get_length_seq_s(result_R.MSAs)+1;
        unsigned long d=(h-1)+(w-1)+1;

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
        signed long seed = std::chrono::system_clock::now().time_since_epoch().count();
        std::default_random_engine generator(seed);
        std::uniform_real_distribution<double> distribution(0.0,1.0);
        //.-----.------.------.//
//	std::cout<<"random generator OFF\n";
//	unsigned seed = 0;
//	std::default_random_engine generator(seed);
//	std::uniform_real_distribution<double> distribution(0.0,1.0);
        //.-----.------.------.-------.------.-------.-------.--------.-------.--------.//

        auto epsilon=DBL_EPSILON;

        //***************************************************************************************
        //***************************************************************************************
        if(local){
            pc0 = compute_pr_gap_local_tree_s(node,
                                              col_gap_Ls,
                                              col_gap_Rs,
                                              pi,
                                              originalAlphabetSize,
                                              alphabet);
        }else{
            pc0 = compute_pr_gap_all_edges_s(node,
                                             col_gap_Ls,
                                             col_gap_Rs,
                                             pi,
                                             originalAlphabetSize,
                                             alphabet);
        }

        //***************************************************************************************
        //***************************************************************************************

        auto** LogM = new double*[2];
        auto** LogX = new double*[2];
        auto** LogY = new double*[2];

        auto** TR = new int*[d];

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
        //double max_lk=-INFINITY;
        double max_lk=-inf;
        //double prev_max_lk=-INFINITY;
        double prev_max_lk=-inf;
        signed long level_max_lk=INT_MIN;
        double val;
        unsigned long m_binary_this;
        unsigned long m_binary_prev;

        double valM;
        double valX;
        double valY;

        int idx;

        unsigned long coordSeq_1;
        unsigned long coordSeq_2;
        unsigned long coordTriangle_this_i;
        unsigned long coordTriangle_this_j;
        unsigned long coordTriangle_prev_i;
        unsigned long coordTriangle_prev_j;

        //int counter;

        double scores=-inf;
        int start_depth;
        unsigned long depth;

        bool flag_exit=false;
        unsigned long last_d=d-1;
        unsigned long size_tr,tr_up_i,tr_up_j,tr_down_i,tr_down_j;
        std::map<std::string,double> lkM;
        std::map<std::string,double> lkX;
        std::map<std::string,double> lkY;

        //counter=0;
        for(unsigned long m=1;m<d;m++){

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
                for(unsigned long i=up_corner_i;i<=bot_corner_i;i++){

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
                            //valM=-INFINITY;
                            valM=-inf;
                        }

                        idx=get_indices_X(coordTriangle_prev_i,coordTriangle_prev_j,up_corner_i,up_corner_j,bot_corner_i,bot_corner_j,m-1,h,w);
                        if(idx>=0){
                            valX=LogX[m_binary_prev][idx];
                        }else{
                            //valX=-INFINITY;
                            valX=-inf;
                        }

                        idx=get_indices_Y(coordTriangle_prev_i,coordTriangle_prev_j,up_corner_i,up_corner_j,bot_corner_i,bot_corner_j,m-1,h,w);
                        if(idx>=0){
                            valY=LogY[m_binary_prev][idx];
                        }else{
                            //valY=-INFINITY;
                            valY=-inf;
                        }

                        if(std::isinf(valM) && std::isinf(valX) && std::isinf(valY)){
                            exit(EXIT_FAILURE);
                        }

                        if(local){
                            val=computeLK_M_local_tree_s_opt(valM,
                                                             valX,
                                                             valY,
                                                             nu,
                                                             node,
                                                             sLs, sRs,
                                                             pi,
                                                             m,
                                                             lkM,
                                                             originalAlphabetSize, alphabet);
                        }else{
                            val=computeLK_M_all_edges_s_opt(valM,
                                                            valX,
                                                            valY,
                                                            nu,
                                                            node,
                                                            sLs, sRs,
                                                            pi,
                                                            m,
                                                            lkM,
                                                            originalAlphabetSize, alphabet);
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
                for(unsigned long i=up_corner_i;i<=bot_corner_i;i++){

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
                            //valM=-INFINITY;
                            valM=-inf;
                        }

                        idx=get_indices_X(coordTriangle_prev_i,coordTriangle_prev_j,up_corner_i,
                                          up_corner_j,bot_corner_i,bot_corner_j,m-1,h,w);
                        if(idx>=0){
                            valX=LogX[m_binary_prev][idx];
                        }else{
                            //valX=-INFINITY;
                            valX=-inf;
                        }

                        idx=get_indices_Y(coordTriangle_prev_i,coordTriangle_prev_j,up_corner_i,
                                          up_corner_j,bot_corner_i,bot_corner_j,m-1,h,w);
                        if(idx>=0){
                            valY=LogY[m_binary_prev][idx];
                        }else{
                            //valY=-INFINITY;
                            valY=-inf;
                        }

                        if(std::isinf(valM) && std::isinf(valX) && std::isinf(valY)){
                            exit(EXIT_FAILURE);
                        }

                        if(local){
                            val=computeLK_X_local_tree_s_opt(valM,
                                                             valX,
                                                             valY,
                                                             nu,
                                                             node,
                                                             sLs, col_gap_Rs,
                                                             pi,
                                                             m,
                                                             lkX,
                                                             originalAlphabetSize, alphabet);
                        }else{
                            val=computeLK_X_all_edges_s_opt(valM,
                                                            valX,
                                                            valY,
                                                            nu,
                                                            node,
                                                            sLs, col_gap_Rs,
                                                            pi,
                                                            m,
                                                            lkX,
                                                            originalAlphabetSize, alphabet);
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
                for(unsigned long i=up_corner_i;i<=bot_corner_i;i++){
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
                            //valM=-INFINITY;
                            valM=-inf;
                        }

                        idx=get_indices_X(coordTriangle_prev_i,coordTriangle_prev_j,up_corner_i,
                                          up_corner_j,bot_corner_i,bot_corner_j,m-1,h,w);
                        if(idx>=0){
                            valX=LogX[m_binary_prev][idx];
                        }else{
                            //valX=-INFINITY;
                            valX=-inf;
                        }

                        idx=get_indices_Y(coordTriangle_prev_i,coordTriangle_prev_j,up_corner_i,
                                          up_corner_j,bot_corner_i,bot_corner_j,m-1,h,w);
                        if(idx>=0){
                            valY=LogY[m_binary_prev][idx];
                        }else{
                            //valY=-INFINITY;
                            valY=-inf;
                        }

                        if(std::isinf(valM) && std::isinf(valX) && std::isinf(valY)){
                            exit(EXIT_FAILURE);
                        }

                        if(local){
                            val=computeLK_Y_local_tree_s_opt(valM,
                                                             valX,
                                                             valY,
                                                             nu,
                                                             node,
                                                             col_gap_Ls, sRs,
                                                             pi,
                                                             m,
                                                             lkY,
                                                             originalAlphabetSize, alphabet);
                        }else{
                            val=computeLK_Y_all_edges_s_opt(valM,
                                                            valX,
                                                            valY,
                                                            nu,
                                                            node,
                                                            col_gap_Ls, sRs,
                                                            pi,
                                                            m,
                                                            lkY,
                                                            originalAlphabetSize, alphabet);
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
            size_tr=(unsigned long)ceil((tr_down_i-tr_up_i+1)*(tr_up_j-tr_down_j+1+1)/2);
            TR[m] = new int[size_tr](); /*TODO: optimize size TR*/
            memset(TR[m],0,size_tr*sizeof(TR[m][0]));
            set_indeces_T(up_corner_i,up_corner_j,bot_corner_i,bot_corner_j,m,h,w);

            if(checkboundary(up_corner_i,up_corner_j,bot_corner_i,bot_corner_j,h,w)){

                lw=0;
                for(unsigned long i=up_corner_i;i<=bot_corner_i;i++){
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
                            //mval=-INFINITY;
                            mval=-inf;
                        }

                        idx=get_indices_X(coordTriangle_this_i,coordTriangle_this_j,up_corner_i,
                                          up_corner_j,bot_corner_i,bot_corner_j,m,h,w);
                        if(idx>=0){
                            xval=LogX[m_binary_this][idx];
                        }else{
                            //xval=-INFINITY;
                            xval=-inf;
                        }

                        idx=get_indices_Y(coordTriangle_this_i,coordTriangle_this_j,up_corner_i,
                                          up_corner_j,bot_corner_i,bot_corner_j,m,h,w);
                        if(idx>=0){
                            yval=LogY[m_binary_this][idx];
                        }else{
                            //yval=-INFINITY;
                            yval=-inf;
                        }

                        //mval=fabs((long double)mval)<epsilon?-INFINITY:mval;
                        //xval=fabs((long double)xval)<epsilon?-INFINITY:xval;
                        //yval=fabs((long double)yval)<epsilon?-INFINITY:yval;
                        mval=fabs((long double)mval)<epsilon?-inf:mval;
                        xval=fabs((long double)xval)<epsilon?-inf:xval;
                        yval=fabs((long double)yval)<epsilon?-inf:yval;

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

                            //prev_max_lk=max_lk;
                            //max_lk=max_of_3;
                            if(max_of_3>scores){
                                scores=max_of_3;
                                level_max_lk=m;
                            }



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
        //depth=start_depth;
        depth=level_max_lk;

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
            int state = TR[lev][idx];
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

        //VLOG(1)<<traceback_path;

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

    void add_sequence_to_alignment(ProgressivePIPResult &result,
                                   VirtualNode *tree,
                                   bpp::SequenceContainer *sequences){

//        unsigned long index=0;
//        bool is_found=false;
//        for(unsigned long i=0;i<alignment->align_dataset.size();i++){
//            if(tree->vnode_name.compare(alignment->align_dataset.at(i)->seq_name)==0){
//                index=i;
//                is_found=true;
//                break;
//            }
//        }
//
//        if(!is_found){
//            perror("ERROR: sequence not found\n");
//            exit(EXIT_FAILURE);
//        }

        std::string seqname = sequences->getSequencesNames().at(tree->vnode_seqid);
        std::string seqdata = sequences->getSequence(seqname).toString();


        result.MSAs.emplace_back(std::make_pair(seqname,seqdata));

    }

    //DP-PIP
    ProgressivePIPResult compute_DP3D_PIP_tree_cross(VirtualNode *node,
                                                     bpp::Tree *tree,
                                                     UtreeBppUtils::treemap *tm,
                                                     Eigen::VectorXd pi,
                                                     double lambda,
                                                     double mu,
                                                     bpp::SequenceContainer *sequences,
                                                     const bpp::Alphabet *alphabet,
                                                     double gamma_rate,
                                                     bool local_tree) {


        ProgressivePIPResult result;

        if (node->isTerminalNode()) {

            add_sequence_to_alignment(result, node, sequences);
            //VLOG(2) << "[PPIP] Processing node " << node->vnode_name;

        } else{

            ProgressivePIPResult result_L = compute_DP3D_PIP_tree_cross(node->getNodeLeft(), tree, tm, pi, lambda, mu, sequences, alphabet, gamma_rate, local_tree);
            ProgressivePIPResult result_R = compute_DP3D_PIP_tree_cross(node->getNodeRight(), tree, tm, pi, lambda, mu, sequences, alphabet, gamma_rate, local_tree);

            result = compute_DP3D_PIP(result_L, result_R, node, tree, tm, pi, lambda, mu, sequences, alphabet, gamma_rate, local_tree);
            //VLOG(2) << "[PPIP] Processing node " << node->vnode_name;

        }

        return result;
    }


    tshlib::Utree *marginalizationOverMSAs(tshlib::TreeSearch *treesearch,
                                           bpp::Alphabet *alpha,
                                           Eigen::VectorXd pi,
                                           double lambda,
                                           double mu,
                                           bpp::SequenceContainer *sequences,
                                           UtreeBppUtils::treemap &tm){

        std::vector<tree_msas> trees_msas;
        std::vector<double> lk_tree;
        int num_trees;
        int num_msas;

        double max_lk;
        int index_max_lk;

        num_msas=10;
        num_trees=10;
        max_lk=-inf;
        index_max_lk=-1;
        for(int i=0;i<num_trees;i++){
            auto utree = new::tshlib::Utree();

            tree_msas tree_msa;

            tree_msa.utree=utree;
            VirtualNode *root;
            bpp::Tree *tree = nullptr;

            double marginal_lk=0.0;
            for(int j=0;j<num_msas;j++){

                progressivePIP::ProgressivePIPResult MSA;
                MSA = progressivePIP::compute_DP3D_PIP_tree_cross(root, tree, &tm, pi, lambda, mu, sequences, alpha, 1.0, false);

                /*
                auto sequences = new bpp::VectorSequenceContainer(alpha);
                for (int i = 0; i < MSA_t.MSAs.size(); i++) {
                    sequences->addSequence(*(new bpp::BasicSequence(MSA_t.MSAs.at(i).first, MSA_t.MSAs.at(i).second, alpha)), true);
                }
                auto sites = new bpp::VectorSiteContainer(*sequences);
                tree_msa.msas.push_back(sites);
                tree_msa.lks.push_back(MSA_t.score);
                */

                marginal_lk+=MSA.score;

            }
            tree_msa.marginal_lk=marginal_lk;

            trees_msas.push_back(tree_msa);

            if(tree_msa.marginal_lk>max_lk){
                max_lk=tree_msa.marginal_lk;
                index_max_lk=i;
            }

        }

        if(index_max_lk>-1){
            return trees_msas.at(index_max_lk).utree;
        }else{
            return NULL;
        }

    }


}
#endif //MINIJATI_PROGRESSIVEPIP_HPP
