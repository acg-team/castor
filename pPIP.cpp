/*******************************************************************************
 * Licensed Materials - Property of Lorenzo Gatti & Massimo Maiolo
 *
 *
 * Copyright (C) 2015-2018 by Lorenzo Gatti & Massimo Maiolo
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

/**
 * @file pPIP.cpp
 * @author Lorenzo Gatti
 * @author Massimo Maiolo
 * @date 19 02 2018
 * @version 1.0
 * @maintainer Lorenzo Gatti
 * @email lg@lorenzogatti.me
 * @maintainer Massimo Maiolo
 * @email massimo.maiolo@zhaw.ch
 * @status Development
 *
 * @brief
 * @details
 * @pre
 * @bug
 * @warning
 *
 * @see For more information visit: 
 */

#include <chrono>
#include <random>

#include "pPIP.hpp"
#include <Bpp/Numeric/Matrix/MatrixTools.h>

#define ERR_STATE (-999)

#define MATCH_STATE 1
#define GAP_X_STATE 2
#define GAP_Y_STATE 3
#define STOP_STATE 4

#define MATCH_CHAR '1'
#define GAP_X_CHAR '2'
#define GAP_Y_CHAR '3'

using namespace bpp;

pPIP::pPIP(bpp::Alphabet *alphabet){

    _alphabet=alphabet;
    _alphabetSize=_alphabet->getSize();
    _extendedAlphabetSize=_alphabet->getSize()+1;

};
void pPIP::init(const Tree *tree,
                UtreeBppUtils::treemap *tm,
                std::vector<tshlib::VirtualNode *> &listNodes,
                const Vdouble &pi,
                double lambda,
                double mu){


    _reserve(listNodes.size());

    _setTree(tree);
    _setLambda(lambda);
    _setMu(mu);
    _setPi(pi);
    _setTau();
    _setNu();
    _setAllIotas(tm,listNodes);
    _setAllBetas(tm,listNodes);
    _computePr(tm,listNodes);

}
void pPIP::_reserve(unsigned long numNodes){

    //_fv.resize(numNodes);
    //_lkxy.resize(numNodes);
    _score.resize(numNodes);
    _iota.resize(numNodes);
    _beta.resize(numNodes);
    _pr.resize(numNodes);
    _seqNames.resize(numNodes);
    _MSA.resize(numNodes);

}
void pPIP::_setTree(const Tree *tree) {
    _tree = new TreeTemplate<Node>(*tree);
}
bool pPIP::is_inside(unsigned long x0,unsigned long y0,unsigned long xf,unsigned long yf,unsigned long xt,unsigned long yt){

    if((xt<x0) || (yt>y0) || (xt>xf) || (yt<yf)){
        return false;
    }

    if( (y0-yt)>(xt-x0) ){
        return false;
    }

    return true;
}
void pPIP::set_indeces_M(unsigned long &up_corner_i,unsigned long &up_corner_j,unsigned long &bot_corner_i,
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
void pPIP::set_indeces_X(unsigned long &up_corner_i,unsigned long &up_corner_j,unsigned long &bot_corner_i,
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
void pPIP::set_indeces_Y(unsigned long &up_corner_i,unsigned long &up_corner_j,unsigned long &bot_corner_i,
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
signed long pPIP::get_indices_M(unsigned long nx,unsigned long ny,unsigned long up_corner_i,unsigned long up_corner_j,
                            unsigned long bot_corner_i,unsigned long bot_corner_j,unsigned long m,unsigned long h,unsigned long w){

    signed long idx;

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
signed long pPIP::get_indices_X(unsigned long nx,unsigned long ny,unsigned long up_corner_i,unsigned long up_corner_j,
                            unsigned long bot_corner_i,unsigned long bot_corner_j,unsigned long m,unsigned long h,unsigned long w){

    signed long idx;

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
signed long pPIP::get_indices_Y(unsigned long nx,unsigned long ny,unsigned long up_corner_i,unsigned long up_corner_j,
                            unsigned long bot_corner_i,unsigned long bot_corner_j,unsigned long m,unsigned long h,unsigned long w){

    signed long idx;

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
void pPIP::set_indeces_T(unsigned long &up_corner_i,unsigned long &up_corner_j,unsigned long &bot_corner_i,
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
void pPIP::reset_corner(unsigned long &up_corner_i,unsigned long &up_corner_j,unsigned long &bot_corner_i,
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
unsigned long pPIP::get_indices_T(unsigned long nx,unsigned long ny,unsigned long up_corner_i,unsigned long up_corner_j,
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
int pPIP::index_of_max(double m, double x, double y,double epsilon,
                 std::default_random_engine &generator,
                 std::uniform_real_distribution<double> &distribution){

    double random_number;

    if(not(std::isinf(m)) & not(std::isinf(x)) & (fabs((m-x))<epsilon)){
        x=m;
    }

    if(not(std::isinf(m)) & not(std::isinf(y)) & (fabs((m-y))<epsilon)){
        y=m;
    }

    if(not(std::isinf(x)) & not(std::isinf(y)) & (fabs((x-y))<epsilon)){
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

}
double pPIP::max_of_three(double a, double b, double c,double epsilon){

    //------------------------------------
    if(fabs(a)<epsilon){
        a=-std::numeric_limits<double>::infinity();
    }
    if(fabs(b)<epsilon){
        b=-std::numeric_limits<double>::infinity();
    }
    if(fabs(c)<epsilon){
        c=-std::numeric_limits<double>::infinity();
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
bool pPIP::checkboundary(unsigned long up_corner_i,unsigned long up_corner_j,unsigned long bot_corner_i,
                   unsigned long bot_corner_j,unsigned long h,unsigned long w){

    if( (up_corner_i  >=0) & (up_corner_i  <h) &\
	   (up_corner_j  >=0) & (up_corner_j  <w) &\
	   (bot_corner_i >=0) & (bot_corner_i <h) &\
	   (bot_corner_j >=0) & (bot_corner_j <w)){
        return true;
    }

    return false;
}
/*
bpp::ColMatrix<double> pPIP::computeLKgapColLocal(bpp::Node *node,
                                   double &val,
                                   double &p0){

    double fv0;
    double pr;
    //Vdouble fv;

    auto sons = node->getSons();
    auto s1 = sons.at(0);
    auto s2 = sons.at(1);
    int s1ID = s1->getId();
    int s2ID = s2->getId();

    //auto fvL=&(_fv.at(s1)[0]);
    //auto fvR=&(_fv.at(s2)[0]);
    //auto fvL=&(_fv.at(s1ID)[0]);
    //auto fvR=&(_fv.at(s2ID)[0]);
    auto fvL=(_fv.at(s1ID).at(0));
    auto fvR=(_fv.at(s2ID).at(0));

    //bpp::RowMatrix<double> *pr1 = &(_pr.at(s1));
    //bpp::RowMatrix<double> *pr2 = &(_pr.at(s2));
    bpp::RowMatrix<double> pr1 = (_pr.at(s1ID));
    bpp::RowMatrix<double> pr2 = (_pr.at(s2ID));

    //auto PrfvL=MatrixBppUtils::matrixVectorProd(*pr1,*fvL);
    //auto PrfvR=MatrixBppUtils::matrixVectorProd(*pr2,*fvR);

    bpp::ColMatrix<double> PrfvL;
    bpp::ColMatrix<double> PrfvR;

    bpp::MatrixTools::mult(pr1,fvL,PrfvL);
    bpp::MatrixTools::mult(pr2,fvR,PrfvR);

//    std::cout<<t2.getNumberOfRows()<<"x"<<t2.getNumberOfColumns()<<std::endl;
//    for(int i=0;i<5;i++){
//        std::cout<<"[t2]"<<t2(i,0)<<std::endl;
//    }


    //fv=MatrixBppUtils::cwiseProd(&PrfvL,&PrfvR);
    bpp::ColMatrix<double> fv;
    bpp::MatrixTools::hadamardMult(PrfvL,PrfvR,fv);

    fv0=MatrixBppUtils::dotProd(fv,_pi);

    //pr=_iota[node]*fv0;
    pr=_iota[node->getId()]*fv0;

    //double pL=(_lkxy[s1]).at(0);
    //double pR=(_lkxy[s2]).at(0);
    double pL=(_lkxy.at(s1ID)).at(0);
    double pR=(_lkxy.at(s2ID)).at(0);

    val=pr+pL+pR;

    p0=val;

    return fv;
}
 */
/*
bpp::ColMatrix<double> pPIP::computeLKmatchLocal(double valM,
                            double valX,
                            double valY,
                            bpp::Node *node,unsigned long col_i, unsigned long col_j,
                            unsigned long m,
                            double &val){


    double pr;
    //Vdouble fv;
    double fv0;

    auto sons = node->getSons();
    auto s1 = sons.at(0);
    auto s2 = sons.at(1);
    int s1ID = s1->getId();
    int s2ID = s2->getId();

    //auto fvL=&(_fv.at(s1)[col_i]);
    //auto fvR=&(_fv.at(s2)[col_j]);
    auto fvL=(_fv.at(s1ID).at(col_i));
    auto fvR=(_fv.at(s2ID).at(col_j));

    //bpp::RowMatrix<double> *pr1 = &(_pr.at(s1));
    //bpp::RowMatrix<double> *pr2 = &(_pr.at(s2));
    bpp::RowMatrix<double> pr1 = (_pr.at(s1ID));
    bpp::RowMatrix<double> pr2 = (_pr.at(s2ID));

    //auto PrfvL=MatrixBppUtils::matrixVectorProd(*pr1,*fvL);
    //auto PrfvR=MatrixBppUtils::matrixVectorProd(*pr2,*fvR);

    bpp::ColMatrix<double> PrfvL;
    bpp::ColMatrix<double> PrfvR;

    bpp::MatrixTools::mult(pr1,fvL,PrfvL);
    bpp::MatrixTools::mult(pr2,fvR,PrfvR);

    bpp::ColMatrix<double> fv;
    bpp::MatrixTools::hadamardMult(PrfvL,PrfvR,fv);

    fv0=MatrixBppUtils::dotProd(fv,_pi);

    //pr=_iota[node]*_beta[node]*fv0;
    pr=_iota.at(node->getId())*_beta.at(node->getId())*fv0;

    val=-log(m)+log(_nu)+log(pr)+max_of_three(valM,valX,valY,DBL_EPSILON);

    return fv;

}
 */
/*
bpp::ColMatrix<double> pPIP::computeLKgapxLocal(double valM,
                                    double valX,
                                    double valY,
                                    bpp::Node *node,
                                    unsigned long col_i,
                                    unsigned long col_j,
                                    unsigned long m,
                                    double &val,
                                    double &lkx){


    double pr;
    //Vdouble fv;
    double fv0;

    auto sons = node->getSons();
    auto s1 = sons.at(0);
    auto s2 = sons.at(1);
    int s1ID = s1->getId();
    int s2ID = s2->getId();

    //auto fvL=&(_fv.at(s1)[col_i]);
    //auto fvR=&(_fv.at(s2)[col_j]);
    auto fvL=(_fv.at(s1ID).at(col_i));
    auto fvR=(_fv.at(s2ID).at(col_j));

    //bpp::RowMatrix<double> *pr1 = &(_pr.at(s1));
    //bpp::RowMatrix<double> *pr2 = &(_pr.at(s2));
    bpp::RowMatrix<double> pr1 = (_pr.at(s1ID));
    bpp::RowMatrix<double> pr2 = (_pr.at(s2ID));

    bpp::ColMatrix<double> PrfvL;
    bpp::ColMatrix<double> PrfvR;

    bpp::MatrixTools::mult(pr1,fvL,PrfvL);
    bpp::MatrixTools::mult(pr2,fvR,PrfvR);

    //auto PrfvL=MatrixBppUtils::matrixVectorProd(*pr1,*fvL);
    //auto PrfvR=MatrixBppUtils::matrixVectorProd(*pr2,*fvR);

    //fv=MatrixBppUtils::cwiseProd(&PrfvL,&PrfvR);

    bpp::ColMatrix<double> fv;
    bpp::MatrixTools::hadamardMult(PrfvL,PrfvR,fv);

    fv0=MatrixBppUtils::dotProd(fv,_pi);

    //pr=_iota[node]*_beta[node]*fv0;
    pr=_iota.at(node->getId())*_beta.at(node->getId())*fv0;

    //double pL=(_lkxy[s1]).at(col_i);
    double pL=(_lkxy.at(s1ID)).at(col_i);

    pr+=pL;

    lkx=log(pr);

    val=-log(m)+log(_nu)+log(pr)+max_of_three(valM,valX,valY,DBL_EPSILON);

    return fv;
}
 */
/*
bpp::ColMatrix<double> pPIP::computeLKgapyLocal(double valM,
                                    double valX,
                                    double valY,
                                    bpp::Node *node,
                                    unsigned long col_i,
                                    unsigned long col_j,
                                    unsigned long m,
                                    double &val,
                                    double &lky){


    double pr;
    //Vdouble fv;
    double fv0;

    auto sons = node->getSons();
    auto s1 = sons.at(0);
    auto s2 = sons.at(1);
    int s1ID = s1->getId();
    int s2ID = s2->getId();

    //auto fvL=&(_fv.at(s1)[col_i]);
    //auto fvR=&(_fv.at(s2)[col_j]);
    auto fvL=(_fv.at(s1ID).at(col_i));
    auto fvR=(_fv.at(s2ID).at(col_j));

    //bpp::RowMatrix<double> *pr1 = &(_pr.at(s1));
    //bpp::RowMatrix<double> *pr2 = &(_pr.at(s2));
    bpp::RowMatrix<double> pr1 = (_pr.at(s1ID));
    bpp::RowMatrix<double> pr2 = (_pr.at(s2ID));

    bpp::ColMatrix<double> PrfvL;
    bpp::ColMatrix<double> PrfvR;

    bpp::MatrixTools::mult(pr1,fvL,PrfvL);
    bpp::MatrixTools::mult(pr2,fvR,PrfvR);
    //auto PrfvL=MatrixBppUtils::matrixVectorProd(*pr1,*fvL);
    //auto PrfvR=MatrixBppUtils::matrixVectorProd(*pr2,*fvR);

    //fv=MatrixBppUtils::cwiseProd(&PrfvL,&PrfvR);

    bpp::ColMatrix<double> fv;
    bpp::MatrixTools::hadamardMult(PrfvL,PrfvR,fv);

    fv0=MatrixBppUtils::dotProd(fv,_pi);

    //pr=_iota[node]*_beta[node]*fv0;
    pr=_iota.at(node->getId())*_beta.at(node->getId())*fv0;

    //double pR=(_lkxy[s2]).at(col_j);
    double pR=(_lkxy.at(s2ID)).at(col_j);

    pr+=pR;

    lky=log(pr);

    val=-log(m)+log(_nu)+log(pr)+max_of_three(valM,valX,valY,DBL_EPSILON);

    return fv;
}
 */
//bool pPIP::checkUniformLen(std::vector<std::pair<std::string,std::string>> &result){
//    unsigned long len;
//
//    if(result.empty()){
//        return false;
//    }
//
//    len=result.at(0).second.length();
//
//    for(unsigned int i=1;i<result.size();i++){
//        if(len!=result.at(i).second.length()){
//            return false;
//        }
//    }
//
//    return true;
//}
//unsigned long pPIP::getMSAlength(std::vector< std::string > &result){
//
//    if(result.empty()){
//        return 0;
//    }
//
////    if(!checkUniformLen(result)){
////        perror("ERROR in get_length_seq: non aligned");
////    }
//
//    return result.size();
//}
std::string pPIP::createGapCol(unsigned long len){

    std::string colMSA (len,'-');

    return colMSA;
}
//std::string pPIP::createMSAcol(std::vector< std::string > &result, unsigned long index){
//    std::string colMSA;
//
//    //for(unsigned int i=0;i<result.size();i++){
//    for(const auto &r : result){
//        //colMSA.append(result.at(i).second,index,1);
//        colMSA.append(r.second,index,1);
//    }
//
//    return colMSA;
//}
//std::vector<std::pair<std::string,std::string>> pPIP::align_seq_left(	std::vector<std::pair<std::string,std::string>> &MSA_in,
//                                                                   std::string &traceback_path){
//
//    unsigned int idx;
//
//    std::vector<std::pair<std::string,std::string>> MSA_out;
//
//    //typedef typename std::vector<std::pair<std::string,std::string>>::iterator vect_iterator;
//    //for (auto iter = MSA_in.begin(); iter != MSA_in.end(); iter++){
//    for (auto&& iter : MSA_in) {
//        //std::pair<std::string,std::string> seq = (*ter);
//        std::pair<std::string,std::string> seq = (iter);
//
//        std::string seq_name=seq.first;
//        std::string seq_not_aligned=seq.second;
//        std::string seq_aligned(traceback_path.size(),'-');
//
//        idx=0;
//        for(unsigned int j=0;j<traceback_path.size();j++){
//
//            if(traceback_path.at(j)==MATCH_CHAR){
//
//                seq_aligned[j]=seq_not_aligned.at(idx);
//                idx++;
//
//            }else if(traceback_path.at(j)==GAP_X_CHAR){
//
//                seq_aligned[j]=seq_not_aligned.at(idx);
//                idx++;
//
//            }else if(traceback_path.at(j)==GAP_Y_CHAR){
//
//                seq_aligned[j]='-';
//
//            }else{
//                perror("ERROR in align_seq_left");
//            }
//        }
//
//        MSA_out.emplace_back(std::make_pair(seq_name,seq_aligned));
//    }
//
//    return MSA_out;
//}
//std::vector< std::string > pPIP::align_seq_left(std::vector< std::string > &MSA_in, std::string &traceback_path){
//
//    unsigned int idx;
//
//    std::vector< std::string > MSA_out;
//
//    for (auto&& seq : MSA_in) {
//        //std::string seq = (iter);
//
//        //std::string seq_name=seq.first;
//        //std::string seq_not_aligned=seq.second;
//        std::string seq_aligned(traceback_path.size(),'-');
//
//        idx=0;
//        for(unsigned int j=0;j<traceback_path.size();j++){
//
//            if(traceback_path.at(j)==MATCH_CHAR){
//
//                seq_aligned[j]=seq_not_aligned.at(idx);
//                idx++;
//
//            }else if(traceback_path.at(j)==GAP_X_CHAR){
//
//                seq_aligned[j]=seq_not_aligned.at(idx);
//                idx++;
//
//            }else if(traceback_path.at(j)==GAP_Y_CHAR){
//
//                seq_aligned[j]='-';
//
//            }else{
//                perror("ERROR in align_seq_left");
//            }
//        }
//
//        MSA_out.emplace_back(std::make_pair(seq_name,seq_aligned));
//    }
//
//    return MSA_out;
//}
//std::vector<std::pair<std::string,std::string>> pPIP::align_seq_right(std::vector<std::pair<std::string,std::string>> &result,
//                                                                std::string &traceback_path){
//    unsigned int idx;
//
//    std::vector<std::pair<std::string,std::string>> MSA_out;
//
//    //for (auto iter = result.begin(); iter != result.end(); iter++){
//    for (auto&& iter : result) {
//        //std::pair<std::string,std::string> seq = (*iter);
//        std::pair<std::string,std::string> seq = (iter);
//
//        std::string seq_name=seq.first;
//        std::string seq_not_aligned=seq.second;
//        std::string seq_aligned(traceback_path.size(),'-');
//
//        idx=0;
//        for(unsigned int j=0;j<traceback_path.size();j++){
//
//            if(traceback_path.at(j)==MATCH_CHAR){
//
//                seq_aligned[j]=seq_not_aligned.at(idx);
//                idx++;
//
//            }else if(traceback_path.at(j)==GAP_X_CHAR){
//
//                seq_aligned[j]='-';
//
//            }else if(traceback_path.at(j)==GAP_Y_CHAR){
//
//                seq_aligned[j]=seq_not_aligned.at(idx);
//                idx++;
//
//            }else{
//                perror("ERROR in align_seq");
//            }
//        }
//
//        MSA_out.emplace_back(std::make_pair(seq_name,seq_aligned));
//    }
//
//    return MSA_out;
//}
//std::vector<std::pair<std::string,std::string>> pPIP::build_MSA(std::string traceback_path,
//                                                                std::vector<std::pair<std::string, std::string>> &MSA_L,
//                                                                std::vector<std::pair<std::string,std::string>> &MSA_R){
//    std::vector<std::pair<std::string,std::string>> MSA;
//    std::vector<std::pair<std::string,std::string>> MSA_L_out;
//    std::vector<std::pair<std::string,std::string>> MSA_R_out;
//
//    MSA_L_out=align_seq_left(MSA_L,traceback_path);
//    MSA_R_out=align_seq_right(MSA_R,traceback_path);
//
//    //for (unsigned int i=0;i<MSA_L_out.size();i++){
//    for (const auto &m : MSA_L_out){
//        //MSA.push_back(MSA_L_out.at(i));
//        MSA.push_back(m);
//    }
//
//    //for (unsigned int i=0;i<MSA_R_out.size();i++){
//    for (const auto &m : MSA_L_out){
//        //MSA.push_back(MSA_R_out.at(i));
//        MSA.push_back(m);
//    }
//
//    return MSA;
//}
void pPIP::build_MSA(bpp::Node *node, std::string traceback_path){


    auto sons = node->getSons();
    auto s1 = sons.at(0);
    auto s2 = sons.at(1);
    unsigned long s1ID = (unsigned long)s1->getId();
    unsigned long s2ID = (unsigned long)s2->getId();

    //std::vector< std::string > *MSA_L=&(_MSA[s1]);
    //std::vector< std::string > *MSA_R=&(_MSA[s2]);
    std::vector< std::string > *MSA_L=&(_MSA.at(s1ID));
    std::vector< std::string > *MSA_R=&(_MSA.at(s2ID));

    unsigned long lenColL=MSA_L->at(0).size();
    unsigned long lenColR=MSA_R->at(0).size();

    std::vector<std::string> MSA;

    int idx_i=0;
    int idx_j=0;
    for(unsigned int j=0;j<traceback_path.size();j++){

        if(traceback_path.at(j)==MATCH_CHAR){

            MSA.push_back(MSA_L->at(idx_i)+MSA_R->at(idx_j));
            idx_i++;
            idx_j++;

        }else if(traceback_path.at(j)==GAP_X_CHAR){

            std::string gapCol(lenColR,GAP_CHAR);
            MSA.push_back(MSA_L->at(idx_i)+gapCol);
            idx_i++;

        }else if(traceback_path.at(j)==GAP_Y_CHAR){

            std::string gapCol(lenColL,GAP_CHAR);
            MSA.push_back(gapCol+MSA_R->at(idx_j));
            idx_j++;

        }else{
            perror("ERROR");
        }
    }

    //_MSA[node]=MSA;
    _MSA.at(node->getId())=MSA;
}
void pPIP::setMSAsequenceNames(bpp::Node *node){

    auto sons = node->getSons();
    int s1ID = sons.at(LEFT)->getId();
    int s2ID = sons.at(RIGHT)->getId();

    std::vector<std::string> seqNames;

    //for(int i=0;i<_seqNames[s1].size();i++){
    for(int i=0;i<_seqNames.at(s1ID).size();i++){
        //seqNames.push_back(_seqNames[s1].at(i));
        seqNames.push_back(_seqNames.at(s1ID).at(i));
    }

    //for(int i=0;i<_seqNames[s2].size();i++){
    for(int i=0;i<_seqNames.at(s2ID).size();i++){
        //seqNames.push_back(_seqNames[s2].at(i));
        seqNames.push_back(_seqNames.at(s2ID).at(i));
    }

    //_seqNames[node]=seqNames;
    _seqNames.at(node->getId())=seqNames;

}
void pPIP::setMSAsequenceNames(bpp::Node *node,std::string seqname){

    std::vector<std::string> seqNames;

    seqNames.push_back(seqname);

    _seqNames.at(node->getId())=seqNames;

}
void pPIP::setMSAleaves(bpp::Node *node,const std::string &MSA){

    /* convert a string into a vector of single char strings */
    std::vector<std::string> msa;
    msa.resize(MSA.size());
    for(int i=0;i<MSA.size();i++){
        std::string s(1, MSA.at(i));
        msa.at(i)=s;
    }

    _MSA.at(node->getId())=msa;

}

/*
void pPIP::setIndicatorFun(bpp::Node *node){

    int idx;

    std::vector<std::string> *MSA = &(_MSA.at(node->getId()));

    std::vector<ColMatrix<double>> indicatorFunctions;

    indicatorFunctions.resize(MSA->size());

    for(unsigned i=0;i<MSA->size();i++){

        bpp::ColMatrix<double> indFun;
        indFun.resize(_extendedAlphabetSize,1);
        bpp::MatrixTools::fill(indFun,0.0);

        idx=_alphabet->charToInt(MSA->at(i));

        indFun(idx,0)=1.0;
        indicatorFunctions.at(i)=indFun;

    }

    _fv.at(node->getId())=indicatorFunctions;

}
 */
/*
void pPIP::setLKxyLeaves(bpp::Node *node){

    //int len=_MSA[node].size();
    int len=_MSA.at(node->getId()).size();

    Vdouble lkxy;

    lkxy.resize(len);
    lkxy.assign(len,0.0);

    //_lkxy[node]=lkxy;
    _lkxy.at(node->getId())=lkxy;

}
 */
void pPIP::_setNu() {

    _nu  = _lambda * (_tau + 1 / _mu);

}
void pPIP::_setLambda(double lambda){
    _lambda=lambda;
}
void pPIP::_setMu(double mu){

    if (fabs(mu) < SMALL_DOUBLE) {
        perror("ERROR: mu is too small");
    }

    _mu=mu;
}
void pPIP::_setPi(const Vdouble &pi){
    _pi.resize(pi.size(),1);
    for(int i=0;i<pi.size();i++){
        _pi(i,0)=pi.at(i);
    }
    //_pi=pi;
}
void pPIP::_setTau() {
    _tau = _tree->getTotalLength();
}
void pPIP::_setAllIotas(UtreeBppUtils::treemap *tm,std::vector<tshlib::VirtualNode *> &listNodes) {
    double T;

    T = _tau + 1/_mu;

    if (fabs(T) < SMALL_DOUBLE) {
        perror("ERROR in set_iota: T too small");
    }

    for (auto &vnode:listNodes) {

        auto node = _tree->getNode(tm->right.at(vnode),false);

        if (!node->hasFather()) {
            //_iota[node]=(1/_mu)/T;
            _iota.at(node->getId())=(1/_mu)/T;
        } else {
            //_iota[node]=node->getDistanceToFather()/T;
            _iota.at(node->getId())=node->getDistanceToFather()/T;
        }
    }
}
void pPIP::_setAllBetas(UtreeBppUtils::treemap *tm,std::vector<tshlib::VirtualNode *> &listNodes) {

    //TreeTemplate<Node> ttree(*tree);

//    if (fabs(mu) < SMALL_DOUBLE) {
//        perror("ERROR in set_iota: mu too small");
//    }

    for (auto &vnode:listNodes) {

        auto node = _tree->getNode(tm->right.at(vnode),false);

        if (!node->hasFather()) {
            //_beta[node]=1.0;
            _beta.at(node->getId())=1.0;
        } else {

            double muT = _mu * node->getDistanceToFather();
            if (fabs(muT) < SMALL_DOUBLE) {
                perror("ERROR mu * T is too small");
            }

            //_beta[node]= (1.0 - exp(-_mu * node->getDistanceToFather())) / (_mu * node->getDistanceToFather());
            _beta.at(node->getId())= (1.0 - exp(-_mu * node->getDistanceToFather())) / muT;
        }

    }


}
void pPIP::_computePr(UtreeBppUtils::treemap *tm,std::vector<tshlib::VirtualNode *> &listNodes) {

    for (auto &vnode:listNodes) {

        auto node = _tree->getNode(tm->right.at(vnode),false);

        if (!node->hasFather()) {

        } else {
            //_pr[node] = MatrixBppUtils::Eigen2Matrix(vnode->getPr());
            _pr.at(node->getId()) = MatrixBppUtils::Eigen2Matrix(vnode->getPr());
        }

    }


//    for (auto &vnode:listNodes) {
//
//        auto node = tree_->getNode(tm->right.at(vnode),false);
//
//        if (!node->hasFather()) {
//
//        } else {
//
//
//            std::cout<<node->getName()<<std::endl;
//
//            bpp::RowMatrix<double> M = _pr[node];
//
//            for(int i=0;i<M.getNumberOfRows();i++){
//                for(int j=0;j<M.getNumberOfRows();j++){
//                    std::cout<<M(i,j)<<" ";
//                }
//                std::cout<<std::endl;
//            }
//            std::cout<<std::endl;
//
//        }
//
//    }



}
//void pPIP::DP3D_PIP(bpp::Node *node,
//                    UtreeBppUtils::treemap *tm,
//                    double gamma_rate,
//                    bool local){
//
//    //TODO: place as argument
//    bool randomSeed = true;
//
//    //TODO: re-implement gamma distribution
//    double lambda_gamma = _lambda * gamma_rate;
//    double mu_gamma = _mu * gamma_rate;
//
//    if(local){
//        //bpp::Tree *subtree = tree_->cloneSubtree(tm->right.at(node));
//        bpp::Tree *subtree = _tree->cloneSubtree(node->getId());
//        _tau = subtree->getTotalLength();
//        _setNu();
//    }else{
//    }
//
//    unsigned long up_corner_i;
//    unsigned long up_corner_j;
//    unsigned long bot_corner_i;
//    unsigned long bot_corner_j;
//    unsigned long lw;
//    unsigned long h,w;
//
//    auto sons = node->getSons();
//
//    auto s1 = sons.at(0);
//    auto s2 = sons.at(1);
//
//    int s1ID = s1->getId();
//    int s2ID = s2->getId();
//
//    int nodeID = node->getId();
//
//    //h = _MSA[s1].size()+1;
//    //w = _MSA[s2].size()+1;
//    h = _MSA.at(s1ID).size()+1;
//    w = _MSA.at(s2ID).size()+1;
//
//    unsigned long d=(h-1)+(w-1)+1;
//
//    double pc0;
//
//    std::string *sLs;
//    std::string *sRs;
//    std::string col_gap_Ls;
//    std::string col_gap_Rs;
//
//    //col_gap_Ls=createGapCol(_MSA[s1].size());
//    //col_gap_Rs=createGapCol(_MSA[s2].size());
//    col_gap_Ls=createGapCol(_MSA.at(s1ID).size());
//    col_gap_Rs=createGapCol(_MSA.at(s2ID).size());
//
//    signed long seed;
//    if(randomSeed){
//        seed = std::chrono::system_clock::now().time_since_epoch().count();
//
//    }else{
//        seed = 0;
//    }
//
//    std::default_random_engine generator(seed);
//    std::uniform_real_distribution<double> distribution(0.0,1.0);
//
//    auto epsilon=DBL_EPSILON;
//
//    unsigned ilk=0;
//    unsigned ifv=0;
//
//
//    (_fv.at(nodeID)).resize(d+1);
//
//    //Vdouble Fvgap;
//    bpp::ColMatrix<double> Fvgap;
//    double p0;
//    if(local){
//        Fvgap = computeLKgapColLocal(node,pc0,p0);
//    }else{
//        /*
//        pc0 = compute_pr_gap_all_edges_s(node,
//                                         col_gap_Ls,
//                                         col_gap_Rs,
//                                         pi,
//                                         originalAlphabetSize,
//                                         alphabet);
//                                         */
//    }
//
//    //(_fv[node]).push_back(Fvgap);
//    //(_fv.at(nodeID)).push_back(Fvgap);
//    (_fv.at(nodeID)).at(ifv)=Fvgap;
//    ifv++;
//
//
//
//    auto** LogM = new double*[2];
//    auto** LogX = new double*[2];
//    auto** LogY = new double*[2];
//
//    auto** TR = new int*[d];
//
//    LogM[0] = new double[int((w*(h+1))/2)];
//    LogX[0] = new double[int((w*(h+1))/2)];
//    LogY[0] = new double[int((w*(h+1))/2)];
//    LogM[1] = new double[int((w*(h+1))/2)];
//    LogX[1] = new double[int((w*(h+1))/2)];
//    LogY[1] = new double[int((w*(h+1))/2)];
//
//    LogM[0][0]=_nu*(pc0-1.0);
//    LogX[0][0]=_nu*(pc0-1.0);
//    LogY[0][0]=_nu*(pc0-1.0);
//
//    TR[0] = new int[1]();
//    TR[0][0]=STOP_STATE;
//
//    double max_of_3;
//    double max_lk=-std::numeric_limits<double>::infinity();
//    double prev_max_lk=-std::numeric_limits<double>::infinity();
//    signed long level_max_lk=INT_MIN;
//    double val;
//    unsigned long m_binary_this;
//    unsigned long m_binary_prev;
//
//    double valM;
//    double valX;
//    double valY;
//
//    signed long idx;
//
//    unsigned long coordSeq_1;
//    unsigned long coordSeq_2;
//    unsigned long coordTriangle_this_i;
//    unsigned long coordTriangle_this_j;
//    unsigned long coordTriangle_prev_i;
//    unsigned long coordTriangle_prev_j;
//
//    double score=-std::numeric_limits<double>::infinity();
//    signed long depth;
//
//    bool flag_exit=false;
//    unsigned long last_d=d-1;
//    unsigned long size_tr,tr_up_i,tr_up_j,tr_down_i,tr_down_j;
//
//    Vdouble lkxy;
//    lkxy.resize(d+1);
//    lkxy.assign(d+1,0.0);
//
//
//    std::cout<<"Aligning at node:"<<node->getName()<<std::endl;
//
//
//    double lkx,lky;
//    for(unsigned long m=1;m<d;m++){
//
//        std::cout<<"m:"<<m<<std::endl;
//
//        //Vdouble Fvmatch;
//        //Vdouble Fvgapx;
//        //Vdouble Fvgapy;
//        bpp::ColMatrix<double> Fvmatch;
//        bpp::ColMatrix<double> Fvgapx;
//        bpp::ColMatrix<double> Fvgapy;
//
//        if(flag_exit){
//            break;
//        }
//
//        m_binary_this=m%2;
//        m_binary_prev=(m+1)%2;
//        //***************************************************************************************
//        //***************************************************************************************
//        set_indeces_M(up_corner_i,up_corner_j,bot_corner_i,bot_corner_j,m,h,w);
//
//        if(checkboundary(up_corner_i,up_corner_j,bot_corner_i,bot_corner_j,h,w)){
//
//            lw=0;
//            for(unsigned long i=up_corner_i;i<=bot_corner_i;i++){
//
//                coordTriangle_this_i=i;
//                coordSeq_1=coordTriangle_this_i-1;
//                coordTriangle_prev_i=coordTriangle_this_i-1;
//
//                //sLs=&(_MSA[s1].at(coordSeq_1));
//                sLs=&(_MSA.at(s1ID).at(coordSeq_1));
//
//                for(int j=0;j<=lw;j++){
//
//                    coordTriangle_this_j=up_corner_j-j;
//                    coordSeq_2=coordTriangle_this_j-1;
//                    coordTriangle_prev_j=coordTriangle_this_j-1;
//
//                    //sRs=&(_MSA[s2].at(coordSeq_2));
//                    sRs=&(_MSA.at(s2ID).at(coordSeq_2));
//
//                    idx=get_indices_M(coordTriangle_prev_i,coordTriangle_prev_j,up_corner_i,up_corner_j,bot_corner_i,bot_corner_j,m-1,h,w);
//                    if(idx>=0){
//                        valM=LogM[m_binary_prev][idx];
//                    }else{
//                        valM=-std::numeric_limits<double>::infinity();
//                    }
//
//                    idx=get_indices_X(coordTriangle_prev_i,coordTriangle_prev_j,up_corner_i,up_corner_j,bot_corner_i,bot_corner_j,m-1,h,w);
//                    if(idx>=0){
//                        valX=LogX[m_binary_prev][idx];
//                    }else{
//                        valX=-std::numeric_limits<double>::infinity();
//                    }
//
//                    idx=get_indices_Y(coordTriangle_prev_i,coordTriangle_prev_j,up_corner_i,up_corner_j,bot_corner_i,bot_corner_j,m-1,h,w);
//                    if(idx>=0){
//                        valY=LogY[m_binary_prev][idx];
//                    }else{
//                        valY=-std::numeric_limits<double>::infinity();
//                    }
//
//                    if(std::isinf(valM) && std::isinf(valX) && std::isinf(valY)){
//                        exit(EXIT_FAILURE);
//                    }
//
//                    std::cout<<"M"<<std::endl;
//
//                    if(local){
//                        Fvmatch=computeLKmatchLocal(valM,
//                                                 valX,
//                                                 valY,
//                                                 node,
//                                                 coordSeq_1,
//                                                 coordSeq_2,
//                                                 m,
//                                                 val);
//                    }else{
//                        /*
//                        val=computeLK_M_all_edges_s_opt(valM,
//                                                        valX,
//                                                        valY,
//                                                        nu,
//                                                        node,
//                                                        sLs, sRs,
//                                                        pi,
//                                                        m,
//                                                        lkM,
//                                                        originalAlphabetSize, alphabet);
//                                                        */
//                    }
//
//                    if(std::isinf(val)){
//                        exit(EXIT_FAILURE);
//                    }
//
//                    if(std::isnan(val)){
//                        exit(EXIT_FAILURE);
//                    }
//
//                    idx=get_indices_M(coordTriangle_this_i,coordTriangle_this_j,up_corner_i,
//                                      up_corner_j,bot_corner_i,bot_corner_j,m,h,w);
//
//                    LogM[m_binary_this][idx]=val;
//                }
//                lw++;
//            }
//
//        }
//        //***************************************************************************************
//        //***************************************************************************************
//        set_indeces_X(up_corner_i,up_corner_j,bot_corner_i,bot_corner_j,m,h,w);
//        tr_down_i=bot_corner_i;
//        tr_down_j=bot_corner_j;
//        if(checkboundary(up_corner_i,up_corner_j,bot_corner_i,bot_corner_j,h,w)){
//
//            lw=0;
//            for(unsigned long i=up_corner_i;i<=bot_corner_i;i++){
//
//                coordTriangle_this_i=i;
//                coordTriangle_prev_i=coordTriangle_this_i-1;
//                coordSeq_1=coordTriangle_this_i-1;
//
//                //sLs=&(_MSA[s1].at(coordSeq_1));
//                sLs=&(_MSA.at(s1ID).at(coordSeq_1));
//
//                for(int j=0;j<=lw;j++){
//
//                    coordTriangle_this_j=up_corner_j-j;
//                    coordTriangle_prev_j=coordTriangle_this_j;
//
//                    idx=get_indices_M(coordTriangle_prev_i,coordTriangle_prev_j,up_corner_i,
//                                      up_corner_j,bot_corner_i,bot_corner_j,m-1,h,w);
//                    if(idx>=0){
//                        valM=LogM[m_binary_prev][idx];
//                    }else{
//                        valM=-std::numeric_limits<double>::infinity();
//                    }
//
//                    idx=get_indices_X(coordTriangle_prev_i,coordTriangle_prev_j,up_corner_i,
//                                      up_corner_j,bot_corner_i,bot_corner_j,m-1,h,w);
//                    if(idx>=0){
//                        valX=LogX[m_binary_prev][idx];
//                    }else{
//                        valX=-std::numeric_limits<double>::infinity();
//                    }
//
//                    idx=get_indices_Y(coordTriangle_prev_i,coordTriangle_prev_j,up_corner_i,
//                                      up_corner_j,bot_corner_i,bot_corner_j,m-1,h,w);
//                    if(idx>=0){
//                        valY=LogY[m_binary_prev][idx];
//                    }else{
//                        valY=-std::numeric_limits<double>::infinity();
//                    }
//
//                    if(std::isinf(valM) && std::isinf(valX) && std::isinf(valY)){
//                        exit(EXIT_FAILURE);
//                    }
//
//                    std::cout<<"X"<<std::endl;
//
//                    if(local){
//                        Fvgapx=computeLKgapxLocal(valM,
//                                                 valX,
//                                                 valY,
//                                                 node,
//                                                 coordSeq_1,
//                                                 coordSeq_2,
//                                                 m,
//                                                val,lkx);
//                    }else{
//
//                        /*val=computeLK_X_all_edges_s_opt(valM,
//                                                        valX,
//                                                        valY,
//                                                        nu,
//                                                        node,
//                                                        sLs, col_gap_Rs,
//                                                        pi,
//                                                        m,
//                                                        lkX,
//                                                        originalAlphabetSize, alphabet);
//                                                        */
//                    }
//
//                    if(std::isinf(val)){
//                        exit(EXIT_FAILURE);
//                    }
//
//                    if(std::isnan(val)){
//                        exit(EXIT_FAILURE);
//                    }
//
//                    idx=get_indices_X(coordTriangle_this_i,coordTriangle_this_j,up_corner_i,
//                                      up_corner_j,bot_corner_i,bot_corner_j,m,h,w);
//
//                    LogX[m_binary_this][idx]=val;
//                }
//                lw++;
//            }
//
//        }
//        //***************************************************************************************
//        //***************************************************************************************
//        set_indeces_Y(up_corner_i,up_corner_j,bot_corner_i,bot_corner_j,m,h,w);
//        tr_up_i=up_corner_i;
//        tr_up_j=up_corner_j;
//        if(checkboundary(up_corner_i,up_corner_j,bot_corner_i,bot_corner_j,h,w)){
//
//            lw=0;
//            for(unsigned long i=up_corner_i;i<=bot_corner_i;i++){
//                coordTriangle_this_i=i;
//                coordTriangle_prev_i=coordTriangle_this_i;
//                for(int j=0;j<=lw;j++){
//
//                    coordTriangle_this_j=up_corner_j-j;
//                    coordTriangle_prev_j=coordTriangle_this_j-1;
//                    coordSeq_2=coordTriangle_this_j-1;
//
//                    //sRs=&(_MSA[s2].at(coordSeq_2));
//                    sRs=&(_MSA.at(s2ID).at(coordSeq_2));
//
//                    idx=get_indices_M(coordTriangle_prev_i,coordTriangle_prev_j,up_corner_i,
//                                      up_corner_j,bot_corner_i,bot_corner_j,m-1,h,w);
//                    if(idx>=0){
//                        valM=LogM[m_binary_prev][idx];
//                    }else{
//                        valM=-std::numeric_limits<double>::infinity();
//                    }
//
//                    idx=get_indices_X(coordTriangle_prev_i,coordTriangle_prev_j,up_corner_i,
//                                      up_corner_j,bot_corner_i,bot_corner_j,m-1,h,w);
//                    if(idx>=0){
//                        valX=LogX[m_binary_prev][idx];
//                    }else{
//                        valX=-std::numeric_limits<double>::infinity();
//                    }
//
//                    idx=get_indices_Y(coordTriangle_prev_i,coordTriangle_prev_j,up_corner_i,
//                                      up_corner_j,bot_corner_i,bot_corner_j,m-1,h,w);
//                    if(idx>=0){
//                        valY=LogY[m_binary_prev][idx];
//                    }else{
//                        valY=-std::numeric_limits<double>::infinity();
//                    }
//
//                    if(std::isinf(valM) && std::isinf(valX) && std::isinf(valY)){
//                        exit(EXIT_FAILURE);
//                    }
//
//                    std::cout<<"Y"<<std::endl;
//
//                    if(local){
//                        Fvgapy=computeLKgapyLocal(valM,
//                                                 valX,
//                                                 valY,
//                                                 node,
//                                                 coordSeq_1,
//                                                 coordSeq_2,
//                                                 m,
//                                                 val,lky);
//                    }else{
//                        /*
//                        val=computeLK_Y_all_edges_s_opt(valM,
//                                                        valX,
//                                                        valY,
//                                                        nu,
//                                                        node,
//                                                        col_gap_Ls, sRs,
//                                                        pi,
//                                                        m,
//                                                        lkY,
//                                                        originalAlphabetSize, alphabet);
//                                                        */
//                    }
//
//                    if(std::isinf(val)){
//                        exit(EXIT_FAILURE);
//                    }
//
//
//                    if(std::isnan(val)){
//                        exit(EXIT_FAILURE);
//                    }
//
//                    idx=get_indices_Y(coordTriangle_this_i,coordTriangle_this_j,up_corner_i,
//                                      up_corner_j,bot_corner_i,bot_corner_j,m,h,w);
//
//                    LogY[m_binary_this][idx]=val;
//                }
//                lw++;
//            }
//
//        }
//        //***************************************************************************************
//        //***************************************************************************************
//
//        std::cout<<"TR"<<std::endl;
//
//        size_tr=(unsigned long)ceil((tr_down_i-tr_up_i+1)*(tr_up_j-tr_down_j+1+1)/2);
//        TR[m] = new int[size_tr](); /*TODO: optimize size TR*/
//        memset(TR[m],0,size_tr*sizeof(TR[m][0]));
//        set_indeces_T(up_corner_i,up_corner_j,bot_corner_i,bot_corner_j,m,h,w);
//
//        if(checkboundary(up_corner_i,up_corner_j,bot_corner_i,bot_corner_j,h,w)){
//
//            lw=0;
//            for(unsigned long i=up_corner_i;i<=bot_corner_i;i++){
//                coordTriangle_this_i=i;
//                for(int j=0;j<=lw;j++){
//                    coordTriangle_this_j=up_corner_j-j;
//
//                    double mval;
//                    double xval;
//                    double yval;
//
//                    idx=get_indices_M(coordTriangle_this_i,coordTriangle_this_j,up_corner_i,
//                                      up_corner_j,bot_corner_i,bot_corner_j,m,h,w);
//                    if(idx>=0){
//                        mval=LogM[m_binary_this][idx];
//                    }else{
//                        mval=-std::numeric_limits<double>::infinity();
//                    }
//
//                    idx=get_indices_X(coordTriangle_this_i,coordTriangle_this_j,up_corner_i,
//                                      up_corner_j,bot_corner_i,bot_corner_j,m,h,w);
//                    if(idx>=0){
//                        xval=LogX[m_binary_this][idx];
//                    }else{
//                        xval=-std::numeric_limits<double>::infinity();
//                    }
//
//                    idx=get_indices_Y(coordTriangle_this_i,coordTriangle_this_j,up_corner_i,
//                                      up_corner_j,bot_corner_i,bot_corner_j,m,h,w);
//                    if(idx>=0){
//                        yval=LogY[m_binary_this][idx];
//                    }else{
//                        yval=-std::numeric_limits<double>::infinity();
//                    }
//
//                    mval=fabs(mval)<epsilon?-std::numeric_limits<double>::infinity():mval;
//                    xval=fabs(xval)<epsilon?-std::numeric_limits<double>::infinity():xval;
//                    yval=fabs(yval)<epsilon?-std::numeric_limits<double>::infinity():yval;
//
//                    //.-----.------.------.-------.------.-------.-------.--------.-------.--------.//
//                    int ttrr;
//
//                    ttrr=index_of_max(mval,xval,yval,epsilon,generator,distribution);
//
//
//                    std::cout<<"TRII:"<<_fv.at(nodeID).size()<<":"<<ifv<<std::endl;
//                    std::cout<<_fv.at(nodeID).size()<<":"<<ifv<<std::endl;
//                    std::cout<<lkxy.size()<<":"<<ilk<<std::endl;
//
//                    switch(ttrr){
//                        case MATCH_STATE:
//                            //_fv[node].push_back(Fvmatch);
//                            //_fv.at(node->getId()).push_back(Fvmatch);
//                            _fv.at(nodeID).at(ifv)=Fvmatch;
//                            ifv++;
//                            break;
//                        case GAP_X_STATE:
//                            //_fv[node].push_back(Fvgapx);
//                            //_fv.at(node->getId()).push_back(Fvgapx);
//                            _fv.at(nodeID).at(ifv)=Fvgapx;
//                            ifv++;
//                            //lkxy.push_back(lkx);
//                            lkxy.at(ilk)=lkx;
//                            ilk++;
//                            break;
//                        case GAP_Y_STATE:
//                            //_fv[node].push_back(Fvgapy);
//                            //_fv.at(nodeID).push_back(Fvgapy);
//                            _fv.at(nodeID).at(ifv)=Fvgapy;
//                            ifv++;
//                            //lkxy.push_back(lky);
//                            lkxy.at(ilk)=lky;
//                            ilk++;
//                            break;
//                        default:
//                            perror("ERROR!!!");
//                            exit(EXIT_FAILURE);
//                    }
//
//
//                    std::cout<<"TRIII"<<std::endl;
//
//
//
//                    idx=get_indices_T(coordTriangle_this_i,coordTriangle_this_j,up_corner_i,
//                                      up_corner_j,bot_corner_i,bot_corner_j,m,h,w);
//
//                    if(TR[m][idx]!=0){
//                        exit(EXIT_FAILURE);
//                    }
//
//                    TR[m][idx]=ttrr;
//
//                    if( (coordTriangle_this_i==(h-1)) & (coordTriangle_this_j==(w-1)) ){
//
//                        max_of_3=max_of_three(mval,xval,yval,epsilon);
//
//                        //TODO: check this part
//                        if(max_of_3>score){
//                            score=max_of_3;
//                            level_max_lk=m;
//                        }
//
//
//
//                    }
//
//                }
//                lw++;
//            }
//        }
//    }
//
//    //TODO:resize
//    //_lkxy[node]=lkxy;
//
//    lkxy.resize(depth+1);
//
//    _lkxy.at(nodeID)=lkxy;
//
//    depth=level_max_lk;
//
//    (_fv.at(nodeID)).resize(depth+1);
//
//    _score=score;
//
//    std::string traceback_path ((unsigned long)depth, ' ');
//    unsigned long id1=h-1;
//    unsigned long id2=w-1;
//    for(long lev=depth;lev>0;lev--){
//        set_indeces_T(up_corner_i,up_corner_j,bot_corner_i,bot_corner_j,(unsigned)lev,h,w);
//        idx=get_indices_T(id1,id2,up_corner_i,up_corner_j,bot_corner_i,bot_corner_j,(unsigned)lev,h,w);
//
//        switch(TR[lev][idx]){
//            case MATCH_STATE:
//                id1=id1-1;
//                id2=id2-1;
//                traceback_path[lev-1]=MATCH_CHAR;
//                break;
//            case GAP_X_STATE:
//                id1=id1-1;
//                traceback_path[lev-1]=GAP_X_CHAR;
//                break;
//            case GAP_Y_STATE:
//                id2=id2-1;
//                traceback_path[lev-1]=GAP_Y_CHAR;
//                break;
//            default:
//                perror("ERROR in alignment_reconstruction !!!");
//                exit(EXIT_FAILURE);
//        }
//    }
//
//    _traceback_path=traceback_path;
//
//    build_MSA(node,traceback_path);
//
//    setMSAsequenceNames(node);
//
//
//    free(LogM[1]);
//    free(LogM[0]);
//    free(LogM);
//
//    free(LogX[1]);
//    free(LogX[0]);
//    free(LogX);
//
//    free(LogY[1]);
//    free(LogY[0]);
//    free(LogY);
//
//    for(long i=last_d;i>=0;i--){
//        free(TR[i]);
//    }
//    free(TR);
//
//}

//void pPIP::PIPAligner(UtreeBppUtils::treemap *tm,
//                      std::vector<tshlib::VirtualNode *> &list_vnode_to_root,
//                      bpp::SequenceContainer *sequences,
//                      double gamma_rate,
//                      bool local) {
//
//    for (auto &vnode:list_vnode_to_root) {
//
//        auto node = _tree->getNode(tm->right.at(vnode),false);
//
//        if(node->isLeaf()){
//
//            std::string seqname = sequences->getSequencesNames().at((unsigned long)vnode->vnode_seqid);
//
//            setMSAsequenceNames(node,seqname);
//
//            setMSAleaves(node,sequences->getSequence(seqname).toString());
//
//            setIndicatorFun(node);
//
//            setLKxyLeaves(node);
//
//        }else{
//
//            DP3D_PIP(node, tm, gamma_rate, local);
//
//        }
//    }
//
//}
bpp::ColMatrix<double> pPIP::fv_observed(std::string &s, unsigned long &idx){
    bpp::ColMatrix<double> fv;
    int ii;
    char ch=s[idx];

    fv.resize(1,_extendedAlphabetSize);
    bpp::MatrixTools::fill(fv,0.0);

    ii=_alphabet->charToInt(&ch);
    //ii=ii<0?_alphabetSize:ii;

    fv(ii,0)=1.0;
    idx++;

    return fv;
}
bpp::ColMatrix<double> pPIP::go_down(bpp::Node *tree,std::string &s, unsigned long &idx){
    //Eigen::VectorXd fv;
    //Eigen::VectorXd fvL;
    //Eigen::VectorXd fvR;
    bpp::ColMatrix<double> fv;
    bpp::ColMatrix<double> fvL;
    bpp::ColMatrix<double> fvR;
    bpp::ColMatrix<double> PrfvL;
    bpp::ColMatrix<double> PrfvR;


    if(tree->isLeaf()){

        fv=fv_observed(s,idx);

    }else{

        auto sons = tree->getSons();
        int sonLeftID = sons.at(LEFT)->getId();
        int sonRightID = sons.at(RIGHT)->getId();

        fvL=go_down(sons.at(LEFT),s,idx);
        fvR=go_down(sons.at(RIGHT),s,idx);

        bpp::MatrixTools::mult(_pr.at(sonLeftID),fvL,PrfvL);
        bpp::MatrixTools::mult(_pr.at(sonRightID),fvR,PrfvR);
        bpp::MatrixTools::hadamardMult(PrfvL,PrfvR,fv);

        //fv=(tree->getNodeLeft()->getPr()*fvL).cwiseProduct(tree->getNodeRight()->getPr()*fvR);

    }

    return fv;
}
void pPIP::allgaps(bpp::Node *tree,std::string &s, unsigned long &idx,bool &flag){

    if(tree->isLeaf()){
        char ch=s[idx];

        idx++;

        if(ch!='-'){
            flag=false;
        }

    }else{

        auto sons = tree->getSons();
        allgaps(sons.at(LEFT),s,idx,flag);
        allgaps(sons.at(RIGHT),s,idx,flag);
    }

}
double pPIP::compute_lk_gap_down(bpp::Node *tree,std::string &s){

    double pr=0;
    double pL=0;
    double pR=0;
    unsigned long idx;
    //Eigen::VectorXd fvL;
    //Eigen::VectorXd fvR;
    //Eigen::VectorXd fv
    bpp::ColMatrix<double> fv;
    bpp::ColMatrix<double> fvL;
    bpp::ColMatrix<double> fvR;
    double fv0;

    auto sons = tree->getSons();

    int nodeID = tree->getId();
    //int sonLeftID = sons.at(LEFT)->getId();
    //int sonRightID = sons.at(RIGHT)->getId();

    if(tree->isLeaf()){
        idx=0;
        fv=go_down(tree,s,idx);
        //fv0=fv.dot(pi);
        fv0=MatrixBppUtils::dotProd(fv,_pi);
        //pr=tree->getIota()-tree->getIota()*tree->getBeta()+tree->getIota()*tree->getBeta()*fv0;
        pr=_iota.at(nodeID)-_iota.at(nodeID)*_beta.at(nodeID)+_iota.at(nodeID)*_beta.at(nodeID)*fv0;
        return pr;
    }else{
        idx=0;
        fv=go_down(tree,s,idx);
        //fv0=fv.dot(pi);
        fv0=MatrixBppUtils::dotProd(fv,_pi);
        //pr=tree->getIota()-tree->getIota()*tree->getBeta()+tree->getIota()*tree->getBeta()*fv0;
        pr=_iota.at(nodeID)-_iota.at(nodeID)*_beta.at(nodeID)+_iota.at(nodeID)*_beta.at(nodeID)*fv0;
        bool flagL=true;
        bool flagR=true;
        idx=0;
        allgaps(sons.at(LEFT),s,idx,flagL);
        unsigned long ixx=idx;
        allgaps(sons.at(RIGHT),s,idx,flagR);
        unsigned long len;

        std::string sL;
        len=ixx;
        sL=s.substr(0,len);
        pL=compute_lk_gap_down(sons.at(LEFT),sL);

        std::string sR;
        sR=s.substr(ixx);
        pR=compute_lk_gap_down(sons.at(RIGHT),sR);
    }

    return pr+pL+pR;
}
double pPIP::compute_lk_down(bpp::Node *tree,std::string &s){

    double pr;
    unsigned long idx;
//    Eigen::VectorXd fvL;
//    Eigen::VectorXd fvR;
//    Eigen::VectorXd fv;
    bpp::ColMatrix<double> fv;
    bpp::ColMatrix<double> fvL;
    bpp::ColMatrix<double> fvR;
    double fv0;

    auto sons = tree->getSons();

    int nodeID = tree->getId();

    if(tree->isLeaf()){

        idx=0;
        fv=go_down(tree,s,idx);
        //fv0=fv.dot(pi);
        fv0=MatrixBppUtils::dotProd(fv,_pi);
        //pr=tree->getIota()*tree->getBeta()*fv0;
        pr=_iota.at(nodeID)*_beta.at(nodeID)*fv0;

        return pr;

    }else{

        idx=0;
        fv=go_down(tree,s,idx);
        //fv0=fv.dot(pi);
        fv0=MatrixBppUtils::dotProd(fv,_pi);
        //pr=tree->getIota()*tree->getBeta()*fv0;
        pr=_iota.at(nodeID)*_beta.at(nodeID)*fv0;

        bool flagL=true;
        bool flagR=true;
        idx=0;
        allgaps(sons.at(LEFT),s,idx,flagL);
        unsigned long ixx=idx;
        allgaps(sons.at(RIGHT),s,idx,flagR);

        unsigned long len;
        if(flagR){
            std::string sL;//=stringFromSequence(s);
            len=ixx;
            sL=s.substr(0,len);
            return pr + compute_lk_down(sons.at(LEFT),sL);
        }

        if(flagL){
            std::string sR;//=stringFromSequence(s);
            sR=s.substr(ixx);
            return pr + compute_lk_down(sons.at(RIGHT),sR);
        }

    }

    return pr;
}
double pPIP::compute_pr_gap_local_tree_s(bpp::Node *tree, std::string &sL, std::string &sR){
    double fv0;
    double pr;
    //Eigen::VectorXd fvL;
    //Eigen::VectorXd fvR;
    //Eigen::VectorXd fv;
    bpp::ColMatrix<double> fvL;
    bpp::ColMatrix<double> fvR;
    bpp::ColMatrix<double> PrfvL;
    bpp::ColMatrix<double> PrfvR;
    bpp::ColMatrix<double> fv;
    unsigned long idx;

    auto sons = tree->getSons();

    int sonLeftID = sons.at(LEFT)->getId();
    int sonRightID = sons.at(RIGHT)->getId();

    idx=0;
    fvL=go_down(sons.at(LEFT),sL,idx);

    idx=0;
    fvR=go_down(sons.at(RIGHT),sR,idx);

    //fv=(_pr.at(sonLeftID)*fvL).cwiseProduct(_pr.at(sonRightID)*fvR);

    bpp::MatrixTools::mult(_pr.at(sonLeftID),fvL,PrfvL);
    bpp::MatrixTools::mult(_pr.at(sonRightID),fvR,PrfvR);
    bpp::MatrixTools::hadamardMult(PrfvL,PrfvR,fv);

    //fv0=fv.dot(_pi);
    fv0=MatrixBppUtils::dotProd(fv,_pi);

    //pr=tree->getIota()*fv0;
    pr = _iota.at(tree->getId()) * fv0;

    double pL,pR;
    pL=compute_lk_gap_down(sons.at(LEFT),sL);
    pR=compute_lk_gap_down(sons.at(RIGHT),sR);

    pr=pr+pL+pR;

    return pr;
}
double pPIP::computeLK_M_local_tree_s_opt(double valM,
                                    double valX,
                                    double valY,
                                    double nu,
                                    bpp::Node *tree,
                                    std::string &sL,
                                    std::string &sR,
                                    unsigned long m,
                                    std::map<std::string,double> &lkM){


    double pr;
    double val;

    std::string s;
//    Eigen::VectorXd fvL;
//    Eigen::VectorXd fvR;
//    Eigen::VectorXd fv;
    bpp::ColMatrix<double> fvL;
    bpp::ColMatrix<double> fvR;
    bpp::ColMatrix<double> PrfvL;
    bpp::ColMatrix<double> PrfvR;
    bpp::ColMatrix<double> fv;
    double fv0;
    unsigned long idx;

    auto sons = tree->getSons();
    int nodeID = tree->getId();
    int sonLeftID = sons.at(LEFT)->getId();
    int sonRightID = sons.at(RIGHT)->getId();

    s.append(sL);
    s.append(sR);

    auto it=lkM.find(s);
    if(it == lkM.end()){
        idx=0;
        fvL=go_down(sons.at(LEFT),sL,idx);

        idx=0;
        fvR=go_down(sons.at(RIGHT),sR,idx);

        //fv=(tree->getNodeLeft()->getPr()*fvL).cwiseProduct(tree->getNodeRight()->getPr()*fvR);
        bpp::MatrixTools::mult(_pr.at(sonLeftID),fvL,PrfvL);
        bpp::MatrixTools::mult(_pr.at(sonRightID),fvR,PrfvR);
        bpp::MatrixTools::hadamardMult(PrfvL,PrfvR,fv);

        //fv0=fv.dot(_pi);
        fv0=MatrixBppUtils::dotProd(fv,_pi);

        //pr=tree->getIota()*tree->getBeta()*fv0;
        pr=_iota.at(nodeID)*_beta.at(nodeID)*fv0;

        pr=log((long double)pr);

        lkM[s]=pr;

    }else{
        pr=it->second;
    }

    val=-log((long double)m)+log((long double)nu)+pr+max_of_three(valM,valX,valY,DBL_EPSILON);

    return val;
}
double pPIP::computeLK_X_local_tree_s_opt(double valM,
                                    double valX,
                                    double valY,
                                    double nu,
                                    bpp::Node *tree,
                                    std::string &sL,
                                    std::string &col_gap_R,
                                    unsigned long m,
                                    std::map<std::string,double> &lkX){


    double pr;
    double val;

    std::string s;
//    Eigen::VectorXd fvL;
//    Eigen::VectorXd fvR;
//    Eigen::VectorXd fv;
    bpp::ColMatrix<double> fvL;
    bpp::ColMatrix<double> fvR;
    bpp::ColMatrix<double> PrfvL;
    bpp::ColMatrix<double> PrfvR;
    bpp::ColMatrix<double> fv;
    unsigned long idx;
    double fv0;

    auto sons = tree->getSons();
    int nodeID = tree->getId();
    int sonLeftID = sons.at(LEFT)->getId();
    int sonRightID = sons.at(RIGHT)->getId();

    s.append(sL);
    s.append(col_gap_R);

    auto it=lkX.find(s);
    if(it == lkX.end()){
        idx=0;
        fvL=go_down(sons.at(LEFT),sL,idx);

        idx=0;
        fvR=go_down(sons.at(RIGHT),col_gap_R,idx);

        //fv=(tree->getNodeLeft()->getPr()*fvL).cwiseProduct(tree->getNodeRight()->getPr()*fvR);
        bpp::MatrixTools::mult(_pr.at(sonLeftID),fvL,PrfvL);
        bpp::MatrixTools::mult(_pr.at(sonRightID),fvR,PrfvR);
        bpp::MatrixTools::hadamardMult(PrfvL,PrfvR,fv);

        //fv0=fv.dot(pi);
        fv0=MatrixBppUtils::dotProd(fv,_pi);

        //pr=tree->getIota()*tree->getBeta()*fv0;
        pr=_iota.at(nodeID)*_beta.at(nodeID)*fv0;

        double pL;

        //pL=compute_lk_down(tree->getNodeLeft(),sL);
        pL=compute_lk_down(sons.at(LEFT),sL);

        pr+=pL;

        pr=log((long double)pr);

        lkX[s]=pr;

    }else{
        pr=it->second;
    }

    val=-log((long double)m)+log((long double)nu)+pr+max_of_three(valM,valX,valY,DBL_EPSILON);

    return val;
}
double pPIP::computeLK_Y_local_tree_s_opt(double valM,
                                    double valX,
                                    double valY,
                                    double nu,
                                    bpp::Node *tree,
                                    std::string &col_gap_L,
                                    std::string &sR,
                                    unsigned long m,
                                    std::map<std::string,double> &lkY){


    double pr;
    double val;

    std::string s;
//    Eigen::VectorXd fvL;
//    Eigen::VectorXd fvR;
//    Eigen::VectorXd fv;
//    double fv0;
//    unsigned long idx;

    bpp::ColMatrix<double> fvL;
    bpp::ColMatrix<double> fvR;
    bpp::ColMatrix<double> PrfvL;
    bpp::ColMatrix<double> PrfvR;
    bpp::ColMatrix<double> fv;
    unsigned long idx;
    double fv0;

    auto sons = tree->getSons();
    int nodeID = tree->getId();
    int sonLeftID = sons.at(LEFT)->getId();
    int sonRightID = sons.at(RIGHT)->getId();

    s.append(col_gap_L);
    s.append(sR);

    auto it=lkY.find(s);
    if(it == lkY.end()){
        idx=0;
        fvL=go_down(sons.at(LEFT),col_gap_L,idx);

        idx=0;
        fvR=go_down(sons.at(RIGHT),sR,idx);

        //fv=(tree->getNodeLeft()->getPr()*fvL).cwiseProduct(tree->getNodeRight()->getPr()*fvR);

        bpp::MatrixTools::mult(_pr.at(sonLeftID),fvL,PrfvL);
        bpp::MatrixTools::mult(_pr.at(sonRightID),fvR,PrfvR);
        bpp::MatrixTools::hadamardMult(PrfvL,PrfvR,fv);

        //fv0=fv.dot(pi);
        fv0=MatrixBppUtils::dotProd(fv,_pi);

        //pr=tree->getIota()*tree->getBeta()*fv0;
        pr=_iota.at(nodeID)*_beta.at(nodeID)*fv0;

        double pR;


        pR=compute_lk_down(sons.at(RIGHT),sR);


        pr+=pR;

        pr=log((long double)pr);

        lkY[s]=pr;

    }else{
        pr=it->second;
    }

    val=-log((long double)m)+log((long double)nu)+pr+max_of_three(valM,valX,valY,DBL_EPSILON);


    return val;
}
void pPIP::DP3D_PIP2(bpp::Node *node,UtreeBppUtils::treemap *tm,double gamma_rate, bool local){

    //TODO: place as argument
    bool randomSeed = true;

    //TODO: re-implement gamma distribution
    double lambda_gamma = _lambda * gamma_rate;
    double mu_gamma = _mu * gamma_rate;

    if(local){
        bpp::Tree *subtree = _tree->cloneSubtree(node->getId());
        _tau = subtree->getTotalLength();
        _setNu();
    }else{
    }

    unsigned long up_corner_i;
    unsigned long up_corner_j;
    unsigned long bot_corner_i;
    unsigned long bot_corner_j;
    unsigned long lw;
    unsigned long h,w;


    //--------------------------------
    auto sons = node->getSons();

    int s1ID = sons.at(LEFT)->getId();
    int s2ID = sons.at(RIGHT)->getId();
    //--------------------------------

    int nodeID = node->getId();

    h = _MSA.at(s1ID).size()+1;
    w = _MSA.at(s2ID).size()+1;

    unsigned long d=(h-1)+(w-1)+1;

    double pc0;

    std::string sLs;
    std::string sRs;
    std::string col_gap_Ls;
    std::string col_gap_Rs;

    col_gap_Ls=createGapCol(_MSA.at(s1ID).size());
    col_gap_Rs=createGapCol(_MSA.at(s2ID).size());

    signed long seed;
    if(randomSeed){
        seed = std::chrono::system_clock::now().time_since_epoch().count();

    }else{
        seed = 0;
    }

    std::default_random_engine generator(seed);
    std::uniform_real_distribution<double> distribution(0.0,1.0);

    auto epsilon=DBL_EPSILON;

    //***************************************************************************************
    //***************************************************************************************
    if(local){
        pc0 = compute_pr_gap_local_tree_s(node, col_gap_Ls, col_gap_Rs);
    }else{
        /*
        pc0 = compute_pr_gap_all_edges_s(node,
                                         col_gap_Ls,
                                         col_gap_Rs,
                                         pi,
                                         originalAlphabetSize,
                                         alphabet);
        */
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

    LogM[0][0]=_nu*(pc0-1.0);
    LogX[0][0]=_nu*(pc0-1.0);
    LogY[0][0]=_nu*(pc0-1.0);

    TR[0] = new int[1]();
    TR[0][0]=STOP_STATE;

    double max_of_3;
    double max_lk=-std::numeric_limits<double>::infinity();
    double prev_max_lk=-std::numeric_limits<double>::infinity();
    signed long level_max_lk=INT_MIN;
    double val;
    unsigned long m_binary_this;
    unsigned long m_binary_prev;

    double valM;
    double valX;
    double valY;

    unsigned long idx;

    unsigned long coordSeq_1;
    unsigned long coordSeq_2;
    unsigned long coordTriangle_this_i;
    unsigned long coordTriangle_this_j;
    unsigned long coordTriangle_prev_i;
    unsigned long coordTriangle_prev_j;

    double score=-std::numeric_limits<double>::infinity();
    int start_depth;
    unsigned long depth;

    bool flag_exit=false;
    unsigned long last_d=d-1;
    unsigned long size_tr,tr_up_i,tr_up_j,tr_down_i,tr_down_j;
    std::map<std::string,double> lkM;
    std::map<std::string,double> lkX;
    std::map<std::string,double> lkY;

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

                sLs=(_MSA.at(s1ID).at(coordSeq_1));

                for(int j=0;j<=lw;j++){

                    coordTriangle_this_j=up_corner_j-j;
                    coordSeq_2=coordTriangle_this_j-1;
                    coordTriangle_prev_j=coordTriangle_this_j-1;

                    sRs=(_MSA.at(s2ID).at(coordSeq_2));

                    idx=get_indices_M(coordTriangle_prev_i,coordTriangle_prev_j,up_corner_i,up_corner_j,bot_corner_i,bot_corner_j,m-1,h,w);
                    if(idx>=0){
                        valM=LogM[m_binary_prev][idx];
                    }else
                        valM=-std::numeric_limits<double>::infinity();
                    }

                    idx=get_indices_X(coordTriangle_prev_i,coordTriangle_prev_j,up_corner_i,up_corner_j,bot_corner_i,bot_corner_j,m-1,h,w);
                    if(idx>=0){
                        valX=LogX[m_binary_prev][idx];
                    }else{
                        valX=-std::numeric_limits<double>::infinity();
                    }

                    idx=get_indices_Y(coordTriangle_prev_i,coordTriangle_prev_j,up_corner_i,up_corner_j,bot_corner_i,bot_corner_j,m-1,h,w);
                    if(idx>=0){
                        valY=LogY[m_binary_prev][idx];
                    }else{
                        valY=-std::numeric_limits<double>::infinity();
                    }

                    if(std::isinf(valM) && std::isinf(valX) && std::isinf(valY)){
                        exit(EXIT_FAILURE);
                    }

                    if(local){
                        val=computeLK_M_local_tree_s_opt(valM,
                                                         valX,
                                                         valY,
                                                         _nu,
                                                         node,
                                                         sLs,
                                                         sRs,
                                                         m,
                                                         lkM);
                    }else{
                        /*
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
                        */
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

                sLs=(_MSA.at(s1ID).at(coordSeq_1));

                for(int j=0;j<=lw;j++){

                    coordTriangle_this_j=up_corner_j-j;
                    coordTriangle_prev_j=coordTriangle_this_j;

                    idx=get_indices_M(coordTriangle_prev_i,coordTriangle_prev_j,up_corner_i,
                                      up_corner_j,bot_corner_i,bot_corner_j,m-1,h,w);
                    if(idx>=0){
                        valM=LogM[m_binary_prev][idx];
                    }else{
                        valM=-std::numeric_limits<double>::infinity();
                    }

                    idx=get_indices_X(coordTriangle_prev_i,coordTriangle_prev_j,up_corner_i,
                                      up_corner_j,bot_corner_i,bot_corner_j,m-1,h,w);
                    if(idx>=0){
                        valX=LogX[m_binary_prev][idx];
                    }else{
                        valX=-std::numeric_limits<double>::infinity();
                    }

                    idx=get_indices_Y(coordTriangle_prev_i,coordTriangle_prev_j,up_corner_i,
                                      up_corner_j,bot_corner_i,bot_corner_j,m-1,h,w);
                    if(idx>=0){
                        valY=LogY[m_binary_prev][idx];
                    }else{
                        valY=-std::numeric_limits<double>::infinity();
                    }

                    if(std::isinf(valM) && std::isinf(valX) && std::isinf(valY)){
                        exit(EXIT_FAILURE);
                    }

                    if(local){
                        val=computeLK_X_local_tree_s_opt(valM,
                                                         valX,
                                                         valY,
                                                         _nu,
                                                         node,
                                                         sLs, col_gap_Rs,
                                                         m,
                                                         lkX);
                    }else{
                        /*
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
                        */
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

                    sRs=(_MSA.at(s2ID).at(coordSeq_2));

                    idx=get_indices_M(coordTriangle_prev_i,coordTriangle_prev_j,up_corner_i,
                                      up_corner_j,bot_corner_i,bot_corner_j,m-1,h,w);
                    if(idx>=0){
                        valM=LogM[m_binary_prev][idx];
                    }else{
                        valM=-std::numeric_limits<double>::infinity();
                    }

                    idx=get_indices_X(coordTriangle_prev_i,coordTriangle_prev_j,up_corner_i,
                                      up_corner_j,bot_corner_i,bot_corner_j,m-1,h,w);
                    if(idx>=0){
                        valX=LogX[m_binary_prev][idx];
                    }else{
                        valX=-std::numeric_limits<double>::infinity();
                    }

                    idx=get_indices_Y(coordTriangle_prev_i,coordTriangle_prev_j,up_corner_i,
                                      up_corner_j,bot_corner_i,bot_corner_j,m-1,h,w);
                    if(idx>=0){
                        valY=LogY[m_binary_prev][idx];
                    }else{
                        valY=-std::numeric_limits<double>::infinity();
                    }

                    if(std::isinf(valM) && std::isinf(valX) && std::isinf(valY)){
                        exit(EXIT_FAILURE);
                    }

                    if(local){
                        val=computeLK_Y_local_tree_s_opt(valM,
                                                         valX,
                                                         valY,
                                                         _nu,
                                                         node,
                                                         col_gap_Ls, sRs,
                                                         m,
                                                         lkY);
                    }else{
                        /*
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
                         */
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
                        mval=-std::numeric_limits<double>::infinity();
                    }

                    idx=get_indices_X(coordTriangle_this_i,coordTriangle_this_j,up_corner_i,
                                      up_corner_j,bot_corner_i,bot_corner_j,m,h,w);
                    if(idx>=0){
                        xval=LogX[m_binary_this][idx];
                    }else{
                        xval=-std::numeric_limits<double>::infinity();
                    }

                    idx=get_indices_Y(coordTriangle_this_i,coordTriangle_this_j,up_corner_i,
                                      up_corner_j,bot_corner_i,bot_corner_j,m,h,w);
                    if(idx>=0){
                        yval=LogY[m_binary_this][idx];
                    }else{
                        yval=-std::numeric_limits<double>::infinity();
                    }

                    mval=fabs((long double)mval)<epsilon?-std::numeric_limits<double>::infinity():mval;
                    xval=fabs((long double)xval)<epsilon?-std::numeric_limits<double>::infinity():xval;
                    yval=fabs((long double)yval)<epsilon?-std::numeric_limits<double>::infinity():yval;

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

                        if(max_of_3>score){
                            score=max_of_3;
                            level_max_lk=m;
                        }



                    }
                    //.-----.------.------.-------.------.-------.-------.--------.-------.--------.//



                }
                lw++;
            }
        }
    }


    depth=level_max_lk;

    //result.score=score;

    _score.at(nodeID)=score;

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

    //result.traceback_path=traceback_path;
    //result.MSAs=build_MSA(traceback_path,result_L.MSAs,result_R.MSAs);
    //result_v.push_back(result);

    _traceback_path=traceback_path;

    build_MSA(node,traceback_path);

    setMSAsequenceNames(node);

    
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
    //return result;
}
void pPIP::PIPAligner2(UtreeBppUtils::treemap *tm,
                      std::vector<tshlib::VirtualNode *> &list_vnode_to_root,
                      bpp::SequenceContainer *sequences,
                      double gamma_rate,
                      bool local) {

    for (auto &vnode:list_vnode_to_root) {

        auto node = _tree->getNode(tm->right.at(vnode),false);

        if(node->isLeaf()){

            std::string seqname = sequences->getSequencesNames().at((unsigned long)vnode->vnode_seqid);

            setMSAsequenceNames(node,seqname);

            setMSAleaves(node,sequences->getSequence(seqname).toString());

            //setIndicatorFun(node);

            //setLKxyLeaves(node);

        }else{

            DP3D_PIP2(node, tm, gamma_rate, local);

        }
    }

}