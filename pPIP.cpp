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
                double mu,
                bool local){


    _reserve(listNodes.size());

    _setTree(tree);
    _setLambda(lambda);
    _setMu(mu);
    _setPi(pi);
    _computePr(tm,listNodes);

    if(!local){
        _setLocalTau(this->getRootNode());
        _setNu();
        _setAllIotas(tm,listNodes);
        _setAllBetas(tm,listNodes);
    }

}
void pPIP::_reserve(unsigned long numNodes){

    _score.resize(numNodes);
    _score.assign(numNodes,-std::numeric_limits<double>::infinity());
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

    if(fabs(a)<epsilon){
        a=-std::numeric_limits<double>::infinity();
    }
    if(fabs(b)<epsilon){
        b=-std::numeric_limits<double>::infinity();
    }
    if(fabs(c)<epsilon){
        c=-std::numeric_limits<double>::infinity();
    }

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
std::string pPIP::createGapCol(unsigned long len){

    std::string colMSA (len,'-');

    return colMSA;
}
void pPIP::build_MSA(bpp::Node *node, std::string traceback_path){


    auto sons = node->getSons();
    auto s1 = sons.at(0);
    auto s2 = sons.at(1);
    unsigned long s1ID = (unsigned long)s1->getId();
    unsigned long s2ID = (unsigned long)s2->getId();

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

    _MSA.at(node->getId())=MSA;
}
void pPIP::setMSAsequenceNames(bpp::Node *node){

    auto sons = node->getSons();
    int s1ID = sons.at(LEFT)->getId();
    int s2ID = sons.at(RIGHT)->getId();

    std::vector<std::string> seqNames;

    for(int i=0;i<_seqNames.at(s1ID).size();i++){
        seqNames.push_back(_seqNames.at(s1ID).at(i));
    }

    for(int i=0;i<_seqNames.at(s2ID).size();i++){
        seqNames.push_back(_seqNames.at(s2ID).at(i));
    }

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
}
void pPIP::_setLocalTau(bpp::Node *node) {
    bpp::Tree *subtree = _tree->cloneSubtree(node->getId());
    _tau = subtree->getTotalLength();
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
            _iota.at(node->getId())=(1/_mu)/T;
        } else {
            _iota.at(node->getId())=node->getDistanceToFather()/T;
        }
    }
}
void pPIP::_setAllIotas(bpp::Node *node,bool local_root) {
    double T;

    T = _tau + 1/_mu;

    if (fabs(T) < SMALL_DOUBLE) {
        perror("ERROR in set_iota: T too small");
    }

    if(local_root){
        _iota.at(node->getId())=(1/_mu)/T;
    }else{
        _iota.at(node->getId())=node->getDistanceToFather()/T;
    }

    if(!node->isLeaf()){
        auto sons = node->getSons();

        _setAllIotas(sons.at(LEFT),false);
        _setAllIotas(sons.at(RIGHT),false);

    }

}
void pPIP::_setAllBetas(bpp::Node *node,bool local_root) {

    if(local_root){

        _beta.at(node->getId())=1.0;

    }else{

        double muT = _mu * node->getDistanceToFather();

        if (fabs(muT) < SMALL_DOUBLE) {
            perror("ERROR mu * T is too small");
        }

        _beta.at(node->getId())= (1.0 - exp(-_mu * node->getDistanceToFather())) / muT;

    }

    if(!node->isLeaf()){
        auto sons = node->getSons();

        _setAllBetas(sons.at(LEFT),false);
        _setAllBetas(sons.at(RIGHT),false);

    }

}
void pPIP::_setAllBetas(UtreeBppUtils::treemap *tm,std::vector<tshlib::VirtualNode *> &listNodes) {

    for (auto &vnode:listNodes) {

        auto node = _tree->getNode(tm->right.at(vnode),false);

        if (!node->hasFather()) {
            _beta.at(node->getId())=1.0;
        } else {

            double muT = _mu * node->getDistanceToFather();
            if (fabs(muT) < SMALL_DOUBLE) {
                perror("ERROR mu * T is too small");
            }

            _beta.at(node->getId())= (1.0 - exp(-_mu * node->getDistanceToFather())) / muT;
        }

    }

}
void pPIP::_computePr(UtreeBppUtils::treemap *tm,std::vector<tshlib::VirtualNode *> &listNodes) {

    for (auto &vnode:listNodes) {

        auto node = _tree->getNode(tm->right.at(vnode),false);

        if (!node->hasFather()) {

        } else {
            _pr.at(node->getId()) = MatrixBppUtils::Eigen2Matrix(vnode->getPr());
        }

    }

}

bpp::ColMatrix<double> pPIP::fv_observed(std::string &s, unsigned long &idx){
    bpp::ColMatrix<double> fv;
    int ii;
    char ch=s[idx];

    fv.resize(_extendedAlphabetSize,1);
    bpp::MatrixTools::fill(fv,0.0);

    ii=_alphabet->charToInt(&ch);
    ii=ii<0?_alphabetSize:ii;

    fv(ii,0)=1.0;
    idx++;

    return fv;
}
bpp::ColMatrix<double> pPIP::go_down(bpp::Node *node,std::string &s, unsigned long &idx){
    bpp::ColMatrix<double> fv;
    bpp::ColMatrix<double> fvL;
    bpp::ColMatrix<double> fvR;
    bpp::ColMatrix<double> PrfvL;
    bpp::ColMatrix<double> PrfvR;


    if(node->isLeaf()){

        fv=fv_observed(s,idx);

    }else{

        auto sons = node->getSons();
        int sonLeftID = sons.at(LEFT)->getId();
        int sonRightID = sons.at(RIGHT)->getId();

        fvL=go_down(sons.at(LEFT),s,idx);
        fvR=go_down(sons.at(RIGHT),s,idx);

        bpp::MatrixTools::mult(_pr.at(sonLeftID),fvL,PrfvL);
        bpp::MatrixTools::mult(_pr.at(sonRightID),fvR,PrfvR);
        bpp::MatrixTools::hadamardMult(PrfvL,PrfvR,fv);

    }

    return fv;
}
void pPIP::allgaps(bpp::Node *node,std::string &s, unsigned long &idx,bool &flag){

    if(node->isLeaf()){
        char ch=s[idx];

        idx++;

        if(ch!='-'){
            flag=false;
        }

    }else{

        auto sons = node->getSons();
        allgaps(sons.at(LEFT),s,idx,flag);
        allgaps(sons.at(RIGHT),s,idx,flag);
    }

}
double pPIP::compute_lk_gap_down(bpp::Node *node,std::string &s){

    double pr=0;
    double pL=0;
    double pR=0;
    unsigned long idx;
    bpp::ColMatrix<double> fv;
    bpp::ColMatrix<double> fvL;
    bpp::ColMatrix<double> fvR;
    double fv0;

    auto sons = node->getSons();

    int nodeID = node->getId();

    if(node->isLeaf()){
        idx=0;
        fv=go_down(node,s,idx);

        fv0=MatrixBppUtils::dotProd(fv,_pi);

        pr=_iota.at(nodeID)-_iota.at(nodeID)*_beta.at(nodeID)+_iota.at(nodeID)*_beta.at(nodeID)*fv0;
        return pr;
    }else{
        idx=0;
        fv=go_down(node,s,idx);

        fv0=MatrixBppUtils::dotProd(fv,_pi);

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
double pPIP::compute_lk_down(bpp::Node *node,std::string &s){

    double pr;
    unsigned long idx;
    bpp::ColMatrix<double> fv;
    bpp::ColMatrix<double> fvL;
    bpp::ColMatrix<double> fvR;
    double fv0;

    auto sons = node->getSons();

    int nodeID = node->getId();

    if(node->isLeaf()){

        idx=0;
        fv=go_down(node,s,idx);
        fv0=MatrixBppUtils::dotProd(fv,_pi);
        pr=_iota.at(nodeID)*_beta.at(nodeID)*fv0;

        return pr;

    }else{

        idx=0;
        fv=go_down(node,s,idx);
        fv0=MatrixBppUtils::dotProd(fv,_pi);
        pr=_iota.at(nodeID)*_beta.at(nodeID)*fv0;

        bool flagL=true;
        bool flagR=true;
        idx=0;
        allgaps(sons.at(LEFT),s,idx,flagL);
        unsigned long ixx=idx;
        allgaps(sons.at(RIGHT),s,idx,flagR);

        unsigned long len;
        if(flagR){
            std::string sL;
            len=ixx;
            sL=s.substr(0,len);
            return pr + compute_lk_down(sons.at(LEFT),sL);
        }

        if(flagL){
            std::string sR;
            sR=s.substr(ixx);
            return pr + compute_lk_down(sons.at(RIGHT),sR);
        }

    }

    return pr;
}
double pPIP::computeLK_GapColumn_local(bpp::Node *node, std::string &sL, std::string &sR){
    double fv0;
    double pr;
    bpp::ColMatrix<double> fvL;
    bpp::ColMatrix<double> fvR;
    bpp::ColMatrix<double> PrfvL;
    bpp::ColMatrix<double> PrfvR;
    bpp::ColMatrix<double> fv;
    unsigned long idx;

    auto sons = node->getSons();

    int sonLeftID = sons.at(LEFT)->getId();
    int sonRightID = sons.at(RIGHT)->getId();

    idx=0;
    fvL=go_down(sons.at(LEFT),sL,idx);

    idx=0;
    fvR=go_down(sons.at(RIGHT),sR,idx);

    bpp::MatrixTools::mult(_pr.at(sonLeftID),fvL,PrfvL);
    bpp::MatrixTools::mult(_pr.at(sonRightID),fvR,PrfvR);
    bpp::MatrixTools::hadamardMult(PrfvL,PrfvR,fv);

    fv0=MatrixBppUtils::dotProd(fv,_pi);

    pr = _iota.at(node->getId()) * fv0;

    double pL,pR;
    pL=compute_lk_gap_down(sons.at(LEFT),sL);
    pR=compute_lk_gap_down(sons.at(RIGHT),sR);

    pr=pr+pL+pR;

    return pr;
}
double pPIP::computeLK_M_local(double valM,
                               double valX,
                               double valY,
                               bpp::Node *node,
                               std::string &sL,
                               std::string &sR,
                               unsigned long m,
                               std::map<std::string, double> &lkM){


    double pr;
    double val;

    std::string s;
    bpp::ColMatrix<double> fvL;
    bpp::ColMatrix<double> fvR;
    bpp::ColMatrix<double> PrfvL;
    bpp::ColMatrix<double> PrfvR;
    bpp::ColMatrix<double> fv;
    double fv0;
    unsigned long idx;

    auto sons = node->getSons();
    int nodeID = node->getId();
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

        bpp::MatrixTools::mult(_pr.at(sonLeftID),fvL,PrfvL);
        bpp::MatrixTools::mult(_pr.at(sonRightID),fvR,PrfvR);
        bpp::MatrixTools::hadamardMult(PrfvL,PrfvR,fv);

        fv0=MatrixBppUtils::dotProd(fv,_pi);

        pr=_iota.at(nodeID)*_beta.at(nodeID)*fv0;

        pr=log((long double)pr);

        lkM[s]=pr;

    }else{
        pr=it->second;
    }

    val=-log((long double)m)+log((long double)_nu)+pr+max_of_three(valM,valX,valY,DBL_EPSILON);

    return val;
}
double pPIP::computeLK_X_local(double valM,
                               double valX,
                               double valY,
                               bpp::Node *node,
                               std::string &sL,
                               std::string &col_gap_R,
                               unsigned long m,
                               std::map<std::string, double> &lkX){


    double pr;
    double val;

    std::string s;
    bpp::ColMatrix<double> fvL;
    bpp::ColMatrix<double> fvR;
    bpp::ColMatrix<double> PrfvL;
    bpp::ColMatrix<double> PrfvR;
    bpp::ColMatrix<double> fv;
    unsigned long idx;
    double fv0;

    auto sons = node->getSons();
    int nodeID = node->getId();
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

        bpp::MatrixTools::mult(_pr.at(sonLeftID),fvL,PrfvL);
        bpp::MatrixTools::mult(_pr.at(sonRightID),fvR,PrfvR);
        bpp::MatrixTools::hadamardMult(PrfvL,PrfvR,fv);

        fv0=MatrixBppUtils::dotProd(fv,_pi);

        pr=_iota.at(nodeID)*_beta.at(nodeID)*fv0;

        double pL;

        pL=compute_lk_down(sons.at(LEFT),sL);

        pr+=pL;

        pr=log((long double)pr);

        lkX[s]=pr;

    }else{
        pr=it->second;
    }

    val=-log((long double)m)+log((long double)_nu)+pr+max_of_three(valM,valX,valY,DBL_EPSILON);

    return val;
}
double pPIP::computeLK_Y_local(double valM,
                               double valX,
                               double valY,
                               bpp::Node *node,
                               std::string &col_gap_L,
                               std::string &sR,
                               unsigned long m,
                               std::map<std::string, double> &lkY){


    double pr;
    double val;

    std::string s;

    bpp::ColMatrix<double> fvL;
    bpp::ColMatrix<double> fvR;
    bpp::ColMatrix<double> PrfvL;
    bpp::ColMatrix<double> PrfvR;
    bpp::ColMatrix<double> fv;
    unsigned long idx;
    double fv0;

    auto sons = node->getSons();
    int nodeID = node->getId();
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

        bpp::MatrixTools::mult(_pr.at(sonLeftID),fvL,PrfvL);
        bpp::MatrixTools::mult(_pr.at(sonRightID),fvR,PrfvR);
        bpp::MatrixTools::hadamardMult(PrfvL,PrfvR,fv);

        fv0=MatrixBppUtils::dotProd(fv,_pi);

        pr=_iota.at(nodeID)*_beta.at(nodeID)*fv0;

        double pR;


        pR=compute_lk_down(sons.at(RIGHT),sR);

        pr+=pR;

        pr=log((long double)pr);

        lkY[s]=pr;

    }else{
        pr=it->second;
    }

    val=-log((long double)m)+log((long double)_nu)+pr+max_of_three(valM,valX,valY,DBL_EPSILON);


    return val;
}
void pPIP::DP3D_PIP(bpp::Node *node,UtreeBppUtils::treemap *tm,double gamma_rate, bool local){

    //TODO: place as argument
    bool randomSeed = true;

    //TODO: re-implement gamma distribution
    double lambda_gamma = _lambda * gamma_rate;
    double mu_gamma = _mu * gamma_rate;

    if(local){
        _setLocalTau(node);
        _setNu();
        _setAllIotas(node,true);
        _setAllBetas(node,true);
    }else{
    }

    unsigned long up_corner_i;
    unsigned long up_corner_j;
    unsigned long bot_corner_i;
    unsigned long bot_corner_j;
    unsigned long lw;
    unsigned long h,w;

    auto sons = node->getSons();

    int s1ID = sons.at(LEFT)->getId();
    int s2ID = sons.at(RIGHT)->getId();

    int nodeID = node->getId();

    h = _MSA.at(s1ID).size()+1;
    w = _MSA.at(s2ID).size()+1;

    unsigned long d=(h-1)+(w-1)+1;

    double pc0;

    std::string sLs;
    std::string sRs;
    std::string col_gap_Ls;
    std::string col_gap_Rs;

    unsigned long numLeavesLeft = _seqNames.at(s1ID).size();
    unsigned long numLeavesRight = _seqNames.at(s2ID).size();

    col_gap_Ls=createGapCol(numLeavesLeft);
    col_gap_Rs=createGapCol(numLeavesRight);

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
        pc0 = computeLK_GapColumn_local(node, col_gap_Ls, col_gap_Rs);
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

    double max_of_3=-std::numeric_limits<double>::infinity();

    signed long level_max_lk=INT_MIN;
    double val;
    unsigned long m_binary_this;
    unsigned long m_binary_prev;

    double valM;
    double valX;
    double valY;

    signed long idx;

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

        if(checkboundary(up_corner_i,up_corner_j,bot_corner_i,bot_corner_j,h,w)) {

            lw = 0;
            for (unsigned long i = up_corner_i; i <= bot_corner_i; i++) {

                coordTriangle_this_i = i;
                coordSeq_1 = coordTriangle_this_i - 1;
                coordTriangle_prev_i = coordTriangle_this_i - 1;

                sLs = (_MSA.at(s1ID).at(coordSeq_1));

                for (int j = 0; j <= lw; j++) {

                    coordTriangle_this_j = up_corner_j - j;
                    coordSeq_2 = coordTriangle_this_j - 1;
                    coordTriangle_prev_j = coordTriangle_this_j - 1;

                    sRs = (_MSA.at(s2ID).at(coordSeq_2));

                    idx = get_indices_M(coordTriangle_prev_i, coordTriangle_prev_j, up_corner_i, up_corner_j,
                                        bot_corner_i, bot_corner_j, m - 1, h, w);
                    if (idx >= 0) {
                        valM = LogM[m_binary_prev][idx];
                    } else {
                        valM = -std::numeric_limits<double>::infinity();
                    }

                    idx = get_indices_X(coordTriangle_prev_i, coordTriangle_prev_j, up_corner_i, up_corner_j,
                                        bot_corner_i, bot_corner_j, m - 1, h, w);
                    if (idx >= 0) {
                        valX = LogX[m_binary_prev][idx];
                    } else {
                        valX = -std::numeric_limits<double>::infinity();
                    }

                    idx = get_indices_Y(coordTriangle_prev_i, coordTriangle_prev_j, up_corner_i, up_corner_j,
                                        bot_corner_i, bot_corner_j, m - 1, h, w);
                    if (idx >= 0) {
                        valY = LogY[m_binary_prev][idx];
                    } else {
                        valY = -std::numeric_limits<double>::infinity();
                    }

                    if (std::isinf(valM) && std::isinf(valX) && std::isinf(valY)) {
                        exit(EXIT_FAILURE);
                    }

                    if (local) {
                        val = computeLK_M_local(valM,
                                                valX,
                                                valY,
                                                node,
                                                sLs,
                                                sRs,
                                                m,
                                                lkM);
                    } else {
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

                    if (std::isinf(val)) {
                        exit(EXIT_FAILURE);
                    }

                    if (std::isnan(val)) {
                        exit(EXIT_FAILURE);
                    }

                    idx = get_indices_M(coordTriangle_this_i, coordTriangle_this_j, up_corner_i,
                                        up_corner_j, bot_corner_i, bot_corner_j, m, h, w);

                    LogM[m_binary_this][idx] = val;
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
                        val= computeLK_X_local(valM,
                                               valX,
                                               valY,
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
                        val= computeLK_Y_local(valM,
                                               valX,
                                               valY,
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

        size_tr=(unsigned long)ceil((tr_down_i-tr_up_i+1)*(tr_up_j-tr_down_j+1+1)/2);

        /*TODO: optimize size TR*/
        TR[m] = new int[size_tr]();
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

                    int ttrr;

                    ttrr=index_of_max(mval,xval,yval,epsilon,generator,distribution);

                    idx=get_indices_T(coordTriangle_this_i,coordTriangle_this_j,up_corner_i,
                                      up_corner_j,bot_corner_i,bot_corner_j,m,h,w);

                    if(TR[m][idx]!=0){
                        exit(EXIT_FAILURE);
                    }

                    TR[m][idx]=ttrr;

                    if( (coordTriangle_this_i==(h-1)) & (coordTriangle_this_j==(w-1)) ){

                        max_of_3=max_of_three(mval,xval,yval,epsilon);

                        if(max_of_3>score){
                            score=max_of_3;
                            level_max_lk=m;
                        }

                    }

                }
                lw++;
            }
        }
    }

    depth=level_max_lk;

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
}
std::vector< std::string > pPIP::getMSA(bpp::Node *node){

    return _MSA.at(node->getId());

}
double pPIP::getScore(bpp::Node *node){
    return _score.at(node->getId());
}
std::vector< std::string > pPIP::getSeqnames(bpp::Node *node) {
    return _seqNames.at(node->getId());
}
bpp::Node * pPIP::getRootNode(){
    return _tree->getRootNode();
}
bpp::Alphabet *pPIP::getAlphabet(){
    return _alphabet;
}
void pPIP::PIPAligner(UtreeBppUtils::treemap *tm,
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

        }else{

            DP3D_PIP(node, tm, gamma_rate, local);

        }
    }

}

bpp::SiteContainer * pPIPUtils::pPIPmsa2Sites(bpp::pPIP *progressivePIP){
    auto MSA = progressivePIP->getMSA(progressivePIP->getRootNode());

    auto sequences = new bpp::VectorSequenceContainer(progressivePIP->getAlphabet());

    auto seqNames = progressivePIP->getSeqnames(progressivePIP->getRootNode());

    int msaLen = MSA.size();

    int numLeaves = seqNames.size();
    for(int j=0;j<numLeaves;j++){
        std::string seqname = seqNames.at(j);
        std::string seqdata;
        seqdata.resize(msaLen);
        for (int i = 0; i < msaLen; i++) {
            seqdata.at(i)=MSA.at(i).at(j);
        }
        sequences->addSequence(*(new bpp::BasicSequence(seqname, seqdata, progressivePIP->getAlphabet())), true);
    }

    return new bpp::VectorSiteContainer(*sequences);
}