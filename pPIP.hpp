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
 * @file pPIP.hpp
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
#ifndef MINIJATI_PPIP_HPP
#define MINIJATI_PPIP_HPP

#include <Bpp/Numeric/VectorTools.h>
#include <Bpp/Phyl/Node.h>
#include <random>
#include "Utree.hpp"
#include "Utilities.hpp"

#define SMALL_DOUBLE 1e-8
#define LEFT 0
#define RIGHT 1

namespace bpp {

    class pPIP {

    public:

        typedef std::string MSAcolumn;
        typedef std::vector<MSAcolumn> MSA;
        typedef std::vector<MSA> MSAensemble;

        pPIP(bpp::Alphabet *alphabet);

        ~pPIP(){};

        void init(const Tree *tree, bpp::SubstitutionModel *smodel,
             UtreeBppUtils::treemap *tm,
             std::vector<tshlib::VirtualNode *> &listNodes, bool local);

        void PIPAligner(UtreeBppUtils::treemap *tm,
                        std::vector<tshlib::VirtualNode *> &list_vnode_to_root,
                        bpp::SequenceContainer *sequences,
                        double gamma_rate,
                        bool local);


        std::vector< std::string > getMSA(bpp::Node *node);
        double getScore(bpp::Node *node);
        std::vector< std::string > getSeqnames(bpp::Node *node);
        bpp::Node *getRootNode();
        bpp::Alphabet *getAlphabet();

    protected:

    private:
        mutable TreeTemplate<Node> *_tree;
        mutable bpp::SubstitutionModel *_substModel;
        std::vector< double > _iota;
        std::vector< double > _beta;
        std::vector< bpp::RowMatrix<double> > _pr;
        std::vector< std::vector< std::string > > _seqNames;
        std::vector< std::vector< std::string > > _MSA;
        std::string _traceback_path;
        std::vector< double > _score;
        double _lambda;
        double _mu;
        double _nu;
        double _tau;

        bpp::ColMatrix<double> _pi;

        bpp::Alphabet *_alphabet;

        long _alphabetSize;

        long _extendedAlphabetSize;

        void _reserve(unsigned long numNodes);

        void _setNu();

        void _setSubstModel(bpp::SubstitutionModel *smodel);

        void _setTree(const Tree *tree);

        void _setLambda(double lambda);

        void _setMu(double mu);

        void _setPi(const Vdouble &pi);

        void _setLocalTau(bpp::Node *node);

        void _setAllIotas(UtreeBppUtils::treemap *tm,std::vector<tshlib::VirtualNode *> &listNodes);

        void _setAllIotas(bpp::Node *node,bool local_root);

        void _setAllBetas(UtreeBppUtils::treemap *tm,std::vector<tshlib::VirtualNode *> &listNodes);

        void _setAllBetas(bpp::Node *node,bool local_root);

        void _computePr(UtreeBppUtils::treemap *tm,std::vector<tshlib::VirtualNode *> &listNodes);

        bool is_inside(unsigned long x0,unsigned long y0,unsigned long xf,unsigned long yf,unsigned long xt,unsigned long yt);

        void set_indeces_M(unsigned long &up_corner_i,unsigned long &up_corner_j,unsigned long &bot_corner_i,
                                 unsigned long &bot_corner_j,unsigned long level,unsigned long h,unsigned long w);

        void set_indeces_X(unsigned long &up_corner_i,unsigned long &up_corner_j,unsigned long &bot_corner_i,
                                 unsigned long &bot_corner_j,unsigned long level,unsigned long h,unsigned long w);

        void set_indeces_Y(unsigned long &up_corner_i,unsigned long &up_corner_j,unsigned long &bot_corner_i,
                                 unsigned long &bot_corner_j,unsigned long level,unsigned long h,unsigned long w);

        signed long get_indices_M(unsigned long nx,unsigned long ny,unsigned long up_corner_i,unsigned long up_corner_j,
                                        unsigned long bot_corner_i,unsigned long bot_corner_j,unsigned long m,unsigned long h,unsigned long w);

        signed long get_indices_X(unsigned long nx,unsigned long ny,unsigned long up_corner_i,unsigned long up_corner_j,
                                        unsigned long bot_corner_i,unsigned long bot_corner_j,unsigned long m,unsigned long h,unsigned long w);

        signed long get_indices_Y(unsigned long nx,unsigned long ny,unsigned long up_corner_i,unsigned long up_corner_j,
                                        unsigned long bot_corner_i,unsigned long bot_corner_j,unsigned long m,unsigned long h,unsigned long w);

        void set_indeces_T(unsigned long &up_corner_i,unsigned long &up_corner_j,unsigned long &bot_corner_i,
                                 unsigned long &bot_corner_j,unsigned long level,unsigned long h,unsigned long w);

        void reset_corner(unsigned long &up_corner_i,unsigned long &up_corner_j,unsigned long &bot_corner_i,
                                unsigned long &bot_corner_j,unsigned long h,unsigned long w);

        unsigned long get_indices_T(unsigned long nx,unsigned long ny,unsigned long up_corner_i,unsigned long up_corner_j,
                                          unsigned long bot_corner_i,unsigned long bot_corner_j,unsigned long m,unsigned long h,unsigned long w);

        int index_of_max(double m, double x, double y,double epsilon,
                               std::default_random_engine &generator,
                               std::uniform_real_distribution<double> &distribution);

        double max_of_three(double a, double b, double c,double epsilon);

        bool checkboundary(unsigned long up_corner_i,unsigned long up_corner_j,unsigned long bot_corner_i,
                                 unsigned long bot_corner_j,unsigned long h,unsigned long w);

        std::string createGapCol(unsigned long len);

        void build_MSA(bpp::Node *node, std::string traceback_path);

        void setMSAsequenceNames(bpp::Node *node);

        void setMSAsequenceNames(bpp::Node *node,std::string seqname);

        void setMSAleaves(bpp::Node *node,const std::string &sequence);

        bpp::ColMatrix<double> fv_observed(MSAcolumn &s, unsigned long &idx);

        bpp::ColMatrix<double> go_down(bpp::Node *node,MSAcolumn &s, unsigned long &idx);

        void allgaps(bpp::Node *node,MSAcolumn &s, unsigned long &idx,bool &flag);

        double compute_lk_gap_down(bpp::Node *node,MSAcolumn &s);

        double computeLK_GapColumn_local(bpp::Node *node, MSAcolumn &sL, MSAcolumn &sR);

        double compute_lk_down(bpp::Node *node,MSAcolumn &s);

        double computeLK_M_local(double valM,
                                 double valX,
                                 double valY,
                                 bpp::Node *node,
                                 std::string &sL,
                                 std::string &sR,
                                 unsigned long m,
                                 std::map<MSAcolumn, double> &lkM);

        double computeLK_X_local(double valM,
                                 double valX,
                                 double valY,
                                 bpp::Node *node,
                                 MSAcolumn &sL,
                                 MSAcolumn &col_gap_R,
                                 unsigned long m,
                                 std::map<MSAcolumn, double> &lkX);

        double computeLK_Y_local(double valM,
                                 double valX,
                                 double valY,
                                 bpp::Node *node,
                                 MSAcolumn &col_gap_L,
                                 MSAcolumn &sR,
                                 unsigned long m,
                                 std::map<MSAcolumn, double> &lkY);

        void DP3D_PIP(bpp::Node *node, UtreeBppUtils::treemap *tm, double gamma_rate, bool local);


        double add_lns_2(double a_ln,double b_ln);

        void DP3D_PIP_SB(bpp::Node *node,UtreeBppUtils::treemap *tm,double gamma_rate, bool local,double temperature,int num_SB);

    };


}

namespace pPIPUtils {
    bpp::SiteContainer *pPIPmsa2Sites(bpp::pPIP *progressivePIP);
}

#endif //MINIJATI_PPIP_HPP
