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

        typedef std::string MSAcolumn_t;
        typedef std::vector<MSAcolumn_t> MSA_t;
        typedef std::vector<MSA_t> MSAensemble_t; //for SB
        typedef std::string TracebackPath_t;
        typedef std::vector<TracebackPath_t> TracebackEnsemble_t; //for SB

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
        std::vector< MSA_t > _MSA; //MSA at each node
        std::vector< MSAensemble_t > _MSAensemble; //MSAensemble at each node
        std::vector<TracebackPath_t> _traceback_path;
        std::vector<TracebackEnsemble_t> _traceback_path_ensemble;
        std::vector< double > _score;
        std::vector<vector< double >> _score_ensemble;
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

        bpp::ColMatrix<double> fv_observed(MSAcolumn_t &s, unsigned long &idx);

        bpp::ColMatrix<double> go_down(bpp::Node *node,MSAcolumn_t &s, unsigned long &idx);

        void allgaps(bpp::Node *node,MSAcolumn_t &s, unsigned long &idx,bool &flag);

        double compute_lk_gap_down(bpp::Node *node,MSAcolumn_t &s);

        double computeLK_GapColumn_local(bpp::Node *node, MSAcolumn_t &sL, MSAcolumn_t &sR);

        double compute_lk_down(bpp::Node *node,MSAcolumn_t &s);

        double computeLK_M_local(double valM,
                                 double valX,
                                 double valY,
                                 bpp::Node *node,
                                 std::string &sL,
                                 std::string &sR,
                                 unsigned long m,
                                 std::map<MSAcolumn_t, double> &lkM);

        double computeLK_X_local(double valM,
                                 double valX,
                                 double valY,
                                 bpp::Node *node,
                                 MSAcolumn_t &sL,
                                 MSAcolumn_t &col_gap_R,
                                 unsigned long m,
                                 std::map<MSAcolumn_t, double> &lkX);

        double computeLK_Y_local(double valM,
                                 double valX,
                                 double valY,
                                 bpp::Node *node,
                                 MSAcolumn_t &col_gap_L,
                                 MSAcolumn_t &sR,
                                 unsigned long m,
                                 std::map<MSAcolumn_t, double> &lkY);

        void DP3D_PIP(bpp::Node *node, UtreeBppUtils::treemap *tm, double gamma_rate, bool local);

        void DP3D_PIP_SB(bpp::Node *node,UtreeBppUtils::treemap *tm,double gamma_rate, bool local,
                         double temperature,int num_SB);

    };


}

namespace pPIPUtils {

    bpp::SiteContainer *pPIPmsa2Sites(bpp::pPIP *progressivePIP);

    double add_lns(double a_ln,double b_ln);

    void max_val_in_column(double ***M,int depth, int height, int width, double &val, int &level);

}

#endif //MINIJATI_PPIP_HPP
