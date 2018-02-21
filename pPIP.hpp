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
#include "Utree.hpp"
#include "Utilities.hpp"

namespace bpp {

    class pPIP {

    public:

        pPIP();

        ~pPIP(){};

        void setLKxyLeaves(bpp::Node *node);

        void setAllIotas(UtreeBppUtils::treemap *tm,std::vector<tshlib::VirtualNode *> list_vnode_to_root,double mu,double tau);
        void setAllBetas(UtreeBppUtils::treemap *tm,std::vector<tshlib::VirtualNode *> &listNodes,double mu);
        void computePr(UtreeBppUtils::treemap *tm,std::vector<tshlib::VirtualNode *> &listNodes);

        void PIPAligner(UtreeBppUtils::treemap *tm,
                              std::vector<tshlib::VirtualNode *> list_vnode_to_root,
                              bpp::SequenceContainer *sequences,
                              Vdouble &pi,
                              double lambda,
                              double mu,
                              const bpp::Alphabet *alphabet,
                              double gamma_rate,
                              bool local);


    protected:
    private:
        mutable TreeTemplate<Node> *tree_;

    public:
        void setTree(const Tree *tree) {
            tree_ = new TreeTemplate<Node>(*tree);
        }

    private:
        mutable std::map<bpp::Node *, bpp::VVdouble > _fv;
        mutable std::map<bpp::Node *, bpp::Vdouble > _lkxy;
        mutable std::map<bpp::Node *, double> _iota;
        mutable std::map<bpp::Node *, double> _beta;
        mutable std::map<bpp::Node *, bpp::RowMatrix<double> > _pr;
        mutable std::map<bpp::Node *, std::vector< std::pair<std::string,std::string> >> _MSA;
        std::string _traceback_path;
        double _score;

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
        Vdouble computeLKgapColLocal(bpp::Node *node,
                                           Vdouble &pi,
                                           double &val,
                                           double &p0);
        Vdouble computeLKmatchLocal(double valM,
                                          double valX,
                                          double valY,
                                          double nu,
                                          bpp::Node *node,
                                    unsigned long col_i, unsigned long col_j,
                                          Vdouble &pi,
                                          unsigned long m,
                                          double &val);
        Vdouble computeLKgapxLocal(double valM,
                                         double valX,
                                         double valY,
                                         double nu,
                                         bpp::Node *node,
                                         unsigned long col_i,
                                         unsigned long col_j,
                                         Vdouble &pi,
                                         unsigned long m,
                                         double &val,
                                         double &lkx);
        Vdouble computeLKgapyLocal(double valM,
                                         double valX,
                                         double valY,
                                         double nu,
                                         bpp::Node *node,
                                         unsigned long col_i,
                                         unsigned long col_j,
                                         Vdouble &pi,
                                         unsigned long m,
                                         double &val,
                                         double &lky);
        bool checkUniformLen(std::vector<std::pair<std::string,std::string>> &result);
        unsigned long getMSAlength(std::vector<std::pair<std::string,std::string>> &result);
        std::string createGapCol(unsigned long len);
        std::string createMSAcol(std::vector<std::pair<std::string,std::string>> &result, unsigned long index);
        std::vector<std::pair<std::string,std::string>> align_seq_left(	std::vector<std::pair<std::string,std::string>> &MSA_in,
                                                                                 std::string &traceback_path);
        std::vector<std::pair<std::string,std::string>> align_seq_right(std::vector<std::pair<std::string,std::string>> &result,
                                                                              std::string &traceback_path);
        std::vector<std::pair<std::string,std::string>> build_MSA(std::string traceback_path,
                                                                        std::vector<std::pair<std::string,std::string>> &MSA_L,
                                                                        std::vector<std::pair<std::string,std::string>> &MSA_R);

        void setIndicatorFun(bpp::Node *node,int extAlphabetSize);

        void DP3D_PIP(bpp::Node *node,
                            UtreeBppUtils::treemap *tm,
                            Vdouble &pi,
                            double lambda,
                            double mu,
                            const bpp::Alphabet *alphabet,
                            double gamma_rate,
                            bool local);



    };


}


#endif //MINIJATI_PPIP_HPP
