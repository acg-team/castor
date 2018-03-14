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

        pPIP(bpp::Alphabet *alphabet);

        ~pPIP(){};

        void init(const Tree *tree,
             UtreeBppUtils::treemap *tm,
             std::vector<tshlib::VirtualNode *> &listNodes,
             const Vdouble &pi,
             double lambda,
             double mu, bool local);

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
        mutable TreeTemplate<Node> *tree_;
        std::vector<double> iota_;
        std::vector<double> beta_;
        std::vector<bpp::RowMatrix<double> > pr_;
        std::vector<std::vector<std::string> > sequenceNames_;
        std::vector<std::vector<std::string> > MSA_;
        std::string tracebackPath_;
        std::vector<double> score_;
        double lambda_;
        double mu_;
        double nu_;
        double tau_;

        bpp::ColMatrix<double> pi_;

        bpp::Alphabet *Alphabet_;

        long alphabetSize_;

        long extendedAlphabetSize_;

        void _reserve(unsigned long numNodes);

        void _setNu();

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

        bool isInside(unsigned long x0,
                      unsigned long y0,
                      unsigned long xf,
                      unsigned long yf,
                      unsigned long xt,
                      unsigned long yt);

        void setIndicesM(unsigned long &up_corner_i,
                         unsigned long &up_corner_j,
                         unsigned long &bot_corner_i,
                         unsigned long &bot_corner_j,
                         unsigned long level,
                         unsigned long h,
                         unsigned long w);

        void setIndicesX(unsigned long &up_corner_i,
                         unsigned long &up_corner_j,
                         unsigned long &bot_corner_i,
                         unsigned long &bot_corner_j,
                         unsigned long level,
                         unsigned long h,
                         unsigned long w);

        void setIndicesY(unsigned long &up_corner_i,
                         unsigned long &up_corner_j,
                         unsigned long &bot_corner_i,
                         unsigned long &bot_corner_j,
                         unsigned long level,
                         unsigned long h,
                         unsigned long w);

        signed long getIndicesM(unsigned long nx,
                                unsigned long ny,
                                unsigned long up_corner_i,
                                unsigned long up_corner_j,
                                unsigned long bot_corner_i,
                                unsigned long bot_corner_j,
                                unsigned long m,
                                unsigned long h,
                                unsigned long w);

        signed long getIndicesX(unsigned long nx,
                                unsigned long ny,
                                unsigned long up_corner_i,
                                unsigned long up_corner_j,
                                unsigned long bot_corner_i,
                                unsigned long bot_corner_j,
                                unsigned long m,
                                unsigned long h,
                                unsigned long w);

        signed long getIndicesY(unsigned long nx,
                                unsigned long ny,
                                unsigned long up_corner_i,
                                unsigned long up_corner_j,
                                unsigned long bot_corner_i,
                                unsigned long bot_corner_j,
                                unsigned long m,
                                unsigned long h,
                                unsigned long w);

        void setIndicesT(unsigned long &up_corner_i,
                         unsigned long &up_corner_j,
                         unsigned long &bot_corner_i,
                         unsigned long &bot_corner_j,
                         unsigned long level,
                         unsigned long h,
                         unsigned long w);

        void resetCorner(unsigned long &up_corner_i,
                         unsigned long &up_corner_j,
                         unsigned long &bot_corner_i,
                         unsigned long &bot_corner_j,
                         unsigned long h,
                         unsigned long w);

        unsigned long getIndicesT(unsigned long nx,
                                  unsigned long ny,
                                  unsigned long up_corner_i,
                                  unsigned long up_corner_j,
                                  unsigned long bot_corner_i,
                                  unsigned long bot_corner_j,
                                  unsigned long m,
                                  unsigned long h,
                                  unsigned long w);

        int getIndexOfMax(double m,
                          double x,
                          double y,
                          double epsilon,
                          std::default_random_engine &generator,
                          std::uniform_real_distribution<double> &distribution);

        double getMaxOfThree(double a,
                             double b,
                             double c,
                             double epsilon);

        bool checkBoundary(unsigned long up_corner_i,
                           unsigned long up_corner_j,
                           unsigned long bot_corner_i,
                           unsigned long bot_corner_j,
                           unsigned long h,
                           unsigned long w);

        std::string createGapColumn(unsigned long len);

        void buildMSA(bpp::Node *node, std::string traceback_path);

        void setMSASequenceNames(bpp::Node *node);

        void setMSASequenceNames(bpp::Node *node, std::string seqname);

        void setMSALeaves(bpp::Node *node, const std::string &MSA);

        bpp::ColMatrix<double> fv_observed(std::string &s, unsigned long &idx);

        bpp::ColMatrix<double> go_down(bpp::Node *node,std::string &s, unsigned long &idx);

        void allgaps(bpp::Node *node,std::string &s, unsigned long &idx,bool &flag);

        double computeLikelihoodGapDown(bpp::Node *node, std::string &s);

        double computeLikelihoodGapColumnLocal(bpp::Node *node, std::string &sL, std::string &sR);

        double computeLikelihoodDown(bpp::Node *node, std::string &s);

        double computeLikelihood_M_local(double valM,
                                         double valX,
                                         double valY,
                                         bpp::Node *node,
                                         std::string &sL,
                                         std::string &sR,
                                         unsigned long m,
                                         std::map<std::string, double> &lkM);

        double computeLikelihood_X_local(double valM,
                                         double valX,
                                         double valY,
                                         bpp::Node *node,
                                         std::string &sL,
                                         std::string &col_gap_R,
                                         unsigned long m,
                                         std::map<std::string, double> &lkX);

        double computeLikelihood_Y_local(double valM,
                                         double valX,
                                         double valY,
                                         bpp::Node *node,
                                         std::string &col_gap_L,
                                         std::string &sR,
                                         unsigned long m,
                                         std::map<std::string, double> &lkY);

        void DP3D_PIP(bpp::Node *node, UtreeBppUtils::treemap *tm, double gamma_rate, bool local);


    };


}

namespace pPIPUtils {
    bpp::SiteContainer *pPIPmsa2Sites(bpp::pPIP *progressivePIP);
}

#endif //MINIJATI_PPIP_HPP
