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
 * @version 1.0.7
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

#include <Bpp/Numeric/VectorTools.h>
#include <Bpp/Phyl/Node.h>
#include <random>
#include <Utree.hpp>
#include <glog/logging.h>

#include "Utilities.hpp"

#include "progressivePIP.hpp"
#include "FactoryPIPnode.hpp"
#include "CompositePIPmsa.hpp"

#ifndef MINIJATI_FACTORYNODERAM_HPP
#define MINIJATI_FACTORYNODERAM_HPP

namespace bpp {

    class nodeRAM : public PIPnode {

    private:

        //***************************************************************************************
        // PRIVATE METHODS
        //***************************************************************************************

        void DP3D_PIP_leaf(); // DP method to align a sequence at a leaf PIPnode
                              // (which reduces to data preparation)

        void DP3D_PIP_node(); // DP method to align 2 MSAs at an internal node

        // get the index and the max value among the three input values (m,x,y)
        bool _index_of_max(double m,            // match value
                           double x,            // gapx value
                           double y,            // gapy value
                           double epsilon,      // small number for the comparison between to numbers
                           std::default_random_engine &generator,   // random number generator (when two or three numbers have the same value)
                           std::uniform_real_distribution<double> &distribution,    // uniform distribution
                           int &index,          // index of max (1: MATCH, 2: GAPX, 3: GAPY)
                           double &val);        // max value between the three (m,x,y)

        // get max value among the three input values (m,x,y)
        double max_of_three(double m,           // match value
                            double x,           // gapx value
                            double y,           // gapy value
                            double epsilon);    // small number for the comparison between to numbers

        void _computeLkLeaf(); // compute the lk at the leaf

        void _computeLkEmptyLeaf(); // compute the lk of an empty column at the leaf

        std::vector<double> _computeLkEmptyNode(); // compute the lk of an empty column at an internal node

        void _compressLK(std::vector<double> &lk_down_not_compressed);  // compress an array of lk values

        // compress an array of fv values and fv_sigma values
        void _compress_Fv(std::vector<std::vector<double>> &fv_sigma_not_compressed,
                          std::vector<vector<bpp::ColMatrix<double> > > &fv_data_not_compressed);

        // compute the MATCH lk
        void _computeLK_M(std::vector<bpp::ColMatrix<double> > &fvL, // fv array of the left child
                                   std::vector<bpp::ColMatrix<double> > &fvR, // fv array of the right child
                                   std::vector<bpp::ColMatrix<double> > &Fv_M_ij, // result of Fv_M_ij = Pr_R * fv_R hadamard_prod Pr_L * fv_L
                                   std::vector<double> &Fv_sigma_M_ij, // result of Fv_M_ij dot pi
                                   double &pr_match, // match probability (stored for the next layer)
                                   double &pr_match_full_path); // full match probability (used at this layer)
                                                                // it encompasses the probability of an insertion along
                                                                // the whole path between the root and this node

        // compute the GAPX lk
        void _computeLK_X(std::vector<bpp::ColMatrix<double> > &fvL, // fv array of the left child
                                   std::vector<bpp::ColMatrix<double> > &fvR, // fv array of the right child
                                   std::vector<bpp::ColMatrix<double> > &Fv_X_ij, // result of Fv_X_ij = Pr_R * fv_R hadamard_prod Pr_L * fv_L
                                   std::vector<double> &Fv_sigma_X_ij, // result of Fv_sigma_X_ij dot pi
                                   double &pr_gapx, // gapx probability (stored for the next layer)
                                   double &pr_gapx_full_path); // full gapx probability (used at this layer)
                                                               // it encompasses the probability of an insertion along
                                                               // the whole path between the root and this node

        // compute the GAPY lk
        void _computeLK_Y(std::vector<bpp::ColMatrix<double> > &fvL, // fv array of the left child
                                   std::vector<bpp::ColMatrix<double> > &fvR, // fv array of the right child
                                   std::vector<bpp::ColMatrix<double> > &Fv_Y_ij, // result of Fv_Y_ij = Pr_R * fv_R hadamard_prod Pr_L * fv_L
                                   std::vector<double> &Fv_sigma_Y_ij, // result of Fv_sigma_Y_ij dot pi
                                   double &pr_gapy, // gapy probability (stored for the next layer)
                                   double &pr_gapy_full_path); // full gapy probability (used at this layer)
                                                               // it encompasses the probability of an insertion along
                                                               // the whole path between the root and this node

        // compute the lk at a given matrix entry extending the previous best lk (valM,valX,valY) together with the
        // actual lk value (log_pr) and the marginal lk of an empty column
        double _computeLK_MXY(double log_phi_gamma,
                              double valM,
                              double valX,
                              double valY,
                              double log_pr);

    public:

        //***************************************************************************************
        // PUBLIC METHODS
        //***************************************************************************************

        // constructor
        nodeRAM(const progressivePIP *pPIP, tshlib::VirtualNode *vnode, bpp::Node *bnode) : PIPnode(pPIP, vnode,
                                                                                                    bnode) {
        }

        void DP3D_PIP(); // DP algorithm to align (leaf/internal node) under the PIP model

    };

}

#endif //MINIJATI_FACTORYNODERAM_HPP
