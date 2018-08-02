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

#ifndef MINIJATI_PIPNODE_HPP
#define MINIJATI_PIPNODE_HPP

#include <Bpp/Numeric/VectorTools.h>
#include <Bpp/Phyl/Node.h>
#include <random>
#include <Utree.hpp>
#include <glog/logging.h>

#include "Utilities.hpp"

#include "progressivePIP.hpp"
#include "CompositePIPmsa.hpp"

namespace bpp {

    class PIPnode{ // pure virtual class

    private:

        //***************************************************************************************
        // PRIVATE FIELDS
        //***************************************************************************************

        int nodeID_;

        //***************************************************************************************
        // PRIVATE METHODS
        //***************************************************************************************

    public:

        //***************************************************************************************
        // PUBLIC FIELDS
        //***************************************************************************************

        // pointers that build a binary tree of PIPnodes
        // each PIPnode has a father and two children (left and right)
        // the parent node of the root is NULL
        // the children nodes of a leaf node are both NULL
        PIPnode *parent;
        PIPnode *childL;
        PIPnode *childR;

        const progressivePIP *progressivePIP_; // pointer to progressivePIP

        //iPIPmsa *MSA_; //contains the MSA

        tshlib::VirtualNode *vnode_; // pointer to vnode
        bpp:: Node *bnode_; // pointer to bnode

        double subTreeLenL_; // left subtree length (subtree rooted at this node)
        double subTreeLenR_; // right subtree length (subtree rooted at this node)

        std::vector<double> iotasNode_; //map of nodeIDs and vector of iotas (1 for each rate (Gamma,...) category
        std::vector<double> betasNode_; //map of nodeIDs and vector of betas (1 for each rate (Gamma,...) category
        std::vector<bpp::RowMatrix<double> > prNode_; // map of NodeIDs of Pr = exp(branchLength * rate * Q), rate under Gamma distribution

        std::vector<double> alphaNode_;

//        std::vector<double> log_lk_down_; //each node a vector of lk
//        std::vector<double> lk_empty_; //each node a vector of lk_empty (for each gamma category)

        //std::vector< std::vector< bpp::ColMatrix<double> > > fv_data_; // [site][catg][fv]
//        std::vector< bpp::ColMatrix<double> > fv_empty_data_; // [catg][fv]
//
//        std::vector< std::vector<double> > fv_sigma_; // [site][catg]
//        std::vector<double>  fv_empty_sigma_; // [catg]

        //double distanceToRoot; // length of the path from this node to root (sum of branch length)

        //***************************************************************************************
        // PUBLIC METHODS
        //***************************************************************************************

        PIPnode(const progressivePIP *pPIP,
                tshlib::VirtualNode *vnode,
                bpp::Node *bnode); // constructor

        virtual ~PIPnode(){};

        int _getId(){ return nodeID_; }; // get the Id of the node

        tshlib::VirtualNode *_getVnode(){ return vnode_; }; // get the tshlib node pointer
        bpp:: Node *_getBnode(){ return bnode_; }; // get the bpp node pointer

        bool _isRootNode(); // true if is the PIPnode root, false otherwise

        bool _isTerminalNode(); // true if is a PIPnode leaf, false otherwise

        //void _setFVemptyLeaf(); // compute the fv array at leaf (indicator array) for an empty column

        //void _setFVsigmaEmptyLeaf(); // compute fv_sigma = fv dot pi for an empty column

        //void _setFVsigmaEmptyNode();

        //void _setFVemptyNode();

        //void _setFVleaf(MSA_t &MSA); // compute the fv array at leaf (indicator array)

        //void _setFVsigmaLeaf(); // compute fv_sigma = fv dot pi

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

        void _getPrFromSubstitutionModel(); // compute exp(br_len * Q)

        void _computeLkEmptyLeaf();

        virtual void _computeAllFvEmptySigmaRec() {};

        virtual void DP3D_PIP_leaf(){}; // align a leaf PIPnode

        virtual void DP3D_PIP_node(){}; // align an internal PIPnode

        virtual void DP3D_PIP() = 0; // pure virtual method
        //***************************************************************************************

    };

}

#endif //MINIJATI_PIPNODE_HPP
