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

    class PIPnode{

    private:

        //***************************************************************************************
        // PRIVATE FIELDS
        //***************************************************************************************

        int nodeID_;

    public:

        //***************************************************************************************
        // PUBLIC FIELDS
        //***************************************************************************************

        PIPnode *parent;
        PIPnode *childL;
        PIPnode *childR;

        const progressivePIP *progressivePIP_; // pointer to progressivePIP

        iPIPmsa *MSA_;

        tshlib::VirtualNode *vnode_; // pointer to vnode
        bpp:: Node *bnode_; // pointer to bnode

        double subTreeLenL_; // left subtree length (subtree rooted at this node)
        double subTreeLenR_; // right subtree length (subtree rooted at this node)

        std::vector<double> iotasNode_; //map of nodeIDs and vector of iotas (1 for each rate (Gamma,...) category
        std::vector<double> betasNode_; //map of nodeIDs and vector of betas (1 for each rate (Gamma,...) category
        std::vector<bpp::RowMatrix<double> > prNode_; // map of NodeIDs of Pr = exp(branchLength * rate * Q), rate under Gamma distribution

        vector<double> log_lk_down_; //each nodeInterface a vector of lk
        vector<double> log_lk_empty_down_; //each nodeInterface a vector of lk_empty (for each gamma category)

        std::vector< std::vector< bpp::ColMatrix<double> > > fv_data_; // [site][catg][fv]
        std::vector< bpp::ColMatrix<double> > fv_empty_data_; // [catg][fv]

        std::vector< std::vector<double> > fv_sigma_; // [site][catg]
        std::vector<double>  fv_empty_sigma_; // [catg]

        //***************************************************************************************
        // PUBLIC METHODS
        //***************************************************************************************
        PIPnode(const progressivePIP *pPIP,
                tshlib::VirtualNode *vnode,
                bpp::Node *bnode);

        int _getId(){ return nodeID_; };

        tshlib::VirtualNode *_getVnode(){ return vnode_; };
        bpp:: Node *_getBnode(){ return bnode_; };

        void _setFVemptyLeaf();
        void _setFVsigmaLeaf();

        void _setFVsigmaEmptyLeaf();
        void _setFVleaf(MSA_t *MSA);

        void _reserve(int numCatg);

        double _computeSubtreeTau();

        void _getPrFromSubstitutionModel();

        //-----------------------------------------------------------
        // TODO: remove from here??? used only by nodeCPU
        bpp::ColMatrix<double> computeFVrec(MSAcolumn_t &s,
                                            int &idx,
                                            int catg);
        void allgaps(std::string &s,
                     int &idx,
                     bool &flag);

        double compute_lk_gap_down(MSAcolumn_t &s,
                                   int catg);

        double _compute_lk_down_rec(int idx,
                                    double lk);

        double _compute_lk_down(MSAcolumn_t &s,
                                int catg);

        std::vector<double> _compute_lk_down();
        //-----------------------------------------------------------

        void DP3D_PIP_leaf();

        void DP3D_PIP_node();

        void DP3D_PIP(){};
        //***************************************************************************************

    };

}

#endif //MINIJATI_PIPNODE_HPP
