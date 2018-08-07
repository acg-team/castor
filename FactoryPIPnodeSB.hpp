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

#ifndef MINIJATI_FACTORYNODESB_HPP
#define MINIJATI_FACTORYNODESB_HPP

namespace bpp {

    class nodeSB : public PIPnode {

    private:

        //***************************************************************************************
        // PRIVATE FIELDS
        //***************************************************************************************

        std::vector<int> subMSAidxL_;
        std::vector<int> subMSAidxR_;

        //***************************************************************************************
        // PRIVATE METHODS
        //***************************************************************************************

        void DP3D_PIP_leaf(); // DP method to align a sequence at a leaf PIPnode
                              // (which reduces to data preparation)

        void DP3D_PIP_node(int msa_idx_L,
                           int msa_idx_R,
                           int &position); // DP method to align 2 MSAs at an internal node

        //std::vector<double> _computeLkEmptyNode();

        void _computeLkEmptyLeaf();

        void _computeLkLeaf();

        void startingLevelSB(std::vector< vector< vector<double> > > &Log3DM,
                                     std::vector< vector< vector<double> > > &Log3DX,
                                     std::vector< vector< vector<double> > > &Log3DY,
                                     double epsilon,
                                     std::default_random_engine &generator,
                                     std::uniform_real_distribution<double> &distribution,
                                     int d,
                                     int h,
                                     int w,
                                     int &lev,
                                     double &val,
                                     int &state);

    public:

        //***************************************************************************************
        // PUBLIC FIELDS
        //***************************************************************************************

//        nodeSB *parent;
//        nodeSB *childL;
//        nodeSB *childR;

        //PIPmsaComp *MSA_; //contains the MSA

        //***************************************************************************************
        // PUBLIC METHODS
        //***************************************************************************************

        // constructor
        nodeSB(const progressivePIP *pPIP, tshlib::VirtualNode *vnode, bpp::Node *bnode) : PIPnode(pPIP,
                                                                                                   vnode,
                                                                                                   bnode) {
            if(bnode->isLeaf()){

                // create a PIPmsaComp object
                MSA_  = new PIPmsaComp(1);

                // create a new PIPmsa
                dynamic_cast<PIPmsaComp *>(MSA_)->pipmsa.at(0) = new PIPmsa();

            }else{

                /*
                // create a PIPmsaComp object
                MSA_  = new PIPmsaComp(progressivePIP_->num_sb_);

                for(int i=0;i<progressivePIP_->num_sb_;i++){
                    // create a new PIPmsa
                    dynamic_cast<PIPmsaComp *>(MSA_)->pipmsa.at(i) = new PIPmsa();
                }

                */

            }

        }

        virtual ~nodeSB() = default;

        void DP3D_PIP();

        //void _computeAllFvEmptySigmaRec();

    };

}

#endif //MINIJATI_FACTORYNODESB_HPP
