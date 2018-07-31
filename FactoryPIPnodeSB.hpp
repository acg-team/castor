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

#ifndef MINIJATI_FACTORYNODESB_HPP
#define MINIJATI_FACTORYNODESB_HPP

namespace bpp {

    class nodeSB : public PIPnode {

    private:

        //***************************************************************************************
        // PUBLIC FIELDS
        //***************************************************************************************

        PIPmsaComp *MSA_;

        //***************************************************************************************
        // PRIVATE METHODS
        //***************************************************************************************

        //void max_val_in_column(double ***M,int depth, int height, int width, double &val, int &level);

        void DP3D_PIP_leaf();

        void DP3D_PIP_node();

    public:

        //***************************************************************************************
        // PUBLIC METHODS
        //***************************************************************************************

        nodeSB(const progressivePIP *pPIP, tshlib::VirtualNode *vnode, bpp::Node *bnode) : PIPnode(pPIP, vnode, bnode) {
        }

        void DP3D_PIP();
    };

}

#endif //MINIJATI_FACTORYNODESB_HPP
