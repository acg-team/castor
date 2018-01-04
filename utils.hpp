/*******************************************************************************
 * Licensed Materials - Property of Lorenzo Gatti & Massimo Maiolo
 *
 *
 * Copyright (C) 2015-2017 by Lorenzo Gatti & Massimo Maiolo
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
 * @file utils.hpp
 * @author Lorenzo Gatti
 * @author Massimo Maiolo
 * @date 21 12 2017
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
#ifndef MINIJATI_UTILS_HPP
#define MINIJATI_UTILS_HPP
#include <Utree.hpp>

#include <Bpp/Phyl/Tree.h>
#include <Bpp/Phyl/TreeTemplate.h>

namespace UtreeBppUtils{
using namespace tshlib;
    void convertUtree(bpp::TreeTemplate<bpp::Node> *in_tree, Utree *out_tree);
    void _traverseTree(Utree *in_tree, VirtualNode *target, bpp::Node *source);


}

namespace MatrixBppUtils{


    Eigen::MatrixXd Matrix2Eigen(const bpp::Matrix<double> &inMatrix);
    Eigen::VectorXd Vector2Eigen(const std::vector<double> &inVector);

}


#endif //MINIJATI_UTILS_HPP
