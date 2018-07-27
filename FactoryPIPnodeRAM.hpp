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

        void DP3D_PIP_leaf();

        void DP3D_PIP_node();

        bool _index_of_max(double m,
                                    double x,
                                    double y,
                                    double epsilon,
                                    std::default_random_engine &generator,
                                    std::uniform_real_distribution<double> &distribution,
                                    int &index,
                                    double &val);

        double max_of_three(double a, double b, double c, double epsilon);

        void _compute_lk_empty_leaf_();

        void _compute_lk_leaf_();

        void _compress_lk_components(std::vector<double> &lk_down_not_compressed,
                                     std::vector<vector<bpp::ColMatrix<double> > > &fv_data_not_compressed);

        std::vector<double> _computeLK_empty(std::vector<bpp::ColMatrix<double> > &fvL,
                                             std::vector<bpp::ColMatrix<double> > &fvR,
                                             std::vector<bpp::ColMatrix<double> > &Fv_gap,
                                             std::vector<double> &fv_empty_sigma_,
                                             std::vector<double> &lk_empty_down_L,
                                             std::vector<double> &lk_empty_down_R);

        double _computeLK_M(std::vector<bpp::ColMatrix<double> > &fvL,
                            std::vector<bpp::ColMatrix<double> > &fvR,
                            std::vector<bpp::ColMatrix<double> > &Fv_M_ij,
                            std::vector<double> &Fv_sigma_M_ij);

        double _computeLK_X(std::vector<bpp::ColMatrix<double> > &fvL,
                            std::vector<bpp::ColMatrix<double> > &fvR,
                            std::vector<bpp::ColMatrix<double> > &Fv_X_ij,
                            std::vector<double> &Fv_sigma_X_ij);

        double _computeLK_Y(std::vector<bpp::ColMatrix<double> > &fvL,
                            std::vector<bpp::ColMatrix<double> > &fvR,
                            std::vector<bpp::ColMatrix<double> > &Fv_Y_ij,
                            std::vector<double> &Fv_sigma_Y_ij);

        double _computeLK_MXY(double log_phi_gamma,
                              double valM,
                              double valX,
                              double valY,
                              double log_pr);

        void _compress_Fv(std::vector<std::vector<double>> &fv_sigma_not_compressed,
                          std::vector<vector<bpp::ColMatrix<double> > > &fv_data_not_compressed);

        void _compute_lk_empty_down_rec(std::vector<double> &lk);

    public:

        //***************************************************************************************
        // PUBLIC METHODS
        //***************************************************************************************

        nodeRAM(const progressivePIP *pPIP, tshlib::VirtualNode *vnode, bpp::Node *bnode) : PIPnode(pPIP, vnode,
                                                                                                    bnode) {
        }

        void DP3D_PIP();

    };

}

#endif //MINIJATI_FACTORYNODERAM_HPP
