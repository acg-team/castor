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
#include <boost/bimap.hpp>
#include <Bpp/Phyl/Tree.h>
#include <Bpp/Phyl/TreeTemplate.h>
#include <Bpp/Phyl/Likelihood/AbstractHomogeneousTreeLikelihood.h>
#include <Bpp/Phyl/Distance/DistanceEstimation.h>
#include <Bpp/Seq/GeneticCode/GeneticCode.h>

//#include "pPIP.hpp"

namespace UtreeBppUtils {
    using namespace tshlib;

    typedef boost::bimap<int, tshlib::VirtualNode *> treemap;
    typedef treemap::value_type nodeassoc;

    //void convertTree_b2u(bpp::TreeTemplate<bpp::Node> *in_tree, Utree *out_tree, treemap &tm);
    void convertTree_b2u(bpp::Tree *in_tree, Utree *out_tree, treemap &tm);

    void _traverseTree_b2u(Utree *in_tree, VirtualNode *target, bpp::Tree *refTree, int nodeId, treemap &tm);

    // void _traverseTree_b2u(Utree *in_tree, VirtualNode *target, bpp::Node *source, treemap &tm);

    bpp::TreeTemplate<bpp::Node> *convertTree_u2b(Utree *in_tree);

    void _traverseTree_u2b(bpp::Node *target, VirtualNode *source);


    void updateTree_b2u(bpp::TreeTemplate<bpp::Node> inBTree, Utree *inUTree, treemap &tm);

    void updateTree_u2b(bpp::Tree *inBTree, Utree *inUTree, treemap &tm);


    void associateNode2Alignment(bpp::SiteContainer *sites, Utree *in_tree);

    void associateNode2Alignment(bpp::SequenceContainer *sequences, tshlib::Utree *in_tree);

    void renameInternalNodes(bpp::Tree *in_tree, std::string prefix = "V");

    std::vector<bpp::Node *> remapNodeLists(std::vector<tshlib::VirtualNode *> &inputList, bpp::TreeTemplate<bpp::Node> *tree, UtreeBppUtils::treemap tm);




}

namespace MatrixBppUtils {


    //Eigen::MatrixXd Matrix2Eigen(const bpp::Matrix<double> &inMatrix);
    //bpp::Matrix<double> Eigen2Matrix(Eigen::MatrixXd &inMatrix);

    //bpp::RowMatrix<double> Eigen2Matrix(const Eigen::MatrixXd &M);

    //Eigen::VectorXd Vector2Eigen(const std::vector<double> &inVector);

    double dotProd(const std::vector<double> *x, const std::vector<double> *y);

    std::vector<double> cwiseProd(std::vector<double> *x, std::vector<double> *y);

    double dotProd(const bpp::ColMatrix<double> &x, const bpp::ColMatrix<double> &y);

    double sumVector(std::vector<double> *x);

    std::vector<double> matrixVectorProd(bpp::RowMatrix<double> &M, std::vector<double> &A);

}

namespace InputUtils {

    bpp::DistanceMatrix *parseDistanceMatrix(std::string filepath);


}

namespace TextUtils {


    std::string appendToFilePath(std::string inputFilePath, std::string string2append);


}

namespace OutputUtils {

    void printParametersLikelihood(bpp::AbstractHomogeneousTreeLikelihood *tl);

    std::string tree2string(bpp::Tree *tree);

    void exportOutput2JSON(bpp::AbstractHomogeneousTreeLikelihood *tl, bpp::SiteContainer *sites, std::string prefix="");

    void exportTreeAnnotations2TSV(bpp::Tree *tree, std::string outputfile);

}

namespace DistanceUtils {

    bpp::DistanceEstimation computeDistanceMethod(std::string seqfilename, bpp::Alphabet *alphabet, bpp::GeneticCode *gCode, std::map<std::string, std::string> &params);

    bpp::TreeTemplate<bpp::Node> *computeDistanceTree(bpp::TransitionModel *model, bpp::DiscreteDistribution *rDist, bpp::DistanceEstimation &distEstimation, std::map<std::string, std::string> &
    params);

}


namespace AlignmentUtils{

    void CheckAlignmentConsistency(bpp::SiteContainer &sites);

}


namespace ComparisonUtils{

    bool areLogicallyEqual(double a, double b);
}

class DistanceMethodsUtils {

private:
    std::string seqfilename_;
    bpp::Alphabet *alphabet_;
    bpp::GeneticCode *gCode_;
    bpp::DiscreteDistribution *rdist_;
    bpp::TransitionModel *model_;
    bpp::DistanceEstimation *distEstimation_;
    bpp::DistanceMatrix *distMatrix_;
    std::map<std::string, std::string> params_;


public:

    DistanceMethodsUtils(std::string seqfilename, bpp::Alphabet *alphabet, bpp::GeneticCode *gCode, std::map<std::string, std::string> &params) :
            seqfilename_(seqfilename),
            alphabet_(alphabet),
            gCode_(gCode),
            params_(params),
            distMatrix_(nullptr),
            distEstimation_(nullptr) {}

    ~DistanceMethodsUtils() {
        delete model_;
        delete distEstimation_;
        delete rdist_;
    };

    void computeDistanceMatrix();

    void setDistanceMatrix();

    bpp::TreeTemplate<bpp::Node> *computeDistanceTree();


};



#endif //MINIJATI_UTILS_HPP
