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
 * @file utils.cpp
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
#include <glog/logging.h>
#include <Bpp/Phyl/Io/Newick.h>
#include <Bpp/Phyl/App/PhylogeneticsApplicationTools.h>
#include <Bpp/Seq/Io/Fasta.h>
#include <Bpp/Phyl/Model/RateDistribution/ConstantRateDistribution.h>
#include <Bpp/Phyl/Distance/PGMA.h>
#include <Bpp/Seq/Container/SiteContainerTools.h>
#include <Bpp/Phyl/Distance/NeighborJoining.h>
#include <Bpp/Phyl/Distance/BioNJ.h>
#include <Bpp/Phyl/OptimizationTools.h>

//#include <elf.h>

#include "Utilities.hpp"
#include "TSHHomogeneousTreeLikelihood.hpp"
#include "Optimizators.hpp"

using namespace bpp;
#include "pPIP.hpp"

/*
void UtreeBppUtils::_traverseTree_b2u(Utree *in_tree, VirtualNode *target, bpp::Node *source, treemap &tm) {

    std::string name;


    if(source->isLeaf()){
        target->vnode_leaf=true;
    }else{
        target->vnode_leaf=false;
    }


    for(auto &bppNode:source->getSons()){

        auto ichild = new VirtualNode();

        // Filling the bidirectional map
        tm.insert(nodeassoc(bppNode, ichild)); // TODO: avoid the map at all costs

        ichild->vnode_id = bppNode->getId();

        if(bppNode->hasName()){
            name = bppNode->getName();
        }else{
            name = "V"+std::to_string(bppNode->getId());
        }

        ichild->vnode_name = name;
        ichild->vnode_branchlength = bppNode->getDistanceToFather();

        if(!bppNode->isLeaf()){

            ichild->vnode_leaf = false;

            target->connectNode(ichild);
            in_tree->addMember(ichild);
            _traverseTree_b2u(in_tree, ichild, bppNode, tm);

        }else{

            // Set all the other directions to null
            ichild->_setNodeLeft(nullptr);
            ichild->_setNodeRight(nullptr);

            // Set the LEAF flag to true
            ichild->vnode_leaf = true;
            in_tree->addMember(ichild);
            target->connectNode(ichild);

        }

    }


}
 */

void UtreeBppUtils::_traverseTree_b2u(Utree *in_tree, VirtualNode *target, bpp::Tree *refTree, int nodeId, treemap &tm) {


    if (refTree->isLeaf(nodeId)) {
        target->vnode_leaf = true;
    } else {
        target->vnode_leaf = false;
    }

    for (auto &sonId:refTree->getSonsId(nodeId)) {
        auto ichild = new VirtualNode();
        // Filling the bidirectional map
        tm.insert(nodeassoc(sonId, ichild));
        ichild->vnode_id = sonId;
        ichild->vnode_name = refTree->getNodeName(sonId);
        ichild->vnode_branchlength = refTree->getDistanceToFather(sonId);

        if (!refTree->isLeaf(sonId)) {

            ichild->vnode_leaf = false;

            target->connectNode(ichild);
            in_tree->addMember(ichild);
            _traverseTree_b2u(in_tree, ichild, refTree, sonId, tm);


        } else {

            // Set all the other directions to null
            ichild->_setNodeLeft(nullptr);
            ichild->_setNodeRight(nullptr);

            // Set the LEAF flag to true
            ichild->vnode_leaf = true;
            in_tree->addMember(ichild);
            target->connectNode(ichild);
        }
    }
}

void UtreeBppUtils::convertTree_b2u(bpp::Tree *in_tree, Utree *out_tree, treemap &tm) {
    int rootId = in_tree->getRootId();

    for (auto &sonId:in_tree->getSonsId(rootId)) {

        auto ichild = new VirtualNode;
        // Filling the bidirectional map
        tm.insert(nodeassoc(sonId, ichild));
        // Filling the bidirectional map

        ichild->vnode_id = sonId;
        ichild->vnode_name = in_tree->getNodeName(sonId);
        ichild->vnode_branchlength = in_tree->getDistanceToFather(sonId);

        _traverseTree_b2u(out_tree, ichild, in_tree, sonId, tm);

        // Add this node as starting point of the tree
        out_tree->addMember(ichild, true);

    }

    // Collapse multiforcating trees to star tree pointing to the same pseudoroot
    // Pick root node at random within the node-vector

    out_tree->startVNodes.at(0)->_setNodeUp(out_tree->startVNodes.at(1));
    out_tree->startVNodes.at(1)->_setNodeUp(out_tree->startVNodes.at(0));

}


/*
void UtreeBppUtils::convertTree_b2u(bpp::TreeTemplate<bpp::Node> *in_tree, Utree *out_tree, treemap &tm) {

    bpp::Node * RootNode = in_tree->getRootNode();

    std::string name;
    // For each node descending the root, create either a new VirtualNode
    for (auto &bppNode:RootNode->getSons()) {

        auto ichild = new VirtualNode;
        // Filling the bidirectional map
        tm.insert(nodeassoc(bppNode, ichild)); // TODO: avoid the map at all costs

        ichild->vnode_id = bppNode->getId();
        if(bppNode->hasName()){
            name = bppNode->getName();
        }else{
            name = "V"+std::to_string(bppNode->getId());
        }
        ichild->vnode_name = name;
        ichild->vnode_branchlength = bppNode->getDistanceToFather();

        _traverseTree_b2u(out_tree, ichild, bppNode, tm);

        // Add this node as starting point of the tree
        out_tree->addMember(ichild, true);


    }

    // Collapse multiforcating trees to star tree pointing to the same pseudoroot
    // Pick root node at random within the node-vector

    out_tree->startVNodes.at(0)->_setNodeUp(out_tree->startVNodes.at(1));
    out_tree->startVNodes.at(1)->_setNodeUp(out_tree->startVNodes.at(0));

}

 */

bpp::TreeTemplate<bpp::Node> *UtreeBppUtils::convertTree_u2b(tshlib::Utree *in_tree) {

    auto *RootNode = new bpp::Node;

    RootNode->setName(in_tree->rootnode->getNodeName());
    RootNode->setDistanceToFather(in_tree->rootnode->vnode_branchlength);

    _traverseTree_u2b(RootNode, in_tree->rootnode->getNodeLeft());
    _traverseTree_u2b(RootNode, in_tree->rootnode->getNodeRight());


    RootNode->setId((int) in_tree->listVNodes.size());

    // Set root node on new bpp tree
    auto *tree = new bpp::TreeTemplate<bpp::Node>();

    tree->setRootNode(RootNode);

    return tree;
}

void UtreeBppUtils::_traverseTree_u2b(bpp::Node *target, tshlib::VirtualNode *source) {

    auto *child = new bpp::Node;

    child->setName(source->getNodeName());
    child->setId(source->vnode_id);
    child->setDistanceToFather(source->vnode_branchlength);

    if (!source->isTerminalNode()) {

        _traverseTree_u2b(child, source->getNodeLeft());

        _traverseTree_u2b(child, source->getNodeRight());

    }

    target->addSon(child);


}

void UtreeBppUtils::associateNode2Alignment(bpp::SiteContainer *sites, tshlib::Utree *in_tree) {

    for (auto &node:in_tree->listVNodes) {

        if (node->isTerminalNode()) {

            std::vector<std::string> seqnames = sites->getSequencesNames();

            for (int i = 0; i < seqnames.size(); i++) {

                if (seqnames.at(i).compare(node->vnode_name) == 0) {
                    node->vnode_seqid = i;
                    break;
                }

            }


        }

    }
}

void UtreeBppUtils::associateNode2Alignment(bpp::SequenceContainer *sequences, tshlib::Utree *in_tree) {

    for (auto &node:in_tree->listVNodes) {

        if (node->isTerminalNode()) {

            std::vector<std::string> seqnames = sequences->getSequencesNames();

            for (int i = 0; i < seqnames.size(); i++) {

                if (seqnames.at(i).compare(node->vnode_name) == 0) {
                    node->vnode_seqid = i;
                    break;
                }

            }


        }

    }
}

void UtreeBppUtils::renameInternalNodes(bpp::Tree *in_tree, std::string prefix) {

    // Rename internal nodes with standard Vxx * where xx is a progressive number
    for (auto &nodeId:in_tree->getNodesId()) {

        if (!in_tree->hasNodeName(nodeId)) {

            std::string stringId;
            std::string stringName;

            stringId = std::to_string(nodeId);
            stringName = prefix + stringId;

            in_tree->setNodeName(nodeId, stringName);

        }
    }

}

std::vector<bpp::Node *> UtreeBppUtils::remapNodeLists(std::vector<tshlib::VirtualNode *> &inputList, bpp::TreeTemplate<bpp::Node> *tree, UtreeBppUtils::treemap tm) {

    std::vector<bpp::Node *> newList;

    for (auto &vnode:inputList) {

        newList.push_back(tree->getNode(tm.right.at(vnode)));
    }

    return newList;
}

void UtreeBppUtils::updateTree_b2u(bpp::TreeTemplate<bpp::Node> inBTree, tshlib::Utree *inUTree, UtreeBppUtils::treemap &tm) {

    //inUTree->addVirtualRootNode();
    std::vector<tshlib::VirtualNode *> nodelist;

    nodelist = inUTree->listVNodes;
    //nodelist.push_back(inUTree->rootnode);
    //bpp::Node *rNode = inBTree.getRootNode();
    //rNode->removeSons();

    std::map<int, bpp::Node *> tempMap;

    // reset inBtree
    for (auto &bnode:inBTree.getNodes()) {

        tempMap.insert(std::pair<int, bpp::Node *>(bnode->getId(), bnode));

        bnode->removeSons();
        bnode->removeFather();

    }


    for (auto &vnode:nodelist) {

        std::cerr << "vnode " << vnode->getNodeName();
        if (!vnode->isTerminalNode()) {

            // get corrisponding sons in inBTree
            //bpp::Node *leftBNode  = inBTree.getNode(tm.right.at(vnode->getNodeLeft()));
            //bpp::Node *rightBNode = inBTree.getNode(tm.right.at(vnode->getNodeRight()));
            bpp::Node *leftBNode = tempMap[tm.right.at(vnode->getNodeLeft())];
            bpp::Node *rightBNode = tempMap[tm.right.at(vnode->getNodeRight())];
            // get corrisponding parent in inBTree
            //bpp::Node *pNode = inBTree.getNode(tm.right.at(vnode));
            bpp::Node *pNode = tempMap[tm.right.at(vnode)];

            // Empty array of sons on the parent node
            //pNode->removeSons();

            //leftBNode->removeFather();
            //rightBNode->removeFather();
            leftBNode->setFather(pNode);
            rightBNode->setFather(pNode);
            //Add new sons
            pNode->setSon(0, leftBNode);
            pNode->setSon(1, rightBNode);

            std::cerr << "\t internal";

        } else {

            std::cerr << "\t leaf";

            //bpp::Node *pNode = inBTree.getNode(tm.right.at(vnode->getNodeUp()));
            //bpp::Node *cNode = inBTree.getNode(tm.right.at(vnode));
            //cNode->removeFather();
            //cNode->setFather(pNode);

        }
        // in case the current vnode is also the pseudo-root
        if (vnode == vnode->getNodeUp()->getNodeUp()) {
            std::cerr << "\tvnode pseudoroot";
            //bpp::Node *leftBNode  = inBTree.getNode(tm.right.at(vnode->getNodeLeft()));
            //bpp::Node *rightBNode = inBTree.getNode(tm.right.at(vnode->getNodeRight()));



            bpp::Node *leftBNode = tempMap[tm.right.at(vnode)];
            bpp::Node *rightBNode = tempMap[tm.right.at(vnode->getNodeUp())];
            // get corrisponding parent in inBTree
            //bpp::Node *pNode = inBTree.getNode(tm.right.at(vnode));
            //bpp::Node *pNode = tempMap[tm.right.at(vnode)].second;


            inBTree.getRootNode()->removeSons();

            //leftBNode->removeFather();
            //rightBNode->removeFather();

            leftBNode->setFather(inBTree.getRootNode());
            rightBNode->setFather(inBTree.getRootNode());

            inBTree.getRootNode()->setSon(0, leftBNode);
            inBTree.getRootNode()->setSon(1, rightBNode);

        }


        std::cerr << "\t done\n";

    }
    //inUTree->removeVirtualRootNode();

}

void UtreeBppUtils::updateTree_u2b(bpp::Tree *inBTree, tshlib::Utree *inUTree, UtreeBppUtils::treemap &tm) {

}

Eigen::MatrixXd MatrixBppUtils::Matrix2Eigen(const bpp::Matrix<double> &inMatrix) {

    size_t rows, cols;

    rows = inMatrix.getNumberOfRows();
    cols = inMatrix.getNumberOfColumns();

    Eigen::MatrixXd m = Eigen::MatrixXd::Zero(rows, cols);

    for (int r = 0; r < rows; r++) {
        for (int c = 0; c < cols; c++) {

            m(r, c) = inMatrix(r, c);

        }

    }

    return m;
}

bpp::RowMatrix<double> MatrixBppUtils::Eigen2Matrix(const Eigen::MatrixXd &M) {

    size_t rows, cols;

    rows = M.rows();
    cols = M.cols();

    bpp::RowMatrix<double> m;
    m.resize(rows, cols);

    for (int r = 0; r < rows; r++) {
        for (int c = 0; c < cols; c++) {

            m(r, c) = M(r, c);

        }

    }

    return m;
}

Eigen::VectorXd MatrixBppUtils::Vector2Eigen(const std::vector<double> &inVector) {

    Eigen::VectorXd vector = Eigen::VectorXd::Zero(inVector.size());

    for (int r = 0; r < inVector.size(); r++) {

        vector(r) = inVector.at(r);

    }

    return vector;
}

/*
bpp::Matrix<double> MatrixBppUtils::Eigen2Matrix(Eigen::MatrixXd &inMatrix) {

    bpp::Matrix<double> outMatrix;
    outMatrix.resize((size_t) inMatrix.cols(), (size_t) inMatrix.cols());

    for(int r=0; r<inMatrix.rows(); r++){
        for (int c=0; c<inMatrix.cols(); c++){

            outMatrix(r,c) = inMatrix(r,c);

        }

    }


    return outMatrix;
}
*/

double MatrixBppUtils::dotProd(const std::vector<double> *x, const std::vector<double> *y) {

    double val;

//    if(x->size() != y->size()){
//        perror("ERROR: MatrixBppUtils::dotProd");
//    }

    val = 0.0;
    for (unsigned long i = 0; i < x->size(); i++) {
        val += (x->at(i) * y->at(i));
    }

    return val;
}

double MatrixBppUtils::dotProd(const bpp::ColMatrix<double> &x, const bpp::ColMatrix<double> &y) {

    double val;

//    if(x->size() != y->size()){
//        perror("ERROR: MatrixBppUtils::dotProd");
//    }

    val = 0.0;
    for (unsigned long i = 0; i < x.getNumberOfRows(); i++) {
        val += (x(i, 0) * y(i, 0));
    }

    return val;
}

std::vector<double> MatrixBppUtils::cwiseProd(std::vector<double> *x, std::vector<double> *y) {

    std::vector<double> val;


//    if(x->size() != y->size()){
//        perror("ERROR: MatrixBppUtils::dotProd");
//    }

    val.resize(x->size());

    for (unsigned long i = 0; i < x->size(); i++) {
        val.at(i) = (x->at(i) * y->at(i));
    }

    return val;
}

double MatrixBppUtils::sumVector(std::vector<double> *x) {

    double val;

    val = 0.0;
    for (unsigned long i = 0; i < x->size(); i++) {
        val += x->at(i);
    }

    return val;
}

std::vector<double> MatrixBppUtils::matrixVectorProd(bpp::RowMatrix<double> &M, std::vector<double> &A) {

    std::vector<double> B;
    B.resize(A.size());

    for (int i = 0; i < A.size(); i++) {
        B[i] = 0;
        for (int j = 0; j < A.size(); j++) {
            B[i] += M(i, j) * A.at(j);

        }
    }

    return B;
}

bpp::DistanceMatrix *InputUtils::parseDistanceMatrix(std::string filepath) {

    std::ifstream inFile;
    inFile.open(filepath);

    int matsize;

    inFile >> matsize;
    std::string character;

    auto outmatrix = new bpp::DistanceMatrix(matsize);
    outmatrix->resize(matsize);
    std::string stringline;

    int x = 0;
    int rownum = 0;
    while (!inFile.eof()) {

        std::getline(inFile, stringline);
        std::stringstream ss(stringline);

        std::string token;
        int y = 0;

        int colnum = 0;
        while (std::getline(ss, token, ' ')) {

            if (!token.empty()) {

                if (y == 0) {
                    std::string seqname = token;
                    (*outmatrix).setName(rownum, seqname);

                } else {
                    double value = std::stod(token);
                    (*outmatrix)(rownum, colnum) = value;
                    (*outmatrix)(colnum, rownum) = value;
                    colnum++;
                }

            }
            //(*outmatrix)(x,y) = std::stod(token);
            y++;


        }
        if (!stringline.empty()) {
            x++;
            rownum++;
        }
    }

    return outmatrix;

}

void OutputUtils::printParametersLikelihood(bpp::AbstractHomogeneousTreeLikelihood *tl) {
    bpp::ParameterList parModel;
    std::ostringstream oss;

    bpp::AbstractHomogeneousTreeLikelihood *ttl;
    if (dynamic_cast<bpp::TSHHomogeneousTreeLikelihood *>( tl )) {
        bpp::TSHHomogeneousTreeLikelihood *flk = dynamic_cast<bpp::TSHHomogeneousTreeLikelihood *>(tl);
        ttl = flk->getLikelihoodFunction();
    } else {
        ttl = tl;
    }


    parModel = ttl->getSubstitutionModelParameters();
    if (parModel.size() > 0) {
        oss << "model=" << ttl->getModel()->getName() << "(";
        for (auto &parameterName:parModel.getParameterNames()) {
            oss << parameterName << "=" << parModel.getParameter(parameterName).getValue() << ",";
        }
        oss << ")";
        LOG(INFO) << oss.str();
    }
    oss.clear();
    oss.str("");


    parModel = ttl->getRateDistributionParameters();
    if (parModel.size() > 0) {
        oss << "rates=" << ttl->getRateDistribution()->getName() << "(";
        for (auto &parameterName:parModel.getParameterNames()) {
            oss << parameterName << "=" << parModel.getParameter(parameterName).getValue() << ",";
        }
        oss << ")";
        LOG(INFO) << oss.str();
    }
    oss.clear();
    oss.str("");

    parModel = ttl->getBranchLengthsParameters();
    if (parModel.size() > 0) {
        oss << "branches=" << ttl->getBranchLengthsParameters().size() << "(";
        for (auto &parameterName:parModel.getParameterNames()) {
            oss << parameterName << "=" << parModel.getParameter(parameterName).getValue() << ",";;
        }
        oss << ")";
        LOG(INFO) << oss.str();
    }
    oss.clear();
    oss.str("");
}

std::string OutputUtils::tree2string(bpp::Tree *tree) {
    bpp::Newick treeWriter;
    bpp::TreeTemplate<bpp::Node> ttree(*tree);
    std::ostringstream oss;
    treeWriter.write(ttree, oss);
    std::string out = oss.str();
    return out;
}

bpp::DistanceEstimation DistanceUtils::computeDistanceMethod(std::string seqfilename, bpp::Alphabet *alphabet, bpp::GeneticCode *gCode, std::map<std::string, std::string> &params) {

    // Create a map containing the required parameters passed by the user
    map<std::string, std::string> parmap;
    parmap = params;
    // Overwrite the model parameter
    parmap["model"] = "JC69";
    // Read the data using the non-extended alphabet
    bpp::Fasta seqReader;
    bpp::SequenceContainer *sequences = seqReader.readAlignment(seqfilename, alphabet);
    bpp::SiteContainer *sites = new bpp::VectorSiteContainer(*sequences);
    bpp::SiteContainerTools::changeGapsToUnknownCharacters(*sites);

    bpp::TransitionModel *model = bpp::PhylogeneticsApplicationTools::getTransitionModel(alphabet, gCode, sites, params);

    bpp::DiscreteDistribution *rDist = 0;
    if (model->getNumberOfStates() > model->getAlphabet()->getSize()) {
        //Markov-modulated Markov model!
        rDist = new ConstantRateDistribution();
    } else {
        rDist = bpp::PhylogeneticsApplicationTools::getRateDistribution(params);
    }

    bpp::DistanceEstimation distEstimation(model, rDist, sites, 1, false);

    delete sites;

    return distEstimation;
}

bpp::TreeTemplate<bpp::Node> *
DistanceUtils::computeDistanceTree(bpp::TransitionModel *model, bpp::DiscreteDistribution *rDist, bpp::DistanceEstimation &distEstimation, std::map<std::string, std::string> &params) {

    std::string method = ApplicationTools::getStringParameter("distance.method", params, "nj");
    bpp::ApplicationTools::displayResult("Initial tree reconstruction method", method);
    bpp::TreeTemplate<Node> *tree;
    bpp::AgglomerativeDistanceMethod *distMethod = 0;
    if (method == "wpgma") {
        PGMA *wpgma = new PGMA(true);
        distMethod = wpgma;
    } else if (method == "upgma") {
        PGMA *upgma = new PGMA(false);
        distMethod = upgma;
    } else if (method == "nj") {
        NeighborJoining *nj = new NeighborJoining();
        nj->outputPositiveLengths(true);
        distMethod = nj;
    } else if (method == "bionj") {
        bpp::BioNJ *bionj = new bpp::BioNJ();
        bionj->outputPositiveLengths(true);
        distMethod = bionj;
    } else throw Exception("Unknown initial tree reconstruction method.");

    string type = ApplicationTools::getStringParameter("distance.optimization.method", params, "init");
    ApplicationTools::displayResult("Model parameters estimation method for initial tree", type);
    if (type == "init") type = OptimizationTools::DISTANCEMETHOD_INIT;
    else if (type == "pairwise") type = OptimizationTools::DISTANCEMETHOD_PAIRWISE;
    else if (type == "iterations") type = OptimizationTools::DISTANCEMETHOD_ITERATIONS;
    else throw Exception("Unknown parameter estimation procedure '" + type + "'.");

    unsigned int optVerbose = ApplicationTools::getParameter<unsigned int>("distance.optimization.verbose", params, 2);

    string mhPath = ApplicationTools::getAFilePath("distance.optimization.message_handler", params, false, false);
    OutputStream *messenger = (mhPath == "none") ? 0 : (mhPath == "std") ? ApplicationTools::message : new StlOutputStream(new ofstream(mhPath.c_str(), ios::out));
    ApplicationTools::displayResult("Message handler for initial tree optimisation", mhPath);

    string prPath = ApplicationTools::getAFilePath("distance.optimization.profiler", params, false, false);
    OutputStream *profiler = (prPath == "none") ? 0 : (prPath == "std") ? ApplicationTools::message : new StlOutputStream(new ofstream(prPath.c_str(), ios::out));

    if (profiler) profiler->setPrecision(20);
    ApplicationTools::displayResult("Profiler for initial tree optimisation", prPath);

    // Should I ignore some parameters?
    ParameterList allParameters = model->getParameters();
    allParameters.addParameters(rDist->getParameters());
    ParameterList parametersToIgnore;
    string paramListDesc = ApplicationTools::getStringParameter("distance.optimization.ignore_parameter", params, "", "", true, false);
    bool ignoreBrLen = false;
    StringTokenizer st(paramListDesc, ",");
    while (st.hasMoreToken()) {
        try {
            string param = st.nextToken();
            if (param == "BrLen")
                ignoreBrLen = true;
            else {
                if (allParameters.hasParameter(param)) {
                    Parameter *p = &allParameters.getParameter(param);
                    parametersToIgnore.addParameter(*p);
                } else ApplicationTools::displayWarning("Parameter '" + param + "' not found.");
            }
        } catch (ParameterNotFoundException &pnfe) {
            ApplicationTools::displayError("Parameter '" + pnfe.getParameter() + "' not found, and so can't be ignored!");
        }
    }

    auto nbEvalMax = ApplicationTools::getParameter<unsigned int>("distance.optimization.max_number_f_eval", params, 1000000);
    ApplicationTools::displayResult("Max # ML evaluations for initial tree optimisation", TextTools::toString(nbEvalMax));

    double tolerance = ApplicationTools::getDoubleParameter("distance.optimization.tolerance", params, .000001);
    ApplicationTools::displayResult("Tolerance for initial tree optimisation", TextTools::toString(tolerance));

    //Here it is:
    ofstream warn("warnings", ios::out);
    ApplicationTools::warning = new StlOutputStreamWrapper(&warn);
    tree = OptimizationTools::buildDistanceTree(distEstimation, *distMethod, parametersToIgnore, !ignoreBrLen, type, tolerance, nbEvalMax, profiler, messenger, optVerbose);
    warn.close();
    delete ApplicationTools::warning;
    ApplicationTools::warning = ApplicationTools::message;

    //Output some parameters:
    if (type == OptimizationTools::DISTANCEMETHOD_ITERATIONS) {
        // Write parameters to screen:
        ParameterList parameters = model->getParameters();
        for (unsigned int i = 0; i < parameters.size(); i++) {
            ApplicationTools::displayResult(parameters[i].getName(), TextTools::toString(parameters[i].getValue()));
        }
        parameters = rDist->getParameters();
        for (unsigned int i = 0; i < parameters.size(); i++) {
            ApplicationTools::displayResult(parameters[i].getName(), TextTools::toString(parameters[i].getValue()));
        }
        // Write parameters to file:
        string parametersFile = ApplicationTools::getAFilePath("distance.output.estimates", params, false, false);
        if (parametersFile != "none") {
            ofstream out(parametersFile.c_str(), ios::out);
            parameters = model->getParameters();
            for (unsigned int i = 0; i < parameters.size(); i++) {
                out << parameters[i].getName() << " = " << parameters[i].getValue() << endl;
            }
            parameters = rDist->getParameters();
            for (unsigned int i = 0; i < parameters.size(); i++) {
                out << parameters[i].getName() << " = " << parameters[i].getValue() << endl;
            }
            out.close();
        }
    }

    delete distMethod;

    return tree;
}

void DistanceMethodsUtils::computeDistanceMatrix() {

    // Create a map containing the required parameters passed by the user
    map<std::string, std::string> parmap;
    parmap = params_;
    // Overwrite the model parameter
    parmap["model"] = "JC69";
    // Read the data using the non-extended alphabet
    bpp::Fasta seqReader;
    bpp::SequenceContainer *sequences = seqReader.readAlignment(seqfilename_, alphabet_);
    bpp::SiteContainer *sites = new bpp::VectorSiteContainer(*sequences);
    bpp::SiteContainerTools::changeGapsToUnknownCharacters(*sites);

    model_ = bpp::PhylogeneticsApplicationTools::getTransitionModel(alphabet_, gCode_, sites, parmap);

    if (model_->getNumberOfStates() > model_->getAlphabet()->getSize()) {
        //Markov-modulated Markov model!
        rdist_ = new ConstantRateDistribution();
    } else {
        rdist_ = bpp::PhylogeneticsApplicationTools::getRateDistribution(params_);
    }

    bpp::DistanceEstimation distEstimation(model_, rdist_, sites, 1, false);
    distEstimation_ = &distEstimation;

    delete sites;

    //return distEstimation;
}


bpp::TreeTemplate<bpp::Node> *DistanceMethodsUtils::computeDistanceTree() {
    std::string method = ApplicationTools::getStringParameter("distance.method", params_, "nj");
    bpp::ApplicationTools::displayResult("Initial tree reconstruction method", method);
    bpp::TreeTemplate<Node> *tree;
    bpp::AgglomerativeDistanceMethod *distMethod = 0;
    if (method == "wpgma") {
        PGMA *wpgma = new PGMA(true);
        distMethod = wpgma;
    } else if (method == "upgma") {
        PGMA *upgma = new PGMA(false);
        distMethod = upgma;
    } else if (method == "nj") {
        NeighborJoining *nj = new NeighborJoining();
        nj->outputPositiveLengths(true);
        distMethod = nj;
    } else if (method == "bionj") {
        bpp::BioNJ *bionj = new bpp::BioNJ();
        bionj->outputPositiveLengths(true);
        distMethod = bionj;
    } else throw Exception("Unknown initial tree reconstruction method.");

    string type = ApplicationTools::getStringParameter("distance.optimization.method", params_, "init");
    ApplicationTools::displayResult("Model parameters estimation method for initial tree", type);
    if (type == "init") type = OptimizationTools::DISTANCEMETHOD_INIT;
    else if (type == "pairwise") type = OptimizationTools::DISTANCEMETHOD_PAIRWISE;
    else if (type == "iterations") type = OptimizationTools::DISTANCEMETHOD_ITERATIONS;
    else throw Exception("Unknown parameter estimation procedure '" + type + "'.");

    unsigned int optVerbose = ApplicationTools::getParameter<unsigned int>("distance.optimization.verbose", params_, 2);

    string mhPath = ApplicationTools::getAFilePath("distance.optimization.message_handler", params_, false, false);
    OutputStream *messenger = (mhPath == "none") ? 0 : (mhPath == "std") ? ApplicationTools::message : new StlOutputStream(new ofstream(mhPath.c_str(), ios::out));
    ApplicationTools::displayResult("Message handler for initial tree optimisation", mhPath);

    string prPath = ApplicationTools::getAFilePath("distance.optimization.profiler", params_, false, false);
    OutputStream *profiler = (prPath == "none") ? 0 : (prPath == "std") ? ApplicationTools::message : new StlOutputStream(new ofstream(prPath.c_str(), ios::out));

    if (profiler) profiler->setPrecision(20);
    ApplicationTools::displayResult("Profiler for initial tree optimisation", prPath);

    // Should I ignore some parameters?
    ParameterList allParameters = model_->getParameters();
    allParameters.addParameters(rdist_->getParameters());
    ParameterList parametersToIgnore;
    string paramListDesc = ApplicationTools::getStringParameter("distance.optimization.ignore_parameter", params_, "", "", true, false);
    bool ignoreBrLen = false;
    StringTokenizer st(paramListDesc, ",");
    while (st.hasMoreToken()) {
        try {
            string param = st.nextToken();
            if (param == "BrLen")
                ignoreBrLen = true;
            else {
                if (allParameters.hasParameter(param)) {
                    Parameter *p = &allParameters.getParameter(param);
                    parametersToIgnore.addParameter(*p);
                } else ApplicationTools::displayWarning("Parameter '" + param + "' not found.");
            }
        } catch (ParameterNotFoundException &pnfe) {
            ApplicationTools::displayError("Parameter '" + pnfe.getParameter() + "' not found, and so can't be ignored!");
        }
    }

    auto nbEvalMax = ApplicationTools::getParameter<unsigned int>("distance.optimization.max_number_f_eval", params_, 1000000);
    ApplicationTools::displayResult("Max # ML evaluations for initial tree optimisation", TextTools::toString(nbEvalMax));

    double tolerance = ApplicationTools::getDoubleParameter("distance.optimization.tolerance", params_, .000001);
    ApplicationTools::displayResult("Tolerance for initial tree optimisation", TextTools::toString(tolerance));


    if (distMatrix_) {

        tree = Optimizators::buildDistanceTreeGenericFromDistanceMatrix(distMatrix_, *distMethod, optVerbose);

    } else {
        //Here it is:
        ofstream warn("warnings", ios::out);
        ApplicationTools::warning = new StlOutputStreamWrapper(&warn);

        tree = Optimizators::buildDistanceTreeGeneric(*distEstimation_, *distMethod, parametersToIgnore, !ignoreBrLen, type, tolerance, nbEvalMax, profiler, messenger, optVerbose);

        warn.close();
        delete ApplicationTools::warning;
        ApplicationTools::warning = ApplicationTools::message;

        //Output some parameters:
        if (type == OptimizationTools::DISTANCEMETHOD_ITERATIONS) {
            // Write parameters to screen:
            ParameterList parameters = model_->getParameters();
            for (unsigned int i = 0; i < parameters.size(); i++) {
                ApplicationTools::displayResult(parameters[i].getName(), TextTools::toString(parameters[i].getValue()));
            }
            parameters = rdist_->getParameters();
            for (unsigned int i = 0; i < parameters.size(); i++) {
                ApplicationTools::displayResult(parameters[i].getName(), TextTools::toString(parameters[i].getValue()));
            }
            // Write parameters to file:
            string parametersFile = ApplicationTools::getAFilePath("distance.output.estimates", params_, false, false);
            if (parametersFile != "none") {
                ofstream out(parametersFile.c_str(), ios::out);
                parameters = model_->getParameters();
                for (unsigned int i = 0; i < parameters.size(); i++) {
                    out << parameters[i].getName() << " = " << parameters[i].getValue() << endl;
                }
                parameters = rdist_->getParameters();
                for (unsigned int i = 0; i < parameters.size(); i++) {
                    out << parameters[i].getName() << " = " << parameters[i].getValue() << endl;
                }
                out.close();
            }
        }
    }

    delete distMethod;

    return tree;
}

void DistanceMethodsUtils::setDistanceMatrix() {

}

