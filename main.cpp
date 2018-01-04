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
 * @file main.hpp
 * @author Lorenzo Gatti
 * @author Massimo Maiolo
 * @date 20 12 2017
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

#include <iostream>
/*
* From Core:
*/
#include <Bpp/Io/OutputStream.h>

/*
* From SeqLib:
*/
#include <Bpp/Seq/Alphabet/DNA.h>
#include <Bpp/Seq/Io/Fasta.h>
#include <Bpp/Seq/Container/SiteContainerTools.h>
#include <Bpp/Seq/Alphabet/AlphabetTools.h>
#include <Bpp/Numeric/Prob/ConstantDistribution.h>
/*
* From PhylLib:
*/
#include <Bpp/Phyl/Model/SubstitutionModel.h>
#include <Bpp/Phyl/Model/Nucleotide/JCnuc.h>
#include <Bpp/Phyl/Model/Nucleotide/K80.h>
#include <Bpp/Phyl/Model/Nucleotide/GTR.h>
#include <Bpp/Phyl/Distance/DistanceMethod.h>
#include <Bpp/Phyl/Distance/DistanceEstimation.h>
#include <Bpp/Phyl/Distance/BioNJ.h>
#include <Bpp/Phyl/Io/Newick.h>
#include <Bpp/Phyl/Tree.h>
#include <Bpp/Phyl/TreeTemplate.h>
#include <Bpp/Phyl/BipartitionList.h>

#include <Eigen/Core>
#include <Eigen/Dense>
#include <unsupported/Eigen/MatrixFunctions>

#include <glog/logging.h>

#include <Alignment.hpp>
#include <Likelihood.hpp>
#include <TreeRearrangment.hpp>

#include "utils.hpp"
#include "PIP.hpp"
#include "cli_parser.hpp"


using namespace tshlib;



int main(int argc, char *argv[]) {

    FLAGS_alsologtostderr = true;
    google::InitGoogleLogging("JATI-minimal");
    gflags::ParseCommandLineFlags(&argc, &argv, true);

    LOG(INFO) << "JATI-minimal v.0.0.1 build: cpp0001 " << std::endl;

    //------------------------------------------------------------------------------------------------------------------
    // LOAD MSA FROM FILE
    // Parse fasta file containing aligned DNA sequences


    auto alignment = new Alignment_DNA;

    LOG(INFO) << "Input file: " << FLAGS_input_sequences;

    bpp::Fasta seqReader;
    bpp::SequenceContainer *sequences = seqReader.readAlignment(FLAGS_input_sequences, &bpp::AlphabetTools::DNA_ALPHABET);
    std::vector<std::string> seqNames = sequences->getSequencesNames();

    for(int i=0; i<seqNames.size(); i++){
        std::string stringSeq = sequences->getSequence(seqNames.at(i)).toString();
        alignment->addSequence(seqNames.at(i),stringSeq);
    }
    alignment->getAlignmentSize();
    alignment->align_num_characters.resize((unsigned long) alignment->align_length);
    alignment->align_alphabetsize += 1; // DNA +1 per PIP
    alignment->countNumberCharactersinColumn();



    bpp::SiteContainer *sites = new bpp::VectorSiteContainer(*sequences);

    size_t num_leaves = sequences->getNumberOfSequences();

    LOG(INFO) << "[Sequences in MSA] Leaves: " << num_leaves;
    std::string testSeq = sequences->getSequence(seqNames.at(0)).toString();
    bpp::Site testSite = nullptr;

    std::vector<int> countCharNumber(sites->getNumberOfSites(), 0);

    for (int i=0; i<sites->getNumberOfSites(); i++){
       testSite = sites->getSite((size_t) i);
        for (int j=0; j<testSite.size(); j++){
            if(testSite.getValue(j)>=0){
                countCharNumber.at(i) += 1;
            }

        }
    }

    alignment->align_num_characters = countCharNumber;

    delete sequences;

    //------------------------------------------------------------------------------------------------------------------
    // INIT ROOTED TREE

    auto * newickReader = new bpp::Newick(false); //No comment allowed!
    bpp::Tree *tree = nullptr;
    try {
        tree = newickReader->read(FLAGS_input_tree); // Tree in file MyTestTree.dnd
        LOG(INFO) << "[Tree parser] Input tree has " << tree->getNumberOfLeaves() << " leaves.";
    } catch (bpp::Exception e) {
        LOG(FATAL) << "[Tree parser] Error when reading tree due to: " << e.message();
    }

    // Convert bpp::tree into thslib::utree
    auto ttTree = bpp::TreeTemplate<bpp::Node>(*tree);
    auto utree = new Utree();
    UtreeBppUtils::convertUtree(&ttTree, utree);
    LOG(INFO) << "[Initial Utree Topology] " << utree->printTreeNewick(true);

    delete tree;
    delete newickReader;


    utree->prepareSetADesCountOnNodes((int) alignment->getAlignmentSize(), alignment->align_alphabetsize);
    UtreeUtils::associateNode2Alignment(alignment, utree);

    // Add the root
    utree->addVirtualRootNode();
    utree->rootnode->initialiseLikelihoodComponents((int) alignment->getAlignmentSize(), alignment->align_alphabetsize);


    //------------------------------------------------------------------------------------------------------------------
    // INIT SubModels + Indel

    // Set the substitution model
    bpp::SubstitutionModel *submodel;

    if(!FLAGS_model_substitution.empty()){

        if(FLAGS_model_substitution == "GTR"){
            submodel = new bpp::GTR(&bpp::AlphabetTools::DNA_ALPHABET);
        }

        if(FLAGS_model_substitution == "JC69") {
            submodel = new bpp::JCnuc(&bpp::AlphabetTools::DNA_ALPHABET);
        }

        if(FLAGS_model_substitution == "K80") {
            double k80_kappa = 0.5;
            submodel = new bpp::K80(&bpp::AlphabetTools::DNA_ALPHABET, k80_kappa);

        }

    }

    Eigen::MatrixXd Q;
    Eigen::VectorXd pi;
    double lambda;
    double mu;

    // Extend the substitution model with PIP
    if(FLAGS_model_indels){

        lambda= 0.2;
        mu=0.1;

        submodel = new PIP_Nuc(&bpp::AlphabetTools::DNA_ALPHABET, lambda, mu, submodel);

        // Fill Q matrix
        Q = MatrixBppUtils::Matrix2Eigen(submodel->getGenerator());
        pi = MatrixBppUtils::Vector2Eigen(submodel->getFrequencies());

        std::cout << Q << std::endl;

    }




    //------------------------------------------------------------------------------------------------------------------
    // INITIAL LIKELIHOOD COMPUTATION

    auto likelihood = new Likelihood();

    likelihood->Init(utree, pi, Q, mu, lambda);

    // COMPUTE LK GIVEN TREE TOPOLOGY AND MSA
    double logLK = 0.0;

    // Traverse the tree in post-order filling a list of node ready for traversal
    std::vector<VirtualNode *> allnodes_postorder;
    likelihood->compileNodeList_postorder(allnodes_postorder, utree->rootnode);

    // Set survival probability to each node in the list
    likelihood->setAllIotas(allnodes_postorder);

    // Set deletion probability to each node in the list
    likelihood->setAllBetas(allnodes_postorder);

    // Set insertion histories on each node of the list
    likelihood->setInsertionHistories(allnodes_postorder,*alignment);

    // set probability matrix -- exponential of substitution rates
    likelihood->computePr(allnodes_postorder, alignment->align_alphabetsize);

    // Initialise likelihood components on the tree
    likelihood->computeFV(allnodes_postorder, *alignment); //TODO: Add weight per column

    // Make a backup of the newly computed components
    likelihood->saveLikelihoodComponents();

    //likelihood->loadParametersOperative();
    //likelihood->unloadParametersOperative();

    // Compute the model likelihood
    logLK = likelihood->computePartialLK_WholeAlignment(allnodes_postorder, *alignment);
    VLOG(1) << "[Tree likelihood] -- full traversal -- " << logLK;

    //----------------------------------------------
    // Remove the root
    utree->removeVirtualRootNode();

    //----------------------------------------------
    // Save tree to file
    //utree->saveTreeOnFile("../data/test.txt");

    //----------------------------------------------
    // Print tree structure on console
    utree->printAllNodesNeighbors();



    //------------------------------------------------------------------------------------------------------------------
    // DEFINE, APPLY & REVERT TREE REARRANGEMENTS
    // Get all the nodes between the radius boundaries and for each of them build the move list

    std::chrono::high_resolution_clock::time_point t1 = std::chrono::high_resolution_clock::now();
    unsigned long total_exec_moves = 0;

    int min_radius = 3;  // Minimum radius for an NNI move is 3 nodes
    int max_radius = utree->getMaxNodeDistance(); // Full tree traversing from any node of the tree

    bool computeMoveLikelihood = true;
    std::vector<VirtualNode *> list_vnode_to_root;

    // Print node description with neighbors
    for (auto &vnode:utree->listVNodes) {
        VLOG(2) << "[utree neighbours] " << vnode->printNeighbours() << std::endl;

        // Initialise a new rearrangement list
        auto rearrangmentList = new TreeRearrangment;
        rearrangmentList->initTreeRearrangment(utree, min_radius, max_radius, true, vnode);

        // Get all the target nodes with distance == radius from the source node
        // excluding the starting node.
        rearrangmentList->defineMoves(false);

        // Print the list of moves for the current P node (source node)
        //rearrangmentList.printMoves();

        VLOG(1) << "[tsh] Strategy " << rearrangmentList->mset_strategy << std::endl;
        VLOG(1) << "[utree rearrangment] Found " << rearrangmentList->getNumberOfMoves() << " possible moves for node " << vnode->vnode_name << std::endl;

        std::string start_col_line, end_col_line;

        // For each potential move computed before, apply it to the tree topology, print the resulting newick tree, and revert it.
        for (unsigned long i = 0; i < rearrangmentList->getNumberOfMoves(); i++) {
            logLK = 0;

            VirtualNode *pnode = rearrangmentList->getSourceNode();
            VirtualNode *qnode = rearrangmentList->getMove(i)->getTargetNode();

            // ------------------------------------
            // Prepare the list of nodes involved in the move
            // TODO: This list belongs specifically to the tree-rearrangement definition -> move to TreeRearrangment class
            list_vnode_to_root.clear();
            list_vnode_to_root = utree->computePathBetweenNodes(pnode, qnode);
            list_vnode_to_root.push_back(utree->rootnode);

            // ------------------------------------
            // Apply the move
            rearrangmentList->applyMove(i);

            // ------------------------------------
            // Print root reachability from every node (includes rotations)
            //utree->_testReachingPseudoRoot();

            // ------------------------------------
            // Print tree on file
            //utree->saveTreeOnFile("../data/test.txt");
            bool isLKImproved = false;
            computeMoveLikelihood = true;
            // ------------------------------------
            if (computeMoveLikelihood) {

                // ------------------------------------
                //utree->printAllNodesNeighbors();
                //testSetAinRootPath(MSA_len, alignment, utree, list_vnode_to_root)

                // ------------------------------------
                // Add the root
                utree->addVirtualRootNode();

                // ------------------------------------
                // Compute the full likelihood from the list of nodes involved in the rearrangment
                allnodes_postorder.clear();
                likelihood->compileNodeList_postorder(allnodes_postorder, utree->rootnode);


                likelihood->recombineAllFv(list_vnode_to_root);
                likelihood->setInsertionHistories(list_vnode_to_root,*alignment);

                logLK = LKFunc::LKRearrangment(*likelihood, allnodes_postorder, *alignment);

                //logLK = likelihood->computePartialLK(list_vnode_to_root, *alignment);

                // ------------------------------------
                // Apply brent -- test only
                //VLOG(2) << "[Tree LK] Before Brent: " << logLK;
                double max_lenght = pi[1] * 1.1;
                double min_lenght = pi[1] * 0.9;

                //isLKImproved = Generic_Brent_Lk(&likelihood->pi[1], min_lenght, max_lenght, SMALL, BRENT_ITMAX,
                //                                LKFunc::LKcore , *likelihood, list_vnode_to_root, *alignment, logLK);

                // ------------------------------------
                // Store likelihood of the move
                rearrangmentList->getMove(i)->move_lk = logLK;
                //VLOG(2) << "[Tree LK] Ater Brent: " << logLK;

                // ------------------------------------
                // Remove virtual root
                utree->removeVirtualRootNode();

            }

            // ------------------------------------
            // Some abbellishments for the console output
            if(rearrangmentList->getMove(i)->move_lk>0){
                start_col_line = "\033[1;34m";
                end_col_line = "\033[0m";
            }else{
                start_col_line = "";
                end_col_line = "";

            }

            // Move exection details
            VLOG(2) << "[apply  move]\t" << rearrangmentList->getMove(i)->move_class << "." << std::setfill('0') << std::setw(3) << i
                    << " | (" << isLKImproved <<") " << start_col_line<< rearrangmentList->getMove(i)->move_lk<<end_col_line << "\t"
                    << " | (" << rearrangmentList->getSourceNode()->vnode_name << "->" << rearrangmentList->getMove(i)->getTargetNode()->vnode_name << ")"
                    << "\t[" << rearrangmentList->getMove(i)->move_radius << "] | "
                    << utree->printTreeNewick(true) << std::endl;

            // ------------------------------------
            // Print root reachability from every node (includes rotations)
            //utree->_testReachingPseudoRoot();

            // ------------------------------------
            // Revert the move, and return to the original tree
            rearrangmentList->revertMove(i);

            // ------------------------------------
            // Print tree on file
            //utree->saveTreeOnFile("../data/test.txt");

            // ------------------------------------
            // Add the root
            utree->addVirtualRootNode();

            computeMoveLikelihood = false;
            // ------------------------------------
            if (FLAGS_lkmove_bothways) {

                // ------------------------------------
                // Compute the full likelihood from the list of nodes involved in the rearrangment
                likelihood->recombineAllFv(list_vnode_to_root);
                likelihood->setInsertionHistories(list_vnode_to_root,*alignment);

                logLK = LKFunc::LKRearrangment(*likelihood, allnodes_postorder, *alignment);

                // ------------------------------------
                // Store likelihood of the move
                rearrangmentList->getMove(i)->move_lk = logLK;

            }else{

                likelihood->restoreLikelihoodComponents();
                //likelihood->setInsertionHistories(allnodes_postorder,*alignment);
                logLK = LKFunc::LKRearrangment(*likelihood, allnodes_postorder, *alignment);
            }
            // ------------------------------------
            // Remove virtual root
            utree->removeVirtualRootNode();
            // ------------------------------------
            // Some abbellishments for the console output
            if(logLK>0){
                start_col_line = "\033[1;34m";
                end_col_line = "\033[0m";
            }else{
                start_col_line = "";
                end_col_line = "";

            }
            // Move exection details
            VLOG(2) << "[revert move]\t" << rearrangmentList->getMove(i)->move_class << "." << std::setfill('0') << std::setw(3) << i
                    << " | (" << isLKImproved <<") " << start_col_line<< logLK <<end_col_line << "\t"
                    << " | (" << rearrangmentList->getMove(i)->getTargetNode()->vnode_name << "->" << rearrangmentList->getSourceNode()->vnode_name << ")"
                    << "\t[" << rearrangmentList->getMove(i)->move_radius << "] | "
                    << utree->printTreeNewick(true) << std::endl;

            // ------------------------------------
            // Print root reachability from every node (includes rotations)
            //utree->_testReachingPseudoRoot();

            // ------------------------------------
            // Count moves performed
            total_exec_moves += rearrangmentList->getNumberOfMoves() * 2;
        }

        // ------------------------------------
        // Clean memory
        delete rearrangmentList;
    }

    // ------------------------------------
    std::chrono::high_resolution_clock::time_point t2 = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::microseconds>(t2 - t1).count();

    VLOG(0) << "Moves applied and reverted: " << total_exec_moves << std::endl;
    VLOG(0) << "Elapsed time: " << duration << " microseconds" << std::endl;
    VLOG(0) << "*** " << (double) duration / total_exec_moves << " microseconds/move *** " << std::endl;

    //treesearchheuristics::testTSH(utree, TreeSearchHeuristics::classic_Mixed);


    exit(0);

}