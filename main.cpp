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
#include <fstream>
/*
* From Core:
*/
#include <Bpp/Io/OutputStream.h>
#include <Bpp/App/BppApplication.h>
#include <Bpp/App/ApplicationTools.h>

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
#include <Bpp/Phyl/App/PhylogeneticsApplicationTools.h>
#include <Bpp/Phyl/Model/Nucleotide/JCnuc.h>
#include <Bpp/Phyl/Model/Nucleotide/K80.h>
#include <Bpp/Phyl/Model/Nucleotide/GTR.h>
#include <Bpp/Phyl/Distance/DistanceMethod.h>
#include <Bpp/Phyl/Distance/DistanceEstimation.h>
#include <Bpp/Phyl/Model/RateDistribution/ConstantRateDistribution.h>
#include <Bpp/Phyl/Likelihood/RHomogeneousTreeLikelihood.h>
#include <Bpp/Phyl/Distance/BioNJ.h>
#include <Bpp/Phyl/Io/Newick.h>

#include <Eigen/Core>
#include <Eigen/Dense>
#include <unsupported/Eigen/MatrixFunctions>

#include <glog/logging.h>
#include <chrono>
#include <Alignment.hpp>
#include <Likelihood.hpp>
#include <TreeRearrangment.hpp>

#include "Version.hpp"
#include "Utilities.hpp"
#include "PIP.hpp"
#include "CommandLineFlags.hpp"
#include "ExtendedAlphabet.hpp"
#include "RHomogeneousTreeLikelihood_PIP.hpp"

#include "progressivePIP.hpp"


using namespace tshlib;


int main(int argc, char *argv[]) {

    FLAGS_alsologtostderr = true;
    gflags::SetUsageMessage("some usage message");
    gflags::SetVersionString(software::build);
    google::InitGoogleLogging(software::name.c_str());
    gflags::ParseCommandLineFlags(&argc, &argv, true);

    LOG(INFO) << software::desc;
    bpp::BppApplication bppml(argc, argv, software::desc);

    //------------------------------------------------------------------------------------------------------------------
    // LOAD MSA FROM FILE
    // Parse fasta file containing aligned DNA sequences


    auto alignment = new Alignment_DNA;

    LOG(INFO) << "[Input alignment file] " << FLAGS_input_sequences;

    NucleicAlphabet *alpha;

    if (FLAGS_model_indels) {
        alpha = new bpp::DNA_EXTENDED();
        //alpha = new bpp::DNA();
    } else {
        alpha = new bpp::DNA();
    }

    bpp::Fasta seqReader;
    bpp::SequenceContainer *sequences;
    bpp::SiteContainer *sites;

    if (FLAGS_alignment) {
        sequences = seqReader.readSequences(FLAGS_input_sequences, alpha);
    } else {
        sequences = seqReader.readAlignment(FLAGS_input_sequences, alpha);
        seqReader.readSequences(FLAGS_input_sequences, alpha);
        std::vector<std::string> seqNames = sequences->getSequencesNames();

        for (int i = 0; i < seqNames.size(); i++) {
            std::string stringSeq = sequences->getSequence(seqNames.at(i)).toString();
            alignment->addSequence(seqNames.at(i), stringSeq);
        }
        alignment->getAlignmentSize();
        alignment->align_num_characters.resize((unsigned long) alignment->align_length);
        alignment->align_alphabetsize += 1; // DNA +1 per PIP
        alignment->countNumberCharactersinColumn();

        sites = new bpp::VectorSiteContainer(*sequences);
        size_t num_leaves = sequences->getNumberOfSequences();


        LOG(INFO) << "[Alignment] Input alignmetn has " << num_leaves << " sequences";

        /*
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
        */

    }




    //------------------------------------------------------------------------------------------------------------------
    // INIT ROOTED TREE

    auto *newickReader = new bpp::Newick(false); //No comment allowed!
    bpp::Tree *tree = nullptr;
    try {
        tree = newickReader->read(FLAGS_input_tree); // Tree in file MyTestTree.dnd
        LOG(INFO) << "[Input tree file] " << FLAGS_input_tree;
        LOG(INFO) << "[Tree parser] Input tree has " << tree->getNumberOfLeaves() << " leaves.";
    } catch (bpp::Exception e) {
        LOG(FATAL) << "[Tree parser] Error when reading tree due to: " << e.message();
    }

    // Convert bpp::tree into thslib::utree
    auto ttTree = bpp::TreeTemplate<bpp::Node>(*tree);
    auto utree = new Utree();
    UtreeBppUtils::treemap tm;
    UtreeBppUtils::convertTree_b2u(&ttTree, utree, tm);
    LOG(INFO) << "[Initial Utree Topology] " << utree->printTreeNewick(true);

    delete tree;
    delete newickReader;

    //------------------------------------------------------------------------------------------------------------------
    // Printing of the bidirectional map [test]
    VLOG(3) << "Bidirectional map size: "<<  tm.size();
    UtreeBppUtils::treemap::left_map& map_view = tm.left;
    for (UtreeBppUtils::treemap::left_map::const_iterator it(map_view.begin()), end(map_view.end()); it != end; ++it) {
        VLOG(3) << (*it).first->getFather()->getId() << " --> " << (*it).second->getNodeUp()->getNodeName();
    }

    //------------------------------------------------------------------------------------------------------------------

    if (!FLAGS_alignment) {
        utree->prepareSetADesCountOnNodes((int) alignment->getAlignmentSize(), alignment->align_alphabetsize);
        UtreeUtils::associateNode2Alignment(alignment, utree);

        // Add the root
        utree->addVirtualRootNode();
        utree->rootnode->initialiseLikelihoodComponents((int) alignment->getAlignmentSize(),
                                                        alignment->align_alphabetsize);
    } else {

        UtreeBppUtils::associateNode2Alignment(sequences, utree);
        utree->addVirtualRootNode();


    }

    //------------------------------------------------------------------------------------------------------------------
    // INIT SubModels + Indel

    // Set the substitution model
    bpp::SubstitutionModel *submodel = nullptr;
    unique_ptr<bpp::GeneticCode> gCode;
    map<std::string, std::string> parmap;
    bpp::TransitionModel *transmodel = nullptr;
    bpp::DiscreteDistribution *rDist = nullptr;
    bpp::TreeLikelihood *tl;

    Eigen::MatrixXd Q;
    Eigen::VectorXd pi;
    double lambda;
    double mu;

    if (!FLAGS_model_substitution.empty()) {
        unique_ptr<GeneticCode> gCode;
        parmap["model"] = FLAGS_model_substitution;

        submodel = bpp::PhylogeneticsApplicationTools::getSubstitutionModel(alpha, gCode.get(), sites, parmap, "", true, false, 0);
        if (!FLAGS_alignment) { if(FLAGS_model_setfreqsfromdata) submodel->setFreqFromData(*sites); }
        //submodel->setFreqFromData(*sites);

        rDist = new bpp::ConstantRateDistribution();
        if (!FLAGS_alignment) { bpp::SiteContainerTools::changeGapsToUnknownCharacters(*sites); }


    }

    // Extend the substitution model with PIP
    if (FLAGS_model_indels) {

        lambda=FLAGS_lambda_PIP;
        mu=FLAGS_mu_PIP;

        //VLOG(1) << alpha->getGapCharacterCode();
        //auto testCS = new CanonicalStateMap(new bpp::DNA(), true);
        //auto testCS_G = new CanonicalStateMap(new bpp::DNA_EXTENDED(), false);
        //VLOG(1) << "size: " << testCS->getAlphabet()->getSize();

        submodel = new PIP_Nuc(alpha, lambda, mu, submodel);
        //submodel->setFreqFromData(*sites);

        // Fill Q matrix
        Q = MatrixBppUtils::Matrix2Eigen(submodel->getGenerator());
        pi = MatrixBppUtils::Vector2Eigen(submodel->getFrequencies());
        //std::cerr << Q << std::endl;
        //std::cerr << pi << std::endl;

    }
    VLOG(1) << "[Substitution model] Number of states: " << (int) submodel->getNumberOfStates();

    bpp::StdStr s1;
    bpp::PhylogeneticsApplicationTools::printParameters(submodel, s1, 1, true);
    VLOG(1) << s1.str();

    //------------------------------------------------------------------------------------------------------------------
    // INITIAL LIKELIHOOD COMPUTATION

    double logLK = 0;
    auto likelihood = new Likelihood();
    std::vector<VirtualNode *> allnodes_postorder;
    progressivePIP::ProgressivePIPResult MSA;
    if (!FLAGS_alignment) {
        tree = UtreeBppUtils::convertTree_u2b(utree);

        bpp::RowMatrix<double> testProb;
        if (!FLAGS_model_indels) {

            transmodel = bpp::PhylogeneticsApplicationTools::getTransitionModel(alpha, gCode.get(), sites, parmap, "", true, false, 0);
            tl = new bpp::RHomogeneousTreeLikelihood(*tree, *sites, transmodel, rDist, false, false, false);

        } else {

            unique_ptr<TransitionModel> test;
            test.reset(submodel);
            transmodel = test.release();

            tl = new bpp::RHomogeneousTreeLikelihood_PIP(*tree, *sites, transmodel, rDist, false, false, false);

        }


        VLOG(1) << "[Transition model] Number of states: " << (int) transmodel->getNumberOfStates();

        tl->initialize();
        logLK = tl->getValue();


        VLOG(1) << "[Tree likelihood] -- full traversal -- (on model " << submodel->getName() << ") = " << -logLK << " \t[BPP METHODS]";


    }


    if (FLAGS_model_indels) {
        //tl = new bpp::RHomogeneousTreeLikelihood_PIP(*tree, *sites, transmodel, rDist, false, false, false);
        likelihood->Init(utree, pi, Q, mu, lambda);

        // COMPUTE LK GIVEN TREE TOPOLOGY AND MSA
        logLK = 0.0;

        // Traverse the tree in post-order filling a list of node ready for traversal
        likelihood->compileNodeList_postorder(allnodes_postorder, utree->rootnode);

        // Set survival probability to each node in the list
        likelihood->setAllIotas(allnodes_postorder);

        // Set deletion probability to each node in the list
        likelihood->setAllBetas(allnodes_postorder);

        // set probability matrix -- exponential of substitution rates
        likelihood->computePr(allnodes_postorder, alignment->align_alphabetsize);

        if (FLAGS_alignment) {

            VLOG(1) << "[ProPIP] starting MSA inference...";

            VirtualNode *root = utree->rootnode;
            MSA = progressivePIP::compute_DP3D_PIP_tree_cross(root, likelihood, sequences, alpha, 1.0, false);
            //sites = new bpp::VectorSiteContainer(*sequences);

            sequences = new bpp::VectorSequenceContainer(alpha);

            for (int i = 0; i < MSA.MSAs.size(); i++) {
                //auto sequence = new bpp::BasicSequence(MSA.MSAs.at(i).first,MSA.MSAs.at(i).second,alpha);
                //sequences->setSequence(MSA.MSAs.at(i).first,*sequence,true);
                sequences->addSequence(*(new bpp::BasicSequence(MSA.MSAs.at(i).first,MSA.MSAs.at(i).second,alpha)), true);
            }

            sites = new bpp::VectorSiteContainer(*sequences);

            delete sequences;

            bpp::Fasta seqWriter;
            seqWriter.writeAlignment(FLAGS_output_msa,*sites,true);

            std::ofstream lkFile;
            lkFile << std::setprecision(18);
            lkFile.open (FLAGS_output_lk);
            lkFile << MSA.score;
            lkFile.close();

            VLOG(1) << "LK ProPIP" << MSA.score;

            VLOG(1) << "[ProPIP] ...done";

            //exit(EXIT_SUCCESS);

            /*
            int dim=20;
            int num_times=10;
            bpp::RowMatrix<double> AA;
            AA.resize(dim,dim);

            Eigen::MatrixXd A;

            A.resize(dim,dim);
            for(int i=0;i<dim;i++){
                for(int j=0;j<dim;j++){
                    double val = (double)(i/10+j/10);
                    A(i,j) = val;
                    AA(i,j) = val;
                }
            }

            bpp::RowMatrix<double> BB;
            BB.resize(dim,dim);
            std::chrono::high_resolution_clock::time_point t1_BPP = std::chrono::high_resolution_clock::now();
            for(int i=0;i<num_times;i++){
                bpp::MatrixTools::exp(AA,BB);
            }
            std::chrono::high_resolution_clock::time_point t2_BPP = std::chrono::high_resolution_clock::now();
            auto durationBPP = std::chrono::duration_cast<std::chrono::microseconds>(t2_BPP - t1_BPP).count();
            VLOG(0) << "Elapsed time (BPP): " << durationBPP << " microseconds" << std::endl;

            Eigen::MatrixXd B;
            B.resize(dim,dim);
            std::chrono::high_resolution_clock::time_point t1_EIG = std::chrono::high_resolution_clock::now();
            for(int i=0;i<num_times;i++){
                B=A.exp();
            }
            std::chrono::high_resolution_clock::time_point t2_EIG = std::chrono::high_resolution_clock::now();
            auto durationEIG = std::chrono::duration_cast<std::chrono::microseconds>(t2_EIG - t1_EIG).count();
            VLOG(0) << "Elapsed time (EIGEN): " << durationEIG << " microseconds" << std::endl;

            //std::cout<<B;
            //std::cout<<std::endl;
            exit(1);
            */

        }


        // Set insertion histories on each node of the list
        likelihood->setInsertionHistories(allnodes_postorder, *alignment);

        // Initialise likelihood components on the tree
        likelihood->computeFV(allnodes_postorder, *alignment); //TODO: Add weight per column

        // Make a backup of the newly computed components
        likelihood->saveLikelihoodComponents();

        //likelihood->loadParametersOperative();
        //likelihood->unloadParametersOperative();

        // Compute the model likelihood
        logLK = likelihood->computePartialLK_WholeAlignment(allnodes_postorder, *alignment);
        VLOG(1) << "[Tree likelihood] -- full traversal -- (on model " << submodel->getName() << ") = " << logLK << " \t[TSHLIB METHODS]";
    }
    //----------------------------------------------
    // Remove the root
    utree->removeVirtualRootNode();

    //----------------------------------------------
    // Save tree to file
    //utree->saveTreeOnFile("../data/test.txt");

    //----------------------------------------------
    // Print tree structure on console
    //utree->printAllNodesNeighbors();

    //------------------------------------------------------------------------------------------------------------------
    // DEFINE, APPLY & REVERT TREE REARRANGEMENTS
    // Get all the nodes between the radius boundaries and for each of them build the move list

    std::chrono::high_resolution_clock::time_point t1 = std::chrono::high_resolution_clock::now();
    unsigned long total_exec_moves = 0;
    int min_radius;
    int max_radius;

    if (FLAGS_optim_topology.find("full-search") != std::string::npos) {

        min_radius = 3;  // Minimum radius for an NNI move is 3 nodes
        max_radius = utree->getMaxNodeDistance(); // Full tree traversing from any node of the tree

    } else if (FLAGS_optim_topology.find("nni-search") != std::string::npos) {

        min_radius = 3;
        max_radius = 3;

    } else if (FLAGS_optim_topology.find("spr-search") != std::string::npos) {

        min_radius = 4;
        max_radius = utree->getMaxNodeDistance();

    } else {

        LOG(FATAL) << "Exiting program without tree search optimisation";
    }


    std::vector<VirtualNode *> list_vnode_to_root;

    // Print node description with neighbors
    for (auto &vnode:utree->listVNodes) {
        //VLOG(2) << "[utree neighbours] " << vnode->printNeighbours() << std::endl;

        // Initialise a new rearrangement list
        auto rearrangmentList = new TreeRearrangment;
        rearrangmentList->initTreeRearrangment(utree, min_radius, max_radius, true, vnode);

        // Get all the target nodes with distance == radius from the source node
        // excluding the starting node.
        rearrangmentList->defineMoves(false);

        // Print the list of moves for the current P node (source node)
        //rearrangmentList.printMoves();

        VLOG(1) << "[utree rearrangment] [" << rearrangmentList->mset_strategy << "] Found "
                << rearrangmentList->getNumberOfMoves() << " candidate moves for node "
                << vnode->vnode_name;

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
            bool computeMoveLikelihood = true;

            // ------------------------------------
            // Add the root
            utree->addVirtualRootNode();
            // ------------------------------------
            if (computeMoveLikelihood && FLAGS_model_indels) {

                // ------------------------------------
                //utree->printAllNodesNeighbors();
                //testSetAinRootPath(MSA_len, alignment, utree, list_vnode_to_root)



                // ------------------------------------
                // Compute the full likelihood from the list of nodes involved in the rearrangment
                allnodes_postorder.clear();
                likelihood->compileNodeList_postorder(allnodes_postorder, utree->rootnode);


                likelihood->recombineAllFv(list_vnode_to_root);
                likelihood->setInsertionHistories(list_vnode_to_root, *alignment);

                logLK = LKFunc::LKRearrangment(*likelihood, list_vnode_to_root, *alignment);



                //std::vector<tshlib::VirtualNode *> &listNodes, UtreeBppUtils::treemap *tm
                //logLK = tl->computeLikelihoodOnTreeRearrangment(list_vnode_to_root, tm);




                //logLK = likelihood->computePartialLK(list_vnode_to_root, *alignment);

                // ------------------------------------
                // Apply brent -- test only
                //VLOG(2) << "[Tree LK] Before Brent: " << logLK;
                //double max_lenght = pi[1] * 1.1;
                //double min_lenght = pi[1] * 0.9;

                //isLKImproved = Generic_Brent_Lk(&likelihood->pi[1], min_lenght, max_lenght, SMALL, BRENT_ITMAX,
                //                                LKFunc::LKcore , *likelihood, list_vnode_to_root, *alignment, logLK);

                // ------------------------------------
                // Store likelihood of the move
                rearrangmentList->getMove(i)->move_lk = logLK;
                //VLOG(2) << "[Tree LK] Ater Brent: " << logLK;



            }

            if (!FLAGS_model_indels) {
                tree = UtreeBppUtils::convertTree_u2b(utree);
                tl = new bpp::RHomogeneousTreeLikelihood(*tree, *sites, transmodel, rDist, false, false, false);
                tl->initialize();
                logLK = tl->getLogLikelihood();
                rearrangmentList->getMove(i)->move_lk = logLK;
            }


            // ------------------------------------
            // Remove virtual root
            utree->removeVirtualRootNode();

            // ------------------------------------
            // Some abbellishments for the console output
            if (rearrangmentList->getMove(i)->move_lk > 0) {
                start_col_line = "\033[1;34m";
                end_col_line = "\033[0m";
            } else {
                start_col_line = "";
                end_col_line = "";

            }

            // Move exection details
            VLOG(2) << "[apply  move]\t" << rearrangmentList->getMove(i)->move_class << "." << std::setfill('0') << std::setw(3) << i
                    << " | (" << isLKImproved << ") " << start_col_line << rearrangmentList->getMove(i)->move_lk << end_col_line << "\t"
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
            if (FLAGS_model_indels) {
                if (FLAGS_lkmove_bothways) {

                    // ------------------------------------
                    // Compute the full likelihood from the list of nodes involved in the rearrangment
                    likelihood->recombineAllFv(list_vnode_to_root);
                    likelihood->setInsertionHistories(list_vnode_to_root, *alignment);

                    logLK = LKFunc::LKRearrangment(*likelihood, list_vnode_to_root, *alignment);

                    // ------------------------------------
                    // Store likelihood of the move
                    rearrangmentList->getMove(i)->move_lk = logLK;

                } else {

                    likelihood->restoreLikelihoodComponents();
                    //likelihood->setInsertionHistories(allnodes_postorder,*alignment);
                    logLK = LKFunc::LKRearrangment(*likelihood, allnodes_postorder, *alignment);
                }

            }
            if (!FLAGS_model_indels) {
                tree = UtreeBppUtils::convertTree_u2b(utree);
                tl = new bpp::RHomogeneousTreeLikelihood(*tree, *sites, transmodel, rDist, false, false, false);
                tl->initialize();
                logLK = tl->getLogLikelihood();
                rearrangmentList->getMove(i)->move_lk = logLK;
            }

            // ------------------------------------
            // Remove virtual root
            utree->removeVirtualRootNode();
            // ------------------------------------
            // Some abbellishments for the console output
            if (logLK > 0) {
                start_col_line = "\033[1;34m";
                end_col_line = "\033[0m";
            } else {
                start_col_line = "";
                end_col_line = "";

            }
            // Move exection details
            VLOG(2) << "[revert move]\t" << rearrangmentList->getMove(i)->move_class << "." << std::setfill('0') << std::setw(3) << i
                    << " | (" << isLKImproved << ") " << start_col_line << logLK << end_col_line << "\t"
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

    delete sequences;
    exit(0);

}