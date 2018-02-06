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
#include <Bpp/Seq/App/SequenceApplicationTools.h>

#include "Version.hpp"
#include "Utilities.hpp"
#include "PIP.hpp"
#include "CommandLineFlags.hpp"
#include "ExtendedAlphabet.hpp"
#include "RHomogeneousTreeLikelihood_PIP.hpp"

#include "progressivePIP.hpp"
#include "JATIApplication.hpp"


using namespace tshlib;

int main(int argc, char *argv[]) {

    FLAGS_alsologtostderr = true;
    gflags::SetUsageMessage("some usage message");
    gflags::SetVersionString(software::build);
    google::InitGoogleLogging(software::name.c_str());
    //gflags::ParseCommandLineFlags(&argc, &argv, true);

    LOG(INFO) << "*****************************************************************************************************************************************";
    LOG(INFO) << "* " << software::desc << "  *";
    LOG(INFO) << "* Authors: Lorenzo Gatti & Massimo Maiolo                                                                                               *";
    LOG(INFO) << "* ------------------------------------------------------------------------------------------------------------------------------------- *";
    LOG(INFO) << "* Based on Bio++ by J. Dutheil, B. Boussau, L. Guéguen, M. Groussin                                                                     *";
    LOG(INFO) << "* Inspired on codonPhyML (Zanetti M. et al.)                                                                                            *";
    LOG(INFO) << "* Inspired on PrographMSA (Szalkowski A. et al.)                                                                                        *";
    LOG(INFO) << "* Implements the Poisson Indel Model (Bouchard-Cote A. et al.)                                                                          *";
    LOG(INFO) << "*****************************************************************************************************************************************";

    try {
        bpp::JATIApplication jatiapp(argc, argv, software::desc);
        jatiapp.startTimer();
        std::string PAR_input_sequences = ApplicationTools::getAFilePath("input_sequences", jatiapp.getParams(), true, true, "", false, "", 1);
        bool PAR_alignment = ApplicationTools::getBooleanParameter("alignment", jatiapp.getParams(), false);
        std::string PAR_input_tree = ApplicationTools::getAFilePath("input_tree", jatiapp.getParams(), false, true, "", false, "", 1);
        std::string PAR_model_substitution = ApplicationTools::getStringParameter("model_substitution", jatiapp.getParams(), "JC69", "", true, true);
        bool PAR_model_setfreqsfromdata = ApplicationTools::getBooleanParameter("model_setfreqsfromdata", jatiapp.getParams(), false);
        bool PAR_model_indels = ApplicationTools::getBooleanParameter("model_indels", jatiapp.getParams(), false);
        std::string PAR_optim_topology = ApplicationTools::getStringParameter("optim_topology", jatiapp.getParams(), "full-search", "", true, true);

        /* ***************************************************
         * Standard workflow
         * [INPUT]
         * 1. tree + alignment => (1.1) just parse everything
         * 2. alignment  => (2.1) parse alignment => (2.2) generate tree using bioNJ
         * 3. sequences  => (3.1) parse sequences => (3.2) generate tree using bioNJ => (3.3) perform alignment
         * (it should not be supported in production)
         * 4. sequences + tree => (4.1) parse sequence => (4.2) parse tree => (4.3) perform alignment
         */
        // Get alphabet from parameters
        bpp::Alphabet *alphabet = bpp::SequenceApplicationTools::getAlphabet(jatiapp.getParams(), "", false, false);
        unique_ptr<GeneticCode> gCode;
        bpp::CodonAlphabet *codonAlphabet = dynamic_cast<bpp::CodonAlphabet *>(alphabet);
        if (codonAlphabet) {
            std::string codeDesc = ApplicationTools::getStringParameter("genetic_code", jatiapp.getParams(), "Standard", "", true, true);
            //ApplicationTools::displayResult("Genetic Code", codeDesc);
            gCode.reset(bpp::SequenceApplicationTools::getGeneticCode(codonAlphabet->getNucleicAlphabet(), codeDesc));
        }

        bpp::Fasta seqReader;
        bpp::SequenceContainer *sequences;
        bpp::SiteContainer *sites;
        NucleicAlphabet *alpha;

        // The alphabet object should be set according to the correct alphabet

        if (PAR_model_indels) {
            alpha = new bpp::DNA_EXTENDED();
        } else {
            alpha = new bpp::DNA();
        }
        //---------------------------------------
        // Input: data
        try {


            LOG(INFO) << "[Input data parser] " << PAR_input_sequences;

            if (PAR_alignment) {
                // If the user requires the computation of an alignment, then the input file is made of unaligned sequences
                sequences = seqReader.readSequences(PAR_input_sequences, alpha);

            } else {

                sequences = seqReader.readAlignment(PAR_input_sequences, alpha);
                seqReader.readSequences(PAR_input_sequences, alpha);
                //std::vector<std::string> seqNames = sequences->getSequencesNames();
                sites = new bpp::VectorSiteContainer(*sequences);

            }

            size_t num_leaves = sequences->getNumberOfSequences();
            LOG(INFO) << "[Input data parser] Sequences " << num_leaves << " sequences";

        } catch (bpp::Exception e) {
            LOG(FATAL) << "[Input data parser] Error when reading sequence file due to: " << e.message();
        }

        //---------------------------------------
        // Input: tree
        bpp::Tree *tree = nullptr;

        if (PAR_input_tree.find("none") == std::string::npos) {
            try {
                auto *newickReader = new bpp::Newick(false); //No comment allowed!
                tree = newickReader->read(PAR_input_tree); // Tree in file MyTestTree.dnd
                LOG(INFO) << "[Input tree file] " << PAR_input_tree;
                LOG(INFO) << "[Tree parser] Input tree has " << tree->getNumberOfLeaves() << " leaves.";
                LOG(INFO) << "[Initial Utree Topology]";
                delete newickReader;

            } catch (bpp::Exception e) {

                LOG(FATAL) << "[Tree parser] Error when reading tree due to: " << e.message();
            }
        } else {

            // Compute bioNJ tree using the JC69 model
            map<std::string, std::string> parmap;
            parmap["model"] = "JC69";

            bpp::SubstitutionModel *submodel_bioNJ = bpp::PhylogeneticsApplicationTools::getSubstitutionModel(alpha, gCode.get(), sites, parmap, "", true, false, 0);
            bpp::DiscreteDistribution *rDist = new bpp::ConstantRateDistribution();
            bpp::SiteContainerTools::changeGapsToUnknownCharacters(*sites);
            bpp::DistanceEstimation distanceMethod(submodel_bioNJ, rDist, sites);
            bpp::DistanceMatrix *distances = distanceMethod.getMatrix();
            bpp::BioNJ bionj(*distances);
            tree = bionj.getTree();

            // Remove all the object used only for computing bioNJ tree
            delete submodel_bioNJ;
            delete rDist;
        }
        // Rename internal nodes with standard Vxx * where xx is a progressive number
        tree->setNodeName(tree->getRootId(), "root");
        UtreeBppUtils::renameInternalNodes(tree);

        // Convert the bpp into utree for tree-search engine
        auto utree = new Utree();
        UtreeBppUtils::treemap tm;
        UtreeBppUtils::convertTree_b2u(tree, utree, tm);
        VLOG(3) << "Bidirectional map size: " << tm.size();
        VLOG(3) << "[Initial Utree Topology] " << utree->printTreeNewick(true);

        UtreeBppUtils::associateNode2Alignment(sequences, utree);
        utree->addVirtualRootNode();

        // Once the tree has the root, then map it as well
        tm.insert(UtreeBppUtils::nodeassoc(tree->getRootId(), utree->rootnode));

        //------------------------------------------------------------------------------------------------------------------
        // INIT SubModels + Indel
        // Set the substitution model
        bpp::SubstitutionModel *submodel = nullptr;
        //unique_ptr<bpp::GeneticCode> gCode;
        map<std::string, std::string> parmap;
        bpp::TransitionModel *transmodel = nullptr;
        bpp::DiscreteDistribution *rDist = nullptr;
        bpp::AbstractHomogeneousTreeLikelihood *tl;

        Eigen::MatrixXd Q;
        Eigen::VectorXd pi;
        double lambda;
        double mu;


        if (!PAR_model_substitution.empty()) {
            //unique_ptr<GeneticCode> gCode;
            parmap["model"] = PAR_model_substitution;

            submodel = bpp::PhylogeneticsApplicationTools::getSubstitutionModel(alpha, gCode.get(), sites, parmap, "", true, false, 0);
            if (!PAR_alignment) { if (PAR_model_setfreqsfromdata) submodel->setFreqFromData(*sites); }
            //submodel->setFreqFromData(*sites);

            rDist = new bpp::ConstantRateDistribution();
            if (!PAR_alignment) { bpp::SiteContainerTools::changeGapsToUnknownCharacters(*sites); }

        }

        // Extend the substitution model with PIP
        if (PAR_model_indels) {

            double PAR_model_pip_lambda = ApplicationTools::getDoubleParameter("model_pip_lambda", jatiapp.getParams(), 0.1, "", true, 1);
            double PAR_model_pip_mu = ApplicationTools::getDoubleParameter("model_pip_mu", jatiapp.getParams(), 0.1, "", true, 1);

            lambda = PAR_model_pip_lambda;
            mu = PAR_model_pip_mu;

            submodel = new PIP_Nuc(alpha, lambda, mu, submodel);

            // Fill Q matrix
            Q = MatrixBppUtils::Matrix2Eigen(submodel->getGenerator());
            pi = MatrixBppUtils::Vector2Eigen(submodel->getFrequencies());

        }
        VLOG(1) << "[Substitution model] Number of states: " << (int) submodel->getNumberOfStates();

        bpp::StdStr s1;
        bpp::PhylogeneticsApplicationTools::printParameters(submodel, s1, 1, true);
        LOG(INFO) << s1.str();

        //------------------------------------------------------------------------------------------------------------------
        // COMPUTE ALIGNMENT USING PROGRESSIVE-PIP
        progressivePIP::ProgressivePIPResult MSA;
        auto likelihood = new tshlib::Likelihood;
        std::vector<tshlib::VirtualNode *> fullTraversalNodes;
        if (PAR_alignment) {

            likelihood->Init(utree, pi, Q, mu, lambda);

            // Traverse the tree in post-order filling a list of node ready for traversal
            likelihood->compileNodeList_postorder(fullTraversalNodes, utree->rootnode);

            // Set survival probability to each node in the list
            likelihood->setAllIotas(fullTraversalNodes);

            // Set deletion probability to each node in the list
            likelihood->setAllBetas(fullTraversalNodes);

            // set probability matrix -- exponential of substitution rates
            likelihood->computePr(fullTraversalNodes, alpha->getSize());


            VLOG(1) << "[ProPIP] starting MSA inference...";

            VirtualNode *root = utree->rootnode;

            MSA = progressivePIP::compute_DP3D_PIP_tree_cross(root, tree, &tm, pi, lambda, mu, sequences, alpha, 1.0, false);
            //sites = new bpp::VectorSiteContainer(*sequences);

            sequences = new bpp::VectorSequenceContainer(alpha);

            for (int i = 0; i < MSA.MSAs.size(); i++) {
                //auto sequence = new bpp::BasicSequence(MSA.MSAs.at(i).first,MSA.MSAs.at(i).second,alpha);
                //sequences->setSequence(MSA.MSAs.at(i).first,*sequence,true);
                sequences->addSequence(*(new bpp::BasicSequence(MSA.MSAs.at(i).first, MSA.MSAs.at(i).second, alpha)), true);
            }

            sites = new bpp::VectorSiteContainer(*sequences);

            delete sequences;

            std::string PAR_output_file_msa = ApplicationTools::getAFilePath("output_file_msa", jatiapp.getParams(), false, false, "", true, "", 1);
            std::string PAR_output_file_lk = ApplicationTools::getAFilePath("output_file_lk", jatiapp.getParams(), false, false, "", true, "", 1);

            bpp::Fasta seqWriter;
            seqWriter.writeAlignment(PAR_output_file_msa, *sites, true);

            std::ofstream lkFile;
            lkFile << std::setprecision(18);
            lkFile.open(PAR_output_file_lk);
            lkFile << MSA.score;
            lkFile.close();

            VLOG(1) << "LK ProPIP" << MSA.score;

            VLOG(1) << "[ProPIP] ...done";
            exit(0);

        }

        //------------------------------------------------------------------------------------------------------------------
        // Initialization likelihood functions

        double logLK = 0;
        //auto likelihood = new Likelihood();
        //std::vector<VirtualNode *> fullTraversalNodes;
        if (!PAR_alignment) {
            //tree = UtreeBppUtils::convertTree_u2b(utree);
            if (!PAR_model_indels) {

                transmodel = bpp::PhylogeneticsApplicationTools::getTransitionModel(alpha, gCode.get(), sites, parmap, "", true, false, 0);
                tl = new bpp::RHomogeneousTreeLikelihood(*tree, *sites, transmodel, rDist, false, false, false);

            } else {

                unique_ptr<TransitionModel> test;
                test.reset(submodel);
                transmodel = test.release();

                tl = new bpp::RHomogeneousTreeLikelihood_PIP(*tree, *sites, transmodel, rDist, &tm, false, false, false);
            }

            VLOG(1) << "[Transition model] Number of states: " << (int) transmodel->getNumberOfStates();

            tl->initialize();
            logLK = tl->getValue();

            LOG(INFO) << "[Tree likelihood] -- full traversal -- (on model " << submodel->getName() << ") = " << -logLK << " \t[BPP METHODS]";
        }

        //----------------------------------------------
        // Remove the root
        utree->removeVirtualRootNode();

        //------------------------------------------------------------------------------------------------------------------
        // DEFINE, APPLY & REVERT TREE REARRANGEMENTS
        // Get all the nodes between the radius boundaries and for each of them build the move list

        std::chrono::high_resolution_clock::time_point t1 = std::chrono::high_resolution_clock::now();
        unsigned long total_exec_moves = 0;
        int min_radius;
        int max_radius;

        if (PAR_optim_topology.find("full-search") != std::string::npos) {

            min_radius = 3;  // Minimum radius for an NNI move is 3 nodes
            max_radius = utree->getMaxNodeDistance(); // Full tree traversing from any node of the tree

        } else if (PAR_optim_topology.find("nni-search") != std::string::npos) {

            min_radius = 3;
            max_radius = 3;

        } else if (PAR_optim_topology.find("spr-search") != std::string::npos) {

            min_radius = 4;
            max_radius = utree->getMaxNodeDistance();

        } else {

            LOG(FATAL) << "Exiting program without tree search optimisation";
        }

        std::string PAR_lkmove = ApplicationTools::getStringParameter("lk_move", jatiapp.getParams(), "bothways", "", true, true);


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
                // Prepare the list of nodes involved in the move (Required here!)
                std::vector<VirtualNode *> listNodesWithinPath;
                listNodesWithinPath = utree->computePathBetweenNodes(pnode, qnode);
                listNodesWithinPath.push_back(utree->rootnode);

                // ------------------------------------
                // Apply the move
                rearrangmentList->applyMove(i);

                // ------------------------------------
                // Print root reachability from every node (includes rotations)
                //utree->_testReachingPseudoRoot();

                // ------------------------------------
                // Print tree on file
                //utree->saveTreeOnFile("../data/test.txt");

                // ------------------------------------
                bool isLKImproved = false;
                bool computeMoveLikelihood = true;

                // ------------------------------------
                // Add the root
                utree->addVirtualRootNode();

                // ------------------------------------
                if (computeMoveLikelihood && PAR_model_indels) {

                    // the dynamic_cast is necessary to access methods which belong to the class itself and not to the parent class
                    // in this case the class is the RHomogeneousTreeLikelihood_PIP, a derived class for PIP likelihood.
                    bpp::RHomogeneousTreeLikelihood_PIP *ttl = dynamic_cast<bpp::RHomogeneousTreeLikelihood_PIP *>(tl);
                    // we use a map to navigate between utree and bpp tree. The map is constant.
                    logLK = ttl->getLogLikelihood(listNodesWithinPath);

                    // ------------------------------------
                    // Store likelihood of the move
                    rearrangmentList->getMove(i)->move_lk = logLK;
                }

                if (!PAR_model_indels) {
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

                // ------------------------------------
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

                // ------------------------------------
                computeMoveLikelihood = false;

                // ------------------------------------
                if (PAR_model_indels) {
                    if (PAR_lkmove.find("bothways") != std::string::npos) {

                        // the dynamic_cast is necessary to access methods which belong to the class itself and not to the parent class
                        // in this case the class is the RHomogeneousTreeLikelihood_PIP, a derived class for PIP likelihood.
                        bpp::RHomogeneousTreeLikelihood_PIP *ttl = dynamic_cast<bpp::RHomogeneousTreeLikelihood_PIP *>(tl);
                        // we use a map to navigate between utree and bpp tree. The map is constant.

                        logLK = ttl->getLogLikelihood(listNodesWithinPath);
                        //VLOG(2) << "BPP::LogLK move revert " <<


                        // ------------------------------------
                        // Store likelihood of the move
                        rearrangmentList->getMove(i)->move_lk = logLK;


                    } else {
                        /*
                        likelihood->restoreLikelihoodComponents();
                        // ------------------------------------
                        // Compute the list of nodes for a full traversal
                        fullTraversalNodes.clear();
                        likelihood->compileNodeList_postorder(fullTraversalNodes, utree->rootnode);
                        //likelihood->setInsertionHistories(allnodes_postorder,*alignment);
                        logLK = LKFunc::LKRearrangment(*likelihood, fullTraversalNodes, *alignment);
                         */
                    }

                }
                if (!PAR_model_indels) {
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
        jatiapp.done();
        exit(0);

    } catch (exception &e) {
        cout << e.what() << endl;
        exit(1);
    }

}