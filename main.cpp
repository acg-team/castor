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
#include "TSHTopologySearch.hpp"
#include "pPIP.hpp"

using namespace tshlib;

int main(int argc, char *argv[]) {

    FLAGS_alsologtostderr = true;
    gflags::SetUsageMessage("some usage message");
    gflags::SetVersionString(software::build);
    google::InitGoogleLogging(software::name.c_str());
    //gflags::ParseCommandLineFlags(&argc, &argv, true);

#pragma omp master

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
        std::string PAR_Alphabet = ApplicationTools::getStringParameter("alphabet", jatiapp.getParams(), "DNA", "", true, true);
        std::string PAR_input_sequences = ApplicationTools::getAFilePath("input_sequences", jatiapp.getParams(), true, true, "", false, "", 1);
        bool PAR_alignment = ApplicationTools::getBooleanParameter("alignment", jatiapp.getParams(), false);
        std::string PAR_input_tree = ApplicationTools::getAFilePath("input_tree", jatiapp.getParams(), false, true, "", false, "", 1);
        std::string PAR_model_substitution = ApplicationTools::getStringParameter("model_substitution", jatiapp.getParams(), "JC69", "", true, true);
        bool PAR_model_setfreqsfromdata = ApplicationTools::getBooleanParameter("model_setfreqsfromdata", jatiapp.getParams(), false);
        bool PAR_model_indels = ApplicationTools::getBooleanParameter("model_indels", jatiapp.getParams(), false);
        std::string PAR_optim_topology_algorithm = ApplicationTools::getStringParameter("optim_topology_algorithm", jatiapp.getParams(), "no-search", "", true, true);
        bool PAR_profile_ppip = ApplicationTools::getBooleanParameter("profile_ppip", jatiapp.getParams(), false);
        std::string PAR_output_file_msa = ApplicationTools::getAFilePath("output_file_msa", jatiapp.getParams(), false, false, "", true, "", 1);
        std::string PAR_output_file_tree = ApplicationTools::getAFilePath("output_file_tree", jatiapp.getParams(), false, false, "", true, "", 1);

        /* ***************************************************
         * Standard workflow
         * [INPUT]
         * 1. tree + alignment => (1.1) just parse everything
         * 2. alignment  => (2.1) parse alignment => (2.2) generate tree using bioNJ
         * 3. sequences  => (3.1) parse sequences => (3.2) generate tree using bioNJ => (3.3) perform alignment
         * (it should not be supported in production)
         * 4. sequences + tree => (4.1) parse sequence => (4.2) parse tree => (4.3) perform alignment
         */


        bpp::Fasta seqReader;
        bpp::SequenceContainer *sequences;
        bpp::SiteContainer *sites;
        bpp::Alphabet *alpha;

        // The alphabet object should be set according to the correct alphabet

        if (PAR_model_indels) {
            if (PAR_Alphabet.find("DNA") != std::string::npos) {
                alpha = new bpp::DNA_EXTENDED();
            } else if (PAR_Alphabet.find("Protein") != std::string::npos) {
                alpha = new bpp::ProteicAlphabet_Extended();
            }
        } else {

            alpha = new bpp::DNA();
        }

        // Get alphabet from parameters
        bpp::Alphabet *alphabet = bpp::SequenceApplicationTools::getAlphabet(jatiapp.getParams(), "", false, false);
        unique_ptr<GeneticCode> gCode;
        bpp::CodonAlphabet *codonAlphabet = dynamic_cast<bpp::CodonAlphabet *>(alphabet);
        if (codonAlphabet) {
            std::string codeDesc = ApplicationTools::getStringParameter("genetic_code", jatiapp.getParams(), "Standard", "", true, true);
            //ApplicationTools::displayResult("Genetic Code", codeDesc);
            gCode.reset(bpp::SequenceApplicationTools::getGeneticCode(codonAlphabet->getNucleicAlphabet(), codeDesc));
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
            LOG(INFO) << "[Input data parser] The input file contains " << num_leaves << " sequences";

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
                delete newickReader;

            } catch (bpp::Exception e) {

                LOG(FATAL) << "[Tree parser] Error when reading tree due to: " << e.message();
            }
        } else {


            bpp::DistanceMatrix *distances;

            if (!PAR_alignment) {
                // Compute bioNJ tree using the JC69 model
                map<std::string, std::string> parmap;
                parmap["model"] = "JC69";

                bpp::SequenceContainer *sequences_bioNJ = seqReader.readAlignment(PAR_input_sequences, alphabet);
                seqReader.readSequences(PAR_input_sequences, alphabet);
                bpp::SiteContainer *sites_bioNJ = new bpp::VectorSiteContainer(*sequences_bioNJ);

                bpp::SubstitutionModel *submodel_bioNJ = bpp::PhylogeneticsApplicationTools::getSubstitutionModel(alphabet, gCode.get(), sites_bioNJ, parmap, "", true, false, 0);
                bpp::DiscreteDistribution *rDist = new bpp::ConstantRateDistribution();
                bpp::SiteContainerTools::changeGapsToUnknownCharacters(*sites_bioNJ);
                bpp::DistanceEstimation distanceMethod(submodel_bioNJ, rDist, sites_bioNJ);
                distances = distanceMethod.getMatrix();
                delete sequences_bioNJ;
                LOG(INFO) << "[BioNJ Pairwise distance matrix] The pairwise distance matrix is computed using JC69";

            } else {

                std::string PAR_distance_matrix = ApplicationTools::getAFilePath("distance_matrix", jatiapp.getParams(), false, true, "", false, "", 0);
                distances = InputUtils::parseDistanceMatrix(PAR_distance_matrix);

                LOG(INFO) << "[BioNJ Pairwise distance matrix] The pairwise distance matrix is computed using LZ compression ";

            }

            bpp::BioNJ bionj(*distances, true, true, false);
            tree = bionj.getTree();

        }
        // Rename internal nodes with standard Vxx * where xx is a progressive number
        tree->setNodeName(tree->getRootId(), "root");
        UtreeBppUtils::renameInternalNodes(tree);


        Newick treeWriter;
        TreeTemplate<Node> ttree(*tree);
        std::ostringstream oss;
        treeWriter.write(ttree, oss);

        LOG(INFO) << "[Initial Tree Topology] " << oss.str();

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

            if (!PAR_alignment) {
                if (PAR_model_setfreqsfromdata) {
                    submodel->setFreqFromData(*sites);
                }
            }

            rDist = new bpp::ConstantRateDistribution();
            if (!PAR_alignment) { bpp::SiteContainerTools::changeGapsToUnknownCharacters(*sites); }

        }

        // Extend the substitution model with PIP
        if (PAR_model_indels) {

            double PAR_model_pip_lambda = ApplicationTools::getDoubleParameter("model_pip_lambda", jatiapp.getParams(), 0.1, "", true, 1);
            double PAR_model_pip_mu = ApplicationTools::getDoubleParameter("model_pip_mu", jatiapp.getParams(), 0.1, "", true, 1);

            lambda = PAR_model_pip_lambda;
            mu = PAR_model_pip_mu;

            if (PAR_Alphabet.find("DNA") != std::string::npos) {
                submodel = new PIP_Nuc(dynamic_cast<NucleicAlphabet *>(alpha), lambda, mu, submodel);
            } else if (PAR_Alphabet.find("Protein") != std::string::npos) {
                submodel = new PIP_AA(dynamic_cast<ProteicAlphabet *>(alpha), lambda, mu, submodel);
            }

            //submodel = new PIP_Nuc(dynamic_cast<NucleicAlphabet *>(alpha), lambda, mu, submodel);

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
            //likelihood->setAllIotas(fullTraversalNodes);

            // Set deletion probability to each node in the list
            //likelihood->setAllBetas(fullTraversalNodes);

            // set probability matrix -- exponential of substitution rates
            likelihood->computePr(fullTraversalNodes, alpha->getSize());


            LOG(INFO) << "[Alignment sequences] Starting MSA inference using Pro-PIP...";

            /*
            VirtualNode *root = utree->rootnode;

            MSA = progressivePIP::compute_DP3D_PIP_tree_cross(root, tree, &tm, pi, lambda, mu, sequences, alpha, 1.0, false);
             */

            //********************************************************************************
            //********************************************************************************
            //********************************************************************************
            //********************************************************************************

            double tau;

            auto progressivePIP=new bpp::pPIP(alphabet);

            progressivePIP->init(tree, &tm, fullTraversalNodes, submodel->getFrequencies(),lambda, mu);

            progressivePIP->PIPAligner(&tm,fullTraversalNodes, sequences, 1.0, true);

            std::cout<<"PIPAligner done...\n";

            //********************************************************************************
            //********************************************************************************
            //********************************************************************************
            //********************************************************************************

            /*
            sequences = new bpp::VectorSequenceContainer(alpha);

            for (int i = 0; i < MSA.MSAs.size(); i++) {
                //auto sequence = new bpp::BasicSequence(MSA.MSAs.at(i).first,MSA.MSAs.at(i).second,alpha);
                //sequences->setSequence(MSA.MSAs.at(i).first,*sequence,true);
                sequences->addSequence(*(new bpp::BasicSequence(MSA.MSAs.at(i).first, MSA.MSAs.at(i).second, alpha)), true);
            }

            sites = new bpp::VectorSiteContainer(*sequences);

            //delete sequences;

            if (PAR_profile_ppip) {
                std::string PAR_output_file_lk = ApplicationTools::getAFilePath("output_file_lk", jatiapp.getParams(), false, false, "", true, "", 1);

                //bpp::Fasta seqWriter;
                //seqWriter.writeAlignment(PAR_output_file_msa, *sites, true);

                std::ofstream lkFile;
                lkFile << std::setprecision(18);
                lkFile.open(PAR_output_file_lk);
                lkFile << MSA.score;
                lkFile.close();
            }

            LOG(INFO) << "[Alignment sequences] MSA inference using Pro-PIP terminated successfully!";
            LOG(INFO) << "[Alignment sequences] Alignment has likelihood: " << MSA.score;

            */

        }




        //--------------------------------------------------------------------------------------------------------------
        // best tree from MSA marginalization
        if(false){
            auto treesearch = new tshlib::TreeSearch;
            Utree *best_tree_from_MSA=progressivePIP::marginalizationOverMSAs(treesearch,alpha,pi,lambda, mu, sequences, tm);
        }
        //--------------------------------------------------------------------------------------------------------------




        //--------------------------------------------------------------------------------------------------------------
        // Initialization likelihood functions

        double logLK = 0;
        //auto likelihood = new Likelihood();
        //std::vector<VirtualNode *> fullTraversalNodes;
        //if (!PAR_alignment) {
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

        LOG(INFO) << "[Intial tree likelihood] (on model " << submodel->getName() << ") = " << -logLK;
        //}

        //----------------------------------------------
        // Optimise parameters automatically
        //if (!PAR_model_indels) {
            tl = dynamic_cast<AbstractHomogeneousTreeLikelihood *>(PhylogeneticsApplicationTools::optimizeParameters(tl, tl->getParameters(),
                                                                                                                     jatiapp.getParams(), "", true, true, 0));
        //}
        //----------------------------------------------
        // Remove the root
        utree->removeVirtualRootNode();

        //------------------------------------------------------------------------------------------------------------------
        // DEFINE, APPLY & REVERT TREE REARRANGEMENTS
        // Get all the nodes between the radius boundaries and for each of them build the move list

        tshlib::TreeRearrangmentOperations treesearch_operations;
        tshlib::TreeSearchHeuristics treesearch_heuristics;

        if (PAR_optim_topology_algorithm.find("greedy") != std::string::npos) {
            treesearch_heuristics = tshlib::TreeSearchHeuristics::greedy;
        } else if (PAR_optim_topology_algorithm.find("hillclimbing") != std::string::npos) {
            treesearch_heuristics = tshlib::TreeSearchHeuristics::hillclimbing;
        } else if (PAR_optim_topology_algorithm.find("no-search") != std::string::npos) {
            treesearch_heuristics = tshlib::TreeSearchHeuristics::nosearch;
        }

        if (treesearch_heuristics != tshlib::TreeSearchHeuristics::nosearch) {
            std::string PAR_lkmove = ApplicationTools::getStringParameter("lk_move", jatiapp.getParams(), "bothways", "", true, true);
            std::string PAR_optim_topology_operations = ApplicationTools::getStringParameter("optim_topology_operations", jatiapp.getParams(), "best-search", "", true, true);
            int PAR_optim_topology_maxcycles = ApplicationTools::getIntParameter("optim_topology_maxcycles", jatiapp.getParams(), 1, "", true, 0);
            int PAR_optim_topology_hillclimbing_startnodes = ApplicationTools::getIntParameter("optim_topology_numnodes", jatiapp.getParams(), 1, "", true, 0);



            if (PAR_optim_topology_operations.find("best-search") != std::string::npos) {
                treesearch_operations = tshlib::TreeRearrangmentOperations::classic_Mixed;
            } else if (PAR_optim_topology_operations.find("nni-search") != std::string::npos) {
                treesearch_operations = tshlib::TreeRearrangmentOperations::classic_NNI;
            } else if (PAR_optim_topology_operations.find("spr-search") != std::string::npos) {
                treesearch_operations = tshlib::TreeRearrangmentOperations::classic_SPR;
            } else {
                LOG(FATAL) << "Exiting program without tree search optimisation";
            }

            auto treesearch = new tshlib::TreeSearch;
            treesearch->setTreeSearchStrategy(treesearch_heuristics, treesearch_operations);
            treesearch->setInitialLikelihoodValue(-logLK);
            treesearch->setScoringMethod(PAR_lkmove);
            treesearch->setStartingNodes(PAR_optim_topology_hillclimbing_startnodes);
            treesearch->setStopCondition(tshlib::TreeSearchStopCondition::iterations, (double) PAR_optim_topology_maxcycles);
            if (PAR_model_indels) {
                treesearch->setModelIndels(true);
                treesearch->setLikelihoodFunc(tl);
            } else {

            }
            logLK = treesearch->performTreeSearch(utree);

        }


        LOG(INFO) << "[TSH Best Topology] (" << -logLK << ") " << utree->printTreeNewick(false);

        if (PAR_output_file_tree.find("none") == std::string::npos) {
            LOG(INFO) << "[Output tree]\t The final topology can be found in " << PAR_output_file_tree;
            std::ofstream file;
            file.open(PAR_output_file_tree);
            file << utree->printTreeNewick(false);
            file.close();
        }
        if (PAR_output_file_msa.find("none") == std::string::npos) {
            LOG(INFO) << "[Output alignment]\t The final alignment can be found in " << PAR_output_file_msa;
            bpp::Fasta seqWriter;
            seqWriter.writeAlignment(PAR_output_file_msa, *sites, true);
        }

        delete sequences;

        jatiapp.done();
        exit(0);

    } catch (exception &e) {
        cout << e.what() << endl;
        exit(1);
    }

}