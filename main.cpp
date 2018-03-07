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
#include <Bpp/Text/KeyvalTools.h>

#include "Version.hpp"
#include "Utilities.hpp"
#include "PIP.hpp"
#include "CommandLineFlags.hpp"
#include "ExtendedAlphabet.hpp"
#include "RHomogeneousTreeLikelihood_PIP.hpp"
#include "TSHHomogeneousTreeLikelihood.hpp"

#include "progressivePIP.hpp"
#include "JATIApplication.hpp"
#include "TSHTopologySearch.hpp"
#include "pPIP.hpp"
#include "Optimizators.hpp"

using namespace tshlib;

int main(int argc, char *argv[]) {

    FLAGS_alsologtostderr = true;
    int OMP_max_avail_threads = 1;
    gflags::SetVersionString(software::build);
    google::InitGoogleLogging(software::name.c_str());
    //gflags::ParseCommandLineFlags(&argc, &argv, true);

    try {
        bpp::JATIApplication jatiapp(argc, argv, software::desc);

        if (argc < 2) {
            jatiapp.help();
            return 0;
        } else {
            jatiapp.banner();
            jatiapp.startTimer();
        };

        bool PAR_alignment = ApplicationTools::getBooleanParameter("alignment", jatiapp.getParams(), false);
        //std::string PAR_input_tree = ApplicationTools::getAFilePath("input_tree", jatiapp.getParams(), false, true, "", false, "", 1);

        std::string PAR_optim_topology_algorithm = ApplicationTools::getStringParameter("optim_topology_algorithm", jatiapp.getParams(), "no-search", "", true, true);
        bool PAR_profile_ppip = ApplicationTools::getBooleanParameter("profile_ppip", jatiapp.getParams(), false);
        std::string PAR_output_file_msa = ApplicationTools::getAFilePath("output_file_msa", jatiapp.getParams(), false, false, "", true, "", 1);
        std::string PAR_output_file_tree = ApplicationTools::getAFilePath("output_file_tree", jatiapp.getParams(), false, false, "", true, "", 1);


        if (OMPENABLED) OMP_max_avail_threads = omp_get_max_threads();
        int PAR_execution_numthreads = ApplicationTools::getIntParameter("exec_numthreads", jatiapp.getParams(), OMP_max_avail_threads, "", true, 0);


        /* ***************************************************
         * Standard workflow
         * [INPUT]
         * 1. tree + alignment => (1.1) just parse everything
         * 2. alignment  => (2.1) parse alignment => (2.2) generate tree using bioNJ
         * 3. sequences  => (3.1) parse sequences => (3.2) generate tree using bioNJ => (3.3) perform alignment
         * (it should not be supported in production)
         * 4. sequences + tree => (4.1) parse sequence => (4.2) parse tree => (4.3) perform alignment
         */

        std::string PAR_model_substitution = ApplicationTools::getStringParameter("model", jatiapp.getParams(), "JC69", "", true, true);
        bool PAR_model_setfreqsfromdata = ApplicationTools::getBooleanParameter("model_setfreqsfromdata", jatiapp.getParams(), false);

        // Split model string description and test if PIP is required
        std::string modelStringName;
        std::map<std::string, std::string> modelMap;
        KeyvalTools::parseProcedure(PAR_model_substitution, modelStringName, modelMap);
        bool PAR_model_indels = (modelStringName == "PIP") ? true : false;
        modelMap["model"];
        //////////////////////////////////////////////
        // ALPHABET
        // The alphabet object contains the not-extended alphabet as requested by the user,
        // while alpha contains the extended version of the same alphabet.
        std::string PAR_Alphabet = ApplicationTools::getStringParameter("alphabet", jatiapp.getParams(), "DNA", "", true, true);
        bpp::Alphabet *alpha;

        if (PAR_model_indels) {
            if (PAR_Alphabet.find("DNA") != std::string::npos) {
                alpha = new bpp::DNA_EXTENDED();
            } else if (PAR_Alphabet.find("Protein") != std::string::npos) {
                alpha = new bpp::ProteicAlphabet_Extended();
            }
        } else {
            if (PAR_Alphabet.find("DNA") != std::string::npos) {
                alpha = new bpp::DNA();
            } else if (PAR_Alphabet.find("Protein") != std::string::npos) {
                alpha = new bpp::ProteicAlphabet();
            }
        }

        bpp::Alphabet *alphabet = bpp::SequenceApplicationTools::getAlphabet(jatiapp.getParams(), "", false, false);
        unique_ptr<GeneticCode> gCode;
        bpp::CodonAlphabet *codonAlphabet = dynamic_cast<bpp::CodonAlphabet *>(alphabet);
        if (codonAlphabet) {
            std::string codeDesc = ApplicationTools::getStringParameter("genetic_code", jatiapp.getParams(), "Standard", "", true, true);
            //ApplicationTools::displayResult("Genetic Code", codeDesc);
            gCode.reset(bpp::SequenceApplicationTools::getGeneticCode(codonAlphabet->getNucleicAlphabet(), codeDesc));
        }

        //////////////////////////////////////////////
        // DATA
        std::string PAR_input_sequences = ApplicationTools::getAFilePath("input.sequence.file", jatiapp.getParams(), true, true, "", false, "", 1);

        bpp::SequenceContainer *sequences = nullptr;
        bpp::SiteContainer *sites = nullptr;

        try {

            LOG(INFO) << "[Input data parser] " << PAR_input_sequences;

            if (PAR_alignment) {
                bpp::Fasta seqReader;
                // If the user requires the computation of an alignment, then the input file is made of unaligned sequences
                sequences = seqReader.readSequences(PAR_input_sequences, alpha);

            } else {

                //sequences = seqReader.readAlignment(PAR_input_sequences, alpha);
                //seqReader.readSequences(PAR_input_sequences, alpha);
                //std::vector<std::string> seqNames = sequences->getSequencesNames();
                //sites = new bpp::VectorSiteContainer(*sequences);

                VectorSiteContainer *allSites = SequenceApplicationTools::getSiteContainer(alpha, jatiapp.getParams());
                sites = SequenceApplicationTools::getSitesToAnalyse(*allSites, jatiapp.getParams(), "", true, !PAR_model_indels, true, 1);
                delete allSites;


            }

            //size_t num_leaves = sequences->getNumberOfSequences();
            LOG(INFO) << "[Input data parser] Number of sequences: " << sites->getNumberOfSequences();
            LOG(INFO) << "[Input data parser] Number of sites: " << sites->getNumberOfSites();

        } catch (bpp::Exception e) {
            LOG(FATAL) << "[Input data parser] Error when reading sequence file due to: " << e.message();
        }

        /////////////////////////////////////////
        // TREE
        bpp::Tree *tree = nullptr;

        string initTreeOpt = ApplicationTools::getStringParameter("init.tree", jatiapp.getParams(), "user", "", false, 1);
        ApplicationTools::displayResult("Input tree", initTreeOpt);

        if (initTreeOpt == "user") {

            tree = PhylogeneticsApplicationTools::getTree(jatiapp.getParams());
            LOG(INFO) << "[Input tree parser] Number of leaves" << tree->getNumberOfLeaves();

        } else if (initTreeOpt == "random") {

            vector<string> names = sites->getSequencesNames();
            tree = TreeTemplateTools::getRandomTree(names);
            tree->setBranchLengths(1.);

        } else if (initTreeOpt == "distance") {

            bpp::DistanceMatrix *distances;

            if (!PAR_alignment) {
                // Compute bioNJ tree using the JC69 model
                map<std::string, std::string> parmap;
                parmap["model"] = "JC69";
                bpp::Fasta seqReader;
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

                std::string PAR_distance_matrix = ApplicationTools::getAFilePath("init.distance.matrix.file", jatiapp.getParams(), false, true, "", false, "", 0);
                distances = InputUtils::parseDistanceMatrix(PAR_distance_matrix);

                LOG(INFO) << "[BioNJ Pairwise distance matrix] The pairwise distance matrix is computed using LZ compression ";

            }

            bpp::BioNJ bionj(*distances, true, true, false);
            tree = bionj.getTree();

        } else throw Exception("Unknown init tree method.");

        // Rename internal nodes with standard Vxx * where xx is a progressive number
        tree->setNodeName(tree->getRootId(), "root");
        UtreeBppUtils::renameInternalNodes(tree);


        // Try to write the current tree to file. This will be overwritten by the optimized tree,
        // but allow to check file existence before running optimization!
        PhylogeneticsApplicationTools::writeTree(*tree, jatiapp.getParams());

        // Setting branch lengths?
        string initBrLenMethod = ApplicationTools::getStringParameter("init.brlen.method", jatiapp.getParams(), "Input", "", true, 1);
        string cmdName;
        map<string, string> cmdArgs;
        KeyvalTools::parseProcedure(initBrLenMethod, cmdName, cmdArgs);
        if (cmdName == "Input") {
            // Is the root has to be moved to the midpoint position along the branch that contains it ? If no, do nothing!
            bool midPointRootBrLengths = ApplicationTools::getBooleanParameter("midpoint_root_branch", cmdArgs, false, "", true, 2);
            if (midPointRootBrLengths)
                TreeTools::constrainedMidPointRooting(*tree);
        } else if (cmdName == "Equal") {
            double value = ApplicationTools::getDoubleParameter("value", cmdArgs, 0.1, "", true, 2);
            if (value <= 0)
                throw Exception("Value for branch length must be superior to 0");
            ApplicationTools::displayResult("Branch lengths set to", value);
            tree->setBranchLengths(value);
        } else if (cmdName == "Clock") {
            TreeTools::convertToClockTree(*tree, tree->getRootId(), true);
        } else if (cmdName == "Grafen") {
            string grafenHeight = ApplicationTools::getStringParameter("height", cmdArgs, "input", "", true, 2);
            double h;
            if (grafenHeight == "input") {
                h = TreeTools::getHeight(*tree, tree->getRootId());
            } else {
                h = TextTools::toDouble(grafenHeight);
                if (h <= 0) throw Exception("Height must be positive in Grafen's method.");
            }
            ApplicationTools::displayResult("Total height", TextTools::toString(h));

            double rho = ApplicationTools::getDoubleParameter("rho", cmdArgs, 1., "", true, 2);
            ApplicationTools::displayResult("Grafen's rho", rho);
            TreeTools::computeBranchLengthsGrafen(*tree, rho);
            double nh = TreeTools::getHeight(*tree, tree->getRootId());
            tree->scaleTree(h / nh);
        } else throw Exception("Method '" + initBrLenMethod + "' unknown for computing branch lengths.");
        ApplicationTools::displayResult("Branch lengths", cmdName);


        /*
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



        }
        */
        LOG(INFO) << "[Initial Tree Topology] " << OutputUtils::tree2string(tree);

        // Convert the bpp into utree for tree-search engine
        auto utree = new Utree();
        UtreeBppUtils::treemap tm;
        UtreeBppUtils::convertTree_b2u(tree, utree, tm);
        UtreeBppUtils::associateNode2Alignment(sites, utree);

        DLOG(INFO) << "Bidirectional map size: " << tm.size();
        DLOG(INFO) << "[Initial Utree Topology] " << utree->printTreeNewick(true);

        utree->addVirtualRootNode();
        // Once the tree has the root, then map it as well
        tm.insert(UtreeBppUtils::nodeassoc(tree->getRootId(), utree->rootnode));


        /////////////////////////
        // MODEL  & LIKELIHOOD

        bpp::SubstitutionModel *smodel = nullptr;
        bpp::TransitionModel *model = nullptr;


        Eigen::MatrixXd Q;
        Eigen::VectorXd pi;
        double lambda;
        double mu;

        if (PAR_model_indels) {
            smodel = bpp::PhylogeneticsApplicationTools::getSubstitutionModel(alpha, gCode.get(), sites, modelMap, "", true, false, 0);

            // Extend the substitution model with PIP
            lambda = (modelMap.find("lambda") == modelMap.end()) ? 0.1 : std::stod(modelMap["lambda"]);
            mu = (modelMap.find("mu") == modelMap.end()) ? 0.2 : std::stod(modelMap["mu"]);
            // Instatiate the corrisponding PIP model given the alphabet
            if (PAR_Alphabet.find("DNA") != std::string::npos) {
                smodel = new PIP_Nuc(dynamic_cast<NucleicAlphabet *>(alpha), lambda, mu, smodel);
            } else if (PAR_Alphabet.find("Protein") != std::string::npos) {
                smodel = new PIP_AA(dynamic_cast<ProteicAlphabet *>(alpha), lambda, mu, smodel);
            } else if (PAR_Alphabet.find("Codon") != std::string::npos) {
                LOG(FATAL) << "Pip model is not implemented for Codon alphabets! :(";
            }
            // Fill Q matrix
            Q = MatrixBppUtils::Matrix2Eigen(smodel->getGenerator());
            pi = MatrixBppUtils::Vector2Eigen(smodel->getFrequencies());

        } else {
            // if the alphabet is not extended, then the gap character is not supported
            if (!PAR_alignment) bpp::SiteContainerTools::changeGapsToUnknownCharacters(*sites);
            smodel = bpp::PhylogeneticsApplicationTools::getSubstitutionModel(alpha, gCode.get(), sites, jatiapp.getParams(), "", true, false, 0);
        }
        //if (!PAR_alignment) { bpp::SiteContainerTools::changeGapsToUnknownCharacters(*sites); }

        DLOG(INFO) << "[Substitution model] Number of states: " << (int) smodel->getNumberOfStates();

        bpp::StdStr s1;
        bpp::PhylogeneticsApplicationTools::printParameters(smodel, s1, 1, true);
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


            LOG(INFO) << "[Alignment sequences] Starting MSA inference using PIPAligner...";

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

            progressivePIP->init(tree, &tm, fullTraversalNodes, smodel->getFrequencies(), lambda, mu);

            progressivePIP->PIPAligner(&tm,fullTraversalNodes, sequences, 1.0, true);

            LOG(INFO) << "[Alignment sequences] Alignment completed succesfully using PIPAligner.";

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

        ///// homogeneous modeling
        // Initialization likelihood functions
        bpp::DiscreteDistribution *rDist = nullptr;
        bpp::AbstractHomogeneousTreeLikelihood *tl;

        // Get transition model from substitution model
        if (!PAR_model_indels) {
            model = bpp::PhylogeneticsApplicationTools::getTransitionModel(alpha, gCode.get(), sites, jatiapp.getParams(), "", true, false, 0);
        } else {
            unique_ptr<TransitionModel> test;
            test.reset(smodel);
            model = test.release();
        }

        // Among site rate variation (ASVR)
        if (model->getNumberOfStates() >= 2 * model->getAlphabet()->getSize()) {
            // Markov-modulated Markov model!
            rDist = new ConstantRateDistribution();
        } else {
            rDist = PhylogeneticsApplicationTools::getRateDistribution(jatiapp.getParams());
        }

        // Initialise likelihood functions
        if (!PAR_model_indels) {
            tl = new bpp::RHomogeneousTreeLikelihood(*tree, *sites, model, rDist, false, false, false);
        } else {
            tl = new bpp::RHomogeneousTreeLikelihood_PIP(*tree, *sites, model, rDist, &tm, false, false, false);
        }

        tl->initialize();

        //Listing parameters
        string paramNameFile = ApplicationTools::getAFilePath("output.parameter_names.file", jatiapp.getParams(), false, false, "", true, "none", 1);
        if (paramNameFile != "none") {
            ApplicationTools::displayResult("List parameters to", paramNameFile);
            ofstream pnfile(paramNameFile.c_str(), ios::out);
            ParameterList pl = tl->getParameters();
            for (size_t i = 0; i < pl.size(); ++i) {
                pnfile << pl[i].getName() << endl;
            }
            pnfile.close();
            cout << "BppML's done." << endl;
            exit(0);
        }

        //Check initial likelihood:
        double logL = tl->getValue();
        if (std::isinf(logL)) {
            // This may be due to null branch lengths, leading to null likelihood!
            ApplicationTools::displayWarning("!!! Warning!!! Initial likelihood is zero.");
            ApplicationTools::displayWarning("!!! This may be due to branch length == 0.");
            ApplicationTools::displayWarning("!!! All null branch lengths will be set to 0.000001.");
            ParameterList pl = tl->getBranchLengthsParameters();
            for (unsigned int i = 0; i < pl.size(); i++) {
                if (pl[i].getValue() < 0.000001) pl[i].setValue(0.000001);
            }
            tl->matchParametersValues(pl);
            logL = tl->getValue();
        }
        ApplicationTools::displayResult("Initial log likelihood", TextTools::toString(-logL, 15));
        if (std::isinf(logL)) {
            ApplicationTools::displayError("!!! Unexpected initial likelihood == 0.");
            if (codonAlphabet) {
                bool f = false;
                size_t s;
                for (size_t i = 0; i < sites->getNumberOfSites(); i++) {
                    if (std::isinf(tl->getLogLikelihoodForASite(i))) {
                        const Site &site = sites->getSite(i);
                        s = site.size();
                        for (size_t j = 0; j < s; j++) {
                            if (gCode->isStop(site.getValue(j))) {
                                (*ApplicationTools::error << "Stop Codon at site " << site.getPosition() << " in sequence " << sites->getSequence(j).getName()).endLine();
                                f = true;
                            }
                        }
                    }
                }
                if (f)
                    exit(-1);
            }
            bool removeSaturated = ApplicationTools::getBooleanParameter("input.sequence.remove_saturated_sites", jatiapp.getParams(), false, "", true, 1);
            if (!removeSaturated) {
                ofstream debug("DEBUG_likelihoods.txt", ios::out);
                for (size_t i = 0; i < sites->getNumberOfSites(); i++) {
                    debug << "Position " << sites->getSite(i).getPosition() << " = " << tl->getLogLikelihoodForASite(i) << endl;
                }
                debug.close();
                ApplicationTools::displayError("!!! Site-specific likelihood have been written in file DEBUG_likelihoods.txt .");
                ApplicationTools::displayError("!!! 0 values (inf in log) may be due to computer overflow, particularily if datasets are big (>~500 sequences).");
                ApplicationTools::displayError("!!! You may want to try input.sequence.remove_saturated_sites = yes to ignore positions with likelihood 0.");
                exit(1);
            } else {
                ApplicationTools::displayBooleanResult("Saturated site removal enabled", true);
                for (size_t i = sites->getNumberOfSites(); i > 0; --i) {
                    if (std::isinf(tl->getLogLikelihoodForASite(i - 1))) {
                        ApplicationTools::displayResult("Ignore saturated site", sites->getSite(i - 1).getPosition());
                        sites->deleteSite(i - 1);
                    }
                }
                ApplicationTools::displayResult("Number of sites retained", sites->getNumberOfSites());
                tl->setData(*sites);
                tl->initialize();
                logL = tl->getValue();
                if (std::isinf(logL)) {
                    throw Exception("Likelihood is still 0 after saturated sites are removed! Looks like a bug...");
                }
                ApplicationTools::displayResult("Initial log likelihood", TextTools::toString(-logL, 15));
            }
        }


        //------------------------------------------------------------------------------------------------------------------
        // OPTIMISE PARAMETERS (numerical + topology) according to user parameters
        // Optimise parameters automatically following standard pipeline
        auto ntl = new bpp::TSHHomogeneousTreeLikelihood(tl, (*tl->getData()), (tl->getModel()), (tl->getRateDistribution()), utree, tm);
        tl = dynamic_cast<AbstractHomogeneousTreeLikelihood *>(Optimizators::optimizeParameters(ntl, ntl->getParameters(), jatiapp.getParams(), "", true, true, 0));
        //tl = dynamic_cast<AbstractHomogeneousTreeLikelihood *>(Optimizators::optimizeParameters(tl, tl->getParameters(), jatiapp.getParams(), "", true, true, 0));

        OutputUtils::printParametersLikelihood(tl);


        LOG(INFO) << "[Final likelihood] " << std::setprecision(12) << -tl->getValue();


        if (PAR_output_file_msa.find("none") == std::string::npos) {
            LOG(INFO) << "[Output alignment]\t The final alignment can be found in " << PAR_output_file_msa;
            bpp::Fasta seqWriter;
            seqWriter.writeAlignment(PAR_output_file_msa, *sites, true);
        }

        delete sequences;

        tree = new TreeTemplate<Node>(tl->getTree());
        PhylogeneticsApplicationTools::writeTree(*tree, jatiapp.getParams());

        // Write parameters to screen:
        ApplicationTools::displayResult("Log likelihood", TextTools::toString(-tl->getValue(), 15));
        ParameterList parameters = tl->getSubstitutionModelParameters();
        for (size_t i = 0; i < parameters.size(); i++) {
            ApplicationTools::displayResult(parameters[i].getName(), TextTools::toString(parameters[i].getValue()));
        }
        parameters = tl->getRateDistributionParameters();
        for (size_t i = 0; i < parameters.size(); i++) {
            ApplicationTools::displayResult(parameters[i].getName(), TextTools::toString(parameters[i].getValue()));
        }

        // Checking convergence:
        PhylogeneticsApplicationTools::checkEstimatedParameters(tl->getParameters());

        // Write parameters to file:
        string parametersFile = ApplicationTools::getAFilePath("output.estimates", jatiapp.getParams(), false, false, "none", 1);
        bool withAlias = ApplicationTools::getBooleanParameter("output.estimates.alias", jatiapp.getParams(), true, "", true, 0);

        ApplicationTools::displayResult("Output estimates to file", parametersFile);


        if (parametersFile != "none") {
            StlOutputStream out(new ofstream(parametersFile.c_str(), ios::out));


            out << "# Log likelihood = ";
            out.setPrecision(20) << (-tl->getValue());
            out.endLine();
            out << "# Number of sites = ";
            out.setPrecision(20) << sites->getNumberOfSites();
            out.endLine();
            out.endLine();
            out << "# Substitution model parameters:";
            out.endLine();

            smodel->matchParametersValues(tl->getParameters());
            PhylogeneticsApplicationTools::printParameters(smodel, out, 1, withAlias);

            out.endLine();
            (out << "# Rate distribution parameters:").endLine();
            rDist->matchParametersValues(tl->getParameters());
            PhylogeneticsApplicationTools::printParameters(rDist, out, withAlias);
        }

        // Compute support measures


        jatiapp.done();
        exit(0);

    } catch (exception &e) {
        cout << e.what() << endl;
        exit(1);
    }

}