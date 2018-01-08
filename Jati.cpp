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
 * @file Jati.cpp
 * @author Lorenzo Gatti
 * @author Massimo Maiolo
 * @date 05 01 2018
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

// From the STL:
#include <iostream>
#include <iomanip>

using namespace std;

// From bpp-core:
#include <Bpp/Numeric/Prob/DiscreteDistribution.h>
#include <Bpp/Numeric/Prob/ConstantDistribution.h>
#include <Bpp/Numeric/DataTable.h>
#include <Bpp/App/BppApplication.h>
#include <Bpp/Text/KeyvalTools.h>

// From bpp-seq:
#include <Bpp/Seq/Alphabet/Alphabet.h>
#include <Bpp/Seq/Container/VectorSiteContainer.h>
#include <Bpp/Seq/Container/SiteContainerTools.h>
#include <Bpp/Seq/SiteTools.h>
#include <Bpp/Seq/App/SequenceApplicationTools.h>

// From bpp-phyl:
#include <Bpp/Phyl/Tree.h>
#include <Bpp/Phyl/Likelihood/RHomogeneousMixedTreeLikelihood.h>
#include <Bpp/Phyl/Likelihood/DRHomogeneousMixedTreeLikelihood.h>
#include <Bpp/Phyl/Likelihood/RNonHomogeneousMixedTreeLikelihood.h>
#include <Bpp/Phyl/Likelihood/DRNonHomogeneousTreeLikelihood.h>
#include <Bpp/Phyl/Likelihood/RASTools.h>
#include <Bpp/Phyl/App/PhylogeneticsApplicationTools.h>
#include <Bpp/Phyl/OptimizationTools.h>
#include <Bpp/Phyl/Model/SubstitutionModelSetTools.h>
#include <Bpp/Phyl/Model/Protein/CoalaCore.h>
#include <Bpp/Phyl/Model/RateDistribution/ConstantRateDistribution.h>
#include <Bpp/Phyl/Model/FrequenciesSet/MvaFrequenciesSet.h>
#include <Bpp/Phyl/Io/Newick.h>

using namespace bpp;

// From this project

#include "Version.hpp"


/******************************************************************************/

void help() {
    (*ApplicationTools::message << "__________________________________________________________________________").endLine();
    (*ApplicationTools::message << "Jati parameter1_name=parameter1_value parameter2_name=parameter2_value").endLine();
    (*ApplicationTools::message << "      ... param=option_file").endLine();
    (*ApplicationTools::message).endLine();
    (*ApplicationTools::message << "  Refer to the Jati Program Suite Manual for a list of available options.").endLine();
    (*ApplicationTools::message << "__________________________________________________________________________").endLine();
}

int main(int args, char **argv) {
    ApplicationTools::displayMessage("**********************************************************************");
    ApplicationTools::displayMessage("* "+ softwarename + " *");
    ApplicationTools::displayMessage("* Authors: Lorenzo Gatti & Massimo Maiolo                            *");
    ApplicationTools::displayMessage("*--------------------------------------------------------------------*");
    ApplicationTools::displayMessage("* Inspired on Bio++ Maximum Likelihood Computation                   *");
    ApplicationTools::displayMessage("* by J. Dutheil, B. Boussau, L. Guéguen, M. Groussin                 *");
    ApplicationTools::displayMessage("**********************************************************************");


    if (args == 1) {
        help();
        return 0;
    }

    try {
        BppApplication jatiapp(args, argv, softwarename);
        jatiapp.startTimer();

        Alphabet *alphabet = SequenceApplicationTools::getAlphabet(jatiapp.getParams(), "", false);
        unique_ptr<GeneticCode> gCode;
        CodonAlphabet *codonAlphabet = dynamic_cast<CodonAlphabet *>(alphabet);
        if (codonAlphabet) {
            string codeDesc = ApplicationTools::getStringParameter("genetic_code", jatiapp.getParams(), "Standard", "", true, true);
            ApplicationTools::displayResult("Genetic Code", codeDesc);

            gCode.reset(SequenceApplicationTools::getGeneticCode(codonAlphabet->getNucleicAlphabet(), codeDesc));
        }

        //////////////////////////////////////////////
        // DATA

        VectorSiteContainer *allSites = SequenceApplicationTools::getSiteContainer(alphabet, jatiapp.getParams());

        VectorSiteContainer *sites = SequenceApplicationTools::getSitesToAnalyse(*allSites, jatiapp.getParams(), "", true, true, true, 1);
        delete allSites;

        ApplicationTools::displayResult("Number of sequences", TextTools::toString(sites->getNumberOfSequences()));
        ApplicationTools::displayResult("Number of sites", TextTools::toString(sites->getNumberOfSites()));


        /////////////////////////////////////////
        // TREE

        // Get the initial tree
        Tree *tree = 0;
        string initTreeOpt = ApplicationTools::getStringParameter("init.tree", jatiapp.getParams(), "user", "", false, 1);
        ApplicationTools::displayResult("Input tree", initTreeOpt);
        if (initTreeOpt == "user") {
            tree = PhylogeneticsApplicationTools::getTree(jatiapp.getParams());
            ApplicationTools::displayResult("Number of leaves", TextTools::toString(tree->getNumberOfLeaves()));
        } else if (initTreeOpt == "random") {
            vector<string> names = sites->getSequencesNames();
            tree = TreeTemplateTools::getRandomTree(names);
            tree->setBranchLengths(1.);
        } else throw Exception("Unknown init tree method.");

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

        string treeWIdPath = ApplicationTools::getAFilePath("output.tree_ids.file", jatiapp.getParams(), false, false, "", true, "none", 1);
        if (treeWIdPath != "none") {
            TreeTemplate<Node> ttree(*tree);
            vector<Node *> nodes = ttree.getNodes();
            for (size_t i = 0; i < nodes.size(); i++) {
                if (nodes[i]->isLeaf())
                    nodes[i]->setName(TextTools::toString(nodes[i]->getId()) + "_" + nodes[i]->getName());
                else
                    nodes[i]->setBranchProperty("NodeId", BppString(TextTools::toString(nodes[i]->getId())));
            }
            Newick treeWriter;
            treeWriter.enableExtendedBootstrapProperty("NodeId");
            ApplicationTools::displayResult("Writing tagged tree to", treeWIdPath);
            treeWriter.write(ttree, treeWIdPath);
            delete tree;
            cout << "BppML's done." << endl;
            exit(0);
        }


        /////////////////////////
        // MODEL  & LIKELIHOOD

        // Check if likelihood

        bool computeLikelihood = ApplicationTools::getBooleanParameter("compute.likelihood", jatiapp.getParams(), true, "", false, 1);
        if (!computeLikelihood) {
            delete alphabet;
            delete sites;
            delete tree;
            cout << "BppML's done. Bye." << endl;
            return 0;
        }


        DiscreteRatesAcrossSitesTreeLikelihood *tl;
        string nhOpt = ApplicationTools::getStringParameter("nonhomogeneous", jatiapp.getParams(), "no", "", true, 1);
        ApplicationTools::displayResult("Heterogeneous model", nhOpt);

        bool checkTree = ApplicationTools::getBooleanParameter("input.tree.check_root", jatiapp.getParams(), true, "", true, 2);
        bool optimizeTopo = ApplicationTools::getBooleanParameter("optimization.topology", jatiapp.getParams(), false, "", true, 1);
        unsigned int nbBS = ApplicationTools::getParameter<unsigned int>("bootstrap.number", jatiapp.getParams(), 0, "", true, 1);

        TransitionModel *model = 0;
        SubstitutionModelSet *modelSet = 0;
        DiscreteDistribution *rDist = 0;

        ////////////
        // If optimize topology

        if (optimizeTopo || nbBS > 0) {
            if (nhOpt != "no")
                throw Exception("Topology estimation with NH model not supported yet, sorry :(");
            model = PhylogeneticsApplicationTools::getTransitionModel(alphabet, gCode.get(), sites, jatiapp.getParams());
            if (model->getName() != "RE08") SiteContainerTools::changeGapsToUnknownCharacters(*sites);
            if (model->getNumberOfStates() >= 2 * model->getAlphabet()->getSize()) {
                // Markov-modulated Markov model!
                rDist = new ConstantRateDistribution();
            } else {
                rDist = PhylogeneticsApplicationTools::getRateDistribution(jatiapp.getParams());
            }
            if (dynamic_cast<MixedSubstitutionModel *>(model) == 0)
                tl = new NNIHomogeneousTreeLikelihood(*tree, *sites, model, rDist, checkTree, true);
            else
                throw Exception("Topology estimation with Mixed model not supported yet, sorry :(");
        }

            //////////////////////
            // If not topology optimization


            ///// homogeneous modeling
        else if (nhOpt == "no") {
            model = PhylogeneticsApplicationTools::getTransitionModel(alphabet, gCode.get(), sites, jatiapp.getParams());
            if (model->getName() != "RE08") SiteContainerTools::changeGapsToUnknownCharacters(*sites);
            if (model->getNumberOfStates() >= 2 * model->getAlphabet()->getSize()) {
                // Markov-modulated Markov model!
                rDist = new ConstantRateDistribution();
            } else {
                rDist = PhylogeneticsApplicationTools::getRateDistribution(jatiapp.getParams());
            }
            string recursion = ApplicationTools::getStringParameter("likelihood.recursion", jatiapp.getParams(), "simple", "", true, 1);
            ApplicationTools::displayResult("Likelihood recursion", recursion);
            if (recursion == "simple") {
                string compression = ApplicationTools::getStringParameter("likelihood.recursion_simple.compression", jatiapp.getParams(), "recursive", "", true, 2);
                ApplicationTools::displayResult("Likelihood data compression", compression);
                if (compression == "simple")
                    if (dynamic_cast<MixedSubstitutionModel *>(model))
                        tl = new RHomogeneousMixedTreeLikelihood(*tree, *sites, model, rDist, checkTree, true, false);
                    else
                        tl = new RHomogeneousTreeLikelihood(*tree, *sites, model, rDist, checkTree, true, false);

                else if (compression == "recursive")
                    if (dynamic_cast<MixedSubstitutionModel *>(model) == 0)
                        tl = new RHomogeneousTreeLikelihood(*tree, *sites, model, rDist, checkTree, true, true);
                    else
                        tl = new RHomogeneousMixedTreeLikelihood(*tree, *sites, model, rDist, checkTree, true, true);

                else throw Exception("Unknown likelihood data compression method: " + compression);
            } else if (recursion == "double") {
                if (dynamic_cast<MixedSubstitutionModel *>(model))
                    tl = new DRHomogeneousMixedTreeLikelihood(*tree, *sites, model, rDist, checkTree);
                else
                    tl = new DRHomogeneousTreeLikelihood(*tree, *sites, model, rDist, checkTree);
            } else throw Exception("Unknown recursion option: " + recursion);
        }


            ///// one per branch modeling
        else if (nhOpt == "one_per_branch") {
            model = PhylogeneticsApplicationTools::getTransitionModel(alphabet, gCode.get(), sites, jatiapp.getParams());
            if (model->getName() != "RE08") SiteContainerTools::changeGapsToUnknownCharacters(*sites);
            if (model->getNumberOfStates() >= 2 * model->getAlphabet()->getSize()) {
                // Markov-modulated Markov model!
                rDist = new ConstantRateDistribution();
            } else {
                rDist = PhylogeneticsApplicationTools::getRateDistribution(jatiapp.getParams());
            }
            vector<double> rateFreqs;
            if (model->getNumberOfStates() != alphabet->getSize()) {
                // Markov-Modulated Markov Model...
                unsigned int n = static_cast<unsigned int>(model->getNumberOfStates() / alphabet->getSize());
                rateFreqs = vector<double>(n, 1. / static_cast<double>(n)); // Equal rates assumed for now, may be changed later (actually, in the most general case,
                // we should assume a rate distribution for the root also!!!
            }

            bool stationarity = ApplicationTools::getBooleanParameter("nonhomogeneous.stationarity", jatiapp.getParams(), false, "", false, 1);
            FrequenciesSet *rootFreqs = 0;
            std::map<std::string, std::string> aliasFreqNames;
            if (!stationarity) {

                rootFreqs = PhylogeneticsApplicationTools::getRootFrequenciesSet(alphabet, gCode.get(), sites, jatiapp.getParams(), aliasFreqNames, rateFreqs);
                stationarity = !rootFreqs;
                string freqDescription = ApplicationTools::getStringParameter("nonhomogeneous.root_freq", jatiapp.getParams(), "", "", true, 1);
                if (freqDescription == "MVAprotein") {
                    if (dynamic_cast<CoalaCore *>(model)) {
                        dynamic_cast<MvaFrequenciesSet *>(rootFreqs)->setModelName("MVAprotein");
                        dynamic_cast<MvaFrequenciesSet *>(rootFreqs)->initSet(dynamic_cast<CoalaCore *>(model));
                    } else
                        throw Exception("The MVAprotein frequencies set at the root can only be used if a COaLA model is used on branches.");
                }
            }
            ApplicationTools::displayBooleanResult("Stationarity assumed", stationarity);

            vector<string> globalParameters = ApplicationTools::getVectorParameter<string>("nonhomogeneous_one_per_branch.shared_parameters", jatiapp.getParams(), ',', "");
            for (size_t i = 0; i < globalParameters.size(); i++)
                ApplicationTools::displayResult("Global parameter", globalParameters[i]);
            modelSet = SubstitutionModelSetTools::createNonHomogeneousModelSet(model, rootFreqs, tree, aliasFreqNames, globalParameters);
            model = 0;

            string recursion = ApplicationTools::getStringParameter("likelihood.recursion", jatiapp.getParams(), "simple", "", true, 1);
            ApplicationTools::displayResult("Likelihood recursion", recursion);
            if (recursion == "simple") {
                if (dynamic_cast<MixedSubstitutionModelSet *>(modelSet) != NULL)
                    tl = new RNonHomogeneousMixedTreeLikelihood(*tree, *sites, dynamic_cast<MixedSubstitutionModelSet *>(modelSet), rDist, true, true);
                else
                    tl = new RNonHomogeneousTreeLikelihood(*tree, *sites, modelSet, rDist, true, true);
            } else if (recursion == "double") {
                if (dynamic_cast<MixedSubstitutionModelSet *>(modelSet))
                    throw Exception("Double recursion with non homogeneous mixed models is not implemented yet.");
                    //            tl = new DRNonHomogeneousMixedTreeLikelihood(*tree, *sites, modelSet, rDist, true);
                else
                    tl = new DRNonHomogeneousTreeLikelihood(*tree, *sites, modelSet, rDist, true);
            } else throw Exception("Unknown recursion option: " + recursion);
        }

            /////// hand made modeling
        else if (nhOpt == "general") {
            modelSet = PhylogeneticsApplicationTools::getSubstitutionModelSet(alphabet, gCode.get(), sites, jatiapp.getParams());

            if (modelSet->getModel(0)->getName() != "RE08") SiteContainerTools::changeGapsToUnknownCharacters(*sites);

            if (modelSet->getNumberOfStates() >= 2 * modelSet->getAlphabet()->getSize()) {
                // Markov-modulated Markov model!
                rDist = new ConstantRateDistribution();
            } else {
                rDist = PhylogeneticsApplicationTools::getRateDistribution(jatiapp.getParams());
            }

            string recursion = ApplicationTools::getStringParameter("likelihood.recursion", jatiapp.getParams(), "simple", "", true, 1);
            ApplicationTools::displayResult("Likelihood recursion", recursion);
            if (recursion == "simple") {
                if (dynamic_cast<MixedSubstitutionModelSet *>(modelSet) != NULL)
                    tl = new RNonHomogeneousMixedTreeLikelihood(*tree, *sites, dynamic_cast<MixedSubstitutionModelSet *>(modelSet), rDist, true, true);
                else
                    tl = new RNonHomogeneousTreeLikelihood(*tree, *sites, modelSet, rDist, true, true);
            } else if (recursion == "double")
                if (dynamic_cast<MixedSubstitutionModelSet *>(modelSet))
                    throw Exception("Double recursion with non homogeneous mixed models is not implemented yet.");
                    //            tl = new DRNonHomogeneousMixedTreeLikelihood(*tree, *sites, modelSet, rDist, true);
                else
                    tl = new DRNonHomogeneousTreeLikelihood(*tree, *sites, modelSet, rDist, true);
            else throw Exception("Unknown recursion option: " + recursion);
        } else throw Exception("Unknown option for nonhomogeneous: " + nhOpt);

        tl->initialize();

        delete tree;

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
        //TODO: Here is the call to the parameter optimiser that is customised according to the selected model class
        ParameterList test = tl->getParameters();
        tl = dynamic_cast<DiscreteRatesAcrossSitesTreeLikelihood *>(PhylogeneticsApplicationTools::optimizeParameters(tl, tl->getParameters(), jatiapp.getParams())), "", true,
                true, 0;

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
            if (modelSet) {
                modelSet->matchParametersValues(tl->getParameters());
                PhylogeneticsApplicationTools::printParameters(modelSet, out, 1, withAlias);
            } else {
                model->matchParametersValues(tl->getParameters());
                PhylogeneticsApplicationTools::printParameters(model, out, 1, withAlias);
            }
            out.endLine();
            (out << "# Rate distribution parameters:").endLine();
            rDist->matchParametersValues(tl->getParameters());
            PhylogeneticsApplicationTools::printParameters(rDist, out, withAlias);
        }

        // Getting posterior rate class distribution:
        DiscreteDistribution *prDist = RASTools::getPosteriorRateDistribution(*tl);
        ApplicationTools::displayMessage("\nPosterior rate distribution for dataset:\n");
        if (ApplicationTools::message) prDist->print(*ApplicationTools::message);
        ApplicationTools::displayMessage("\n");
        delete prDist;

        // Write infos to file:
        string infosFile = ApplicationTools::getAFilePath("output.infos", jatiapp.getParams(), false, false);
        if (infosFile != "none") {
            ApplicationTools::displayResult("Alignment information logfile", infosFile);
            ofstream out(infosFile.c_str(), ios::out);

            // Get the rate class with maximum posterior probability:
            vector<size_t> classes = tl->getRateClassWithMaxPostProbOfEachSite();

            // Get the posterior rate, i.e. rate averaged over all posterior probabilities:
            Vdouble rates = tl->getPosteriorRateOfEachSite();

            vector<string> colNames;
            colNames.push_back("Sites");
            colNames.push_back("is.complete");
            colNames.push_back("is.constant");
            colNames.push_back("lnL");
            colNames.push_back("rc");
            colNames.push_back("pr");
            vector<string> row(6);
            DataTable *infos = new DataTable(colNames);

            for (unsigned int i = 0; i < sites->getNumberOfSites(); i++) {
                double lnL = tl->getLogLikelihoodForASite(i);
                const Site *currentSite = &sites->getSite(i);
                int currentSitePosition = currentSite->getPosition();
                string isCompl = "NA";
                string isConst = "NA";
                try { isCompl = (SiteTools::isComplete(*currentSite) ? "1" : "0"); }
                catch (EmptySiteException &ex) {}
                try { isConst = (SiteTools::isConstant(*currentSite) ? "1" : "0"); }
                catch (EmptySiteException &ex) {}
                row[0] = (string("[" + TextTools::toString(currentSitePosition) + "]"));
                row[1] = isCompl;
                row[2] = isConst;
                row[3] = TextTools::toString(lnL);
                row[4] = TextTools::toString(classes[i]);
                row[5] = TextTools::toString(rates[i]);
                infos->addRow(row);
            }

            DataTable::write(*infos, out, "\t");

            delete infos;
        }


        // Bootstrap:
        string optimizeClock = ApplicationTools::getStringParameter("optimization.clock", jatiapp.getParams(), "None", "", true, 1);
        if (nbBS > 0 && optimizeClock != "None") {
            ApplicationTools::displayError("Bootstrap is not supported with clock trees.");
        }
        if (nbBS > 0 && optimizeClock == "None") {
            ApplicationTools::displayResult("Number of bootstrap samples", TextTools::toString(nbBS));
            bool approx = ApplicationTools::getBooleanParameter("bootstrap.approximate", jatiapp.getParams(), true, "", true, 2);
            ApplicationTools::displayBooleanResult("Use approximate bootstrap", approx);
            bool bootstrapVerbose = ApplicationTools::getBooleanParameter("bootstrap.verbose", jatiapp.getParams(), false, "", true, 2);

            const Tree *initTree = tree;
            if (!bootstrapVerbose) jatiapp.getParam("optimization.verbose") = "0";
            jatiapp.getParam("optimization.profiler") = "none";
            jatiapp.getParam("optimization.messageHandler") = "none";
            if (!optimizeTopo) {
                jatiapp.getParam("optimization.topology") = "yes";
                tl = dynamic_cast<NNIHomogeneousTreeLikelihood *>(
                        PhylogeneticsApplicationTools::optimizeParameters(tl, tl->getParameters(), jatiapp.getParams(), "", true, false));
                initTree = &tl->getTree();
            }

            string bsTreesPath = ApplicationTools::getAFilePath("bootstrap.output.file", jatiapp.getParams(), false, false);
            ofstream *out = 0;
            if (bsTreesPath != "none") {
                ApplicationTools::displayResult("Bootstrap trees stored in file", bsTreesPath);
                out = new ofstream(bsTreesPath.c_str(), ios::out);
            }
            Newick newick;
            ParameterList paramsToIgnore = tl->getSubstitutionModelParameters();
            paramsToIgnore.addParameters(tl->getRateDistributionParameters());

            ApplicationTools::displayTask("Bootstrapping", true);
            vector<Tree *> bsTrees(nbBS);
            for (unsigned int i = 0; i < nbBS; i++) {
                ApplicationTools::displayGauge(i, nbBS - 1, '=');
                VectorSiteContainer *sample = SiteContainerTools::bootstrapSites(*sites);
                if (!approx) {
                    model->setFreqFromData(*sample);
                }

                if (dynamic_cast<MixedSubstitutionModel *>(model) != NULL)
                    throw Exception("Bootstrap estimation with Mixed model not supported yet, sorry :(");

                NNIHomogeneousTreeLikelihood *tlRep = new NNIHomogeneousTreeLikelihood(*initTree, *sample, model, rDist, true, false);
                tlRep->initialize();
                ParameterList parametersRep = tlRep->getParameters();
                if (approx) {
                    parametersRep.deleteParameters(paramsToIgnore.getParameterNames());
                }
                tlRep = dynamic_cast<NNIHomogeneousTreeLikelihood *>(
                        PhylogeneticsApplicationTools::optimizeParameters(tlRep, parametersRep, jatiapp.getParams(), "", true, false));
                bsTrees[i] = new TreeTemplate<Node>(tlRep->getTree());
                if (out && i == 0) newick.write(*bsTrees[i], bsTreesPath, true);
                if (out && i > 0) newick.write(*bsTrees[i], bsTreesPath, false);
                delete tlRep;
                delete sample;
            }
            if (out) out->close();
            if (out) delete out;
            ApplicationTools::displayTaskDone();


            ApplicationTools::displayTask("Compute bootstrap values");
            TreeTools::computeBootstrapValues(*tree, bsTrees);
            ApplicationTools::displayTaskDone();
            for (unsigned int i = 0; i < nbBS; i++) {
                delete bsTrees[i];
            }

            // Write resulting tree:
            PhylogeneticsApplicationTools::writeTree(*tree, jatiapp.getParams());
        }


        delete alphabet;
        delete sites;
        if (model) delete model;
        if (modelSet) delete modelSet;
        delete rDist;
        delete tl;
        delete tree;
        jatiapp.done();
    }
    catch (exception &e) {
        cout << e.what() << endl;
        return 1;
    }

    return 0;
}

