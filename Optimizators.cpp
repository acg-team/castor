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
 * @file Optimizators.cpp
 * @author Lorenzo Gatti
 * @author Massimo Maiolo
 * @date 24 02 2018
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

// From bpp-core
#include <Bpp/Io/BppODiscreteDistributionFormat.h>
#include <Bpp/Text/KeyvalTools.h>
#include <Bpp/Numeric/AutoParameter.h>
#include <Bpp/Numeric/Function/DownhillSimplexMethod.h>
#include <Bpp/Numeric/Function/PowellMultiDimensions.h>

// From bpp-seq:
#include <Bpp/Phyl/Likelihood/TreeLikelihood.h>
#include <Bpp/Phyl/OptimizationTools.h>
#include <TreeRearrangment.hpp>
#include "Optimizators.hpp"
#include "Utilities.hpp"
#include "TSHHomogeneousTreeLikelihood.hpp"
#include "TSHTopologySearch.hpp"

using namespace bpp;

namespace bpp {
    TreeLikelihood *Optimizators::optimizeParameters(
            TreeLikelihood *tl,
            const ParameterList &parameters,
            std::map<std::string, std::string> &params,
            const std::string &suffix,
            bool suffixIsOptional,
            bool verbose,
            int warn)
    throw(Exception) {
        std::string optimization = ApplicationTools::getStringParameter("optimization", params, "FullD(derivatives=Newton)", suffix, suffixIsOptional, warn);
        if (optimization == "None")
            return tl;
        std::string optName;
        std::map<std::string, std::string> optArgs;
        KeyvalTools::parseProcedure(optimization, optName, optArgs);

        unsigned int optVerbose = ApplicationTools::getParameter<unsigned int>("optimization.verbose", params, 2, suffix, suffixIsOptional, warn + 1);

        std::string mhPath = ApplicationTools::getAFilePath("optimization.message_handler", params, false, false, suffix, suffixIsOptional, "none", warn + 1);
        OutputStream *messageHandler =
                static_cast<OutputStream *>((mhPath == "none") ? 0 :
                                            (mhPath == "std") ? ApplicationTools::message :
                                            new StlOutputStream(new std::ofstream(mhPath.c_str(), std::ios::out)));
        //ApplicationTools::displayResult("Message handler", mhPath);
        LOG(INFO) << "[Parameter optimization]\tMessage handler: " << mhPath;

        std::string prPath = ApplicationTools::getAFilePath("optimization.profiler", params, false, false, suffix, suffixIsOptional, "none", warn + 1);
        OutputStream *profiler =
                static_cast<OutputStream *>((prPath == "none") ? 0 :
                                            (prPath == "std") ? ApplicationTools::message :
                                            new StlOutputStream(new std::ofstream(prPath.c_str(), std::ios::out)));
        if (profiler)
            profiler->setPrecision(20);
        //ApplicationTools::displayResult("Profiler", prPath);
        LOG(INFO) << "[Parameter optimization]\tOptimizator profiler: " << prPath;

        bool scaleFirst = ApplicationTools::getBooleanParameter("optimization.scale_first", params, false, suffix, suffixIsOptional, warn + 1);
        if (scaleFirst) {
            // We scale the tree before optimizing each branch length separately:
            //ApplicationTools::displayMessage("Scaling the tree before optimizing each branch length separately.");
            LOG(INFO) << "[Parameter optimization]\tScaling the tree before optimizing each branch length separately. ";

            double tolerance = ApplicationTools::getDoubleParameter("optimization.scale_first.tolerance", params, .0001, suffix, suffixIsOptional, warn + 1);
            //ApplicationTools::displayResult("Scaling tolerance", TextTools::toString(tolerance));
            LOG(INFO) << "[Parameter optimization]\tScaling tolerance:  " << tolerance;

            unsigned int nbEvalMax = ApplicationTools::getParameter<unsigned int>("optimization.scale_first.max_number_f_eval", params, 1000000, suffix, suffixIsOptional,
                                                                                  warn + 1);
            //ApplicationTools::displayResult("Scaling max # f eval", TextTools::toString(nbEvalMax));
            LOG(INFO) << "[Parameter optimization]\tScaling max # f eval:  " << nbEvalMax;

            OptimizationTools::optimizeTreeScale(tl, tolerance, nbEvalMax, messageHandler, profiler);
            //ApplicationTools::displayResult("New tree likelihood", -tl->getValue());
            LOG(INFO) << "[Parameter optimization]\tNew tree likelihood:  " << -tl->getValue();
        }

        // Should I ignore some parameters?
        ParameterList parametersToEstimate = parameters;
        vector<string> parNames = parametersToEstimate.getParameterNames();

        if (params.find("optimization.ignore_parameter") != params.end())
            throw Exception("optimization.ignore_parameter is deprecated, use optimization.ignore_parameters instead!");
        string paramListDesc = ApplicationTools::getStringParameter("optimization.ignore_parameters", params, "", suffix, suffixIsOptional, warn + 1);
        StringTokenizer st(paramListDesc, ",");
        while (st.hasMoreToken()) {
            try {
                string param = st.nextToken();
                if (param == "BrLen") {
                    vector<string> vs = tl->getBranchLengthsParameters().getParameterNames();
                    parametersToEstimate.deleteParameters(vs);
                    //if (verbose)
                    //    ApplicationTools::displayResult("Parameter ignored", string("Branch lengths"));
                    LOG(INFO) << "[Parameter optimization]\tParameter ignored: Branch lengths";
                } else if (param == "Ancient") {
                    NonHomogeneousTreeLikelihood *nhtl = dynamic_cast<NonHomogeneousTreeLikelihood *>(tl);
                    if (!nhtl)
                        //ApplicationTools::displayWarning("The 'Ancient' parameters do not exist in homogeneous models, and will be ignored.");
                        LOG(WARNING) << "[Parameter optimization]\tThe 'Ancient' parameters do not exist in homogeneous models, and will be ignored.";

                    else {
                        vector<string> vs = nhtl->getRootFrequenciesParameters().getParameterNames();
                        parametersToEstimate.deleteParameters(vs);
                    }
                    //if (verbose)
                    //    ApplicationTools::displayResult("Parameter ignored", string("Root frequencies"));
                    LOG(INFO) << "[Parameter optimization]\tParameter ignored: Root frequencies";

                } else if (param == "Model") {
                    vector<string> vs;
                    vector<string> vs1 = tl->getSubstitutionModelParameters().getParameterNames();
                    NonHomogeneousTreeLikelihood *nhtl = dynamic_cast<NonHomogeneousTreeLikelihood *>(tl);
                    if (nhtl != NULL) {
                        vector<string> vs2 = nhtl->getRootFrequenciesParameters().getParameterNames();
                        VectorTools::diff(vs1, vs2, vs);
                    } else
                        vs = vs1;

                    parametersToEstimate.deleteParameters(vs);
                    //if (verbose)
                    //    ApplicationTools::displayResult("Parameter ignored", string("Model"));
                    LOG(INFO) << "[Parameter optimization]\tParameter ignored: Model";

                } else if (param.find("*") != string::npos) {
                    vector<string> vs = ApplicationTools::matchingParameters(param, parNames);

                    for (vector<string>::iterator it = vs.begin(); it != vs.end(); it++) {
                        parametersToEstimate.deleteParameter(*it);
                        //if (verbose)
                        //    ApplicationTools::displayResult("Parameter ignored", *it);
                        LOG(INFO) << "[Parameter optimization]\tParameter ignored: " << *it;
                    }
                } else {
                    parametersToEstimate.deleteParameter(param);
                    //if (verbose)
                    //    ApplicationTools::displayResult("Parameter ignored", param);
                    LOG(INFO) << "[Parameter optimization]\tParameter ignored: " << param;

                }
            }
            catch (ParameterNotFoundException &pnfe) {
                //ApplicationTools::displayWarning("Parameter '" + pnfe.getParameter() + "' not found, and so can't be ignored!");
                LOG(WARNING) << "[Parameter optimization]\tParameter  " << pnfe.getParameter() << "' not found, and so can't be ignored!";

            }
        }

        // Should I constrain some parameters?
        vector<string> parToEstNames = parametersToEstimate.getParameterNames();

        if (params.find("optimization.constrain_parameter") != params.end())
            throw Exception("optimization.constrain_parameter is deprecated, use optimization.constrain_parameters instead!");
        paramListDesc = ApplicationTools::getStringParameter("optimization.constrain_parameters", params, "", suffix, suffixIsOptional, warn + 1);

        string constraint = "";
        string pc, param = "";

        StringTokenizer st2(paramListDesc, ",");
        while (st2.hasMoreToken()) {
            try {
                pc = st2.nextToken();
                std::string::size_type index = pc.find("=");
                if (index == std::string::npos)
                    throw Exception("PhylogeneticsApplicationTools::optimizeParamaters. Bad constrain syntax, should contain `=' symbol: " + pc);
                param = pc.substr(0, index);
                constraint = pc.substr(index + 1);
                IntervalConstraint ic(constraint);

                vector<string> parNames2;

                if (param == "BrLen")
                    parNames2 = tl->getBranchLengthsParameters().getParameterNames();
                else if (param == "Ancient") {
                    NonHomogeneousTreeLikelihood *nhtl = dynamic_cast<NonHomogeneousTreeLikelihood *>(tl);
                    if (!nhtl)
                        //ApplicationTools::displayWarning("The 'Ancient' parameters do not exist in homogeneous models, and will be ignored.");
                        LOG(WARNING) << "[Parameter optimization]\tThe 'Ancient' parameters do not exist in homogeneous models, and will be ignored.";

                    else {
                        parNames2 = nhtl->getRootFrequenciesParameters().getParameterNames();
                        //ApplicationTools::displayResult("Parameter ignored", string("Root frequencies"));
                        LOG(INFO) << "[Parameter optimization]\tParameter ignored: Root frequencies";
                    }
                } else if (param == "Model") {
                    vector<string> vs1 = tl->getSubstitutionModelParameters().getParameterNames();
                    NonHomogeneousTreeLikelihood *nhtl = dynamic_cast<NonHomogeneousTreeLikelihood *>(tl);
                    if (nhtl != NULL) {
                        vector<string> vs2 = nhtl->getRootFrequenciesParameters().getParameterNames();
                        VectorTools::diff(vs1, vs2, parNames2);
                    } else
                        parNames2 = vs1;
                } else if (param.find("*") != std::string::npos)
                    parNames2 = ApplicationTools::matchingParameters(param, parToEstNames);
                else
                    parNames2.push_back(param);


                for (size_t i = 0; i < parNames2.size(); i++) {
                    Parameter &par = parametersToEstimate.getParameter(parNames2[i]);
                    if (par.hasConstraint()) {
                        par.setConstraint(ic & (*par.getConstraint()), true);
                        if (par.getConstraint()->isEmpty())
                            throw Exception("Empty interval for parameter " + parNames[i] + par.getConstraint()->getDescription());
                    } else
                        par.setConstraint(ic.clone(), true);

                    //if (verbose)
                    //ApplicationTools::displayResult("Parameter constrained " + par.getName(), par.getConstraint()->getDescription());
                    LOG(INFO) << "[Parameter optimization]\tParameter constrained: " << par.getName() << par.getConstraint()->getDescription();
                }
            }
            catch (ParameterNotFoundException &pnfe) {
                //ApplicationTools::displayWarning("Parameter '" + pnfe.getParameter() + "' not found, and so can't be constrained!");
                LOG(WARNING) << "[Parameter optimization]\t Parameter '" << pnfe.getParameter() << "' not found, and so can't be constrained!";
            }
            catch (ConstraintException &pnfe) {
                throw Exception("Parameter '" + param + "' does not fit the constraint " + constraint);
            }
        }


        // /////
        // / optimization options

        unsigned int nbEvalMax = ApplicationTools::getParameter<unsigned int>("optimization.max_number_f_eval", params, 1000000, suffix, suffixIsOptional, warn + 1);
        //if (verbose)
        //    ApplicationTools::displayResult("Max # ML evaluations", TextTools::toString(nbEvalMax));
        LOG(INFO) << "[Parameter optimization]\tMax # ML evaluations: " << nbEvalMax;


        double tolerance = ApplicationTools::getDoubleParameter("optimization.tolerance", params, .000001, suffix, suffixIsOptional, warn + 1);
        //if (verbose)
        //    ApplicationTools::displayResult("Tolerance", TextTools::toString(tolerance));
        LOG(INFO) << "[Parameter optimization]\tTolerance: " << tolerance;


        // Backing up or restoring?
        unique_ptr<BackupListener> backupListener;
        string backupFile = ApplicationTools::getAFilePath("optimization.backup.file", params, false, false, suffix, suffixIsOptional, "none", warn + 1);
        if (backupFile != "none") {
            //ApplicationTools::displayResult("Parameters will be backup to", backupFile);
            LOG(INFO) << "[Parameter optimization]\tParameters will be backup to: " << backupFile;

            backupListener.reset(new BackupListener(backupFile));
            if (FileTools::fileExists(backupFile)) {
                //ApplicationTools::displayMessage("A backup file was found! Try to restore parameters from previous run...");
                LOG(INFO) << "[Parameter optimization]\tA backup file was found! Try to restore parameters from previous run... ";

                ifstream bck(backupFile.c_str(), ios::in);
                vector<string> lines = FileTools::putStreamIntoVectorOfStrings(bck);
                double fval = TextTools::toDouble(lines[0].substr(5));
                ParameterList pl = tl->getParameters();
                for (size_t l = 1; l < lines.size(); ++l) {
                    if (!TextTools::isEmpty(lines[l])) {
                        StringTokenizer stp(lines[l], "=");
                        if (stp.numberOfRemainingTokens() != 2) {
                            cerr << "Corrupted backup file!!!" << endl;
                            cerr << "at line " << l << ": " << lines[l] << endl;
                        }
                        string pname = stp.nextToken();
                        string pvalue = stp.nextToken();
                        size_t p = pl.whichParameterHasName(pname);
                        pl.setParameter(p, AutoParameter(pl[p]));
                        pl[p].setValue(TextTools::toDouble(pvalue));
                    }
                }
                bck.close();
                tl->setParameters(pl);
                if (abs(tl->getValue() - fval) > 0.000001)
                    //ApplicationTools::displayWarning("Warning, incorrect likelihood value after restoring from backup file.");
                    LOG(WARNING) << "[Parameter optimization]\tWarning, incorrect likelihood value after restoring from backup file. ";
                //ApplicationTools::displayResult("Restoring log-likelihood", -fval);
                LOG(INFO) << "[Parameter optimization]\tRestoring log-likelihood: " << -fval;

            }
        }

        // There it goes...
        bool optimizeTopo = ApplicationTools::getBooleanParameter("optimization.topology", params, false, suffix, suffixIsOptional, warn + 1);
        //if (verbose)
        //    ApplicationTools::displayResult("Optimize topology", optimizeTopo ? "yes" : "no");

        LOG(INFO) << "[Parameter optimization]\tOptimize topology: " << (optimizeTopo ? "yes" : "no");

        string nniMethod = ApplicationTools::getStringParameter("optimization.topology.algorithm_nni.method", params, "phyml", suffix, suffixIsOptional, warn + 1);
        string nniAlgo;
        if (nniMethod == "fast") {
            nniAlgo = NNITopologySearch::FAST;
        } else if (nniMethod == "better") {
            nniAlgo = NNITopologySearch::BETTER;
        } else if (nniMethod == "phyml") {
            nniAlgo = NNITopologySearch::PHYML;
        } else
            throw Exception("Unknown NNI algorithm: '" + nniMethod + "'.");


        string order = ApplicationTools::getStringParameter("derivatives", optArgs, "Newton", "", true, warn + 1);
        string optMethodDeriv;
        if (order == "Gradient") {
            optMethodDeriv = OptimizationTools::OPTIMIZATION_GRADIENT;
        } else if (order == "Newton") {
            optMethodDeriv = OptimizationTools::OPTIMIZATION_NEWTON;
        } else if (order == "BFGS") {
            optMethodDeriv = OptimizationTools::OPTIMIZATION_BFGS;
        } else
            throw Exception("Unknown derivatives algorithm: '" + order + "'.");

        //if (verbose)
        //    ApplicationTools::displayResult("Optimization method", optName);
        //if (verbose)
        //    ApplicationTools::displayResult("Algorithm used for derivable parameters", order);

        LOG(INFO) << "[Parameter optimization]\tOptimization method: " << optName;
        LOG(INFO) << "[Parameter optimization]\tAlgorithm used for derivable parameters: " << order;


        // See if we should reparametrize:
        bool reparam = ApplicationTools::getBooleanParameter("optimization.reparametrization", params, false, suffix, suffixIsOptional, warn + 1);
        //if (verbose)
        //    ApplicationTools::displayResult("Reparametrization", (reparam ? "yes" : "no"));

        LOG(INFO) << "[Parameter optimization]\tReparametrization: " << (reparam ? "yes" : "no");


        // See if we should use a molecular clock constraint:
        string clock = ApplicationTools::getStringParameter("optimization.clock", params, "None", suffix, suffixIsOptional, warn + 1);
        if (clock != "None" && clock != "Global")
            throw Exception("Molecular clock option not recognized, should be one of 'Global' or 'None'.");
        bool useClock = (clock == "Global");
        if (useClock && optimizeTopo)
            throw Exception("PhylogeneticsApplicationTools::optimizeParameters. Cannot optimize topology with a molecular clock.");

        //if (verbose)
        //    ApplicationTools::displayResult("Molecular clock", clock);

        LOG(INFO) << "[Parameter optimization]\tMolecular clock: " << clock;


        unsigned int n = 0;
        if ((optName == "D-Brent") || (optName == "D-BFGS")) {
            // Uses Newton-Brent method or Newton-BFGS method
            string optMethodModel;
            if (optName == "D-Brent")
                optMethodModel = OptimizationTools::OPTIMIZATION_BRENT;
            else
                optMethodModel = OptimizationTools::OPTIMIZATION_BFGS;

            unsigned int nstep = ApplicationTools::getParameter<unsigned int>("nstep", optArgs, 1, "", true, warn + 1);

            if (optimizeTopo) {
                tshlib::TreeRearrangmentOperations treesearch_operations;
                tshlib::TreeSearchHeuristics treesearch_heuristics;
                auto ttl = dynamic_cast<AbstractHomogeneousTreeLikelihood *>(tl);
                auto ntl = new bpp::TSHHomogeneousTreeLikelihood(ttl, (*ttl->getData()), (ttl->getModel()), (ttl->getRateDistribution()));

                /*
                bool optNumFirst = ApplicationTools::getBooleanParameter("optimization.topology.numfirst", params, true, suffix, suffixIsOptional, warn + 1);
                unsigned int topoNbStep = ApplicationTools::getParameter<unsigned int>("optimization.topology.nstep", params, 1, suffix, suffixIsOptional, warn + 1);
                double tolBefore = ApplicationTools::getDoubleParameter("optimization.topology.tolerance.before", params, 100, suffix, suffixIsOptional, warn + 1);
                double tolDuring = ApplicationTools::getDoubleParameter("optimization.topology.tolerance.during", params, 100, suffix, suffixIsOptional, warn + 1);
                tl = OptimizationTools::optimizeTreeNNI(
                        dynamic_cast<NNIHomogeneousTreeLikelihood *>(tl), parametersToEstimate,
                        optNumFirst, tolBefore, tolDuring, nbEvalMax, topoNbStep, messageHandler, profiler,
                        reparam, optVerbose, optMethodDeriv, nstep, nniAlgo);
                */
                std::string PAR_optim_topology_algorithm = ApplicationTools::getStringParameter("optim_topology_algorithm", params, "no-search", "", true, true);

                if (PAR_optim_topology_algorithm.find("greedy") != std::string::npos) {
                    treesearch_heuristics = tshlib::TreeSearchHeuristics::greedy;
                } else if (PAR_optim_topology_algorithm.find("hillclimbing") != std::string::npos) {
                    treesearch_heuristics = tshlib::TreeSearchHeuristics::hillclimbing;
                } else if (PAR_optim_topology_algorithm.find("no-search") != std::string::npos) {
                    treesearch_heuristics = tshlib::TreeSearchHeuristics::nosearch;
                }

                if (treesearch_heuristics != tshlib::TreeSearchHeuristics::nosearch) {
                    std::string PAR_lkmove = ApplicationTools::getStringParameter("lk_move", params, "bothways", "", true, true);
                    std::string PAR_optim_topology_operations = ApplicationTools::getStringParameter("optim_topology_operations", params, "best-search", "", true, true);
                    int PAR_optim_topology_maxcycles = ApplicationTools::getIntParameter("optim_topology_maxcycles", params, 1, "", true, 0);
                    int PAR_optim_topology_hillclimbing_startnodes = ApplicationTools::getIntParameter("optim_topology_numnodes", params, 1, "", true, 0);

                    if (PAR_optim_topology_operations.find("best-search") != std::string::npos) {
                        treesearch_operations = tshlib::TreeRearrangmentOperations::classic_Mixed;
                    } else if (PAR_optim_topology_operations.find("nni-search") != std::string::npos) {
                        treesearch_operations = tshlib::TreeRearrangmentOperations::classic_NNI;
                    } else if (PAR_optim_topology_operations.find("spr-search") != std::string::npos) {
                        treesearch_operations = tshlib::TreeRearrangmentOperations::classic_SPR;
                    } else {}

                    auto treesearch = new tshlib::TreeSearch;
                    /*
                    treesearch->setTreeSearchStrategy(treesearch_heuristics, treesearch_operations);
                    treesearch->setInitialLikelihoodValue(-ntl->getLikelihoodFunction()->getValue());
                    treesearch->setScoringMethod(PAR_lkmove);
                    treesearch->setStartingNodes(PAR_optim_topology_hillclimbing_startnodes);
                    treesearch->setTreemap(tm);
                    treesearch->setStopCondition(tshlib::TreeSearchStopCondition::iterations, (double) PAR_optim_topology_maxcycles);
                    if (PAR_model_indels) {
                        treesearch->setModelIndels(true);
                        treesearch->setLikelihoodFunc(ntl);
                    } else {
                        treesearch->setModelIndels(false);

                    }
                    logLK = treesearch->performTreeSearch(utree);
                    */
                }

                //OutputUtils::printParametersLikelihood(tl);


            }

            //if (verbose && nstep > 1)
            //    ApplicationTools::displayResult("# of precision steps", TextTools::toString(nstep));
            LOG(INFO) << "[Parameter optimization]\t# of precision steps: " << nstep;

            parametersToEstimate.matchParametersValues(tl->getParameters());
            n = OptimizationTools::optimizeNumericalParameters(
                    dynamic_cast<DiscreteRatesAcrossSitesTreeLikelihood *>(tl), parametersToEstimate,
                    backupListener.get(), nstep, tolerance, nbEvalMax, messageHandler, profiler, reparam, optVerbose, optMethodDeriv, optMethodModel);
        } else if (optName == "FullD") {
            // Uses Newton-raphson algorithm with numerical derivatives when required.

            if (optimizeTopo) {
                bool optNumFirst = ApplicationTools::getBooleanParameter("optimization.topology.numfirst", params, true, suffix, suffixIsOptional, warn + 1);
                unsigned int topoNbStep = ApplicationTools::getParameter<unsigned int>("optimization.topology.nstep", params, 1, suffix, suffixIsOptional, warn + 1);
                double tolBefore = ApplicationTools::getDoubleParameter("optimization.topology.tolerance.before", params, 100, suffix, suffixIsOptional, warn + 1);
                double tolDuring = ApplicationTools::getDoubleParameter("optimization.topology.tolerance.during", params, 100, suffix, suffixIsOptional, warn + 1);
                tl = OptimizationTools::optimizeTreeNNI2(
                        dynamic_cast<NNIHomogeneousTreeLikelihood *>(tl), parametersToEstimate,
                        optNumFirst, tolBefore, tolDuring, nbEvalMax, topoNbStep, messageHandler, profiler,
                        reparam, optVerbose, optMethodDeriv, nniAlgo);
            }

            parametersToEstimate.matchParametersValues(tl->getParameters());
            n = OptimizationTools::optimizeNumericalParameters2(
                    dynamic_cast<DiscreteRatesAcrossSitesTreeLikelihood *>(tl), parametersToEstimate,
                    backupListener.get(), tolerance, nbEvalMax, messageHandler, profiler, reparam, useClock, optVerbose, optMethodDeriv);
        } else
            throw Exception("Unknown optimization method: " + optName);

        string finalMethod = ApplicationTools::getStringParameter("optimization.final", params, "none", suffix, suffixIsOptional, warn + 1);
        Optimizer *finalOptimizer = 0;
        if (finalMethod == "none") {}
        else if (finalMethod == "simplex") {
            finalOptimizer = new DownhillSimplexMethod(tl);
        } else if (finalMethod == "powell") {
            finalOptimizer = new PowellMultiDimensions(tl);
        } else
            throw Exception("Unknown final optimization method: " + finalMethod);

        if (finalOptimizer) {
            parametersToEstimate.matchParametersValues(tl->getParameters());
            //if (verbose)
            //    ApplicationTools::displayResult("Final optimization step", finalMethod);
            LOG(INFO) << "[Parameter optimization]\tFinal optimization step: " << finalMethod;

            finalOptimizer->setProfiler(profiler);
            finalOptimizer->setMessageHandler(messageHandler);
            finalOptimizer->setMaximumNumberOfEvaluations(nbEvalMax);
            finalOptimizer->getStopCondition()->setTolerance(tolerance);
            finalOptimizer->setVerbose(verbose);
            finalOptimizer->setConstraintPolicy(AutoParameter::CONSTRAINTS_AUTO);
            finalOptimizer->init(parametersToEstimate);
            finalOptimizer->optimize();
            n += finalOptimizer->getNumberOfEvaluations();
            delete finalOptimizer;
        }

        //if (verbose)
        //    ApplicationTools::displayResult("Performed", TextTools::toString(n) + " function evaluations.");
        LOG(INFO) << "[Parameter optimization]\tPerformed " << n << " function evaluations.";

        if (backupFile != "none") {
            string bf = backupFile + ".def";
            rename(backupFile.c_str(), bf.c_str());
        }
        return tl;
    }

}
