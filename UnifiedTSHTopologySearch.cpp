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
 * @file TSHTopologySearch.cpp
 * @author Lorenzo Gatti
 * @author Massimo Maiolo
 * @date 05 02 2018
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
 * @see For more information visit: https://bitbucket.org/acg-team/minijati/wiki/Home
 */
#include <random>
#include <glog/logging.h>

#include <Bpp/Phyl/Likelihood/RHomogeneousTreeLikelihood.h>
#include <Bpp/Phyl/Io/Newick.h>

#include "UnifiedTSHTopologySearch.hpp"


using namespace bpp;

tshlib::TreeRearrangment *tshlib::TreeSearch::defineMoves() {

    // Initialise a new rearrangement list
    auto candidateMoveSet = new tshlib::TreeRearrangment;
    candidateMoveSet->setTreeTopology(utree_);
    candidateMoveSet->setTreeCoverage(tshRearrangementCoverage);
    candidateMoveSet->setStrategy(tshSearchHeuristic);
    candidateMoveSet->initialize();

    std::vector<tshlib::VirtualNode *> pickedNodes((int)search_startingnodes);

    // instantiating the method to define the number of initial nodes on which computing the candidate moves
    switch (tshStartingNodeMethod) {

        case tshlib::StartingNodeHeuristics::greedy:

            // Generate candidate list of possible moves given the tree topology and the rearrangement operation type
            for (auto &node:utree_->listVNodes) {

                // Print node description with neighbors
                //DVLOG(2) << "[utree neighbours] " << vnode->printNeighbours() << std::endl;

                candidateMoveSet->setSourceNode(node);

                // Get all the target nodes with distance == radius from the source node
                // excluding the starting node (false)
                // do not duplicate moves in the list
                candidateMoveSet->defineMoves(false, false);

            }
            break;

        case tshlib::StartingNodeHeuristics::hillclimbing:

            RandomTools::getSample(utree_->listVNodes,pickedNodes);

            // Generate candidate list of possible moves given the node topology and the rearrangement operation type
            for(auto &node:pickedNodes){

                // Print node description with neighbors
                //DVLOG(2) << "[utree neighbours] " << vnode->printNeighbours() << std::endl;

                candidateMoveSet->setSourceNode(node);

                // Get all the target nodes with distance == radius from the source node
                // excluding the starting node (false)
                // do not duplicate moves in the list
                candidateMoveSet->defineMoves(false, false);

            }
            break;

        case tshlib::StartingNodeHeuristics::particle_swarm:
        case tshlib::StartingNodeHeuristics::undef:
            LOG(FATAL) << "[TSH Optimisation] Strarting node heuristic not defined. Execution aborted!";
            break;
    }

    // Print the list of moves for the current P node (source node)
    //rearrangmentList.printMoves();

    return candidateMoveSet;
}


double tshlib::TreeSearch::executeTreeSearch() {
    double newScore = 0;
    setTreeSearchStatus(false);
    // Store initial likelihood value
    setInitialLikelihoodValue(likelihoodFunc->getLogLikelihood());
    //tshinitScore =  likelihoodFunc->getLogLikelihood();
    tshcycleScore = tshinitScore;

    // Remove root on the tree topology
    utree_->removeVirtualRootNode();

   DLOG(INFO) << "[TSH Optimisation] Initial topology: " << utree_->printTreeNewick(true);
   DLOG(INFO) << "[TSH Optimisation] Initial likelihood: " << TextTools::toString(tshinitScore, 15);
   DLOG(INFO) << "[TSH Optimisation] algorithm: " << getStartingNodeHeuristicDescription()
              << ", moves: " << getRearrangmentCoverageDescription() << ", "
              << "starting nodes: " + std::to_string(search_startingnodes);

    // execute tree search according to settings
    newScore = iterate();

    // Add root on the tree topology
    utree_->addVirtualRootNode();

    return newScore;
}


void tshlib::TreeSearch::testMoves(tshlib::TreeRearrangment *candidateMoves) {

    try {
        bool status = false;
        for (unsigned long i = 0; i < candidateMoves->getNumberOfMoves(); i++) {
            ApplicationTools::displayGauge(i + 1, candidateMoves->getNumberOfMoves(), '>', std::string("node " + candidateMoves->getMove(i)->getSourceNode()->vnode_name));

            std::vector<tshlib::VirtualNode *> listNodesWithinPath, updatedNodesWithinPath;
            double moveLogLK = 0;

            // ------------------------------------
            // Prepare the list of nodes involved in the move (Required here!)
            listNodesWithinPath = utree_->computePathBetweenNodes(candidateMoves->getMove(i)->getSourceNode(), candidateMoves->getMove(i)->getTargetNode());
            updatedNodesWithinPath = candidateMoves->updatePathBetweenNodes(i, listNodesWithinPath);

            // ------------------------------------
            // Log move details
           DVLOG(1) << "Move [" << candidateMoves->getMove(i)->getSourceNode()->getNodeName() << "->" << candidateMoves->getMove(i)->getTargetNode()->getNodeName() << "]\t"
                      << candidateMoves->getMove(i)->getClass() << "\t(" << candidateMoves->getMove(i)->getRadius() << ")" << "\tdirection: " << candidateMoves->getMove(i)->getDirection();

            // ------------------------------------
            // Apply the move
            status = candidateMoves->applyMove(i);

            if(status) {

               DVLOG(1) << "[TSH Cycle - Topology] [A] MOVE#" << candidateMoves->getMove(i)->getUID() << " | " << utree_->printTreeNewick(true);

                // ------------------------------------
                // Print root reachability from every node (includes rotations)
                // inputTree->_testReachingPseudoRoot();

                // ------------------------------------
                // Print tree on file
                //inputTree->saveTreeOnFile("../data/test.txt");

                // ------------------------------------
                // Add root node to the list of nodes involved in the tree-rearrangement
                updatedNodesWithinPath.push_back(utree_->rootnode);
                listNodesWithinPath.push_back(utree_->rootnode);

                // ------------------------------------
                // Recompute the likelihood
                if (dynamic_cast<UnifiedTSHomogeneousTreeLikelihood_PIP *>(likelihoodFunc)) {
                    moveLogLK = dynamic_cast<UnifiedTSHomogeneousTreeLikelihood_PIP *>(likelihoodFunc)->updateLikelihoodOnTreeRearrangement(updatedNodesWithinPath);
                } else {
                    moveLogLK = dynamic_cast<UnifiedTSHomogeneousTreeLikelihood *>(likelihoodFunc)->updateLikelihoodOnTreeRearrangement(updatedNodesWithinPath);
                }


                LOG_IF(FATAL, std::isinf(moveLogLK)) << "llk[Move] value is -inf for [MOVE " << candidateMoves->getMove(i)->getUID() << "]" <<
                                                     debugStackTraceMove(candidateMoves->getMove(i), utree_,
                                                                         listNodesWithinPath,
                                                                         updatedNodesWithinPath,
                                                                         tshinitScore,
                                                                         moveLogLK,
                                                                         0);

                // ------------------------------------
                // Store likelihood of the move
                candidateMoves->getMove(i)->setScore(moveLogLK);

                // ------------------------------------
                // Display status of the rearrangement (deprecated)
                //candidateMoves->displayRearrangmentStatus(i, true);

                // ------------------------------------
                // Revert the move, and return to the original tree
                candidateMoves->revertMove(i);
               DVLOG(1) << "[TSH Cycle - Topology] [R] MOVE#" << candidateMoves->getMove(i)->getUID() << " | " << utree_->printTreeNewick(true);

                // ------------------------------------
                // Recompute the likelihood after reverting the move
                double moveLogLK_return = -std::numeric_limits<double>::infinity();

                if (dynamic_cast<UnifiedTSHomogeneousTreeLikelihood_PIP *>(likelihoodFunc)) {
                    moveLogLK_return = dynamic_cast<UnifiedTSHomogeneousTreeLikelihood_PIP *>(likelihoodFunc)->updateLikelihoodOnTreeRearrangement(listNodesWithinPath);
                } else {
                    moveLogLK_return = dynamic_cast<UnifiedTSHomogeneousTreeLikelihood *>(likelihoodFunc)->updateLikelihoodOnTreeRearrangement(listNodesWithinPath);
                }

                // ------------------------------------
                // Display status of the rearrangement (deprecated)
                //candidateMoves->displayRearrangmentStatus(i, true);

                LOG_IF(FATAL, std::isinf(moveLogLK_return)) << "llk[Return] value is -inf for [MOVE " << candidateMoves->getMove(i)->getUID() << "]" <<
                                                            debugStackTraceMove(candidateMoves->getMove(i),
                                                                                utree_,
                                                                                listNodesWithinPath,
                                                                                updatedNodesWithinPath,
                                                                                tshinitScore,
                                                                                moveLogLK,
                                                                                moveLogLK_return);

                LOG_IF(FATAL, !ComparisonUtils::areLogicallyEqual(moveLogLK_return, tshinitScore)) << "Error in evaluating likelihood [MOVE " << candidateMoves->getMove(i)->getUID() << "]" <<
                                                                                                   debugStackTraceMove(candidateMoves->getMove(i),
                                                                                                                       utree_,
                                                                                                                       listNodesWithinPath,
                                                                                                                       updatedNodesWithinPath,
                                                                                                                       tshinitScore,
                                                                                                                       moveLogLK,
                                                                                                                       moveLogLK_return);

                // ------------------------------------
                // Count moves performed
                performed_moves += candidateMoves->getNumberOfMoves() * 2;
            }
        }

    } catch (const std::overflow_error &e) {

        LOG(ERROR) << e.what() ;

    } catch (const std::runtime_error &e) {

        LOG(ERROR) << e.what() ;

    } catch (const std::exception &e) {

        LOG(ERROR) << e.what() ;
    }

    ApplicationTools::displayMessage("");

}


double tshlib::TreeSearch::iterate() {
    std::chrono::high_resolution_clock::time_point t1 = std::chrono::high_resolution_clock::now();

    // Condition handler
    //double cycle_no = 0;
    double c = std::abs(tshinitScore);

    while (toleranceValue < c) {

        // ------------------------------------
        // 1. Define moves according to tree-search criteria
        tshlib::TreeRearrangment *candidateMoves = defineMoves();

        std::ostringstream taskDescription;
        taskDescription << "Tree-search cycle #" << std::setfill('0') << std::setw(3) << performed_cycles + 1
                        << " | Testing " << candidateMoves->getNumberOfMoves() << " tree rearrangements.";
        ApplicationTools::displayMessage(taskDescription.str());
       DVLOG(1) << "[TSH Optimisation] Cycle " << std::setfill('0') << std::setw(3) << performed_cycles + 1 << " - lk= " << TextTools::toString(tshcycleScore, 15);

        // ------------------------------------
        // 2. Test and record likelihood of each and every candidate move
        testMoves(candidateMoves);

        // ------------------------------------
        // 3. Select the best move in the list and store it
        Move *bestMove = candidateMoves->selectBestMove(tshcycleScore);

        // number of cycles
        performed_cycles = performed_cycles+1;
        // ------------------------------------
        // Print winning move
        if (bestMove) {

            setTreeSearchStatus(true);

            DVLOG(1) << "[TSH Cycle - Topology] (cycle" << std::setfill('0') << std::setw(3) << performed_cycles + 1 << ")\t"
                      << "lk = " << std::setprecision(12) << bestMove->getScore() << " | "
                      << bestMove->getClass() << "." << std::setfill('0') << std::setw(3) << bestMove->getUID()
                      << " [" << bestMove->getSourceNode()->getNodeName() << " -> " << bestMove->getTargetNode()->getNodeName() << "]";


            taskDescription.str(std::string());
            taskDescription << "Tree-search cycle #" << std::setfill('0') << std::setw(3) << performed_cycles + 1 << " | ["
                            << bestMove->getClass() << "." << std::setfill('0') << std::setw(3) << bestMove->getUID()
                            << "] New likelihood: " << std::setprecision(12) << bestMove->getScore();
            ApplicationTools::displayMessage(taskDescription.str());


            std::vector < tshlib::VirtualNode * > listNodesWithinPath, updatedNodesWithinPath;
            listNodesWithinPath = utree_->computePathBetweenNodes(bestMove->getSourceNode(), bestMove->getTargetNode());
            updatedNodesWithinPath = candidateMoves->updatePathBetweenNodes(bestMove->getUID(), listNodesWithinPath);
            updatedNodesWithinPath.push_back(utree_->rootnode);

            // ------------------------------------
            // Commit final move on the topology
            candidateMoves->commitMove(bestMove->getUID());

            DVLOG(1) << "utree after commit " << utree_->printTreeNewick(true);

            ApplicationTools::displayTask("Optimising " + TextTools::toString(updatedNodesWithinPath.size()) + " branches");

            if (dynamic_cast<UnifiedTSHomogeneousTreeLikelihood_PIP *>(likelihoodFunc)) {
                dynamic_cast<UnifiedTSHomogeneousTreeLikelihood_PIP *>(likelihoodFunc)->topologyChangeSuccessful(updatedNodesWithinPath);
            } else {
                dynamic_cast<UnifiedTSHomogeneousTreeLikelihood *>(likelihoodFunc)->topologyChangeSuccessful(updatedNodesWithinPath);
            }

            ApplicationTools::displayTaskDone();
            ApplicationTools::displayMessage("");

            DVLOG(1) << "bpp after commit " << OutputUtils::TreeTools::writeTree2String(likelihoodFunc->getTree().clone());

            tshcycleScore = -likelihoodFunc->getValue();

            // Update tolerance
            c = std::abs(tshinitScore) - std::abs(tshcycleScore);
            tshinitScore = tshcycleScore;

            // ------------------------------------
            // Clean memory
            delete candidateMoves;

        } else {

           DLOG(INFO) << "[TSH Cycle] No further likelihood improvement after " << performed_cycles << " cycles and " << performed_moves << " performed moves. Exit loop.";

            // ------------------------------------
            // Clean memory
            delete candidateMoves;

            break;
        }

        if(performed_cycles == maxTSCycles){
           DLOG(INFO) << "[TSH Cycle] Reached max number of tree-search cycles after " << performed_cycles << " cycles";
            break;
        }



    }


    // ------------------------------------
    std::chrono::high_resolution_clock::time_point t2 = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::microseconds>(t2 - t1).count();

   DVLOG(1) << "[TSH Cycle] Moves applied and reverted: " << performed_moves << std::endl;
   DVLOG(1) << "[TSH Cycle] Elapsed time: " << duration << " microseconds" << std::endl;
   DVLOG(1) << "[TSH Cycle] *** " << (double) duration / performed_moves << " microseconds/move *** " << std::endl;

    return -tshcycleScore;
}


std::string tshlib::TreeSearch::debugStackTraceMove(Move *move, Utree *tree,
                                                           std::vector < tshlib::VirtualNode * > listNodesInvolved,
                                                           std::vector < tshlib::VirtualNode * > updatedList,
                                                           double initLK, double moveLK, double returnLK){

    std::ostringstream nodepath_lni;
    for (auto &node:listNodesInvolved) {
        nodepath_lni << node->getNodeName() << ">";
    }
    std::ostringstream nodepath_uni;
    for (auto &node:updatedList) {
        nodepath_uni << node->getNodeName() << ">";
    }


    std::ostringstream stm ;
    stm << std::endl << "*** Stack trace [MOVE " << move->getClass() << "." << move->getUID() <<"] (" <<
        move->getSourceNode()->getNodeName() <<" -> "<<move->getTargetNode()->getNodeName()<<") - " << move->getDirection() <<" ***"<< std::endl;
    stm << "    @        [ReferenNodeList]  " << nodepath_lni.str() << std::endl;
    stm << "    @        [UpdatedNodeList]  " << nodepath_uni.str()  << std::endl;
    stm << "    @        [LLK.Initial]      " << TextTools::toString(initLK, 15)  << std::endl;
    stm << "    @        [LLK.Return]       " << TextTools::toString(returnLK, 15)  << std::endl;
    stm << "    @        [LLK.Move]         " << TextTools::toString(moveLK, 15) << std::endl;

    //stm << "    @ [Tree] " << tree->printTreeNewick(true);

    return stm.str();

}
