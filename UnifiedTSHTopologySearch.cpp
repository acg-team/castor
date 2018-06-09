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
#include <random>
#include <glog/logging.h>

#include <Bpp/Phyl/Likelihood/RHomogeneousTreeLikelihood.h>
#include <Bpp/Phyl/Io/Newick.h>

#include "UnifiedTSHTopologySearch.hpp"


using namespace bpp;

tshlib::TreeRearrangment *tshlib::TreeSearch::defineCandidateMoves() {

    int min_radius;
    int max_radius;

    // Define the radius for pruning and regrafting the input tree.
    switch (tshOperations) {

        case tshlib::TreeRearrangmentOperations::classic_NNI:
            min_radius = 3;
            max_radius = 3;
            break;

        case tshlib::TreeRearrangmentOperations::classic_SPR:
            min_radius = 4;
            max_radius = utree_->getMaxNodeDistance() / 2;
            break;

        case tshlib::TreeRearrangmentOperations::classic_TBR:
            min_radius = 5;
            max_radius = utree_->getMaxNodeDistance() / 2;
            break;

        case tshlib::TreeRearrangmentOperations::classic_Mixed:
            min_radius = 3;  // Minimum radius for an NNI move is 3 nodes
            max_radius = utree_->getMaxNodeDistance(); // Full tree traversing from any nodeInterface of the tree
            break;

    }

    // Initialise a new rearrangement list
    auto candidateMoveSet = new tshlib::TreeRearrangment;
    candidateMoveSet->setTreeTopology(utree_);
    candidateMoveSet->setMinRadius(min_radius);
    candidateMoveSet->setMaxRadius(max_radius);


    switch (tshStrategy) {
        case tshlib::TreeSearchHeuristics::greedy:
            // Generate candidate list of possible moves given the tree topology and the rearrangement operation type
            for (auto &node:utree_->listVNodes) {
                // Print nodeInterface description with neighbors
                //VLOG(2) << "[utree neighbours] " << vnode->printNeighbours() << std::endl;

                candidateMoveSet->setSourceNode(node);
                // Get all the target nodes with distance == radius from the source nodeInterface
                // excluding the starting nodeInterface (false)
                // do not duplicate moves in the list
                candidateMoveSet->defineMoves(false, false);

            }
            break;

        case tshlib::TreeSearchHeuristics::hillclimbing:

            // fill array with [min value, max_value] number of nodes in the tree
            //std::vector<int> node_ids(utree_->listVNodes.size());
            //std::iota(node_ids.begin(), node_ids.end(), 0);

            // shuffle the array
            //std::random_device rd;
            //std::mt19937 e{rd()};
            //std::shuffle(node_ids.begin(), node_ids.end(), e);

            std::vector<tshlib::VirtualNode *> pickedNodes((int)search_startingnodes);
            RandomTools::getSample(utree_->listVNodes,pickedNodes);

            // Generate candidate list of possible moves given the nodeInterface topology and the rearrangement operation type
            //for (int i = 0; i < search_startingnodes; i++) {
            for(auto &node:pickedNodes){
                //VirtualNode *nodeInterface = utree_->listVNodes.at(node_ids.at(i));

                // Print nodeInterface description with neighbors
                //VLOG(2) << "[utree neighbours] " << vnode->printNeighbours() << std::endl;

                candidateMoveSet->setSourceNode(node);
                // Get all the target nodes with distance == radius from the source nodeInterface
                // excluding the starting nodeInterface (false)
                // do not duplicate moves in the list
                candidateMoveSet->defineMoves(false, false);

            }

            break;
    }

    // Print the list of moves for the current P nodeInterface (source nodeInterface)
    //rearrangmentList.printMoves();

    return candidateMoveSet;
}


double tshlib::TreeSearch::performTreeSearch() {
    double newScore = 0;

    // Store initial likelihood value
    tshinitScore = likelihoodFunc->getLogLikelihood();
    tshcycleScore = tshinitScore;


    utree_->removeVirtualRootNode();
    std::string stopcondition;
    switch (stopConditionMethod) {
        case tshlib::TreeSearchStopCondition::convergence:
            stopcondition = "convergence";
            break;
        case tshlib::TreeSearchStopCondition::iterations:
            stopcondition = "max iterations";
            break;
    }
    std::string startingnodes = "";
    if (search_startingnodes > 0)
        startingnodes = "starting nodes: " + std::to_string(search_startingnodes);

    std::string operations;
    switch (tshOperations) {
        case tshlib::TreeRearrangmentOperations::classic_NNI:
            operations = "classic_NNI";
            break;
        case tshlib::TreeRearrangmentOperations::classic_SPR:
            operations = "classic_SPR";
            break;
        case tshlib::TreeRearrangmentOperations::classic_TBR:
            operations = "classic_TBR";
            break;
        case tshlib::TreeRearrangmentOperations::classic_Mixed:
            operations = "classic_Best";
            break;

    }

    LOG(INFO) << "[TSH Optimisation] Initial topology: " << utree_->printTreeNewick(true);
    LOG(INFO) << "[TSH Optimisation] Initial likelihood: " << TextTools::toString(tshinitScore, 15);


    // define the high-level strategy to evaluate the tree rearrangement operations
    std::string strategy;
    switch (tshStrategy) {
        case tshlib::TreeSearchHeuristics::greedy:
            strategy = "Greedy";
            LOG(INFO) << "[TSH Optimisation] algorithm: " << strategy << ", moves: " << operations << ", " << startingnodes << ", stop condition: " << stopcondition;
            newScore = greedy();
            break;

        case tshlib::TreeSearchHeuristics::hillclimbing:
            strategy = "Hillclimbing";
            LOG(INFO) << "[TSH Optimisation] algorithm: " << strategy << ", moves: " << operations << ", " << startingnodes << ", stop condition: " << stopcondition;
            newScore = hillclimbing();
            break;

        case tshlib::TreeSearchHeuristics::particle_swarm:
            strategy = "Particle Swarm";
            LOG(INFO) << "[TSH Optimisation] algorithm: " << strategy << ", moves: " << operations << ", " << startingnodes << ", stop condition: " << stopcondition;
            newScore = particleswarming();
            break;

        case tshlib::TreeSearchHeuristics::nosearch:
            LOG(INFO) << "[TSH Optimisation] No tree search optimisation!";
            break;
        default:
            LOG(FATAL) << "[TSH Optimisation] Tree-search heuristic non implemented. Execution aborted!";
            break;
    }
    utree_->addVirtualRootNode();
    return newScore;
}


void tshlib::TreeSearch::testCandidateMoves(tshlib::TreeRearrangment *candidateMoves) {


    for (unsigned long i = 0; i < candidateMoves->getNumberOfMoves(); i++) {
        ApplicationTools::displayGauge(i + 1, candidateMoves->getNumberOfMoves(), '>', std::string("nodeInterface " + candidateMoves->getMove(i)->getSourceNode()->vnode_name));

        std::vector < tshlib::VirtualNode * > listNodesWithinPath, updatedNodesWithinPath;
        double moveLogLK = 0;

        // ------------------------------------
        // Prepare the list of nodes involved in the move (Required here!)
        //std::vector<tshlib::VirtualNode *> listNodesWithinPath = candidateMoves->getNodesInMovePath(i);
        //listNodesWithinPath = inputTree->computePathBetweenNodes(pnode, qnode);
        listNodesWithinPath = utree_->computePathBetweenNodes(candidateMoves->getMove(i)->getSourceNode(), candidateMoves->getMove(i)->getTargetNode());
        updatedNodesWithinPath = candidateMoves->updatePathBetweenNodes(i, listNodesWithinPath);

        //listNodesWithinPath = candidateMoves->getNodesInMovePath(i);

        // ------------------------------------
        // Log move details
        LOG(INFO) << "Move [" << candidateMoves->getMove(i)->getSourceNode()->getNodeName() << "->" << candidateMoves->getMove(i)->getTargetNode()->getNodeName() << "]\t"
                  << candidateMoves->getMove(i)->move_class << "\t(" << candidateMoves->getMove(i)->move_radius << ")" << "\tdirection: " << candidateMoves->getMove(i)->getMoveDirection();

        // ------------------------------------
        // Apply the move
        candidateMoves->applyMove(i);
        VLOG(1) << "[TSH Cycle - Topology] [A] MOVE#" << candidateMoves->getMove(i)->move_id << " | " << utree_->printTreeNewick(true);

        // ------------------------------------
        // Print root reachability from every nodeInterface (includes rotations)
        // inputTree->_testReachingPseudoRoot();

        // ------------------------------------
        // Print tree on file
        //inputTree->saveTreeOnFile("../data/test.txt");
        updatedNodesWithinPath.push_back(utree_->rootnode);
        listNodesWithinPath.push_back(utree_->rootnode);

        //moveLogLK = likelihoodFunc->updateLikelihood(updatedNodesWithinPath);

        if (dynamic_cast<UnifiedTSHomogeneousTreeLikelihood_PIP *>(likelihoodFunc)) {
            moveLogLK = dynamic_cast<UnifiedTSHomogeneousTreeLikelihood_PIP *>(likelihoodFunc)->updateLikelihoodOnTreeRearrangement(updatedNodesWithinPath);
        } else {
            moveLogLK = dynamic_cast<UnifiedTSHomogeneousTreeLikelihood *>(likelihoodFunc)->updateLikelihoodOnTreeRearrangement(updatedNodesWithinPath);
        }

        if (std::isinf(moveLogLK)) {
            std::ostringstream nodepath;
            for (auto &node:updatedNodesWithinPath) {
                nodepath << node->getNodeName() << ">";
            }
            LOG_IF(ERROR, std::isinf(moveLogLK)) << "llk[Move] is -inf for move " << candidateMoves->getMove(i)->move_id << " nodeInterface-path:" << nodepath.str();
        }
        // ------------------------------------
        // Store likelihood of the move
        candidateMoves->getMove(i)->move_lk = moveLogLK;

        candidateMoves->displayRearrangmentStatus(i, true);


        // ------------------------------------
        // Revert the move, and return to the original tree
        candidateMoves->revertMove(i);
        VLOG(1) << "[TSH Cycle - Topology] [R] MOVE#" << candidateMoves->getMove(i)->move_id << " | " << utree_->printTreeNewick(true);

        //listNodesWithinPath = candidateMoves->getNodesInMovePath(i);

        // ------------------------------------
        double moveLogLK_return = -std::numeric_limits<double>::infinity();
        if (scoringMethod.find("bothways") != std::string::npos) {
            //moveLogLK = likelihoodFunc->updateLikelihood(updatedNodesWithinPath);

            if (dynamic_cast<UnifiedTSHomogeneousTreeLikelihood_PIP *>(likelihoodFunc)) {
                moveLogLK_return = dynamic_cast<UnifiedTSHomogeneousTreeLikelihood_PIP *>(likelihoodFunc)->updateLikelihoodOnTreeRearrangement(listNodesWithinPath);
            } else {
                moveLogLK_return = dynamic_cast<UnifiedTSHomogeneousTreeLikelihood *>(likelihoodFunc)->updateLikelihoodOnTreeRearrangement(listNodesWithinPath);
            }
        }

        if (std::isinf(moveLogLK_return)) {
            std::ostringstream nodepath;
            for (auto &node:listNodesWithinPath) {
                nodepath << node->getNodeName() << ">";
            }
            LOG_IF(ERROR, std::isinf(moveLogLK_return)) << "llk[Return] value is -inf for move " << candidateMoves->getMove(i)->move_id << " nodeInterface-path:" << nodepath.str();
        }

        //candidateMoves->displayRearrangmentStatus(i, true);
        LOG_IF(ERROR, !ComparisonUtils::areLogicallyEqual(moveLogLK_return, tshinitScore)) << "Error in evaluating likelihood during TS @move " << candidateMoves->getMove(i)->move_id <<
                                                        "\tllk[Initial]=" << TextTools::toString(tshinitScore, 15) <<
                                                        "\tllk[Return]=" << TextTools::toString(moveLogLK_return, 15) <<
                                                        "\tllk[Move]=" << TextTools::toString(moveLogLK, 15);

        // ------------------------------------
        // Count moves performed
        performed_moves += candidateMoves->getNumberOfMoves() * 2;

    }

    ApplicationTools::displayMessage("\n");

}


double tshlib::TreeSearch::greedy() {
    std::chrono::high_resolution_clock::time_point t1 = std::chrono::high_resolution_clock::now();

    //tshcycleScore = tshinitScore;

    // Condition handler
    double cycle_no = 0;
    double c;
    switch (stopConditionMethod) {
        case tshlib::TreeSearchStopCondition::iterations:
            c = 0;
            break;
        case tshlib::TreeSearchStopCondition::convergence:
            c = tshinitScore;
            break;
    }

    while (c < stopConditionValue) {

        // Define moves according to tree-search criteria
        tshlib::TreeRearrangment *candidateMoves = defineCandidateMoves();
        std::ostringstream taskDescription;
        taskDescription << "Tree-search cycle #" << std::setfill('0') << std::setw(3) << cycle_no + 1
                        << " | Testing " << candidateMoves->getNumberOfMoves() << " tree rearrangements.";
        ApplicationTools::displayMessage(taskDescription.str());

        LOG(INFO) << "[TSH Optimisation] Cycle " << std::setfill('0') << std::setw(3) << cycle_no + 1 << " - lk= " << TextTools::toString(tshcycleScore, 15);


        // Test and record likelihood of each and every candidate move
        testCandidateMoves(candidateMoves);

        // Select the best move in the list and store it
        Move *bestMove = candidateMoves->selectBestMove(tshcycleScore);

        // ------------------------------------
        // Print winning move
        if (bestMove) {

            LOG(INFO) << "[TSH Cycle - Topology] (cycle" << std::setfill('0') << std::setw(3) << cycle_no + 1 << ")\t"
                      << "lk = " << std::setprecision(12) << bestMove->getLikelihood() << " | "
                      << bestMove->move_class << "." << std::setfill('0') << std::setw(3) << bestMove->move_id
                      << " [" << bestMove->getSourceNode()->getNodeName() << " -> " << bestMove->getTargetNode()->getNodeName() << "]";


            taskDescription.str(std::string());
            taskDescription << "Tree-search cycle #" << std::setfill('0') << std::setw(3) << cycle_no + 1 << " | ["
                            << bestMove->move_class << "." << std::setfill('0') << std::setw(3) << bestMove->move_id
                            << "] New likelihood: " << std::setprecision(12) << bestMove->getLikelihood();
            ApplicationTools::displayMessage(taskDescription.str());


            std::vector < tshlib::VirtualNode * > listNodesWithinPath, updatedNodesWithinPath;
            listNodesWithinPath = utree_->computePathBetweenNodes(bestMove->getSourceNode(), bestMove->getTargetNode());
            updatedNodesWithinPath = candidateMoves->updatePathBetweenNodes(bestMove->move_id, listNodesWithinPath);
            updatedNodesWithinPath.push_back(utree_->rootnode);

            // Commit final move on the topology
            candidateMoves->commitMove(bestMove->move_id);
            DVLOG(1) << "utree after commit " << utree_->printTreeNewick(true);

            ApplicationTools::displayTask("Optimising " + TextTools::toString(updatedNodesWithinPath.size()) + " branches");
            //likelihoodFunc->topologyChange(updatedNodesWithinPath, utree);

            if (dynamic_cast<UnifiedTSHomogeneousTreeLikelihood_PIP *>(likelihoodFunc)) {
                dynamic_cast<UnifiedTSHomogeneousTreeLikelihood_PIP *>(likelihoodFunc)->topologyChangeSuccessful(updatedNodesWithinPath);
            } else {
                dynamic_cast<UnifiedTSHomogeneousTreeLikelihood *>(likelihoodFunc)->topologyChangeSuccessful(updatedNodesWithinPath);
            }

            ApplicationTools::displayTaskDone();

            bpp::Newick treeWriter;
            //bpp::TreeTemplate<Node> ttree(likelihoodFunc->getLikelihoodFunction()->getTree());
            bpp::TreeTemplate <Node> ttree(likelihoodFunc->getTree());
            std::ostringstream oss;
            treeWriter.write(ttree, oss);

            DVLOG(1) << "bpp after commit " << oss.str();

            tshcycleScore = -likelihoodFunc->getValue();//bestMove->getLikelihood();
            tshinitScore = tshcycleScore;
            cycle_no++;

            // ------------------------------------
            // Clean memory
            delete candidateMoves;

        } else {

            LOG(INFO) << "[TSH Cycle] No further likelihood improvement after " << cycle_no << " cycles. Exit loop. ";
            ApplicationTools::displayWarning("No further likelihood improvement. A cycle with a greedier approach will be performed to overcome the local optimum.");
            // ------------------------------------
            // Clean memory
            delete candidateMoves;

            break;
        }

        switch (stopConditionMethod) {
            case tshlib::TreeSearchStopCondition::iterations:
                c++;
                break;
            case tshlib::TreeSearchStopCondition::convergence:
                c = c - tshcycleScore;
                break;
        }
    }

    if (c < stopConditionValue) {

        LOG(INFO) << "[TSH Cycle] Reached Stop-Condition boundary at " << c;
    }

    // ------------------------------------
    std::chrono::high_resolution_clock::time_point t2 = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::microseconds>(t2 - t1).count();

    VLOG(1) << "[TSH Cycle] Moves applied and reverted: " << performed_moves << std::endl;
    VLOG(1) << "[TSH Cycle] Elapsed time: " << duration << " microseconds" << std::endl;
    VLOG(1) << "[TSH Cycle] *** " << (double) duration / performed_moves << " microseconds/move *** " << std::endl;

    return -tshcycleScore;
}


double tshlib::TreeSearch::hillclimbing() {
    return greedy();
}

double tshlib::TreeSearch::particleswarming() {
    return 0;
}
