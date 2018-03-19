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
#include <glog/logging.h>
#include <Bpp/Phyl/Likelihood/RHomogeneousTreeLikelihood.h>
#include <random>
#include <Bpp/Phyl/Io/Newick.h>
#include "TSHTopologySearch.hpp"
#include "Utilities.hpp"
#include "RHomogeneousTreeLikelihood_PIP.hpp"


using namespace bpp;


void TSHTopologySearch::notifyAllPerformed(const TopologyChangeEvent &event) {
    searchableTree_->topologyChangePerformed(event);
    for (size_t i = 0; i < topoListeners_.size(); i++) {
        topoListeners_[i]->topologyChangePerformed(event);
    }
}

void TSHTopologySearch::notifyAllTested(const TopologyChangeEvent &event) {
    searchableTree_->topologyChangeTested(event);
    for (size_t i = 0; i < topoListeners_.size(); i++) {
        topoListeners_[i]->topologyChangeTested(event);
    }
}

void TSHTopologySearch::notifyAllSuccessful(const TopologyChangeEvent &event) {
    searchableTree_->topologyChangeSuccessful(event);
    for (size_t i = 0; i < topoListeners_.size(); i++) {
        topoListeners_[i]->topologyChangeSuccessful(event);
    }
}

void TSHTopologySearch::search() throw(Exception) {

}


tshlib::TreeRearrangment *tshlib::TreeSearch::defineCandidateMoves(tshlib::Utree *inputTree) {

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
            max_radius = inputTree->getMaxNodeDistance() / 2;
            break;

        case tshlib::TreeRearrangmentOperations::classic_TBR:
            min_radius = 5;
            max_radius = inputTree->getMaxNodeDistance() / 2;
            break;

        case tshlib::TreeRearrangmentOperations::classic_Mixed:
            min_radius = 3;  // Minimum radius for an NNI move is 3 nodes
            max_radius = inputTree->getMaxNodeDistance(); // Full tree traversing from any node of the tree
            break;

    }

    // Initialise a new rearrangement list
    auto candidateMoveSet = new tshlib::TreeRearrangment;
    candidateMoveSet->setTreeTopology(inputTree);
    candidateMoveSet->setMinRadius(min_radius);
    candidateMoveSet->setMaxRadius(max_radius);


    switch (tshStrategy) {
        case tshlib::TreeSearchHeuristics::greedy:
            // Generate candidate list of possible moves given the tree topology and the rearrangement operation type
            for (auto &node:inputTree->listVNodes) {
                // Print node description with neighbors
                //VLOG(2) << "[utree neighbours] " << vnode->printNeighbours() << std::endl;

                candidateMoveSet->setSourceNode(node);
                // Get all the target nodes with distance == radius from the source node
                // excluding the starting node (false)
                // do not duplicate moves in the list
                candidateMoveSet->defineMoves(false, false);

            }
            break;

        case tshlib::TreeSearchHeuristics::hillclimbing:

            // fill array with [min value, max_value] number of nodes in the tree
            std::vector<int> node_ids(inputTree->listVNodes.size());
            std::iota(node_ids.begin(), node_ids.end(), 0);

            // shuffle the array
            std::random_device rd;
            std::mt19937 e{rd()};
            std::shuffle(node_ids.begin(), node_ids.end(), e);

            // Generate candidate list of possible moves given the node topology and the rearrangement operation type
            for (int i = 0; i < search_startingnodes; i++) {

                VirtualNode *node = inputTree->listVNodes.at(node_ids.at(i));

                // Print node description with neighbors
                //VLOG(2) << "[utree neighbours] " << vnode->printNeighbours() << std::endl;

                candidateMoveSet->setSourceNode(node);
                // Get all the target nodes with distance == radius from the source node
                // excluding the starting node (false)
                // do not duplicate moves in the list
                candidateMoveSet->defineMoves(false, false);

            }

            break;
    }

    // Print the list of moves for the current P node (source node)
    //rearrangmentList.printMoves();

    return candidateMoveSet;
}


double tshlib::TreeSearch::performTreeSearch(tshlib::Utree *inputTree) {
    double newScore = 0;

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

    // define the high-level strategy to evaluate the tree rearrangement operations
    std::string strategy;
    switch (tshStrategy) {
        case tshlib::TreeSearchHeuristics::greedy:
            strategy = "Greedy";
            LOG(INFO) << "[TSH Optimisation] algorithm: " << strategy << ", moves: " << operations << ", " << startingnodes << ", stop condition: " << stopcondition;
            newScore = greedy(inputTree);
            break;

        case tshlib::TreeSearchHeuristics::hillclimbing:
            strategy = "Hillclimbing";
            LOG(INFO) << "[TSH Optimisation] algorithm: " << strategy << ", moves: " << operations << ", " << startingnodes << ", stop condition: " << stopcondition;
            newScore = hillclimbing(inputTree);
            break;

        case tshlib::TreeSearchHeuristics::particle_swarm:
            strategy = "Particle Swarm";
            LOG(INFO) << "[TSH Optimisation] algorithm: " << strategy << ", moves: " << operations << ", " << startingnodes << ", stop condition: " << stopcondition;
            newScore = particleswarming(inputTree);
            break;

        case tshlib::TreeSearchHeuristics::nosearch:
            LOG(INFO) << "[Tree search] No tree search optimisation!";
            break;
        default:
            LOG(FATAL) << "Tree-search heuristic non implemented. Execution aborted!";
            break;
    }

    return newScore;
}


void tshlib::TreeSearch::testCandidateMoves(tshlib::TreeRearrangment *candidateMoves, tshlib::Utree *inputTree) {


    for (unsigned long i = 0; i < candidateMoves->getNumberOfMoves(); i++) {
        ApplicationTools::displayGauge(i + 1, candidateMoves->getNumberOfMoves(), '>', std::string("node " + candidateMoves->getMove(i)->getSourceNode()->vnode_name));


        double moveLogLK = 0;

        // Get source node and target node references
        VirtualNode *pnode = candidateMoves->getMove(i)->getSourceNode();
        VirtualNode *qnode = candidateMoves->getMove(i)->getTargetNode();

        // ------------------------------------
        // Prepare the list of nodes involved in the move (Required here!)
        std::vector<tshlib::VirtualNode *> listNodesWithinPath;
        listNodesWithinPath = inputTree->computePathBetweenNodes(pnode, qnode);
        listNodesWithinPath.push_back(inputTree->rootnode);

        // ------------------------------------
        // Log move details
        LOG(INFO) << "Move [" << candidateMoves->getMove(i)->getSourceNode()->getNodeName() << "->" << candidateMoves->getMove(i)->getTargetNode()->getNodeName() << "]\t"
                  << candidateMoves->getMove(i)->move_class << "\t(" << candidateMoves->getMove(i)->move_radius << ")" << "\tdirection: " << candidateMoves->getMove(i)->getMoveDirection();

        // ------------------------------------
        // Apply the move
        candidateMoves->applyMove(i);

        // ------------------------------------
        // Print root reachability from every node (includes rotations)
        //utree->_testReachingPseudoRoot();

        // ------------------------------------
        // Print tree on file
        //utree->saveTreeOnFile("../data/test.txt");

        // ------------------------------------
        // Add the root
        //inputTree->addVirtualRootNode();

        moveLogLK = likelihoodFunc->updateLikelihood(listNodesWithinPath);

        if (isinf(moveLogLK)) {
            std::ostringstream nodepath;
            for (auto &node:listNodesWithinPath) {
                nodepath << node->getNodeName() << ">";
            }
            LOG_IF(WARNING, isinf(moveLogLK)) << "Likelihood value is -inf for move " << candidateMoves->getMove(i)->move_name << " node-path:" << nodepath.str();
        }
        // ------------------------------------
        // Store likelihood of the move
        candidateMoves->getMove(i)->move_lk = moveLogLK;

        // ------------------------------------
        // Remove virtual root
        //inputTree->removeVirtualRootNode();

        candidateMoves->displayRearrangmentStatus(i, true);

        // ------------------------------------
        // Revert the move, and return to the original tree
        candidateMoves->revertMove(i);

        // ------------------------------------
        // Add the root
        //inputTree->addVirtualRootNode();

        // ------------------------------------
        if (scoringMethod.find("bothways") != std::string::npos) {

            moveLogLK = likelihoodFunc->updateLikelihood(listNodesWithinPath);

        }

        // ------------------------------------
        // Remove virtual root
        //inputTree->removeVirtualRootNode();

        //candidateMoves->displayRearrangmentStatus(i, true);
        DLOG(INFO) << "return: " << moveLogLK;

        // ------------------------------------
        // Count moves performed
        performed_moves += candidateMoves->getNumberOfMoves() * 2;

    }

    ApplicationTools::displayMessage("\n");

}


double tshlib::TreeSearch::greedy(tshlib::Utree *inputTree) {
    std::chrono::high_resolution_clock::time_point t1 = std::chrono::high_resolution_clock::now();

    double currentBestLk = initialLikelihoodValue;
    Utree *utree = inputTree;

    // Condition handler
    double cycle_no = 0;
    double c;
    switch (stopConditionMethod) {
        case tshlib::TreeSearchStopCondition::iterations:
            c = 0;
            break;
        case tshlib::TreeSearchStopCondition::convergence:
            c = initialLikelihoodValue;
            break;
    }

    while (c < stopConditionValue) {

        // Define moves according to tree-search criteria
        tshlib::TreeRearrangment *candidateMoves = defineCandidateMoves(utree);
        std::ostringstream taskDescription;
        taskDescription << "Tree-search cycle #" << std::setfill('0') << std::setw(3) << cycle_no + 1
                        << " | Testing " << candidateMoves->getNumberOfMoves() << " tree rearrangements.";
        ApplicationTools::displayMessage(taskDescription.str());

        // Test and record likelihood of each and every candidate move
        testCandidateMoves(candidateMoves, utree);

        // Select the best move in the list and store it
        Move *bestMove = candidateMoves->selectBestMove(currentBestLk);

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


            std::vector<tshlib::VirtualNode *> listNodesWithinPath;
            listNodesWithinPath = utree->computePathBetweenNodes(bestMove->getSourceNode(), bestMove->getTargetNode());
            listNodesWithinPath.push_back(utree->rootnode);

            // Commit final move on the topology
            candidateMoves->commitMove(bestMove->move_id);
            DVLOG(1) << "utree after commit " << utree->printTreeNewick(true);

            likelihoodFunc->topologyChange(listNodesWithinPath, utree);


            bpp::Newick treeWriter;
            bpp::TreeTemplate<Node> ttree(likelihoodFunc->getLikelihoodFunction()->getTree());
            std::ostringstream oss;
            treeWriter.write(ttree, oss);

            DVLOG(1) << "bpp after commit " << oss.str();


            currentBestLk = bestMove->getLikelihood();
            cycle_no++;

            // ------------------------------------
            // Clean memory
            delete candidateMoves;

        } else {

            LOG(INFO) << "[TSH Cycle] No further likelihood improvement after " << cycle_no << " cycles. Try greedy approach to overcome the local optimum.";
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
                c = c - currentBestLk;
                break;
        }
    }

    if (c < stopConditionValue) {

        LOG(WARNING) << "[TSH Cycle] Reached Stop-Condition boundary at " << c;
    }

    // ------------------------------------
    std::chrono::high_resolution_clock::time_point t2 = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::microseconds>(t2 - t1).count();

    VLOG(1) << "[TSH Cycle] Moves applied and reverted: " << performed_moves << std::endl;
    VLOG(1) << "[TSH Cycle] Elapsed time: " << duration << " microseconds" << std::endl;
    VLOG(1) << "[TSH Cycle] *** " << (double) duration / performed_moves << " microseconds/move *** " << std::endl;

    return -currentBestLk;
}


double tshlib::TreeSearch::hillclimbing(tshlib::Utree *inputTree) {
    return greedy(inputTree);
}

double tshlib::TreeSearch::particleswarming(tshlib::Utree *inputTree) {
    return 0;
}
