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
 * @file pPIP.hpp
 * @author Lorenzo Gatti
 * @author Massimo Maiolo
 * @date 19 02 2018
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

#ifndef MINIJATI_PROGRESSIVEPIP_HPP
#define MINIJATI_PROGRESSIVEPIP_HPP

#include <Bpp/Numeric/VectorTools.h>
#include <Bpp/Phyl/Node.h>
#include <random>
#include <Utree.hpp>

#include "Utilities.hpp"

#define SMALL_DOUBLE 1e-8
#define LEFT 0
#define RIGHT 1

namespace bpp {

    typedef std::string MSAcolumn_t; // MSA column type
    typedef std::vector<MSAcolumn_t> MSA_t; // MSA as vector of columns

    class progressivePIP; // forward declaration

    //*****************************************************
    // interface
    class PIPcomponent {

    public:
        virtual int getId() = 0; // pure virtual
        virtual void PIPalignNode() = 0; // pure virtual
    };
    //*****************************************************
    // concrete class
    class PIPnode : public PIPcomponent{

    private:
        //***************************************************************************************
        // FIELDS
        //***************************************************************************************
        int nodeID_;
        tshlib::VirtualNode *vnode_;
        bpp:: Node *bnode_;
        const progressivePIP *progressivePIP_;

        std::vector<double> iotasNode_; //map of nodeIDs and vector of iotas (1 for each rate (Gamma,...) category
        std::vector<double> betasNode_; //map of nodeIDs and vector of betas (1 for each rate (Gamma,...) category
        std::vector<bpp::RowMatrix<double> > prNode_; // map of NodeIDs of Pr = exp(branchLength * rate * Q), rate under Gamma distribution

        std::vector<std::string>  seqNames_; // vector[nodeId] of sequence names (MSAs seq. names at each internal nodeInterface) nodeInterface
        std::vector<MSA_t>  MSA_; // vector[nodeId] MSA at each nodeInterface
        std::vector<vector<int> >  traceback_path_; // vector[nodeId] of traceback paths (1 at each internal nodeInterface)
        std::vector< vector< vector<int> > > traceback_map_;

        std::vector<double>  score_; // vector[msas] of likelihood score

        double tau_;
        std::vector<double> nu_;

        vector<double> log_lk_down_; //each nodeInterface a vector of lk
        vector<double> log_lk_empty_down_; //each nodeInterface a vector of lk_empty (for each gamma category)

        std::vector< vector< vector< bpp::ColMatrix<double> > > > fv_data_; // [nodeInterface][site][catg][fv]
        vector< bpp::ColMatrix<double> > fv_empty_data_; // [nodeInterface][catg][fv]
        std::vector< vector<int> > map_compressed_seqs_; // [nodeInterface][idx]
        std::vector< vector<int> > rev_map_compressed_seqs_; // [nodeInterface][idx]

        std::vector< std::vector< std::vector<double> > > fv_sigma_; // [nodeInterface][site][catg]
        std::vector<double>  fv_empty_sigma_; // [catg]
        //***************************************************************************************
        // METHODS
        //***************************************************************************************

        //***************************************************************************************
    public:
        //***************************************************************************************
        // FIELDS
        //***************************************************************************************

        //***************************************************************************************
        // METHODS
        //***************************************************************************************
        PIPnode(const progressivePIP *pPIP,tshlib::VirtualNode *vnode,bpp::Node *bnode);
        int getId(){ return nodeID_; }; // concrete
        tshlib::VirtualNode *getVnode(){ return vnode_; };
        bpp:: Node *getBnode(){ return bnode_; };
        std::vector< std::vector<std::string> > getMSA();
        void _compressMSA(int idx_sb);
        void _setFVemptyLeaf();
        void _setFVsigmaLeaf();
        void _setFVsigmaEmptyLeaf();
        void _setFVleaf();
        void _setMSAsequenceNames();
        void _setMSAleaves();
        void _setTracebackPathleaves();
        void _reserve(int numCatg);
        void _computeLocalTau();
        void _computeLocalNu(int numCatg);
        void _getPrFromSubstutionModel();

        void PIPalignNode();
        //***************************************************************************************
    };
    //*****************************************************

    //*************************
    class msa {
        double score;
    public:
        double getScore(){ return score; };
    };
    //*************************
    class CompositeMSA : public msa {
        vector<msa*> elems;
    public:

        CompositeMSA(int n){
            elems.resize(n);
        }

        void Add(msa* elem,int idx) {
            elems.at(idx) = elem;
        }

        ~CompositeMSA() {
            for (vector<msa*>::const_iterator iter = elems.begin(); iter != elems.end(); ++iter) {
                delete *iter;
            }
        }
    };
    //*************************
    class nodeInterface {
        int val;
        int val1;
    public:
        virtual void DP3D_PIP() = 0;

        int getVal(){
            return val;
        };

        void setVal(int x){
            val = x;
        };

        CompositeMSA* clientCompositeMSA;

    };
    //*************************
    class nodeCPU : public nodeInterface {
        int val;
        int val1;
    public:

        nodeCPU(int x){ val=x;}

        int getVal(){
            return val;
        };

        void setVal(int x){
            val = x;
        };

        void DP3D_PIP() {
            std::cout<<"align nodeCPU: "<<val<<"\n";
            setVal(getVal()*20);
            clientCompositeMSA = new CompositeMSA(10);
        }
    };
    //*************************
    class nodeRAM : public nodeInterface {
        int val;
        int val1;
    public:
        nodeRAM(int x){ val=x;}

        int getVal(){
            return val;
        };

        void setVal(int x){
            val = x;
        };

        void DP3D_PIP() {
            std::cout<<"align nodeRAM: "<<val<<"\n";
            setVal(getVal()*2);
            clientCompositeMSA = new CompositeMSA(5);
        }
    };
    //*************************
    class CompositeInterface {
    public:
        virtual void Align(void) = 0;
        virtual void Add(nodeInterface* elem,int idx){}
    };
    //*************************
    class Composite : public CompositeInterface {
        vector<nodeInterface*> elems;
    public:

        Composite(int n){
            elems.resize(n);
        }

        void Align(void) {
            for (vector<nodeInterface*>::const_iterator iter = elems.begin(); iter != elems.end(); ++iter) {
                (*iter)->DP3D_PIP();
            }
        }

        void Add(nodeInterface* elem,int idx) {
            elems.at(idx) = elem;
        }

        ~Composite() {
            for (vector<nodeInterface*>::const_iterator iter = elems.begin(); iter != elems.end(); ++iter) {
                delete *iter;
            }
        }
    };
    //*************************
    enum genre_e{CPU,RAM};

    class nodeFactory {

    public:
        /* Factory Method */
        nodeInterface *getNode(genre_e genre,int id) {
            switch(genre) {
                case CPU:
                    return new nodeCPU(id);
                    break;
                case RAM:
                    return new nodeRAM(id);
                    break;
                default:
                    return NULL;
            }

        }
    };
    //*****************************************************
    //*****************************************************
    //*****************************************************
    //*****************************************************
    class CompositePIPaligner : public PIPcomponent{

    private:

    public:

        //***************************************************************************************
        // FIELDS
        //***************************************************************************************
        std::vector<PIPcomponent *> pip_nodes_;
        //***************************************************************************************
        // METHODS
        //***************************************************************************************
        CompositePIPaligner(int numNodes);
        void addPIPcomponent(PIPcomponent *pip_node);
        int getId(); // needed but void????
        void PIPalign();
        void PIPalignNode();
        //~CompositePIPaligner(){};
    };
    //*******************************************************************************************
    class progressivePIP {

    public:
        //***************************************************************************************
        // METHODS
        //***************************************************************************************
        progressivePIP(tshlib::Utree *utree,              // tshlib:: tree
             bpp::Tree *tree,                   // bpp::tree
             bpp::SubstitutionModel *smodel,    // extended substitution model
             UtreeBppUtils::treemap &inTreeMap, // bpp::Node * <-> tshlib::VirtualNode *
             bpp::SequenceContainer *sequences, // un-aligned input sequences
             bpp::DiscreteDistribution *rDist,  // distribution for rate variation among sites
             long seed);                         // seed for the random numbers generation


        ~progressivePIP(){};

        void initializePIP(std::vector<tshlib::VirtualNode *> &list_vnode_to_root,
                           int num_sb);

        bpp::Node *getRootNode();

        const Alphabet *getAlphabet() const;
        //***************************************************************************************
        // FIELDS
        //***************************************************************************************
        std::vector<double> lambda_; // vector[rate] of lambda rate with Gamma distribution
        std::vector<double> mu_; // vector[rate] of mu rate with Gamma distribution
        bpp::DiscreteDistribution *rDist_; // distribution for rate variation among sites
        bpp::SubstitutionModel *substModel_;
        bpp::SequenceContainer *sequences_; // un-aligned input sequences
        const bpp::Alphabet *alphabet_; // extended alphabet (alphabet U {'-'})
        long alphabetSize_; // original alphabet size
        long extendedAlphabetSize_; // extended alphabet size
        bpp::ColMatrix<double> pi_; // steady state base frequencies
        //***************************************************************************************
    protected:

    private:

        //***************************************************************************************
        // FIELDS
        //***************************************************************************************
        tshlib::Utree *utree_; // tshlib:: tree
        bpp::TreeTemplate<bpp::Node> *tree_; // bpp::tree
        long seed_; //jatiapp seed for the random numbers generation
        mutable UtreeBppUtils::treemap treemap_; // bpp::Node * <-> tshlib::VirtualNode *
        double lambda0_; // original lambda (no Gamma distribution)
        double mu0_; // original mu (no Gamma distribution)
        //std::vector<double> nu_; // vector[rate] of nu (normalizing constant) with Gamma distribution
        double tau_; // total tree length
        CompositePIPaligner *compositePIPaligner_;


        //CompositePIPalignerNEW *compositePIPalignerNEW_;
        //***************************************************************************************
        // METHODS
        //***************************************************************************************
        void _setTree(const Tree *tree);
        void _reserve(int numCatg);
        void _setLambda(double lambda);
        void _setMu(double mu);
        void _setPi(const Vdouble &pi);
        //***************************************************************************************

    };

}

namespace progressivePIPutils {

    // convert MSA PIP into sites container
    bpp::SiteContainer *pPIPmsa2Sites(bpp::progressivePIP *progressivePIP,bpp::PIPnode *pip_node,int idx_sb);
    // sum of logs
    double add_lns(double a_ln,double b_ln);

    void max_val_in_column(double ***M,int depth, int height, int width, double &val, int &level);

    std::vector<std::string> siteContainer_2_sequence_vector(std::vector<bpp::MSAcolumn_t> &MSA);

    std::vector<int> reverse_map(std::vector<int> &m);

}

#endif //MINIJATI_PROGRESSIVEPIP_HPP
