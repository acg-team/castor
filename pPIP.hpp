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
#ifndef MINIJATI_PPIP_HPP
#define MINIJATI_PPIP_HPP

#include <Bpp/Numeric/VectorTools.h>
#include <Bpp/Phyl/Node.h>
#include <random>
#include <Utree.hpp>

#include "Utilities.hpp"

#define SMALL_DOUBLE 1e-8
#define LEFT 0
#define RIGHT 1

namespace bpp {

    struct max_val_str
    {
        double val;
        int index;
    };


    class pPIP {

    public:

        // TODO: is it really necessary to redefine std::string? and vectors of strings?
        typedef std::string MSAcolumn_t;           // MSA column type
        typedef std::vector<MSAcolumn_t> MSA_t;    // MSA as vector of columns
        typedef std::vector<MSA_t> MSAensemble_t;  // for SB: set (vector) of MSAs
        typedef std::string TracebackPath_t;       // traceback path type
        typedef std::vector<TracebackPath_t> TracebackEnsemble_t; //for SB: set (vector) of traceback paths

        pPIP(tshlib::Utree *utree,              // tshlib:: tree
             bpp::Tree *tree,                   // bpp::tree
             bpp::SubstitutionModel *smodel,    // extended substitution model
             UtreeBppUtils::treemap &inTreeMap, // bpp::Node * <-> tshlib::VirtualNode *
             bpp::SequenceContainer *sequences, // un-aligned input sequences
             bpp::DiscreteDistribution *rDist,  // distribution for rate variation among sites
             long seed);                        // seed for the random numbers generation

        ~pPIP(){};

        void PIPAligner(std::vector<tshlib::VirtualNode *> &list_vnode_to_root,
                        bool local,
                        bool flag_RAM,
                        bool flag_map,
                        bool flag_pattern,
                        bool flag_fv);

        std::vector< std::string > getMSA(bpp::Node *node);
        double getScore(bpp::Node *node);
        std::vector< std::string > getSeqnames(bpp::Node *node);
        bpp::Node *getRootNode();

        const Alphabet *getAlphabet() const;

        void setSubstModel(bpp::SubstitutionModel *smodel);

        void setTree(const Tree *tree);

        void setFVleaf(bpp::Node *node);

        void set_lk_leaf(bpp::Node *node);

        void set_lk_empty_leaf(bpp::Node *node);

        void setRevMapComprSeqsleaf(bpp::Node *node);

    protected:

    private:

        tshlib::Utree *utree_;                   // tshlib:: tree
        bpp::TreeTemplate<bpp::Node> *tree_;     // bpp::tree
        bpp::SubstitutionModel *substModel_;
        long seed_;                              //jatiapp seed for the random numbers generation

    private:
        // extended substitution model
        mutable UtreeBppUtils::treemap treemap_; // bpp::Node * <-> tshlib::VirtualNode *
        bpp::SequenceContainer *sequences_;      // un-aligned input sequences
        bpp::DiscreteDistribution *rDist_;       // distribution for rate variation among sites

        std::map<unsigned long, std::vector<double>> iotasNode_; //map of nodeIDs and vector of iotas (1 for each rate (Gamma,...) category
        std::map<unsigned long, std::vector<double>> betasNode_; //map of nodeIDs and vector of betas (1 for each rate (Gamma,...) category
        std::map<unsigned long, std::vector<bpp::RowMatrix<double> > >prNode_; // map of NodeIDs of Pr = exp(branchLength * rate * Q), rate under Gamma distribution
        std::vector<std::vector<std::string> > seqNames_;          // vector[nodeId] of sequence names (MSAs seq. names at each internal node) node
        std::vector<MSA_t> MSA_;                                   // vector[nodeId] MSA at each node
        std::vector<MSAensemble_t> MSAensemble_;                   // MSA ensemble at each node (for SB)
        std::vector<TracebackPath_t> traceback_path_;              // vector[nodeId] of traceback paths (1 at each internal node)
        std::vector<TracebackEnsemble_t> traceback_path_ensemble_; // traceback path ensemble at each internal node (for SB)
        std::vector<double> score_;                                // vector[nodeId] of likelihood score
        std::vector<vector<double >> score_ensemble_;              // set of likelihoods at each internal node (for SB)
        double lambda0;                                            // original lambda (no Gamma distribution)
        double mu0;                                                // original mu (no Gamma distribution)
        std::vector<double> lambda_;                               // vector[rate] of lambda rate with Gamma distribution
        std::vector<double> mu_;                                   // vector[rate] of mu rate with Gamma distribution
        std::vector<double> nu_;                                   // vector[rate] of nu (normalizing constant) with Gamma distribution
        double tau_;                                               // total tree length

        std::vector< vector<double> > lk_down_;                      //each node a vector of lk
        std::vector< vector<double> > lk_empty_down_;                //each node a vector of lk_empty (for each gamma category)
        std::vector< vector< vector< bpp::ColMatrix<double> > > > fv_data_; // [node][column][catg][fv]
        std::vector< vector< bpp::ColMatrix<double> > > fv_empty_data_; // [node][catg][fv]
        std::vector< vector<int> > rev_map_compressed_seqs_; // [node][idx]

        bpp::ColMatrix<double> pi_;                                // steady state base frequencies

        const bpp::Alphabet *alphabet_;                            // extended alphabet (alphabet U {'-'}

        long alphabetSize_;                                        // original alphabet size

        long extendedAlphabetSize_;                                // extended alphabet size

        void _reserve(std::vector<tshlib::VirtualNode *> &nodeList);

        void _setNu();

        void _setTree(const Tree *tree);

        void _setLambda(double lambda);

        void _setMu(double mu);

        void _setPi(const Vdouble &pi);

        double _setTauRecursive(tshlib::VirtualNode *vnode);

        void _setTau(tshlib::VirtualNode *vnode);

        void _setAllIotas(bpp::Node *node,bool local_root);

        void _setAllBetas(bpp::Node *node,bool local_root);

        void _getPrFromSubstutionModel(std::vector<tshlib::VirtualNode *> &listNodes);

        bool is_inside(unsigned long x0,
                       unsigned long y0,
                       unsigned long xf,
                       unsigned long yf,
                       unsigned long xt,
                       unsigned long yt);

        void set_indeces_M(unsigned long &up_corner_i,
                           unsigned long &up_corner_j,
                           unsigned long &bot_corner_i,
                           unsigned long &bot_corner_j,
                           unsigned long level,
                           unsigned long h,
                           unsigned long w);

        void set_indeces_X(unsigned long &up_corner_i,
                           unsigned long &up_corner_j,
                           unsigned long &bot_corner_i,
                           unsigned long &bot_corner_j,
                           unsigned long level,
                           unsigned long h,
                           unsigned long w);

        void set_indeces_Y(unsigned long &up_corner_i,
                           unsigned long &up_corner_j,
                           unsigned long &bot_corner_i,
                           unsigned long &bot_corner_j,
                           unsigned long level,
                           unsigned long h,
                           unsigned long w);

        signed long get_indices_M(unsigned long nx,
                                  unsigned long ny,
                                  unsigned long up_corner_i,
                                  unsigned long up_corner_j,
                                  unsigned long bot_corner_i,
                                  unsigned long bot_corner_j,
                                  unsigned long m,
                                  unsigned long h,
                                  unsigned long w);

        signed long get_indices_X(unsigned long nx,
                                  unsigned long ny,
                                  unsigned long up_corner_i,
                                  unsigned long up_corner_j,
                                  unsigned long bot_corner_i,
                                  unsigned long bot_corner_j,
                                  unsigned long m,
                                  unsigned long h,
                                  unsigned long w);

        signed long get_indices_Y(unsigned long nx,
                                  unsigned long ny,
                                  unsigned long up_corner_i,
                                  unsigned long up_corner_j,
                                  unsigned long bot_corner_i,
                                  unsigned long bot_corner_j,
                                  unsigned long m,
                                  unsigned long h,
                                  unsigned long w);

        void set_indeces_T(unsigned long &up_corner_i,
                           unsigned long &up_corner_j,
                           unsigned long &bot_corner_i,
                           unsigned long &bot_corner_j,
                           unsigned long level,
                           unsigned long h,
                           unsigned long w);

        void reset_corner(unsigned long &up_corner_i,
                          unsigned long &up_corner_j,
                          unsigned long &bot_corner_i,
                          unsigned long &bot_corner_j,
                          unsigned long h,
                          unsigned long w);

        unsigned long get_indices_T(unsigned long nx,
                                    unsigned long ny,
                                    unsigned long up_corner_i,
                                    unsigned long up_corner_j,
                                    unsigned long bot_corner_i,
                                    unsigned long bot_corner_j,
                                    unsigned long m,
                                    unsigned long h,
                                    unsigned long w);

        bool index_of_max(double m,
                         double x,
                         double y,
                         double epsilon,
                         std::default_random_engine &generator,
                         std::uniform_real_distribution<double> &distribution,
                         max_val_str &max_val,
                         bool flag_RAM);

        double max_of_three(double a,
                            double b,
                            double c,
                            double epsilon,
                            bool flag_RAM);

        bool checkboundary(unsigned long up_corner_i,
                           unsigned long up_corner_j,
                           unsigned long bot_corner_i,
                           unsigned long bot_corner_j,
                           unsigned long h,
                           unsigned long w);

        std::string createGapCol(unsigned long len);

        void build_MSA(bpp::Node *node, std::string traceback_path);

        void setMSAsequenceNames(bpp::Node *node);

        void setMSAsequenceNames(bpp::Node *node,std::string seqname);

        void setMSAleaves(bpp::Node *node,const std::string &sequence);

        bpp::ColMatrix<double> fv_observed(MSAcolumn_t &s, unsigned long &idx);

        bpp::ColMatrix<double> computeFVrec(bpp::Node *node, MSAcolumn_t &s, unsigned long &idx, int catg);

        void allgaps(bpp::Node *node,MSAcolumn_t &s, unsigned long &idx,bool &flag);

        double compute_lk_gap_down(bpp::Node *node,MSAcolumn_t &s,int catg);

        std::vector<double> computeLK_GapColumn_local(bpp::Node *node,
                                                      MSAcolumn_t &sL,
                                                      MSAcolumn_t &sR,
                                                      bool flag_RAM);

        double compute_lk_down(bpp::Node *node,MSAcolumn_t &s,int catg);

        double computeLK_MXY_local(double NU,
                                       double valM,
                                       double valX,
                                       double valY,
                                       double log_pr,
                                       unsigned long m);

        double computeLK_M_local(double NU,
                                 double valM,
                                 double valX,
                                 double valY,
                                 bpp::Node *node,
                                 std::string &sL,
                                 std::string &sR,
                                 unsigned long m,
                                 std::map<MSAcolumn_t, double> &lkM,
                                 std::vector< std::vector<double> > &lkM_pattern,
                                 bool flag_map,
                                 bool flag_RAM,
                                 bool flag_pattern);

        double computeLK_M_local(int nodeID,
                                 int sonLeftID,
                                 int sonRightID,
                                 MSAcolumn_t &sL,
                                 MSAcolumn_t &sR,
                                 std::vector< bpp::ColMatrix<double> > &fvL,
                                 std::vector< bpp::ColMatrix<double> > &fvR,
                                 std::vector< bpp::ColMatrix<double> > &Fv_M_ij);

        double computeLK_X_local(double NU,
                                 double valM,
                                 double valX,
                                 double valY,
                                 bpp::Node *node,
                                 MSAcolumn_t &sL,
                                 MSAcolumn_t &col_gap_R,
                                 unsigned long m,
                                 std::map<MSAcolumn_t, double> &lkX,
                                 std::vector< std::vector<double> > &lkX_pattern,
                                 bool flag_map,
                                 bool flag_RAM,
                                 bool flag_pattern);

        double computeLK_X_local(int nodeID,
                                 int sonLeftID,
                                 int sonRightID,
                                 MSAcolumn_t &sL,
                                 MSAcolumn_t &col_gap_R,
                                 std::vector< bpp::ColMatrix<double> > &fvL,
                                 std::vector< bpp::ColMatrix<double> > &fvR,
                                 std::vector< bpp::ColMatrix<double> > &Fv_X_ij);

        double computeLK_Y_local(double NU,
                                 double valM,
                                 double valX,
                                 double valY,
                                 bpp::Node *node,
                                 MSAcolumn_t &col_gap_L,
                                 MSAcolumn_t &sR,
                                 unsigned long m,
                                 std::map<MSAcolumn_t, double> &lkY,
                                 std::vector< std::vector<double> > &lkY_pattern,
                                 bool flag_map,
                                 bool flag_RAM,
                                 bool flag_pattern);


        double computeLK_Y_local(int nodeID,
                                 int sonLeftID,
                                 int sonRightID,
                                 MSAcolumn_t &col_gap_L,
                                 MSAcolumn_t &sR,
                                 std::vector< bpp::ColMatrix<double> > &fvL,
                                 std::vector< bpp::ColMatrix<double> > &fvR,
                                 std::vector< bpp::ColMatrix<double> > &Fv_Y_ij);

        void DP3D_PIP(bpp::Node *node, bool local,bool flag_map);

        void DP3D_PIP_RAM(bpp::Node *node,
                          bool local,
                          bool flag_map,
                          bool flag_pattern);

        void DP3D_PIP_RAM_FAST(bpp::Node *node);

        void DP3D_PIP_SB(bpp::Node *node,UtreeBppUtils::treemap *tm,double gamma_rate, bool local,
                         double temperature,int num_SB);

    };

}

namespace pPIPUtils {

    // convert MSA PIP into sites container
    bpp::SiteContainer *pPIPmsa2Sites(bpp::pPIP *progressivePIP);

    // sum of logs
    double add_lns(double a_ln,double b_ln);

    void max_val_in_column(double ***M,int depth, int height, int width, double &val, int &level);

    std::vector<std::string> siteContainer_2_sequence_vector(std::vector<bpp::pPIP::MSAcolumn_t> &MSA);

    std::vector<int> reverse_map(std::vector<int> &m);

}

#endif //MINIJATI_PPIP_HPP
