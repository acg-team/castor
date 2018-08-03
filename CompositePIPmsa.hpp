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
 * @see For more information visit:
 */

#ifndef MINIJATI_COMPOSITEPIPMSA_HPP
#define MINIJATI_COMPOSITEPIPMSA_HPP
#include <string>
#include <vector>
#include <Bpp/Numeric/Matrix/Matrix.h>
#include <Bpp/Seq/Sequence.h>
#include <Bpp/Seq/Container/SiteContainer.h>

namespace bpp {

    //*******************************************************************************************

    typedef std::string MSAcolumn_t; // MSA column type
    typedef std::vector<MSAcolumn_t> MSA_t; // MSA as vector of columns

    //*******************************************************************************************
    class PIPmsa {

    private:

    public:

        //***************************************************************************************
        // PUBLIC FIELDS
        //***************************************************************************************

        double score_;

        std::vector<std::string> seqNames_; // vector of strings (sequence names)

        MSA_t msa_;

        std::vector<int> traceback_path_;

        std::vector<int> traceback_map_;
//        std::vector<int> traceback_mapL_;
//        std::vector<int> traceback_mapR_;

        std::vector<int> map_compressed_seqs_;

        std::vector<int> rev_map_compressed_seqs_;

        std::vector<std::vector<bpp::ColMatrix<double>>> fv_data_; // [site][catg][fv]

        std::vector<std::vector<double> > fv_sigma_; // [site][catg]

        std::vector<bpp::ColMatrix<double> > fv_empty_data_; // [catg][fv]

        std::vector<double> fv_empty_sigma_; // [catg]

        std::vector<double> log_lk_down_; //each node a vector of lk

        std::vector<double> lk_empty_; //each node a vector of lk_empty (for each gamma category)

        //***************************************************************************************
        // PUBLIC METHODS
        //***************************************************************************************

        PIPmsa() { score_ = -std::numeric_limits<double>::infinity(); };

        void _setSeqNameLeaf(std::string &seqName);

        void _setMSAleaf(const bpp::Sequence *sequence);

        void _setSeqNameNode(std::vector<std::string> &seqNamesL,
                             std::vector<std::string> &seqNamesR);

        int getNumSites();

        void _setTracebackPathleaves();

        std::vector<std::string> getSequenceNames() { return seqNames_; };

        void _setFVleaf(int numCatg, const bpp::Alphabet *alphabet);

        void _setFVsigmaLeaf(int lenComprSeqs,
                             int numCatg,
                             const bpp::ColMatrix<double> &pi);

        void _setFVemptyNode(int numCatg,
                             PIPmsa *childL,
                             PIPmsa *childR,
                             std::vector<bpp::RowMatrix<double> > &PrL,
                             std::vector<bpp::RowMatrix<double> > &PrR);

        void _setFVsigmaEmptyLeaf(int numCatg);

        void _setFVsigmaEmptyNode(int numCatg,
                                  PIPmsa *childL,
                                  PIPmsa *childR,
                                  double bL,
                                  double bR,
                                  const std::vector<double> &mu);

        void _setFVemptyLeaf(int numCatg, const bpp::Alphabet *alphabet);

        void _compress_Fv(std::vector<std::vector<double>> &fv_sigma_not_compressed,
                          std::vector<std::vector<bpp::ColMatrix<double> > > &fv_data_not_compressed);

    };

    //*******************************************************************************************
    class iPIPmsa {

    public:

        //***************************************************************************************
        // PUBLIC METHODS
        //***************************************************************************************

        iPIPmsa() {}; // constructor

        virtual ~iPIPmsa() {}; // destructor

//        void _setMSAleaves(bpp::Sequence &sequence,
//                           std::string &seqName);

        virtual void add(PIPmsa *msa) {};

//        virtual void _setFVemptyLeaf(int numCatg,
//                                     const bpp::Alphabet *alphabet){};

        virtual int getCompressedMSAlength() {};

        virtual MSA_t *_getMSA() {};

        virtual std::vector<std::string> *_getseqNames() {};

        virtual double _getScore() {};

    };

    //*******************************************************************************************
    class PIPmsaSingle : public iPIPmsa {

    private:

        //***************************************************************************************
        // PRIVATE FIELDS
        //***************************************************************************************

    public:

        //***************************************************************************************
        // PUBLIC FIELDS
        //***************************************************************************************

        PIPmsa *pipmsa;

//        int subMSAidxL_;
//        int subMSAidxR_;

        //***************************************************************************************
        // PUBLIC METHODS
        //***************************************************************************************

        PIPmsaSingle() {
//            subMSAidxL_ = 0;
//            subMSAidxR_ = 0;
        }

        void _compressMSA(const bpp::Alphabet *alphabet);

        void _build_MSA(MSA_t &msaL, MSA_t &msaR);

        int getMSAlength();

        MSA_t *_getMSA(int idx=0) { return &(pipmsa->msa_); }

        std::vector<std::string> *_getseqNames(int idx=0) { return &(pipmsa->seqNames_); };

        double _getScore(int idx=0) { return pipmsa->score_; }

        void add(PIPmsa *x);

        void _compress_lk_components(std::vector<double> &lk_down_not_compressed,
                                     std::vector<std::vector<bpp::ColMatrix<double> > > &fv_data_not_compressed);

        void _setFVsigmaEmptyNode(int numCatg,
                                  PIPmsa *childL,
                                  PIPmsa *childR,
                                  double bL,
                                  double bR,
                                  const std::vector<double> &mu);

        int getCompressedMSAlength(int idx=0);

        ~PIPmsaSingle() {
            delete pipmsa;
        }
    };

    //*******************************************************************************************
    class PIPmsaComp : public iPIPmsa {

    private:

        //***************************************************************************************
        // PRIVATE FIELDS
        //***************************************************************************************

    public:

        //***************************************************************************************
        // PUBLIC FIELDS
        //***************************************************************************************

        std::vector<PIPmsa *> pipmsa;

//        std::vector<int> subMSAidxL_;
//        std::vector<int> subMSAidxR_;

        //***************************************************************************************
        // PUBLIC METHODS
        //***************************************************************************************

        PIPmsaComp(int size) { pipmsa.resize(size); }

        void _compressMSA(const bpp::Alphabet *alphabet, int idx_sb);

        void _build_MSA(MSA_t &msaL, MSA_t &msaR, int idx_sb);

        int getMSAlength(int idx);

        int getCompressedMSAlength(int idx);

        MSA_t *_getMSA(int idx) { return &(pipmsa.at(idx)->msa_); }

        std::vector<std::string> *_getseqNames(int idx) { return &(pipmsa.at(idx)->seqNames_); };

        double _getScore(int idx) { return pipmsa.at(idx)->score_; }

        void add(PIPmsa *x);

        void _compress_lk_components(std::vector<double> &lk_down_not_compressed,
                                     std::vector<std::vector<bpp::ColMatrix<double> > > &fv_data_not_compressed,
                                     int idx_sb);

//        void _setFVsigmaEmptyNode(int numCatg,
//                                    PIPmsa *childL,
//                                    PIPmsa *childR,
//                                    double bL,
//                                    double bR,
//                                    const std::vector<double> &mu);

//        int getCompressedMSAlength();

        ~PIPmsaComp() {
            for (std::vector<PIPmsa *>::const_iterator iter = pipmsa.begin(); iter != pipmsa.end(); ++iter) {
                delete *iter;
            }
        }
    };
    //*******************************************************************************************

}
//***********************************************************************************************
//***********************************************************************************************
//***********************************************************************************************
namespace compositePIPmsaUtils {

    std::vector<std::string> siteContainer2sequenceVector(std::vector<bpp::MSAcolumn_t> &MSA);

    std::vector<int> reverse_map(std::vector<int> &m);

    bpp::SiteContainer *pPIPmsa2Sites(const bpp::Alphabet *alphabet,
                                      std::vector<std::string> &seqNames,
                                      std::vector<std::string> &MSA);

}
//***********************************************************************************************
//***********************************************************************************************
//***********************************************************************************************

#endif //MINIJATI_COMPOSITEPIPMSA_HPP
