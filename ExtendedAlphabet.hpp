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
 * @file ExtendedAlphabet.hpp
 * @author Lorenzo Gatti
 * @author Massimo Maiolo
 * @date 11 01 2018
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
#ifndef MINIJATI_EXTENDEDALPHABET_HPP
#define MINIJATI_EXTENDEDALPHABET_HPP


#include <Bpp/Seq/Alphabet/LetterAlphabet.h>
#include <Bpp/Seq/Alphabet/NucleicAlphabetState.h>
#include <Bpp/Seq/Alphabet/NucleicAlphabet.h>
#include <Bpp/Seq/Alphabet/ProteicAlphabetState.h>


namespace bpp {

    /**
     * @brief This alphabet is used to deal with DNA_EXTENDED sequences.
     *
     * It supports all 4 nucleotides (A, T, G and C) with their standard denomination.
     * Gaps are coded by '-', unresolved characters are coded by 'X, N, O, 0 or ?'.
     * Extensive support for generic characters (e.g. 'P', 'Y', etc.) is provided.
     */
    class DNA_EXTENDED : public NucleicAlphabet {
    public:
        /**
         * @param exclamationMarkCountsAsGap If yes, '!' characters are replaced by gaps.
         * Otherwise, they are counted as unknown characters.
         */
        DNA_EXTENDED(bool exclamationMarkCountsAsGap = false);

        DNA_EXTENDED(const DNA_EXTENDED &bia) : NucleicAlphabet(bia) {}

        DNA_EXTENDED &operator=(const DNA_EXTENDED &bia) {
            NucleicAlphabet::operator=(bia);
            return *this;
        }

        DNA_EXTENDED *clone() const {
            return new DNA_EXTENDED(*this);
        }

        virtual ~DNA_EXTENDED() {}

    public:
        std::vector<int> getAlias(int state) const throw(BadIntException);

        std::vector<std::string> getAlias(const std::string &state) const throw(BadCharException);

        int getGeneric(const std::vector<int> &states) const throw(BadIntException);

        std::string getGeneric(const std::vector<std::string> &states) const throw(BadCharException);

        std::string getAlphabetType() const { return "DNA_EXTENDED"; }

        // return 4 : A, C, G, T (or U)
        unsigned int getSize() const { return 5; }

        // return 15 : gap isn't included, generic unresolved bases (N, X, ?, O, 0) count for one
        unsigned int getNumberOfTypes() const { return 16; }

        int getUnknownCharacterCode() const { return 15; }

        bool isUnresolved(int state) const { return state > 3; }

        bool isUnresolved(const std::string &state) const { return charToInt(state) > 3; }

        int getGapCharacterCode() const { return 4; }
    };


    /**
     * @brief This alphabet is used to deal with proteins.
     *
     * It supports all 20 amino-acids with their standard denomination.
     * Gaps are coded by '-', unresolved characters are coded by 'X'.
     */

    class ProteicAlphabet_Extended : public LetterAlphabet {
        /**
         * @name Overloaded methods from AbstractAlphabet
         * @{
         */
    public:
        const ProteicAlphabetState &getState(const std::string &letter) const
        throw(BadCharException) {
            return dynamic_cast<const ProteicAlphabetState &>(
                    AbstractAlphabet::getState(letter)
            );
        }

        const ProteicAlphabetState &getState(int num) const
        throw(BadIntException) {
            return dynamic_cast<const ProteicAlphabetState &>(
                    AbstractAlphabet::getState(num)
            );
        }

    protected:

        const ProteicAlphabetState &getStateAt(size_t pos) const
        throw(IndexOutOfBoundsException) {
            return dynamic_cast<const ProteicAlphabetState &>(
                    AbstractAlphabet::getStateAt(pos)
            );
        }

        ProteicAlphabetState &getStateAt(size_t pos)
        throw(IndexOutOfBoundsException) {
            return dynamic_cast<ProteicAlphabetState &>(
                    AbstractAlphabet::getStateAt(pos)
            );
        }

        /** @} */
    public:
        ProteicAlphabet_Extended();

        ProteicAlphabet_Extended(const ProteicAlphabet_Extended &bia) : LetterAlphabet(bia) {}

        ProteicAlphabet_Extended &operator=(const ProteicAlphabet_Extended &bia) {
            LetterAlphabet::operator=(bia);
            return *this;
        }

        ProteicAlphabet_Extended *clone() const {
            return new ProteicAlphabet_Extended(*this);
        }


        virtual ~ProteicAlphabet_Extended() {}


    public:
        unsigned int getSize() const { return 21; }

        unsigned int getNumberOfTypes() const { return 23; }

        int getUnknownCharacterCode() const { return 22; }

        std::vector<int> getAlias(int state) const throw(BadIntException);

        std::vector<std::string> getAlias(const std::string &state) const throw(BadCharException);

        int getGeneric(const std::vector<int> &states) const throw(BadIntException);

        std::string getGeneric(const std::vector<std::string> &states) const throw(BadCharException);

        bool isUnresolved(int state) const { return state > 19; }

        bool isUnresolved(const std::string &state) const { return charToInt(state) > 19; }

        std::string getAlphabetType() const { return "Proteic"; }

    public:

        /**
         * @name Specific methods
         *
         * @{
         */

        /**
         * @brief Get the abbreviation (3 letter code) for a state coded as char.
         *
         * @param aa Char description of the amino-acid to analyse.
         */
        std::string getAbbr(const std::string &aa) const throw(AlphabetException);

        /**
         * @brief Get the abbreviation (3 letter code) for a state coded as int.
         *
         * @param aa Int description of the amino-acid to analyse.
         */
        std::string getAbbr(int aa) const throw(AlphabetException);
        /** @} */

    };


}

#endif //MINIJATI_EXTENDEDALPHABET_HPP
