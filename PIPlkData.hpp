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
 * @file PIPlkData.hpp
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

#ifndef MINIJATI_PIPLKDATA_HPP
#define MINIJATI_PIPLKDATA_HPP

struct LKdata {

    int d_;
    int h_;
    int h_compr_;
    int w_;
    int w_compr_;
    int numCatg_;

    std::vector<vector<vector<double> > > Log3DM;   // DP sparse matrix for MATCH case (only 2 layer are needed)
    std::vector<vector<vector<double> > > Log3DX;   // DP sparse matrix for GAPX case (only 2 layer are needed)
    std::vector<vector<vector<double> > > Log3DY;   // DP sparse matrix for GAPY case (only 2 layer are needed)

    std::vector<vector<vector<int> > > TR;   // DP matrix for tracebacking

    std::vector<vector<double> > Log2DM;
    std::vector<double> Log2DX;
    std::vector<double> Log2DY;
    std::vector<vector<double> > Log2DM_fp;
    std::vector<double> Log2DX_fp;
    std::vector<double> Log2DY_fp;
    std::vector<vector<vector<bpp::ColMatrix<double> > > > Fv_M;
    std::vector<vector<bpp::ColMatrix<double> > > Fv_X;
    std::vector<vector<bpp::ColMatrix<double> > > Fv_Y;
    std::vector<vector<vector<double> > > Fv_sigma_M;
    std::vector<vector<double> > Fv_sigma_X;
    std::vector<vector<double> > Fv_sigma_Y;

    LKdata(int d,int h,int h_compr,int w,int w_compr,int numCatg,bool sparse) {

        d_ = d;
        h_ = h;
        h_compr_ = h_compr;
        w_ = w;
        w_compr_ = w_compr;
        numCatg_ = numCatg;

        int i,j,k;
        double min_inf = -std::numeric_limits<double>::infinity(); // -inf
        //*************************
        int depth;
        if(sparse) {
            depth=2;
        }else{
            depth=d;
        }
        Log3DM.resize(depth);
        Log3DX.resize(depth);
        Log3DY.resize(depth);
        // allocate memory for the 2 layers
        for (k = 0; k < depth; k++) {
            Log3DM[k].resize(h);
            Log3DX[k].resize(h);
            Log3DY[k].resize(h);
            for (i = 0; i < h; i++) {
                Log3DM[k][i].resize(w, min_inf);
                Log3DX[k][i].resize(w, min_inf);
                Log3DY[k][i].resize(w, min_inf);
            }
        }
        //*************************
        Log2DM.resize(h_compr);
        Log2DX.resize(h_compr);
        Log2DY.resize(w_compr);

        Log2DM_fp.resize(h_compr);
        Log2DX_fp.resize(h_compr);
        Log2DY_fp.resize(w_compr);

        Fv_M.resize(h_compr);
        Fv_X.resize(h_compr);
        Fv_Y.resize(w_compr);

        for(i = 0; i < h_compr; i++){
            Log2DM[i].resize(w_compr);
            Log2DM_fp[i].resize(w_compr);
            Fv_M[i].resize(w_compr);
            for(j = 0; j < w_compr; j++){
                Fv_M[i][j].resize(numCatg);
            }
        }
        for(i = 0; i < h_compr; i++){
            Fv_X[i].resize(numCatg);
        }
        for(j = 0; j < w_compr; j++){
            Fv_Y[j].resize(numCatg);
        }

        Fv_sigma_M.resize(h_compr);
        Fv_sigma_X.resize(h_compr);
        Fv_sigma_Y.resize(w_compr);

        for(i = 0; i < h_compr; i++){
            Fv_sigma_M[i].resize(w_compr);
            for(j = 0; j < w_compr; j++){
                Fv_sigma_M[i][j].resize(numCatg);
            }
        }
        for(i = 0; i < h_compr; i++){
            Fv_sigma_X[i].resize(numCatg);
        }
        for(j = 0; j < w_compr; j++){
            Fv_sigma_Y[j].resize(numCatg);
        }


        if(sparse){
            TR.resize(d);
            TR[0].resize(1);
            for (i = 1; i < d; i++) {
                TR[i].resize(h);
                for(j = 0; j < h; j++){
                    TR[i][j].resize(w,0);
                }
            }
        }

    }

};

#endif //MINIJATI_PIPLKDATA_HPP
