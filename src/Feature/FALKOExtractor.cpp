/**
 * FALKOLib - Fast Adaptive Laser Keypoint Orientation-invariant
 * Copyright (C) 2016 Fabjan Kallasi and Dario Lodi Rizzini.
 *
 * This file is part of FALKOLib.
 *
 * FALKOLib is free software: you can redistribute it and/or modify
 * it under the terms of the GNU Lesser General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 * 
 * FALKOLib is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public License
 * along with FALKOLib.  If not, see <http://www.gnu.org/licenses/>.
 */
#include <falkolib/Feature/FALKOExtractor.h>
#include <iostream>

using namespace std;

namespace falkolib {

    FALKOExtractor::FALKOExtractor() {
        minScoreTh = 50;
        minExtractionRange = 0;
        maxExtractionRange = 30;
        subbeam = true;
        NMSRadius = 0.1;
        neighA = 0.1;
        neighB = 0.07;
        neighMinPoint = 2;
        bRatio = 2.5;
        gridSectors = 16;
    }


    void FALKOExtractor::extract(const LaserScan& scan, std::vector<FALKO>& keypoints) {
        extract(scan, keypoints, false);
    }

    void FALKOExtractor::extract(const LaserScan& scan, std::vector<FALKO>& keypoints, bool normalize_scores) {
        const int numBeams = scan.getNumBeams();
        vector<double> scores(numBeams, -10);
        vector<double> areas(numBeams, -10);
        vector<double> base_len(numBeams, -10);
        vector<double> base_height(numBeams, -10);

            vector<pair<double,double>> right(numBeams, pair<double,double>(0,0));
        vector<pair<double,double>> left(numBeams, pair<double,double>(0,0));
        vector<Point2d> neigh;
        vector<double> radius(numBeams);
        int neighSize;
        int neighSizeL;
        int neighSizeR;
        int midIndex;
        double triangleBLength;
        double triangleHLength;
        vector<double> thetaCorner(numBeams);
        int scoreL;
        int scoreR;
        int scoreMax = 0;
        vector<int> peaks;
        vector<int> neighCircularIndexesL;
        vector<int> neighCircularIndexesR;
        int av_neigh_L = 0;
        int av_neigh_R = 0;
        int out_of_range = 0;
        int no_neighbours = 0;
        int trio_len = 0;
        for (int ind = 0; ind < numBeams; ++ind) {
            if (scan.ranges[ind] < minExtractionRange || scan.ranges[ind] > maxExtractionRange) {
                scores[ind] = -10;
                out_of_range++;
                continue;
            }
            neigh.clear();
            radius[ind] = getNeighRadius(scan.ranges[ind]);
            scan.getNeighPoints(ind, radius[ind], neigh, midIndex);
            neighSize = neigh.size();

            neighSizeL = midIndex;
            neighSizeR = neighSize - midIndex - 1;
            av_neigh_R += neighSizeR;
            av_neigh_L += neighSizeL;
            if (neighSizeL < neighMinPoint || neighSizeR < neighMinPoint) {
                scores[ind] = -10;
                no_neighbours++;
                continue;
            }

            triangleBLength = pointsDistance(neigh.front(), neigh.back());
            double area = signedTriangleArea(neigh[midIndex], neigh.front(), neigh.back());
            triangleHLength = std::abs(area) / triangleBLength;
            areas[ind] = area;
            base_len[ind] = triangleBLength;
            base_height[ind] = triangleHLength;
            right[ind] = pair<double,double>(neigh.front()[0], neigh.front()[1]);
            left[ind] = pair<double,double>(neigh.back()[0], neigh.back()[1]);

            //std::cout<<ind<< " length "<<triangleBLength<<" height "<<triangleHLength<<std::endl;
            if (triangleBLength < (radius[ind] / bRatio) || triangleHLength < (radius[ind] / bRatio)) {
                scores[ind] = -10;
                trio_len++;
                continue;
            }

            thetaCorner[ind] = getCornerOrientation(neigh, midIndex);

            neighCircularIndexesL.resize(neighSizeL, 0);
            neighCircularIndexesR.resize(neighSizeR, 0);

            for (int i = 0; i < neighSizeL; ++i) {
                neighCircularIndexesL[i] = getCircularSectorIndex(neigh[i], neigh[midIndex], thetaCorner[ind]);
            }
            for (int i = 0; i < neighSizeR; ++i) {
                neighCircularIndexesR[i] = getCircularSectorIndex(neigh[midIndex + i + 1], neigh[midIndex], thetaCorner[ind]);
            }

            scoreL = 0;
            scoreR = 0;
            uint64_t cL(0), cR(0);
            for (int i = midIndex - 1; i >= 0; --i) {
                for (int j = i; j >= 0; --j) {
                    cL++;
                    scoreL += circularSectorDistance(neighCircularIndexesL[i], neighCircularIndexesL[j], gridSectors);
                }
            }

            for (int i = midIndex + 1; i < neighSize; ++i) {
                for (int j = i; j < neighSize; ++j) {
                    cR++;
                    scoreR += circularSectorDistance(neighCircularIndexesR[i - midIndex - 1], neighCircularIndexesR[j - midIndex - 1], gridSectors);
                }
            }

            if(normalize_scores) {
                scoreL = (cL>0) ? scoreL / cL : 1000;
                scoreR = (cR>0) ? scoreR / cR : 1000;
            }

            scores[ind] = scoreL + scoreR;

            if (scores[ind] > scoreMax) {
                scoreMax = scores[ind];
            }
        }


        av_neigh_L = av_neigh_L / (numBeams-out_of_range);
        av_neigh_R = av_neigh_R / (numBeams-out_of_range);
        std::cout<< "Average Neighs "<<av_neigh_L<<", "<<av_neigh_R<<std::endl;
        std::cout<<"Beams: "<<numBeams<<" OutRange: "<<out_of_range<<" NoNeighs: "<<no_neighbours 
            <<" TriangleFactor: "<<trio_len<<" Remaining: "<<(numBeams-out_of_range-no_neighbours-trio_len)<<std::endl;

        for (int ind = 0; ind < numBeams; ++ind) {
            if (scores[ind] < 0) {
                scores[ind] = scoreMax;
            }
            scores[ind] = scoreMax - scores[ind];
        }

        NMSKeypoint(scores, scan, 0, numBeams, NMSRadius, (scoreMax * minScoreTh / 100.0), peaks);

        for (int i = 0; i < peaks.size(); ++i) {
            FALKO kp;
            kp.index = peaks[i];
            kp.orientation = thetaCorner[peaks[i]];
            kp.radius = radius[peaks[i]];
            kp.score = scores[peaks[i]];
            
            kp.area = areas[peaks[i]];
            std::cout<<"Peak "<< peaks[i] <<" area: "<< areas[peaks[i]]<<std::endl;
            kp.base_len = base_len[peaks[i]];
            kp.base_height = base_height[peaks[i]];
            kp.right_x = right[peaks[i]].first;
            kp.right_y = right[peaks[i]].second;
            kp.left_x = left[peaks[i]].first;
            kp.left_y = left[peaks[i]].second;
            
            std::cout<<"Keypoint " << scan.points[kp.index][0]<<", "<< scan.points[kp.index][1]<< "right "<<kp.right_x << ", "<<kp.right_y<< " left "<< kp.left_x<<", "<< kp.left_y << std::endl;
            double d1 = std::hypot(scan.points[kp.index][0] - kp.right_x, scan.points[kp.index][1]-kp.right_y);
            double d2 = std::hypot(scan.points[kp.index][0] - kp.left_x, scan.points[kp.index][1]-kp.left_y);
            double d3 = std::hypot(kp.left_x - kp.right_x, kp.left_y-kp.right_y);
            std::cout<<"D1 "<<d1<<" D2 "<<d2<<" D3 "<<d3<<std::endl;
            if (subbeam) {
                subBeamCorner(scan, peaks[i], radius[peaks[i]], kp.point);
            } else {
                kp.point = scan.points[peaks[i]];
            }

            keypoints.push_back(std::move(kp));
        }

    }

    int FALKOExtractor::circularSectorDistance(int a1, int a2, int res) {
        const int r2 = res / 2;
        return std::abs(((a1 - a2) + r2) % res - r2);
    }

    int FALKOExtractor::getCircularSectorIndex(const Point2d& p, const Point2d& pmid, double theta) {
        return (int) floor((angleBetweenPoints(p, pmid) + M_PI - theta) / (2.0 * M_PI / gridSectors));
    }

    double FALKOExtractor::getNeighRadius(double rho) {
        double radius = neighA * std::exp(neighB * rho);
        if (radius >= rho) {
            return rho * 0.8;
        }
        return radius;
    }

    double FALKOExtractor::getCornerOrientation(const std::vector<Point2d>& neigh, int midIndex) {
        Point2d oriL(0, 0);
        Point2d oriR(0, 0);
        int size = neigh.size();
        for (int i = 0; i < size; ++i) {
            if (i < midIndex) {
                oriL += (neigh[i] - neigh[midIndex]); // / (s[i] - s[mid]).norm();
            } else if (i > midIndex) {
                oriR += (neigh[i] - neigh[midIndex]); // / (s[i] - s[mid]).norm();
            }
        }
        oriL /= midIndex;
        oriR /= (size - (midIndex + 1));
        Point2d ori = oriL + oriR;
        double theta = atan2(ori(1), ori(0));
    }

    void FALKOExtractor::NMSKeypoint(const std::vector<double>& scores, const LaserScan& scan, unsigned int ibeg, unsigned int iend, double radius, int minval, std::vector<int>& peaks) {
        int i;
        double imax;
        unsigned j, jbeg, jend;
        double jmax;
        std::vector<int> candidates;
        peaks.clear();

        i = ibeg;
        imax = ibeg;
        while (i < iend) {
            int win;
            if (radius >= scan.ranges[i]) {
                win = std::floor(std::asin(0.8) / scan.getAngleInc());
            } else {
                win = std::floor(std::asin(radius / scan.ranges[i]) / scan.getAngleInc());
            }

            jmax = i;

            if (imax + win < i) {
                jbeg = (i >= win ? i - win : 0);
                imax = i;
            } else {
                jbeg = i;
            }

            jend = (i + win + 1 < iend ? i + win + 1 : iend);
            for (j = jbeg; j < jend; ++j) {
                if (scores[j] > scores[jmax]) {
                    jmax = j;
                }
            }
            imax = (scores[jmax] >= scores[imax] ? jmax : imax);

            if (i == imax && scores[i] > minval) {
                candidates.push_back(i);
            }
            if (jmax > i) i = jmax;
            else ++i;
        }

        int i1 = 0;
        int i2 = 0;
        int counter = 0;

        for (i1 = 0, i2 = 0; i1 < candidates.size(); ++i1) {
            if (scores[candidates[i2]] == scores[candidates[i1]]) {
                counter++;
                if (2 * abs(i2 - i1) > counter) {
                    ++i2;
                }
            } else {
                peaks.push_back(candidates[i2]);
                i2 = i1;
                counter = 0;
            }
        }
        if (i2 != candidates.size()) {
            peaks.push_back(candidates[i2]);
        }
    }

    void FALKOExtractor::subBeamCorner(const LaserScan& scan, int index, double radius, Point2d& p) {
        vector<Point2d> neigh;
        int midIndex;
        scan.getNeighPoints(index, radius, neigh, midIndex);
        int neighSize = neigh.size();

        Eigen::Vector3d leftLine;
        Eigen::Vector3d rightLine;

        std::vector<Point2d> leftSide;
        for (int i = 0; i <= midIndex; ++i) {
            leftSide.push_back(neigh[i]);
        }

        std::vector<Point2d> rightSide;
        for (int i = midIndex; i < neighSize; ++i) {
            rightSide.push_back(neigh[i]);
        }

        generateLine(leftSide, leftLine);
        generateLine(rightSide, rightLine);

        double A[4];
        double b[2];
        double x[2];
        A[0] = leftLine[0];
        A[1] = leftLine[1];
        A[2] = rightLine[0];
        A[3] = rightLine[1];

        b[0] = -leftLine[2];
        b[1] = -rightLine[2];

        bool valid = solveSystem2x2(A, b, x);

        if (not valid) {
            p = neigh[midIndex];
            return;
        }

        Point2d temp(x[0], x[1]);

        if (pointsDistance(neigh[midIndex], temp) < 0.20) {
            p = temp;
        } else {
            p = neigh[midIndex];
        }
    }

    void FALKOExtractor::generateLine(const std::vector<Point2d>& points, Eigen::Vector3d& model) {
        double sxx = 0.0;
        double syy = 0.0;
        double sxy = 0.0;
        double sx = 0.0;
        double sy = 0.0;
        int num = 0;
        for (unsigned int i = 0; i < points.size(); ++i) {
            sxx += points[i](0) * points[i](0);
            syy += points[i](1) * points[i](1);
            sxy += points[i](0) * points[i](1);
            sx += points[i](0);
            sy += points[i](1);
            ++num;
        }

        double msxx = (sxx - sx * sx / num) / num;
        double msyy = (syy - sy * sy / num) / num;
        double msxy = (sxy - sx * sy / num) / num;
        double b = 2.0 * msxy;
        double a = msxx;
        double c = msyy;
        double theta = 0.5 * (atan2(b, a - c) + M_PI);
        theta = atan2(sin(theta), cos(theta));
        double rho = (sx * cos(theta) + sy * sin(theta)) / num;

        if (rho < 0) {
            theta = atan2(sin(theta + M_PI), cos(theta + M_PI));
            rho = -rho;
        }
        model(0) = cos(theta);
        model(1) = sin(theta);
        model(2) = -rho;
    }

    bool FALKOExtractor::solveSystem2x2(double* A, double* b, double* x) {
        double det = A[0] * A[3] - A[1] * A[2];
        if (std::abs(det) < 0.0000001)
            return false;
        x[0] = (b[0] * A[3] - A[1] * b[1]) / det;
        x[1] = (A[0] * b[1] - b[0] * A[2]) / det;
        return true;
    }
}