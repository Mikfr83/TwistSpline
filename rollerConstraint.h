/*
MIT License

Copyright (c) 2018 Blur Studio

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.
*/

#pragma once
#include <maya/MPxNode.h>
#include <maya/MDataBlock.h>
#include <maya/MMatrix.h>
#include <maya/MPoint.h>
#include <maya/MQuaternion.h>

class rollerConstraint : public MPxNode {
public:
	rollerConstraint() {}
	static void* creator() {
		return new rollerConstraint();
	}
	virtual MStatus compute(const MPlug& plug, MDataBlock& data);
	static	MStatus initialize();

    void rollerConstraint::getSplines(
        MDataBlock& data,
        std::vector<TwistSplineT*> &splines,
        std::vector<double> &weights,
        std::vector<double> &splineLens,
        std::vector<double> &endParams
    );

    void rollerConstraint::normalizeWeights(std::vector<double> &weights);

    void rollerConstraint::loadSplineParams(
        MDataBlock& data,
        std::vector<double> &params,
        std::vector<double> &pMins,
        std::vector<double> &pMaxs,
        std::vector<bool> &pUseMins,
        std::vector<bool> &pUseMaxs,
        std::vector<MMatrix> &invParMats
    );

    void rollerConstraint::evaluateSplines(
        double doNorm,
        double normVal,
        double gOffset,
        double gSpread,
        bool useCycle,
        bool useGlobalMin,
        double minGlobalParam,
        bool useGlobalMax,
        double maxGlobalParam,

        const std::vector<TwistSplineT*> &splines,
        const std::vector<double> &weights,
        const std::vector<double> &splineLens,
        const std::vector<double> &endParams,
        const std::vector<double> &params,
        const std::vector<double> &pMins,
        const std::vector<double> &pMaxs,
        const std::vector<bool> &pUseMins,
        const std::vector<bool> &pUseMaxs,

        std::vector<std::vector<MPoint>> &trans,
        std::vector<std::vector<MPoint>> &scales,
        std::vector<std::vector<MQuaternion>> &quats,
        std::vector<std::vector<double>> &twists
    );

    void rollerConstraint::blendSplines(
        const std::vector<TwistSplineT*> &splines,
        const std::vector<double> &weights,
        const std::vector<double> &params,
        const std::vector<std::vector<MPoint>> &trans,
        const std::vector<std::vector<MPoint>> &scales,
        const std::vector<std::vector<MQuaternion>> &quats,
        const std::vector<std::vector<double>> &twists,

        std::vector<MPoint> &otrans,
        std::vector<MPoint> &oscales,
        std::vector<MQuaternion> &oquats,
        std::vector<double> &otwists
    );


public:
	// inputs
	static MObject aRotateOrder;
	static MObject aGlobalOffset;
	static MObject aGlobalSpread;
	static MObject aUseCycle;
	static MObject aNormalize;
	static MObject aNormValue;
	static MObject aUseGlobalMin;
	static MObject aMinGlobalParam;
	static MObject aUseGlobalMax;
	static MObject aMaxGlobalParam;

	// inputs
	static MObject aInputSplines;
		static MObject aSpline;
		static MObject aSplineLength;
		static MObject aEndParam;
		static MObject aWeight;

	static MObject aParams;
		static MObject aParam;

		static MObject aUseMin;
		static MObject aMinParam;
		static MObject aUseMax;
		static MObject aMaxParam;

		static MObject aRolling;
		static MObject aRollingRadius;

		static MObject aParentInverseMatrix;

	// output
	static MObject aOutputs;
		static MObject aOutMat;
		static MObject aTranslate;
		static MObject aTranslateX;
		static MObject aTranslateY;
		static MObject aTranslateZ;
		static MObject aRotate;
		static MObject aRotateX;
		static MObject aRotateY;
		static MObject aRotateZ;
		static MObject aScale;
		static MObject aScaleX;
		static MObject aScaleY;
		static MObject aScaleZ;

	static MTypeId id;
};

