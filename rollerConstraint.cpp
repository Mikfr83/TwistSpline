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

#include "rollerConstraint.h"

#include <maya/MFnMatrixAttribute.h>
#include <maya/MFnTypedAttribute.h>
#include <maya/MFnNumericAttribute.h>
#include <maya/MFnUnitAttribute.h>
#include <maya/MFnEnumAttribute.h>
#include <maya/MFnCompoundAttribute.h>
#include <maya/MRampAttribute.h>
#include <maya/MObjectArray.h>
#include <maya/MPxData.h>
#include <maya/MArrayDataBuilder.h>
#include <maya/MEulerRotation.h>
#include <maya/MQuaternion.h>

#include <maya/MMatrix.h>
#include <maya/MTransformationMatrix.h>
#include <maya/MVector.h>
#include <maya/MPoint.h>
#include <maya/MTypes.h>

#include <vector>
#include "twistSpline.h"
#include "twistSplineData.h"


#define ROTATE_ORDER_XYZ        0
#define ROTATE_ORDER_YZX        1
#define ROTATE_ORDER_ZXY        2
#define ROTATE_ORDER_XZY        3
#define ROTATE_ORDER_YXZ        4
#define ROTATE_ORDER_ZYX        5

#define CHECKSTAT(m) if (!status) {status.perror(m); return status;};

MTypeId     rollerConstraint::id(0x001226FC);

MObject     rollerConstraint::aRotateOrder;
MObject     rollerConstraint::aGlobalOffset;
MObject     rollerConstraint::aGlobalSpread;
MObject     rollerConstraint::aUseCycle;
MObject     rollerConstraint::aNormalize;
MObject     rollerConstraint::aNormValue;

MObject     rollerConstraint::aUseGlobalMin;
MObject     rollerConstraint::aMinGlobalParam;
MObject     rollerConstraint::aUseGlobalMax;
MObject     rollerConstraint::aMaxGlobalParam;

// inputs
MObject     rollerConstraint::aInputSplines;
    MObject     rollerConstraint::aSpline;
    MObject     rollerConstraint::aSplineLength;
    MObject     rollerConstraint::aEndParam;
    MObject     rollerConstraint::aWeight;


MObject     rollerConstraint::aParams;
    MObject     rollerConstraint::aParam;
    MObject     rollerConstraint::aUseMin;
    MObject     rollerConstraint::aMinParam;
    MObject     rollerConstraint::aUseMax;
    MObject     rollerConstraint::aMaxParam;

    MObject     rollerConstraint::aRolling;
    MObject     rollerConstraint::aRollingRadius;

    MObject     rollerConstraint::aParentInverseMatrix;


// output
MObject     rollerConstraint::aOutputs;
    MObject     rollerConstraint::aTranslate;
    MObject     rollerConstraint::aTranslateX;
    MObject     rollerConstraint::aTranslateY;
    MObject     rollerConstraint::aTranslateZ;
    MObject     rollerConstraint::aRotate;
    MObject     rollerConstraint::aRotateX;
    MObject     rollerConstraint::aRotateY;
    MObject     rollerConstraint::aRotateZ;
    MObject     rollerConstraint::aScale;
    MObject     rollerConstraint::aScaleX;
    MObject     rollerConstraint::aScaleY;
    MObject     rollerConstraint::aScaleZ;

MStatus rollerConstraint::initialize() {
    MStatus status;
    MFnTypedAttribute tAttr;
    MFnNumericAttribute nAttr;
    MFnMatrixAttribute mAttr;
    MFnUnitAttribute uAttr;
    MFnEnumAttribute eAttr;
    MFnCompoundAttribute cAttr;
    MRampAttribute rAttr;

    // Single inputs
    aRotateOrder = eAttr.create("rotateOrder", "ro", ROTATE_ORDER_XYZ, &status);
    CHECKSTAT("aRotateOrder");
    eAttr.setKeyable(true);
    status = eAttr.addField("xyz", ROTATE_ORDER_XYZ);
    CHECKSTAT("aRotateOrder");
    status = eAttr.addField("yzx", ROTATE_ORDER_YZX);
    CHECKSTAT("aRotateOrder");
    status = eAttr.addField("zxy", ROTATE_ORDER_ZXY);
    CHECKSTAT("aRotateOrder");
    status = eAttr.addField("xzy", ROTATE_ORDER_XZY);
    CHECKSTAT("aRotateOrder");
    status = eAttr.addField("yxz", ROTATE_ORDER_YXZ);
    CHECKSTAT("aRotateOrder");
    status = eAttr.addField("zyx", ROTATE_ORDER_ZYX);
    CHECKSTAT("aRotateOrder");
    status = addAttribute(aRotateOrder);
    CHECKSTAT("aRotateOrder");

    aGlobalOffset = nAttr.create("globalOffset", "go", MFnNumericData::kDouble, 0.0, &status);
    CHECKSTAT("aGlobalOffset");
    nAttr.setKeyable(true);
    status = addAttribute(aGlobalOffset);
    CHECKSTAT("aGlobalOffset");
    
    aGlobalSpread = nAttr.create("globalSpread", "gs", MFnNumericData::kDouble, 1.0, &status);
    CHECKSTAT("aGlobalSpread");
    nAttr.setKeyable(true);
    status = addAttribute(aGlobalSpread);
    CHECKSTAT("aGlobalSpread");

    aUseCycle = nAttr.create("useCycle", "uc", MFnNumericData::kBoolean, false, &status);
    CHECKSTAT("aUseCycle");
    nAttr.setKeyable(true);
    status = addAttribute(aUseCycle);
    CHECKSTAT("aUseCycle");


    aNormalize = nAttr.create("normalize", "n", MFnNumericData::kDouble, 1.0, &status);
    CHECKSTAT("aNormalize");
    nAttr.setMin(0.0);
    nAttr.setMax(1.0);
    nAttr.setKeyable(true);
    status = addAttribute(aNormalize);
    CHECKSTAT("aNormalize");

    aNormValue = nAttr.create("normValue", "nv", MFnNumericData::kDouble, 1.0, &status);
    CHECKSTAT("aNormValue");
    nAttr.setMin(0.0);
    nAttr.setKeyable(true);
    status = addAttribute(aNormValue);
    CHECKSTAT("aNormValue");




    aUseGlobalMin = nAttr.create("useGlobalMin", "ugn", MFnNumericData::kBoolean, false, &status);
    CHECKSTAT("aUseGlobalMin");
    nAttr.setKeyable(true);
    status = addAttribute(aUseGlobalMin);
    CHECKSTAT("aUseGlobalMin");

    aMinGlobalParam = nAttr.create("minGlobalParam", "ngp", MFnNumericData::kDouble, 0.0, &status);
    CHECKSTAT("aMinGlobalParam");
    nAttr.setKeyable(true);
    status = addAttribute(aMinGlobalParam);
    CHECKSTAT("aMinGlobalParam");

    aUseGlobalMax = nAttr.create("useGlobalMax", "ugx", MFnNumericData::kBoolean, false, &status);
    CHECKSTAT("aUseGlobalMax");
    nAttr.setKeyable(true);
    status = addAttribute(aUseGlobalMax);
    CHECKSTAT("aUseGlobalMax");

    aMaxGlobalParam = nAttr.create("maxGlobalParam", "xgp", MFnNumericData::kDouble, 1.0, &status);
    CHECKSTAT("aMaxGlobalParam");
    nAttr.setKeyable(true);
    status = addAttribute(aMaxGlobalParam);
    CHECKSTAT("aMaxGlobalParam");


    // Spline input array
    aSpline = tAttr.create("spline", "s", TwistSplineData::id, MObject::kNullObj, &status);
    tAttr.setHidden(true);
    CHECKSTAT("aSpline");

    aSplineLength = uAttr.create("splineLength", "sl", MFnUnitAttribute::kDistance, 0.0, &status);
    CHECKSTAT("aSplineLength");
    uAttr.setKeyable(true);

    aEndParam = nAttr.create("endParam", "ep", MFnNumericData::kDouble, 1.0, &status);
    CHECKSTAT("aEndParam");
    nAttr.setKeyable(true);

    aWeight = nAttr.create("weight", "w", MFnNumericData::kDouble, 1.0, &status);
    CHECKSTAT("aWeight");
    nAttr.setKeyable(true);

    aInputSplines = cAttr.create("inputSplines", "is", &status);
    CHECKSTAT("aInputSplines");
    cAttr.setArray(true);
    cAttr.addChild(aSpline);
    cAttr.addChild(aSplineLength);
    cAttr.addChild(aEndParam);
    cAttr.addChild(aWeight);

    status = addAttribute(aInputSplines);
    CHECKSTAT("aInputSplines");

    // input param array
    aParam = nAttr.create("param", "p", MFnNumericData::kDouble, 0.0, &status);
    CHECKSTAT("aParam");
    nAttr.setKeyable(true);

    aUseMin = nAttr.create("useMin", "un", MFnNumericData::kBoolean, false, &status);
    CHECKSTAT("aUseMin");
    nAttr.setKeyable(true);

    aMinParam = nAttr.create("minParam", "np", MFnNumericData::kDouble, 0.0, &status);
    CHECKSTAT("aMinParam");
    nAttr.setKeyable(true);

    aUseMax = nAttr.create("useMax", "ux", MFnNumericData::kBoolean, false, &status);
    CHECKSTAT("aUseMax");
    nAttr.setKeyable(true);

    aMaxParam = nAttr.create("maxParam", "xp", MFnNumericData::kDouble, 1.0, &status);
    CHECKSTAT("aMaxParam");
    nAttr.setKeyable(true);

    aRolling = nAttr.create("rolling", "ro", MFnNumericData::kBoolean, false, &status);
    CHECKSTAT("aRolling");

    aRollingRadius = rAttr.createCurveRamp("rollingRadius", "rr", &status);
    CHECKSTAT("aRollingRadius");

    aParentInverseMatrix = mAttr.create("parentInverseMatrix", "pim", MFnMatrixAttribute::kDouble, &status);
    CHECKSTAT("aParentInverseMatrix");
    mAttr.setHidden(true);

    aParams = cAttr.create("params", "ps", &status);
    cAttr.setArray(true);
    status = cAttr.addChild(aParam);
    CHECKSTAT("child aParam");
    status = cAttr.addChild(aUseMin);
    CHECKSTAT("child aUseMin");
    status = cAttr.addChild(aMinParam);
    CHECKSTAT("child aMinParam");
    status = cAttr.addChild(aUseMax);
    CHECKSTAT("child aUseMax");
    status = cAttr.addChild(aMaxParam);
    CHECKSTAT("child aMaxParam");
    status = cAttr.addChild(aParentInverseMatrix);
    CHECKSTAT("child aParentInverseMatrix");
    status = addAttribute(aParams);
    CHECKSTAT("aParams");

    // Output: Matrices
    aTranslateX = uAttr.create("translateX", "tx", MFnUnitAttribute::kDistance, 0.0, &status);
    CHECKSTAT("aTranslateX");
    uAttr.setWritable(false);
    uAttr.setStorable(false);

    aTranslateY = uAttr.create("translateY", "ty", MFnUnitAttribute::kDistance, 0.0, &status);
    CHECKSTAT("aTranslateY");
    uAttr.setWritable(false);
    uAttr.setStorable(false);

    aTranslateZ = uAttr.create("translateZ", "tz", MFnUnitAttribute::kDistance, 0.0, &status);
    CHECKSTAT("aTranslateZ");
    uAttr.setWritable(false);
    uAttr.setStorable(false);

    aTranslate = nAttr.create("translate", "t", aTranslateX, aTranslateY, aTranslateZ, &status);
    nAttr.setHidden(true);
    CHECKSTAT("aTranslate");
    nAttr.setWritable(false);
    nAttr.setStorable(false);

    // Output Attributes
    aRotateX = uAttr.create("rotateX", "rotx", MFnUnitAttribute::kAngle, 0.0, &status);
    CHECKSTAT("aRotateX");
    uAttr.setWritable(false);
    uAttr.setStorable(false);

    aRotateY = uAttr.create("rotateY", "roty", MFnUnitAttribute::kAngle, 0.0, &status);
    CHECKSTAT("aRotateY");
    uAttr.setWritable(false);
    uAttr.setStorable(false);

    aRotateZ = uAttr.create("rotateZ", "rotz", MFnUnitAttribute::kAngle, 0.0, &status);
    CHECKSTAT("aRotateZ");
    uAttr.setWritable(false);
    uAttr.setStorable(false);

    aRotate = nAttr.create("rotate", "rot", aRotateX, aRotateY, aRotateZ, &status);
    nAttr.setHidden(true);
    CHECKSTAT("aRotate");
    nAttr.setWritable(false);

    // Output Attributes
    aScaleX = nAttr.create("scaleX", "sclx", MFnNumericData::kDouble, 0.0, &status);
    CHECKSTAT("aScaleX");
    nAttr.setWritable(false);
    nAttr.setStorable(false);

    aScaleY = nAttr.create("scaleY", "scly", MFnNumericData::kDouble, 0.0, &status);
    CHECKSTAT("aScaleY");
    nAttr.setWritable(false);
    nAttr.setStorable(false);

    aScaleZ = nAttr.create("scaleZ", "sclz", MFnNumericData::kDouble, 0.0, &status);
    CHECKSTAT("aScaleZ");
    nAttr.setWritable(false);
    nAttr.setStorable(false);

    aScale = nAttr.create("scale", "scl", aScaleX, aScaleY, aScaleZ, &status);
    nAttr.setHidden(true);
    CHECKSTAT("aScale");
    nAttr.setWritable(false);
    nAttr.setStorable(false);

    aOutputs = cAttr.create("outputs", "out", &status);
    CHECKSTAT("aOutputs");
    cAttr.setHidden(true);
    cAttr.setArray(true);
    cAttr.setUsesArrayDataBuilder(true);
    cAttr.addChild(aTranslate);
    cAttr.addChild(aRotate);
    cAttr.addChild(aScale);

    status = addAttribute(aOutputs);
    CHECKSTAT("aOutputs");

    attributeAffects(aRotateOrder, aRotate);
    attributeAffects(aRotateOrder, aRotateX);
    attributeAffects(aRotateOrder, aRotateY);
    attributeAffects(aRotateOrder, aRotateZ);

    std::vector<MObject *> iobjs = {
        &aInputSplines, &aSpline, &aSplineLength, &aEndParam, &aWeight, &aParams, &aParam,
        &aParentInverseMatrix, &aGlobalOffset, &aUseCycle, &aGlobalSpread, &aNormalize,
        &aNormValue, &aUseMin, &aMinParam, &aUseMax, &aMaxParam,
        &aUseGlobalMin, &aMinGlobalParam, &aUseGlobalMax, &aMaxGlobalParam
    };

    std::vector<MObject *> oobjs = {
        &aOutputs, &aTranslate, &aRotate, &aScale, &aTranslateX, &aRotateX, &aScaleX,
        &aTranslateY, &aRotateY, &aScaleY, &aTranslateZ, &aRotateZ, &aScaleZ
    };

    for (auto &ii : iobjs) {
        for (auto &oo : oobjs) {
            attributeAffects(*ii, *oo);
        }
    }

    return MS::kSuccess;
}


void getRoll(TwistSplineT* spline, std::vector<double> radii) {
    // for each point in the lut get the normal and binormal
    // The normal is the axis of rotation
    // the binormal is the direction of the ground

    // Assuming a constant radius for the duration of the step
    // we can directly compute the radian angle as stepDist/radius

    // What I will need to do is come up with a way to read the entirety
    // of a curve... maybe it's a baking step? Maybe I can steal from harmonics?
    // By the time we get here, I want a radius per LUT
   


}




void rollerConstraint::getSplines(
        MDataBlock& data,
        std::vector<TwistSplineT*> &splines,
        std::vector<double> &weights,
        std::vector<double> &splineLens,
        std::vector<double> &endParams

    ) {
    unsigned ecount, possibleMax;

    // First, get the splines that need computing, and their weight.
    MArrayDataHandle inSpAH = data.inputArrayValue(aInputSplines);
    ecount = inSpAH.elementCount();
    inSpAH.jumpToArrayElement(ecount - 1);
    inSpAH.jumpToArrayElement(0);

    for (unsigned i = 0; i < ecount; ++i) {
        auto inSpGroup = inSpAH.inputValue();
        MDataHandle inSpwH = inSpGroup.child(aWeight);

        double inSpw = inSpwH.asDouble();
        if (inSpw <= 0.0){
            inSpAH.next();
            continue;
        }

        MDataHandle inSpH = inSpGroup.child(aSpline);
        MPxData* pd = inSpH.asPluginData();
        if (pd == nullptr) {
            inSpAH.next();
            continue;
        }

        MDataHandle inSpLenH = inSpGroup.child(aSplineLength);
        double inSpLen = inSpLenH.asDouble();

        MDataHandle endParamH = inSpGroup.child(aEndParam);
        double endParam = endParamH.asDouble();

        auto inSplineData = (TwistSplineData *)pd;
        TwistSplineT *spline = inSplineData->getSpline();
        if (spline == nullptr) {
            inSpAH.next();
            continue;
        }

        splines.push_back(spline);
        weights.push_back(inSpw);
        splineLens.push_back(inSpLen);
        endParams.push_back(endParam);
        inSpAH.next();
    }
}


void rollerConstraint::loadSplineParams(
        MDataBlock& data,
        std::vector<double> &params,
        std::vector<double> &pMins,
        std::vector<double> &pMaxs,
        std::vector<bool> &pUseMins,
        std::vector<bool> &pUseMaxs,
        std::vector<MMatrix> &invParMats
    ){
    // Get the params
    MStatus status;

    MArrayDataHandle inPAH = data.inputArrayValue(aParams, &status);

    unsigned ecount, possibleMax;
    ecount = inPAH.elementCount();
    inPAH.jumpToArrayElement(ecount - 1);
    possibleMax = inPAH.elementIndex();
    params.resize(possibleMax + 1);
    pMins.resize(possibleMax + 1);
    pUseMins.resize(possibleMax + 1);
    pMaxs.resize(possibleMax + 1);
    pUseMaxs.resize(possibleMax + 1);

    invParMats.resize(possibleMax + 1);
    inPAH.jumpToArrayElement(0);

    for (unsigned i = 0; i < ecount; ++i) {
        unsigned realIndex = inPAH.elementIndex();
        if (realIndex > possibleMax) {
            params.resize(realIndex+1);
            pMins.resize(realIndex+1);
            pMaxs.resize(realIndex+1);
            pUseMins.resize(realIndex+1);
            pUseMaxs.resize(realIndex+1);
            invParMats.resize(realIndex+1);
            possibleMax = realIndex;
        }
        auto inPGrp = inPAH.inputValue(&status);

        MDataHandle inPH = inPGrp.child(aParam);
        double inParam = inPH.asDouble();
        params[realIndex] = inParam;

        MDataHandle inUN = inPGrp.child(aUseMin);
        pUseMins[realIndex] = inUN.asBool();
        MDataHandle inN = inPGrp.child(aMinParam);
        pMins[realIndex] = inN.asDouble();

        MDataHandle inUX = inPGrp.child(aUseMax);
        pUseMaxs[realIndex] = inUX.asBool();
        MDataHandle inX = inPGrp.child(aMaxParam);
        pMaxs[realIndex] = inX.asDouble();

        MDataHandle inPimH = inPGrp.child(aParentInverseMatrix);
        MMatrix inPim = inPimH.asMatrix();
        invParMats[realIndex] = inPim;

        inPAH.next();
    }
}


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
    ){
    // Loop through the splines to get the transform groups
    // There is no single way to interpolate between orientations
    // I'm just going to do the simple thing (for now) and N-LERP
    // between the quaternion outputs in order

    bool twisted = splines.size() == 1;

    for (size_t sIdx = 0; sIdx<splines.size(); ++sIdx) {
        TwistSplineT *spline = splines[sIdx];
        if (spline == nullptr) continue;

        double mp = endParams[sIdx];
        double mrmp = splineLens[sIdx];
        if (mrmp == 0.0){
            const auto &lp = spline->getLockPositions();
            mp = lp[lp.size() - 1];

            const auto &rmp = spline->getRemap();
            mrmp = rmp[rmp.size() - 1];
        }

        std::vector<MPoint> ttrans, tscales;
        std::vector<MQuaternion> tquats;
        std::vector<double> ttwists;

        ttrans.reserve(params.size());
        tscales.reserve(params.size());
        tquats.reserve(params.size());
        ttwists.reserve(params.size());

        for (size_t pIdx = 0; pIdx<params.size(); ++pIdx) {
            // Calculate the transforms
            MVector tan, norm, binorm;
            MPoint tran, scale;

            double twist, p = params[pIdx];
            p *= gSpread;
            p += gOffset;

            if (useGlobalMin)   { p = std::fmax(minGlobalParam, p); }
            if (pUseMins[pIdx]) { p = std::fmax(pMins[pIdx],    p); }
            if (useGlobalMax)   { p = std::fmin(maxGlobalParam, p); }
            if (pUseMaxs[pIdx]) { p = std::fmin(pMaxs[pIdx],    p); }

            if (doNorm > 0.0) {
                p = (doNorm) * (p * mp / normVal) + (1.0 - doNorm) * (p);
            }
            if (useCycle) {
                if (p >= 0.0)
                    p = std::fmod(p, mrmp);
                else
                    p = mrmp - std::fmod(-p, mrmp);
            }

            spline->matrixAtParam(p, tan, norm, binorm, tran, scale, twist, twisted);
            ttrans.push_back(std::move(tran));
            tscales.push_back(std::move(scale));
            ttwists.push_back(std::move(twist));

            // fill the matrix
            double mat[4][4];
            mat[0][0] = tan[0]; mat[1][0] = norm[0]; mat[2][0] = binorm[0]; mat[3][0] = 0.0;
            mat[0][1] = tan[1]; mat[1][1] = norm[1]; mat[2][1] = binorm[1]; mat[3][1] = 0.0;
            mat[0][2] = tan[2]; mat[1][2] = norm[2]; mat[2][2] = binorm[2]; mat[3][2] = 0.0;
            mat[0][3] = 0.0;    mat[1][3] = 0.0;     mat[2][3] = 0.0;       mat[3][3] = 1.0;
            MQuaternion q;
            q = MMatrix(mat);
            tquats.push_back(std::move(q));
        }

        trans.push_back(ttrans);
        scales.push_back(tscales);
        quats.push_back(tquats);
        twists.push_back(ttwists);
    }
}

void rollerConstraint::normalizeWeights(std::vector<double> &weights){
    // normalize the weights
    double s = 0.0;
    for (double &w : weights) s += w;
    for (double &w : weights) w /= s;
}

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
    ){

    if (splines.size() > 1) {
        // weight the twists trans and scales
        for (size_t pIdx = 0; pIdx < params.size(); ++pIdx) {
            MPoint tran, scale;
            double twist = 0.0;
            for (size_t sIdx = 0; sIdx < splines.size(); ++sIdx) {
                double w = weights[sIdx];
                tran += trans[sIdx][pIdx] * w;
                scale += scales[sIdx][pIdx] * w;
                twist += twists[sIdx][pIdx] * w;
            }
            otrans.push_back(std::move(tran));
            oscales.push_back(std::move(scale));
            otwists.push_back(std::move(twist));
        }

        // Weight the quaternions
        for (size_t pIdx = 0; pIdx < params.size(); ++pIdx) {
            MQuaternion prev = quats[0][pIdx].normal();
            for (size_t sIdx = 1; sIdx < splines.size(); ++sIdx) {
                MQuaternion cur = quats[sIdx][pIdx];
                double pw = weights[sIdx - 1]; // preWeight
                double cw = weights[sIdx]; // currentWeight
                double s = pw + cw;
                prev = slerp(prev, cur, cw / s);
            }
            MQuaternion xt(otwists[pIdx], MVector(1.0, 0.0, 0.0));
            prev = xt * prev;
            oquats.push_back(std::move(prev));
        }
    }
    else if (splines.size() == 1) {
        // don't weight anything if we don't have to
        otrans = trans[0];
        oscales = scales[0];
        oquats = quats[0];
    }

}




MStatus rollerConstraint::compute(const MPlug& plug, MDataBlock& data) {
    if (plug == aOutputs || 
        plug == aScale || plug == aScaleX || plug == aScaleY || plug == aScaleZ ||
        plug == aRotate || plug == aRotateX || plug == aRotateY || plug == aRotateZ ||
        plug == aTranslate || plug == aTranslateX || plug == aTranslateY || plug == aTranslateZ
        ) {
        // I can optimize things a lot more if we do everything at once
        // So, whatever plug is asked for, just compute it all

        // Get the spline data
        std::vector<TwistSplineT*> splines;
        std::vector<double> weights;
        std::vector<double> splineLens;
        std::vector<double> endParams;
        getSplines(data, splines, weights, splineLens, endParams);
        normalizeWeights(weights);

        // Get the per-rider params
        std::vector<double> params, pMins, pMaxs;
        std::vector<bool> pUseMins, pUseMaxs;
        std::vector<MMatrix> invParMats;
        loadSplineParams(data, params, pMins, pMaxs, pUseMins, pUseMaxs, invParMats);

        // Get global constraint data
        MDataHandle gOffsetH = data.inputValue(aGlobalOffset);
        MDataHandle gSpreadH = data.inputValue(aGlobalSpread);
        MDataHandle useCycleH = data.inputValue(aUseCycle);

        MDataHandle useGlobalMinH = data.inputValue(aUseGlobalMin);
        MDataHandle minGlobalParamH = data.inputValue(aMinGlobalParam);
        MDataHandle useGlobalMaxH = data.inputValue(aUseGlobalMax);
        MDataHandle maxGlobalParamH = data.inputValue(aMaxGlobalParam);

        MDataHandle normalizeH = data.inputValue(aNormalize);
        MDataHandle normValueH = data.inputValue(aNormValue);
        double doNorm = normalizeH.asDouble();
        double normVal = normValueH.asDouble();
        double gOffset = gOffsetH.asDouble();
        double gSpread = gSpreadH.asDouble();
        bool useCycle = useCycleH.asBool();

        bool useGlobalMin = useGlobalMinH.asBool();
        double minGlobalParam = minGlobalParamH.asDouble();
        bool useGlobalMax = useGlobalMaxH.asBool();
        double maxGlobalParam = maxGlobalParamH.asDouble();

        // Evaluate all the data to S/R/T's per rider, per connected spline
        std::vector<std::vector<MPoint>> trans, scales;
        std::vector<std::vector<MQuaternion>> quats;
        std::vector<std::vector<double>> twists;
        evaluateSplines(
            doNorm, normVal, gOffset, gSpread, useCycle, useGlobalMin, minGlobalParam, useGlobalMax, maxGlobalParam,
            splines, weights, splineLens, endParams, params, pMins, pMaxs, pUseMins, pUseMaxs,
            trans, scales, quats, twists
        );

        // Blend all the spline data together to get one output
        std::vector<MPoint> otrans, oscales;
        std::vector<MQuaternion> oquats;
        std::vector<double> otwists;
        blendSplines(splines, weights, params, trans, scales, quats, twists,
            otrans, oscales, oquats, otwists
        );

        // Finally, set the output data
        MDataHandle hOrder = data.inputValue(aRotateOrder);
        short order = hOrder.asShort();

        // Now we can set all the outputs
        MArrayDataHandle outHandle = data.outputArrayValue(aOutputs);

        MArrayDataBuilder builder = outHandle.builder();
        for (size_t pIdx = 0; pIdx < params.size(); ++pIdx) {

            MDataHandle outH = builder.addElement(pIdx);
            MDataHandle tranH = outH.child(aTranslate);
            MDataHandle rotH = outH.child(aRotate);
            MDataHandle sclH = outH.child(aScale);

            MMatrix invPar = invParMats[pIdx];
            MPoint tran = otrans[pIdx] * invPar;
            MMatrix qmat = oquats[pIdx].asMatrix() * invPar;
            MPoint &scale = oscales[pIdx];

            // Ugh, the only way to convert to quat to euler *with a given rot order*
            // is expanding to matrix and decomposing
            MEulerRotation meu = MEulerRotation::decompose(qmat, (MEulerRotation::RotationOrder)order);

            tranH.set3Double(tran[0], tran[1], tran[2]);
            rotH.set3Double(meu.x, meu.y, meu.z);
            sclH.set3Double(scale[0], scale[1], scale[2]);
        }
        outHandle.set(builder);
        outHandle.setAllClean();

    }
    else {
        return MS::kUnknownParameter;
    }

    return MS::kSuccess;
}

