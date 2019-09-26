#ifndef TMVACUTS_H
#define TMVACUTS_H

namespace tmvaCuts
{
    int const totalNumberOfEvents = 120e6; // Run14 dataset
    int   const nPtBins = 5;
    float const PtBins[nPtBins+1] = {0., 1., 2., 3., 5., 10};

    float const minPt = 0.15;

    float const dcaV0ToPv = 0.05;
    float const decayLength = 0.0005; //0.0005
    float const cosTheta = 0.7;
    float const dcaDaughters = 0.02;
    float const kDca = 0.002;
    float const pDca = 0.002;

}
#endif