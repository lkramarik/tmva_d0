#ifndef TMVACUTS_H
#define TMVACUTS_H

namespace tmvaCuts {
    int const totalNumberOfEvents = 120e6; // Run16 dataset
    int   const nPtBins = 5;
    float const PtBins[nPtBins+1] = {0., 1., 2., 3., 5., 10};
    float const minPt = 0.15;
    float const dcaV0ToPv[nPtBins] = {1,1,1,1,1};
    float const decayLength[nPtBins] = {0.0050, 0.0050, 0.0050, 0.0050, 0.0050};
    float const cosTheta[nPtBins] = {0.5, 0.95, 0.95, 0.95, 0.95};
    float const dcaDaughters[nPtBins] = {0.0100, 0.0100, 0.0100, 0.0100, 0.0100};
    float const kDca[nPtBins] = {0.0050, 0.0050, 0.0050, 0.0050, 0.0050};
    float const pDca[nPtBins] = {0.0050, 0.0050, 0.0050, 0.0050, 0.0050};
}
#endif