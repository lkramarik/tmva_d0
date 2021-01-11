#ifndef analysisSetup_H
#define analysisSetup_H

namespace analysisSetup
{
    int const nPtBins = 3;
    float const ptBinMin[] = {1,2,3};
    float const ptBinMax[] = {2,3,5};

    float const nTrees[] = {100,150,400};
    float const maxDepth[] = {3,3,3};
    float const bdtCut[] = {0.7552,0.64516,0.53154};
}
#endif