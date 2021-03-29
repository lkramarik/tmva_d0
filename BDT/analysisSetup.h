#ifndef analysisSetup_H
#define analysisSetup_H

namespace analysisSetup
{
    int const nPtBins = 3;
    Double_t  ptBinMin[] = {1,2,3};
    Double_t  ptBinMax[] = {2,3,5};

    Double_t  nTrees[] = {150,150,400};
    Double_t  maxDepth[] = {3,3,3};
    Double_t  bdtCut[] = {0.70913, 0.676477, 0.480328};
}
#endif