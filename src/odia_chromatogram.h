// odia_chromatogram.h — the compact XIC store now lives in OpenMS, where its consumer is.
//
// The class was developed and self-tested here, then moved to
// OpenMS/ANALYSIS/OPENSWATH/CompactChromatogramStore.h so that ChromatogramExtractor and
// OpenSwathWorkflow can use it: OpenMS cannot include a header from this tool, and duplicating the
// class would let the two copies drift. The header has no OpenMS dependencies, so
// odia_chromatogram_test still compiles standalone against it.
#ifndef ODIA_CHROMATOGRAM_H
#define ODIA_CHROMATOGRAM_H

#include <OpenMS/ANALYSIS/OPENSWATH/CompactChromatogramStore.h>

namespace odia
{
using ChromatogramStore = OpenMS::CompactChromatogramStore;
}

#endif // ODIA_CHROMATOGRAM_H
