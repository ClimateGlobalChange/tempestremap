#include "GridElements.h"
#include "Exception.h"
#include "Announce.h"

#include <vector>

///////////////////////////////////////////////////////////////////////////////

extern "C"
int GenerateConnectivityData(
	const Mesh & meshIn,
	std::vector< std::set<int> > & vecConnectivity
) {

	// Number of elements
	int nElements = meshIn.faces.size();

	vecConnectivity.resize(nElements);

	EdgeMapConstIterator iter = meshIn.edgemap.begin();
	for (; iter != meshIn.edgemap.end(); iter++) {

		if ((iter->second[0] != InvalidFace) &&
			(iter->second[1] != InvalidFace)
			) {
			if ((iter->second[0] < 0) || (iter->second[0] >= nElements)) {
				_EXCEPTION1("Face index (%i) out of range",
							iter->second[0]);
			}
			if ((iter->second[1] < 0) || (iter->second[1] >= nElements)) {
				_EXCEPTION1("Face index (%i) out of range",
							iter->second[1]);
			}
			
			vecConnectivity[iter->second[0]].insert(iter->second[1]+1);
			vecConnectivity[iter->second[1]].insert(iter->second[0]+1);
		}
	}

	return 0;
}

///////////////////////////////////////////////////////////////////////////////

