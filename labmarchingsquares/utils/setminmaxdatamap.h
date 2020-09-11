#pragma once

#include <labmarchingsquares/labmarchingsquaresmoduledefine.h>
#include <inviwo/core/datastructures/volume/volume.h>

namespace inviwo
{
///Sets the value and data range of the data map.
void IVW_MODULE_LABMARCHINGSQUARES_API SetMinMaxForInviwoDataMap(const char* pData, const size_t NumOfBytes, std::shared_ptr<inviwo::Volume> pVolume);

};
