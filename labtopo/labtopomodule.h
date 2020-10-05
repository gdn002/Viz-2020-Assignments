/*********************************************************************
 *  Author  : Anke Friederici
 *
 *  Project : KTH Inviwo Modules
 *
 *  License : Follows the Inviwo BSD license model
 *********************************************************************
 */

#pragma once

#include <labtopo/labtopomoduledefine.h>
#include <inviwo/core/common/inviwomodule.h>

namespace inviwo
{

class IVW_MODULE_LABTOPO_API LabTopoModule : public InviwoModule {
public:
    LabTopoModule(InviwoApplication* app);
    virtual ~LabTopoModule() = default;
};

} // namespace
