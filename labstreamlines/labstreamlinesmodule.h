/*********************************************************************
 *  Author  : Himangshu Saikia
 *  Init    : Tuesday, September 19, 2017 - 15:07:48
 *
 *  Project : KTH Inviwo Modules
 *
 *  License : Follows the Inviwo BSD license model
 *********************************************************************
 */

#pragma once

#include <labstreamlines/labstreamlinesmoduledefine.h>
#include <inviwo/core/common/inviwomodule.h>

namespace inviwo
{

class IVW_MODULE_LABSTREAMLINES_API LabStreamlinesModule : public InviwoModule
{
public:
    LabStreamlinesModule(InviwoApplication* app);
    virtual ~LabStreamlinesModule() = default;
};

} // namespace
