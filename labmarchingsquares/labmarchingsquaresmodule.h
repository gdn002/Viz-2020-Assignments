/*********************************************************************
 *  Author  : Himangshu Saikia
 *  Init    : Monday, September 11, 2017 - 12:57:56
 *
 *  Project : KTH Inviwo Modules
 *
 *  License : Follows the Inviwo BSD license model
 *********************************************************************
 */

#pragma once

#include <labmarchingsquares/labmarchingsquaresmoduledefine.h>
#include <inviwo/core/common/inviwomodule.h>

namespace inviwo
{

class IVW_MODULE_LABMARCHINGSQUARES_API LabMarchingSquaresModule : public InviwoModule {
public:
    LabMarchingSquaresModule(InviwoApplication* app);
    virtual ~LabMarchingSquaresModule() = default;
};

} // namespace
