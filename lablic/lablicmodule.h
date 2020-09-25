/*********************************************************************
 *  Author  : Himangshu Saikia
 *  Init    : Friday, September 29, 2017 - 15:50:00
 *
 *  Project : KTH Inviwo Modules
 *
 *  License : Follows the Inviwo BSD license model
 *********************************************************************
 */

#pragma once

#include <lablic/lablicmoduledefine.h>
#include <inviwo/core/common/inviwomodule.h>

namespace inviwo
{

class IVW_MODULE_LABLIC_API LabLICModule : public InviwoModule
{
public:
    LabLICModule(InviwoApplication* app);
    virtual ~LabLICModule() = default;
};

} // namespace
