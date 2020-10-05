/*********************************************************************
 *  Author  : Anke Friederici
 *
 *  Project : KTH Inviwo Modules
 *
 *  License : Follows the Inviwo BSD license model
 *********************************************************************
 */

#include <labtopo/labtopomodule.h>
#include <labtopo/topology.h>

namespace inviwo
{

LabTopoModule::LabTopoModule(InviwoApplication* app) : InviwoModule(app, "LabTopo") {
    registerProcessor<Topology>();
}

} // namespace
