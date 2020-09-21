/*********************************************************************
 *  Author  : Himangshu Saikia, Wiebke Koepp
 *  Init    : Tuesday, September 19, 2017 - 15:07:48
 *
 *  Project : KTH Inviwo Modules
 *
 *  License : Follows the Inviwo BSD license model
 *********************************************************************
 */

#include <labstreamlines/labstreamlinesmodule.h>
#include <labstreamlines/eulerrk4comparison.h>
#include <labstreamlines/streamlineintegrator.h>
#include <labstreamlines/utils/vectorfieldgenerator2d.h>
#include <modules/opengl/shader/shadermanager.h>

namespace inviwo
{

LabStreamlinesModule::LabStreamlinesModule(InviwoApplication* app) : InviwoModule(app, "LabStreamlines")
{
    // Add a directory to the search path of the Shadermanager
    ShaderManager::getPtr()->addShaderSearchPath(getPath(ModulePath::GLSL));

    // Register processors
    registerProcessor<EulerRK4Comparison>();
    registerProcessor<StreamlineIntegrator>();
    registerProcessor<VectorFieldGenerator>();
}

} // namespace
