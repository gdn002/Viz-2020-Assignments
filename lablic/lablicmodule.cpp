/*********************************************************************
 *  Author  : Himangshu Saikia
 *  Init    : Friday, September 29, 2017 - 15:50:00
 *
 *  Project : KTH Inviwo Modules
 *
 *  License : Follows the Inviwo BSD license model
 *********************************************************************
 */

#include <lablic/lablicmodule.h>
#include <lablic/licprocessor.h>
#include <lablic/noisetexturegenerator.h>

namespace inviwo
{

LabLICModule::LabLICModule(InviwoApplication* app) : InviwoModule(app, "LabLIC")
{
	registerProcessor<LICProcessor>();
	registerProcessor<NoiseTextureGenerator>();
}

} // namespace
