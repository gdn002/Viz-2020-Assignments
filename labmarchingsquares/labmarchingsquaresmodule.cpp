/*********************************************************************
 *  Author  : Himangshu Saikia
 *  Init    : Monday, September 11, 2017 - 12:57:56
 *
 *  Project : KTH Inviwo Modules
 *
 *  License : Follows the Inviwo BSD license model
 *********************************************************************
 */

#include <labmarchingsquares/labmarchingsquaresmodule.h>
#include <labmarchingsquares/marchingsquares.h>

namespace inviwo {

LabMarchingSquaresModule::LabMarchingSquaresModule(InviwoApplication* app)
    : InviwoModule(app, "LabMarchingSquares")

{
    registerProcessor<MarchingSquares>();
}

}  // namespace inviwo
