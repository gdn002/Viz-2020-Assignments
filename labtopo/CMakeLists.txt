#--------------------------------------------------------------------
# Inviwo LabTopo Module
ivw_module(LabTopo)

#--------------------------------------------------------------------
# Add header files
set(HEADER_FILES
    topology.h
    utils/gradients.h
    utils/tnt/jama_eig.h
    utils/tnt/tnt_array1d.h
    utils/tnt/tnt_array2d.h
    utils/tnt/tnt_i_refvec.h
    utils/tnt/tnt_math_utils.h
)

#--------------------------------------------------------------------
# Add source files
set(SOURCE_FILES
    topology.cpp
)
ivw_group("Sources" ${SOURCE_FILES} ${HEADER_FILES})

#--------------------------------------------------------------------
# Create module
ivw_create_module(${SOURCE_FILES} ${HEADER_FILES} ${SHADER_FILES})

