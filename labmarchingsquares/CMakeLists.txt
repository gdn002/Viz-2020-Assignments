#--------------------------------------------------------------------
# Inviwo LabMarchingSquares Module
ivw_module(LabMarchingSquares)

#--------------------------------------------------------------------
# Add header files
set(HEADER_FILES
    #${CMAKE_CURRENT_SOURCE_DIR}/processor.h
    ${CMAKE_CURRENT_SOURCE_DIR}/marchingsquares.h
    ${CMAKE_CURRENT_SOURCE_DIR}/utils/setminmaxdatamap.h
)
#~ ivw_group("Header Files" ${HEADER_FILES})

#--------------------------------------------------------------------
# Add source files
set(SOURCE_FILES
    #${CMAKE_CURRENT_SOURCE_DIR}/processor.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/marchingsquares.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/utils/setminmaxdatamap.cpp
)
ivw_group("Sources" ${SOURCE_FILES} ${HEADER_FILES})


#--------------------------------------------------------------------
# Add shaders
set(SHADER_FILES
    #${CMAKE_CURRENT_SOURCE_DIR}/glsl/processor.frag
)
ivw_group("Shader Files" ${SHADER_FILES})


#--------------------------------------------------------------------
# Add Unittests
set(TEST_FILES
	#${CMAKE_CURRENT_SOURCE_DIR}/tests/unittests/module-test.cpp
)
ivw_add_unittest(${TEST_FILES})

#--------------------------------------------------------------------
# Create module
ivw_create_module(${SOURCE_FILES} ${HEADER_FILES} ${SHADER_FILES})

#--------------------------------------------------------------------
# Add shader directory to pack
# ivw_add_to_module_pack(${CMAKE_CURRENT_SOURCE_DIR}/glsl)
