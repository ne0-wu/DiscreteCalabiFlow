cmake_minimum_required(VERSION 3.25)

set(CMAKE_CXX_STANDARD 11)

project(DiscreteCalabiFlow)

find_package(Eigen3 REQUIRED)
find_package(OpenMesh REQUIRED)

add_compile_definitions(_USE_MATH_DEFINES)

add_executable(DiscreteCalabiFlow
	main.cpp
	Mesh.h
	CalabiFlow.h
	CalabiFlow.cpp
	stb_image_write.h
	DrawMesh.h
)

target_link_libraries(DiscreteCalabiFlow OpenMeshCore OpenMeshTools)

add_custom_command(TARGET DiscreteCalabiFlow POST_BUILD
    COMMAND ${CMAKE_COMMAND} -E copy_directory
        ${CMAKE_SOURCE_DIR}/data
        $<TARGET_FILE_DIR:DiscreteCalabiFlow>/data
    COMMENT "Copying directory to output directory"
)