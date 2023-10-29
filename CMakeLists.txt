cmake_minimum_required(VERSION 3.25.0)
set(CMAKE_CXX_STANDARD 17)
project(photon LANGUAGES CXX C)

set(SOURCES
    core/photon.cpp
)

if(WIN32)
    set(PLATFORM
        platforms/win32/win32_photon.cpp
    )
    
	set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -incremental:no")
	
	
	
    add_executable("win32_photon" ${PLATFORM})
	add_library(${PROJECT_NAME} SHARED ${SOURCES})

    	target_link_libraries(${PROJECT_NAME} kernel32.lib user32.lib)
	target_link_libraries("win32_photon" kernel32.lib user32.lib)
	target_link_options(${PROJECT_NAME} PRIVATE "-EXPORT:raytrace")
	set_property(TARGET ${PROJECT_NAME} PROPERTY VS_DEBUGGER_WORKING_DIRECTORY "${PROJECT_SOURCE_DIR}/data")
	set_property(TARGET "win32_photon" PROPERTY VS_DEBUGGER_WORKING_DIRECTORY "${PROJECT_SOURCE_DIR}/data")
    	set_property(DIRECTORY ${PROJECT_SOURCE_DIR} PROPERTY VS_STARTUP_PROJECT "win32_photon")
endif(WIN32)