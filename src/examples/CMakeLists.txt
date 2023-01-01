cmake_minimum_required(VERSION 3.0.0)

# Find includes in corresponding build directories
set(CMAKE_INCLUDE_CURRENT_DIR ON)

if(EXISTS "${CMAKE_SOURCE_DIR}/.git")
	execute_process(
		COMMAND git rev-parse HEAD
		WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
		OUTPUT_VARIABLE GIT_BRANCH
		OUTPUT_STRIP_TRAILING_WHITESPACE
	)
else()
	set(GIT_BRANCH "Unknown")
endif()
add_definitions(-DGIT_BRANCH=${GIT_BRANCH})
#set build dir for your library
set(LIBRARY_OUTPUT_PATH ${CMAKE_BUILD_DIR}/lib)
message( STATUS "LIBRARY_OUTPUT_PATH: " ${LIBRARY_OUTPUT_PATH} )

add_executable(pardis pardis.cpp)
target_link_libraries(pardis ppam jmatrix)

add_executable(parpam parpam.cpp)
target_link_libraries(parpam ppam jmatrix)

add_executable(parsil parsil.cpp)
target_link_libraries(parsil ppam jmatrix)

add_executable(tdvalue tdvalue.cpp)
target_link_libraries(tdvalue ppam jmatrix)

#If your app, links to an external lib -ie Boost
#target_link_libraries( jmatrixlib ${Boost_LIBRARIES} )

# Install
install(TARGETS pardis DESTINATION bin)
install(CODE "execute_process(COMMAND ${CMAKE_COMMAND} -E create_symlink ${CMAKE_INSTALL_PREFIX}/bin/pardis ${CMAKE_INSTALL_PREFIX}/bin/pardisd)")
install(TARGETS parpam DESTINATION bin)
install(CODE "execute_process(COMMAND ${CMAKE_COMMAND} -E create_symlink ${CMAKE_INSTALL_PREFIX}/bin/parpam ${CMAKE_INSTALL_PREFIX}/bin/parpamd)")
install(TARGETS parsil DESTINATION bin)
install(CODE "execute_process(COMMAND ${CMAKE_COMMAND} -E create_symlink ${CMAKE_INSTALL_PREFIX}/bin/parsil ${CMAKE_INSTALL_PREFIX}/bin/parsild)")
install(TARGETS tdvalue DESTINATION bin)

