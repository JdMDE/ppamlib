set(HEADERS_INSTALL_PATH ${CMAKE_INSTALL_PREFIX}/include/parallelpam)
message(STATUS "HEADERS_INSTALL_PATH : " ${HEADERS_INSTALL_PATH})

INSTALL (
    DIRECTORY ${CMAKE_SOURCE_DIR}/src/headers/
    DESTINATION include/parallelpam
    FILES_MATCHING PATTERN "*.h*")

#If you place header files in this directory and want them installed
#rememder to add a target for them too.
