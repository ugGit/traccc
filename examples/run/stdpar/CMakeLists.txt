# TRACCC library, part of the ACTS project (R&D line)
#
# (c) 2021-2022 CERN for the benefit of the ACTS project
#
# Mozilla Public License Version 2.0

cmake_minimum_required(VERSION 3.18)

#traccc_add_executable( stdpar_par_example "stdpar_par_example.cpp"
#   LINK_LIBRARIES vecmem::core traccc::core traccc::io 
#   traccc::performance traccc::options)

# Set up the executable.
set( exe_name "traccc_stdpar_par_example" )
add_executable( ${exe_name} "stdpar_par_example.cpp" )

set (STDPAR_TARGET "gpu")
target_compile_options(${exe_name} PRIVATE -stdpar=${STDPAR_TARGET})
target_link_options(${exe_name} PRIVATE -stdpar=${STDPAR_TARGET})

# target_link_libraries( ${exe_name} PRIVATE vecmem::core traccc::core traccc::io traccc::performance traccc::options )

set_target_properties( ${exe_name} PROPERTIES EXPORT_NAME stdpar_par_example )
add_executable( traccc::${name} ALIAS ${exe_name} )

# Set up the installation of the executable.
install( TARGETS ${exe_name}
   EXPORT traccc-exports
   RUNTIME DESTINATION "${CMAKE_INSTALL_BINDIR}" )
