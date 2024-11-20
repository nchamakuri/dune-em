set(noinst_PROGRAMS electromechanical)

#add_dgf_executable(electromechanical electromechanical.cc GRIDDIM 2 WORLDDIM 2 GRIDTYPE YASPGRID)
add_dgf_executable(electromechanical electromechanical.cc)
add_dune_mpi_flags(electromechanical )
dune_symlink_to_source_files(FILES paramL2space)
target_link_libraries(electromechanical PUBLIC ${DUNE_LIBS})
target_link_libraries(electromechanical PUBLIC "/usr/lib64/libboost_system.so.1.66.0")
target_link_libraries(electromechanical PUBLIC "/usr/lib64/libboost_filesystem.so.1.66.0")
target_link_libraries(electromechanical PUBLIC "/usr/lib64/libhdf5.so")
#target_include_directories(electromechanical PUBLIC "/opt/ohpc/pub/libs/gnu8/openmpi3/boost/1.70.0/include/")
#target_include_directories(electromechanical PUBLIC "/opt/ohpc/pub/libs/gnu8/hdf5/1.10.5/include/")
#target_link_libraries(electromechanical PUBLIC "/opt/ohpc/pub/libs/gnu8/openmpi3/boost/1.70.0/lib/libboost_filesystem.so")
#target_link_libraries(electromechanical PUBLIC "/opt/ohpc/pub/libs/gnu8/openmpi3/boost/1.70.0/lib/libboost_system.so")
#target_link_libraries(electromechanical PUBLIC "/opt/ohpc/pub/libs/gnu8/hdf5/1.10.5/lib/libhdf5.so")

set(CLEANFILES  dgfparser.log q1*.vtu *.vtu *.png)

#============================================================
# pass most important options when "make distcheck" is used
set(DISTCHECK_CONFIGURE_FLAGS  --with-dune-common=${DUNE_COMMON_ROOT} --with-dune-grid=${DUNE_GRID_ROOT} --with-dune-istl=${DUNE_ISTL_ROOT} --with-dune-localfunctions=${DUNE_LOCALFUNCTIONS_ROOT} --with-dune-pdelab=${DUNE_PDELAB_ROOT}  CXX="${CXX}" CC="${CC}")
