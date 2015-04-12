macro(precompile_headers target_ private_includes_ global_includes_ definitions_
    compile_flags_ headers relative_path_)
  message(STATUS "sdd=================================")
  message(STATUS "${includes_}")
  message(STATUS "${definitions_}")
  message(STATUS "${compile_flags_}")
  message(STATUS "${headers}")
  message(STATUS "sdd=================================")
  # Get the compiler flags for this build type
  string(TOUPPER "CMAKE_CXX_FLAGS_${CMAKE_BUILD_TYPE}" flags_var_for_the_build_name)
  set(compile_flags "${${flags_var_for_the_build_name}}")
  separate_arguments(compile_flags)
  set(temp_1 "${CMAKE_CXX_FLAGS} -x c++-header")
  separate_arguments(temp_1)
  list(APPEND compile_flags ${temp_1})
  list(APPEND compile_flags ${compile_flags_})

  foreach(item ${private_includes_})
    list(APPEND compile_flags "-I${item}")
  endforeach()

  foreach(item ${global_includes_})
    list(APPEND compile_flags "-isystem${item}")
  endforeach()

  # Get the list of all build-independent preprocessor definitions
  get_directory_property(defines_global COMPILE_DEFINITIONS)
  list(APPEND defines ${defines_global})
  list(APPEND defines ${definitions_})

  # Get the list of all build-dependent preprocessor definitions
  string(TOUPPER "COMPILE_DEFINITIONS_${CMAKE_BUILD_TYPE}" defines_for_build_name)
  get_directory_property(defines_build ${defines_for_build_name})
  list(APPEND defines ${defines_build})

  # Add the "-D" prefix to all of them
  foreach(item ${defines})
    list(APPEND all_define_flags "${item}")
  endforeach()

  list(APPEND compile_flags ${all_define_flags})

  set(pch_include_pathes)
  foreach(header ${headers})
    file(RELATIVE_PATH output_direcotry_path "${relative_path_}" ${header})
    set(output_direcotry_path
      "${CMAKE_CURRENT_BINARY_DIR}/pch/${output_direcotry_path}")
    list(APPEND pch_include_pathes "${output_direcotry_path}.gch")
    target_compile_options(${target_} PRIVATE -include ${output_direcotry_path}.gch)
    get_filename_component(output_direcotry_path ${output_direcotry_path} DIRECTORY)

    file(MAKE_DIRECTORY ${output_direcotry_path})

    get_filename_component(file_name ${header} NAME_WE)
    set(output_file_path "${output_direcotry_path}/${file_name}.hpp.gch")
    # Finally, build the precompiled header.
    # We don't add the build command to add_custom_target
    # because that would force a PCH rebuild even when
    # the ${header_name}.h file hasn't changed. We add it to
    # a special add_custom_command to work around this problem.
    add_custom_target(${file_name}_HPP_GCH ALL DEPENDS ${output_file_path})
    add_custom_command(OUTPUT ${output_file_path}
      COMMAND ${CMAKE_CXX_COMPILER} ${compile_flags} ${header} -o ${output_file_path}
      MAIN_DEPENDENCY ${header}
      WORKING_DIRECTORY ${output_direcotry_path}
      VERBATIM)
    add_dependencies(${target_} ${file_name}_HPP_GCH)
  endforeach()
  list(REMOVE_DUPLICATES pch_include_pathes)

  message(STATUS "sdd=================================")
  message(STATUS "${pch_include_pathes}")
  message(STATUS "sdd=================================")
  # target_include_directories(${target_} PRIVATE ${pch_include_pathes})
endmacro()
