# Include {{{
# ============================================================================
if (NOT EIGEN_INCLUDE_DIRS)
  set(EIGEN_INCLUDE_DIRS $ENV{EIGEN_INCLUDE_DIRS})
endif()

if (EIGEN_INCLUDE_DIRS)
  find_path(Eigen_INCLUDE_DIR
    NAMES Eigen/Eigen
    PATHS ${EIGEN_INCLUDE_DIRS}
    NO_DEFAULT_PATH
    PATH_SUFFIXES eigen3
    )
endif()

if (NOT Eigen_INCLUDE_DIR)
  find_path(Eigen_INCLUDE_DIR
    NAMES Eigen/Eigen
    PATH_SUFFIXES eigen3
    )
endif()

if (Eigen_INCLUDE_DIR)
  set(Eigen_INCLUDE_DIRS ${Eigen_INCLUDE_DIR})
endif()
# ============================================================================
# }}} Include

include(FindPackageHandleStandardArgs)

find_package_handle_standard_args(
  "Eigen3"
  REQUIRED_VARS
  Eigen_INCLUDE_DIRS
  )

mark_as_advanced(
  Eigen_INCLUDE_DIR
  )
