include_directories("${PROJECT_SOURCE_DIR}")
add_catch_test(single_particle LIBRARIES optilib ${library_dependencies})
add_catch_test(translation_addition LIBRARIES optilib ${library_dependencies})
