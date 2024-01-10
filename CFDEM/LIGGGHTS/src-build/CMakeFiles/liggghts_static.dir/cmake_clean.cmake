file(REMOVE_RECURSE
  "libliggghts.a"
  "libliggghts.pdb"
)

# Per-language clean rules from dependency scanning.
foreach(lang CXX)
  include(CMakeFiles/liggghts_static.dir/cmake_clean_${lang}.cmake OPTIONAL)
endforeach()
