file(REMOVE_RECURSE
  "libliggghts.pdb"
  "libliggghts.so"
)

# Per-language clean rules from dependency scanning.
foreach(lang CXX)
  include(CMakeFiles/liggghts_shared.dir/cmake_clean_${lang}.cmake OPTIONAL)
endforeach()
