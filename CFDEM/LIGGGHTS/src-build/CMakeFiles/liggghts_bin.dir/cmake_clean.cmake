file(REMOVE_RECURSE
  "liggghts"
  "liggghts.pdb"
)

# Per-language clean rules from dependency scanning.
foreach(lang CXX)
  include(CMakeFiles/liggghts_bin.dir/cmake_clean_${lang}.cmake OPTIONAL)
endforeach()
