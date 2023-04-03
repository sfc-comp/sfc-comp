cmake_minimum_required(VERSION 3.21)

project(SFCComp
  VERSION 0.1.0
  LANGUAGES CXX)

set(PRODUCT_NAME "SFC Compress")
configure_file(
  ${CMAKE_SOURCE_DIR}/src/version.h.in
  ${CMAKE_SOURCE_DIR}/include/version.h)

set(lib_sfc_comp_shared sfccomp)
set(lib_sfc_comp_static sfccomp_static)

if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE Release)
endif()

# set(CMAKE_VERBOSE_MAKEFILE ON)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

if(CMAKE_COMPILER_IS_GNUCXX)
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -Wextra -pedantic")
  # Debug: -g
  set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -O2")
  # Release: -O3 -DNDEBUG
  set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -s")

  if(${CMAKE_SYSTEM_NAME} MATCHES "Windows")
    set(CMAKE_EXE_LINKER_FLAGS "-static-libgcc -static-libstdc++")
  endif()
else()
  message(FATAL_ERROR "Not a supported compiler: ${CMAKE_CXX_COMPILER}")
endif()

set(sfc_comp_src
  src/encode.cpp
  src/huffman.cpp
  src/io.cpp
  src/utility.cpp

  src/bahamut_lagoon_comp.cpp
  src/bokujou_monogatari_comp.cpp
  src/chrono_trigger_comp.cpp
  src/diet_comp.cpp
  src/dokapon_comp.cpp
  src/doom_comp.cpp
  src/dq12_comp.cpp
  src/dq6_comp.cpp
  src/dragon_knight_4_comp.cpp
  src/estpolis_biography_comp.cpp
  src/famicom_tantei_club_part_ii_comp.cpp
  src/fe3_comp.cpp
  src/fe4_comp.cpp
  src/ff5_comp.cpp
  src/ff6_comp.cpp
  src/ffusa_comp.cpp
  src/front_mission_comp_2.cpp
  src/hal_comp.cpp
  src/hanjuku_hero_comp.cpp
  src/koei_comp.cpp
  src/konami_comp.cpp
  src/live_a_live_comp_1.cpp
  src/madara2_comp.cpp
  src/marios_super_picross_comp.cpp
  src/marvelous_comp.cpp
  src/papuwa_comp.cpp
  src/pokemon_gold_comp.cpp
  src/popful_mail_comp.cpp
  src/rareware_comp.cpp
  src/rayearth_comp.cpp
  src/rob_northen_comp_1.cpp
  src/rob_northen_comp_2.cpp
  src/rs3_comp_1.cpp
  src/sailor_moon_comp_1.cpp
  src/sansara_naga2_comp.cpp
  src/sd_gundam_gnext_comp.cpp
  src/sd_gundam_gx_comp.cpp
  src/sd_gundam_x_comp.cpp
  src/seiken_densetsu_2_comp.cpp
  src/shin_megami_tensei2_comp.cpp
  src/slap_stick_comp.cpp
  src/soul_and_sword_comp.cpp
  src/super_donkey_kong_comp.cpp
  src/super_mario_rpg_comp.cpp
  src/super_robot_wars_comp.cpp
  src/syndicate_comp.cpp
  src/tactics_ogre_comp_1.cpp
  src/tails_of_phantasia_comp.cpp
  src/vortex_comp.cpp
  src/wizardry5_comp_1.cpp
  src/wizardry5_comp_2.cpp
  src/zelda_comp.cpp
)

include_directories("${CMAKE_SOURCE_DIR}/include")

add_library(objlibsfccomp OBJECT ${sfc_comp_src})
set_property(TARGET objlibsfccomp PROPERTY POSITION_INDEPENDENT_CODE 1)

add_library(${lib_sfc_comp_static} STATIC $<TARGET_OBJECTS:objlibsfccomp>)
add_library(${lib_sfc_comp_shared} SHARED $<TARGET_OBJECTS:objlibsfccomp>)

add_executable(bench bench.cpp)
target_link_libraries(bench ${lib_sfc_comp_static})

set(comp_ext_list
  bahamut_lagoon_comp comp
  bokujou_monogatari_comp comp
  chrono_trigger_comp comp
  diet_comp comp
  dokapon_comp comp
  doom_comp_1 comp
  doom_comp_2 comp
  dq12_comp comp
  dq6_comp comp
  dragon_knight_4_comp_1 comp
  estpolis_biography_comp comp
  famicom_tantei_club_part_ii_comp comp
  fe3_comp comp
  fe4_comp comp
  ff5_comp comp
  ff6_comp comp
  ffusa_comp comp
  front_mission_comp_2 comp
  hal_comp comp
  hanjuku_hero_comp comp
  live_a_live_comp_1 comp
  koei_comp comp
  konami_comp_1 comp
  konami_comp_2 comp
  madara2_comp comp
  marios_super_picross_comp comp
  marvelous_comp comp
  papuwa_comp comp
  pokemon_gold_comp comp
  popful_mail_comp comp
  rareware_comp comp
  rayearth_comp comp
  rob_northen_comp_1 comp
  rob_northen_comp_2 comp
  rs3_comp_1 comp
  sailor_moon_comp_1 comp
  sansara_naga2_comp comp
  sd_gundam_gnext_comp comp
  sd_gundam_gx_comp comp
  sd_gundam_x_comp comp
  seiken_densetsu_2_comp comp
  shin_megami_tensei2_comp comp
  slap_stick_comp comp
  slayers_comp comp
  soul_and_sword_comp comp
  super_donkey_kong_comp comp
  super_mario_rpg_comp comp
  super_robot_wars_comp comp
  syndicate_comp comp
  tactics_ogre_comp_1 comp
  tactics_ogre_comp_2 comp
  tales_of_phantasia_comp comp
  vortex_comp comp
  wizardry5_comp_1 comp
  wizardry5_comp_2 comp
  zelda_comp_1 comp
  zelda_comp_2 comp
)

list(LENGTH comp_ext_list _list_len)
math(EXPR _max_index "${_list_len} / 2 - 1")

foreach(_index RANGE ${_max_index})
  math(EXPR _comp_i "${_index} * 2 + 0")
  math(EXPR _ext_i "${_index} * 2 + 1")
  list(GET comp_ext_list ${_comp_i} _comp_func)
  list(GET comp_ext_list ${_ext_i} _comp_ext)

  add_executable(${_comp_func} main.cpp)
  target_link_libraries(${_comp_func} ${lib_sfc_comp_static})
  target_compile_definitions(${_comp_func} PRIVATE
    COMP_FUNC=${_comp_func}
    COMP_EXT=.${_comp_ext}
  )
endforeach()