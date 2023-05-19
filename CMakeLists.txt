cmake_minimum_required(VERSION 3.21)

project(SFCComp
  VERSION 0.3.0
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
    set(CMAKE_EXE_LINKER_FLAGS "-static")
  endif()
else()
  message(FATAL_ERROR "Not a supported compiler: ${CMAKE_CXX_COMPILER}")
endif()

set(sfc_comp_src
  src/encode.cpp
  src/huffman.cpp
  src/image.cpp
  src/io.cpp
  src/utility.cpp

  src/action_pachio_comp.cpp
  src/addams_family_comp.cpp
  src/asameshimae_nyanko_comp.cpp
  src/bokujou_monogatari_comp.cpp
  src/cannon_fodder_comp.cpp
  src/cb_chara_wars_comp.cpp
  src/chrono_trigger_comp.cpp
  src/danzarb_comp.cpp
  src/der_langrisser_comp.cpp
  src/derby_stallion_2_comp.cpp
  src/dokapon_comp.cpp
  src/doom_comp.cpp
  src/doraemon_comp.cpp
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
  src/flintstones_comp.cpp
  src/front_mission_comp_2.cpp
  src/gokinjo_boukentai_comp.cpp
  src/gun_hazard_comp.cpp
  src/hal_comp.cpp
  src/hanjuku_hero_comp.cpp
  src/heberekes_popoon_comp.cpp
  src/kamen_rider_sd_comp.cpp
  src/kiki_kaikai_comp.cpp
  src/knights_of_the_round_comp.cpp
  src/koei_comp.cpp
  src/konami_comp.cpp
  src/legend_comp.cpp
  src/lennus_2_comp.cpp
  src/live_a_live_comp_1.cpp
  src/madara2_comp.cpp
  src/mahoujin_guru_guru_comp.cpp
  src/marios_super_picross_comp.cpp
  src/marvelous_comp.cpp
  src/mujintou_monogatari_comp.cpp
  src/nba_jam_comp.cpp
  src/oscar_comp.cpp
  src/pac_in_time_comp.cpp
  src/papuwa_comp.cpp
  src/picross_np_comp.cpp
  src/pokemon_gold_comp.cpp
  src/popful_mail_comp.cpp
  src/power_piggs_comp.cpp
  src/rareware_comp.cpp
  src/rayearth_comp.cpp
  src/rob_northen_comp_1.cpp
  src/rob_northen_comp_2.cpp
  src/royal_conquest_comp.cpp
  src/rs3_comp_1.cpp
  src/sailor_moon_comp_1.cpp
  src/sansara_naga2_comp.cpp
  src/sd_gundam_gnext_comp.cpp
  src/sd_gundam_gx_comp.cpp
  src/sd_gundam_x_comp.cpp
  src/seiken_densetsu_2_comp.cpp
  src/shin_megami_tensei2_comp.cpp
  src/sky_mission_comp.cpp
  src/slap_stick_comp.cpp
  src/smash_tv_comp.cpp
  src/soccer_kid_comp.cpp
  src/soul_and_sword_comp.cpp
  src/stargate_comp.cpp
  src/super_4wd_the_baja_comp.cpp
  src/super_bomberman_5_comp.cpp
  src/super_donkey_kong_comp.cpp
  src/super_loopz_comp.cpp
  src/super_mario_rpg_comp.cpp
  src/super_robot_wars_comp.cpp
  src/super_soukoban_comp.cpp
  src/syndicate_comp.cpp
  src/tactics_ogre_comp_1.cpp
  src/tails_of_phantasia_comp.cpp
  src/tamolympic_comp.cpp
  src/tenchi_wo_kurau_comp.cpp
  src/vortex_comp.cpp
  src/wizardry5_comp_1.cpp
  src/wizardry5_comp_2.cpp
  src/wizardry6_comp.cpp
  src/zelda_comp.cpp
)

include_directories("${CMAKE_SOURCE_DIR}/include")

add_library(objlibsfccomp OBJECT ${sfc_comp_src})
set_property(TARGET objlibsfccomp PROPERTY POSITION_INDEPENDENT_CODE 1)

add_library(${lib_sfc_comp_shared} SHARED $<TARGET_OBJECTS:objlibsfccomp>)
add_library(${lib_sfc_comp_static} STATIC $<TARGET_OBJECTS:objlibsfccomp>)

set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)

add_executable(bench bench.cpp)
target_link_libraries(bench ${lib_sfc_comp_static})

set(comp_ext_list
  action_pachio_comp comp
  addams_family_comp comp
  asameshimae_nyanko_comp comp
  asameshimae_nyanko_4bpp_comp comp
  assault_suits_valken_comp comp
  bahamut_lagoon_comp comp
  bahamut_lagoon_comp_fast comp
  battle_cross_comp comp
  battletech_comp comp
  bokujou_monogatari_comp comp
  bounty_sword_comp comp
  brandish_comp comp
  cannon_fodder_comp comp
  cb_chara_wars_comp comp
  chrono_trigger_comp comp
  chrono_trigger_comp_fast comp
  danzarb_comp comp
  dekitate_high_school_comp_1 comp
  dekitate_high_school_comp_2 comp
  der_langrisser_comp comp
  derby_stallion_2_comp comp
  diet_comp comp
  dokapon_comp comp
  doom_comp_1 comp
  doom_comp_2 comp
  doraemon_comp comp
  dq12_comp comp
  dq5_comp_2 comp
  dq6_comp comp
  dragon_knight_4_comp comp
  dragon_knight_4_4bpp_comp comp
  estpolis_biography_comp comp
  famicom_tantei_club_part_ii_comp comp
  fe3_comp comp
  fe4_comp comp
  ff5_comp comp
  ff6_comp comp
  ffusa_comp comp
  final_stretch_comp comp
  flintstones_comp comp
  front_mission_comp_2 comp
  gionbana_comp comp
  gokinjo_boukentai_comp comp
  gun_hazard_comp comp
  hal_comp comp
  hanjuku_hero_comp comp
  heberekes_popoon_comp comp
  ihatovo_monogatari_comp comp
  jurassic_park_comp comp
  kamen_rider_sd_comp comp
  keiba_eight_special_2_comp comp
  keirin_ou_comp comp
  kiki_kaikai_comp comp
  knights_of_the_round_comp comp
  koei_comp comp
  konami_comp_1 comp
  konami_comp_2 comp
  konami_comp_2_r comp
  legend_comp comp
  lemmings_comp comp
  lennus_2_comp comp
  live_a_live_comp_1 comp
  love_quest_comp comp
  madara2_comp comp
  mahoujin_guru_guru_comp comp
  maka_maka_comp comp
  marios_super_picross_comp comp
  marvelous_comp comp
  mujintou_monogatari_comp comp
  nba_jam_comp comp
  odekake_lester_comp comp
  olivias_mystery_comp comp
  oscar_comp comp
  pac_in_time_comp comp
  papuwa_comp comp
  picross_np_comp comp
  pokemon_gold_comp comp
  popful_mail_comp comp
  power_piggs_comp comp
  rareware_comp comp
  rayearth_comp comp
  riddick_bowe_boxing_comp comp
  rob_northen_comp_1 comp
  rob_northen_comp_2 comp
  royal_conquest_comp comp
  rs3_comp_1 comp
  sailor_moon_comp_1 comp
  sansara_naga2_comp comp
  sd_gundam_gnext_comp comp
  sd_gundam_gx_comp comp
  sd_gundam_x_comp comp
  seiken_densetsu_2_comp comp
  shadowrun_comp comp
  shima_kousaku_comp comp
  shin_megami_tensei2_comp comp
  sky_mission_comp comp
  slap_stick_comp comp
  slayers_comp comp
  smash_tv_comp comp
  smurfs_comp comp
  sotsugyou_bangai_hen_comp comp
  soccer_kid_comp comp
  soul_and_sword_comp comp
  spirou_comp comp
  stargate_comp comp
  super_4wd_the_baja_comp comp
  super_bomberman_5_comp comp
  super_donkey_kong_comp comp
  super_dunk_star_comp comp
  super_jinsei_game_comp comp
  super_loopz_comp comp
  super_mario_rpg_comp comp
  super_robot_wars_comp comp
  super_soukoban_comp comp
  syndicate_comp comp
  tactics_ogre_comp_1 comp
  tactics_ogre_comp_2 comp
  tales_of_phantasia_comp comp
  tamolympic_comp comp
  tenchi_souzou_comp comp
  tenchi_wo_kurau_comp comp
  time_cop_comp comp
  vortex_comp comp
  wild_guns_comp comp
  wizardry5_comp_1 comp
  wizardry5_comp_2 comp
  wizardry6_comp comp
  yatterman_comp comp
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
