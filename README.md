# SFC Compress

A collection of implementation of compressions mainly used in SFC/SNES Games.

Please see [GameList.md](doc/GameList.md) for supported Games.

Compatibility with [Lunar Compress](https://fusoya.eludevisibility.org/lc/index.html) is described at [LunarCompress.md](doc/LunarCompress.md).

## Build Requirements

- gcc >= 10.1 (due to `-std=c++20`)
- CMake >= 3.21

## Build

```bash
$ mkdir build
$ cd build
$ cmake ..
$ make
```

## Programs

### *_comp (e.g. fe4_comp)

This program will output compressed files for the given input files.

#### Usage

```bash
$ ./*_comp [file1] [file2] ...
```

**Note**: This tool is supposed to be used by dragging and dropping files.

#### Sample Output

```text
$ ./fe4_comp fe3_1.4bpp
SFC Compress version 0.1.0

   fe3_1.4bpp
-> fe3_1.4bpp.comp
     elapsed time: 0.053 sec
uncompressed size: 8000h Byte(s)
  compressed size: 56FEh Byte(s)
compression ratio: 67.96%

Press Enter to exit.
```

### bench

This program will perform every compression for each input file in `<input-dir>`.

#### Usage

```bash
$ ./bench <input-dir>
```

#### Sample Output

| Compression | 2bytes.bin | fe3_1.4bpp | fe3_2.4bpp | ff6.4bpp | ff6.map | lal.event | sample.4bpp | sdk2_1.map | Total Size | Running Time | Hash |
| :---- | :---- | :---- | :---- | :---- | :---- | :---- | :---- | :---- | :---- | :---- | :---- |
| **uncompressed**                         |     2 | 32768 | 32768 | 32768 | 32768 | 17435 | 28672 | 29568 | 206749 | ------ | -------- |
| addams_family_comp                       |     4 | 23437 | 13216 | 31666 |  8851 | 13549 | 25773 |  8601 | 125097 | 0.0687 | FDCD9EB3 |
| bahamut_lagoon_comp                      |     9 | 24430 | 14176 | 32718 |  9694 | 14430 | 27163 |  6648 | 129268 | 0.0742 | 3CEAFFE1 |
| battletech_comp                          |     9 | 22857 | 12765 | 30165 |  8331 | 13438 | 25260 |  5736 | 118561 | 0.3214 | 3CFD6EC6 |
| bokujou_monogatari_comp                  |     7 | 25004 | 14051 | 33168 |  9453 | 14849 | 27809 |  7225 | 131566 | 0.0573 | 4D976D4D |
| cb_chara_wars_comp                       |     6 | 24041 | 13211 | 32279 |  8618 | 14150 | 26690 |  6252 | 125247 | 0.2330 | 8947D6F0 |
| chrono_trigger_comp                      |     9 | 24430 | 14056 | 32718 |  9458 | 14430 | 27163 |  6648 | 128912 | 0.1054 | 0DD1A27C |
| danzarb_comp                             |     5 | 23809 | 13786 | 31534 |  9304 | 13988 | 26485 |  5758 | 124669 | 0.0485 | D5B00688 |
| diet_comp                                |    22 | 22245 | 12582 | 30687 |  8203 | 12620 | 24631 |  6489 | 117479 | 0.1354 | 29EF3CCD |
| dokapon_comp                             |     6 | 23983 | 13940 | 31828 |  9434 | 14033 | 26017 |  6689 | 125930 | 0.4478 | 6C1C92FB |
| doom_comp_1                              |     9 | 23696 | 13550 | 30207 |  8950 | 13706 | 25671 |  6797 | 122586 | 0.3600 | E6A4C014 |
| doom_comp_2                              |     7 | 24139 | 14304 | 30205 |  9915 | 14370 | 25669 |  6795 | 125404 | 0.2046 | A0182B5B |
| doraemon_comp                            |    34 | 20552 | 11139 | 28438 |  7478 | 11961 | 23191 |  5122 | 107915 | 0.2196 | 00B96132 |
| dq12_comp                                |     5 | 24423 | 14169 | 32710 |  9693 | 14424 | 27153 |  6643 | 129220 | 0.0492 | 9E3F4EF0 |
| dq6_comp                                 |     3 | 25688 | 14401 | 33638 |  9875 | 15169 | 28400 |  8489 | 135663 | 0.0608 | EF916FB5 |
| dragon_knight_4_comp_1                   |     4 | 24678 | 13941 | 32770 |  9433 | 14593 | 27388 |  7217 | 130024 | 0.0518 | 2CA7FBC4 |
| estpolis_biography_comp                  |     5 | 23478 | 13366 | 31256 |  8797 | 13507 | 26323 |  5962 | 122694 | 0.0665 | CC9E43C3 |
| famicom_tantei_club_part_ii_comp         |     5 | 24426 | 14169 | 32713 |  9693 | 14424 | 27158 |  6643 | 129231 | 0.0488 | 23FCFD3F |
| fe3_comp                                 |     4 | 23551 | 13412 | 30818 |  9716 | 14144 | 26259 |  6101 | 124005 | 0.3126 | 2635D057 |
| fe4_comp                                 |     4 | 22270 | 12763 | 29833 |  9019 | 13588 | 24654 |  5436 | 117567 | 0.1449 | 69FDF671 |
| ff5_comp                                 |     5 | 25020 | 14050 | 33168 |  9457 | 14847 | 27812 |  7239 | 131598 | 0.0538 | F94B1331 |
| ff6_comp                                 |     5 | 25002 | 14049 | 33166 |  9451 | 14847 | 27807 |  7223 | 131550 | 0.0543 | 1C4BF2B6 |
| ffusa_comp                               |     6 | 25181 | 15386 | 32282 | 12510 | 14844 | 27388 | 11305 | 138902 | 0.0727 | 555BBE73 |
| front_mission_comp_2                     |     5 | 25900 | 14694 | 34101 |  9511 | 15546 | 28894 |  6641 | 135292 | 0.0996 | FA516D58 |
| hal_comp                                 |     4 | 23589 | 13547 | 30810 | 10181 | 14677 | 26442 |  6556 | 125806 | 0.3090 | 2BA7727B |
| hanjuku_hero_comp                        |    11 | 24418 | 14151 | 32712 |  9692 | 14424 | 27127 |  6649 | 129184 | 0.0560 | 8EC736EB |
| jurassic_park_comp                       |     4 | 23438 | 13216 | 31666 |  8852 | 13550 | 25774 |  8602 | 125102 | 0.0607 | 34541B52 |
| koei_comp                                |     5 | 21895 | 12134 | 30632 |  8096 | 12369 | 24554 |  6225 | 115910 | 0.1770 | E344CF3E |
| konami_comp_1                            |     5 | 23914 | 14232 | 30474 | 10451 | 14671 | 25860 |  7630 | 127237 | 0.0667 | 3A7D7CB0 |
| konami_comp_2                            |     5 | 23906 | 14139 | 30474 | 10451 | 14671 | 25860 |  7605 | 127111 | 0.0675 | DBF893B8 |
| legend_comp                              |    11 | 23165 | 12952 | 31118 |  8635 | 13522 | 25467 |  6422 | 121292 | 0.0537 | 96E3F84C |
| live_a_live_comp_1                       |     5 | 24823 | 15027 | 30716 | 11685 | 14135 | 26503 |  7707 | 130601 | 0.1329 | 3141C171 |
| madara2_comp                             |     4 | 22708 | 13517 | 29768 | 10594 | 14364 | 24851 |  7593 | 123399 | 0.1026 | 3381A7EF |
| marios_super_picross_comp                |     5 | 24423 | 14169 | 32710 |  9693 | 14424 | 27153 |  6643 | 129220 | 0.0485 | 43D4AA0A |
| marvelous_comp                           |     4 | 24252 | 14148 | 31052 | 10233 | 14682 | 26820 |  6560 | 127751 | 1.2515 | 51CFB1CD |
| papuwa_comp                              |     5 | 23761 | 13277 | 30810 |  8715 | 13980 | 26296 |  6503 | 123347 | 0.4384 | 7776B6CB |
| pokemon_gold_comp                        |     4 | 22979 | 13032 | 30593 | 10041 | 14177 | 26125 |  6058 | 123009 | 0.4900 | 07FC9A7F |
| popful_mail_comp                         |     6 | 24033 | 13745 | 30810 |  9053 | 14278 | 26388 |  6501 | 124814 | 0.6011 | 29C3E45A |
| power_piggs_comp                         |     5 | 24198 | 14396 | 32558 | 10013 | 14407 | 27093 |  6373 | 129043 | 0.0477 | 84DF3AE0 |
| rareware_comp                            |    41 | 22224 | 12737 | 30608 |  8709 | 13202 | 25221 |  5721 | 118463 | 1.8188 | B91B7DA6 |
| rayearth_comp                            |     6 | 25154 | 14655 | 31247 | 10547 | 14686 | 27108 |  8860 | 132263 | 0.0914 | DF56CDD8 |
| riddick_bowe_boxing_comp                 |     7 | 22377 | 12464 | 31061 |  8552 | 12750 | 25050 |  7634 | 119895 | 0.1200 | 6C9B46B7 |
| rob_northen_comp_1                       |    26 | 21295 | 11784 | 29422 |  7675 | 11970 | 24046 |  5164 | 111382 | 1.3621 | 29D648E8 |
| rob_northen_comp_2                       |    23 | 22773 | 12667 | 30191 |  8110 | 12896 | 25150 |  6517 | 118327 | 0.1164 | F93802D2 |
| rs3_comp_1                               |     5 | 24429 | 14172 | 32714 |  9693 | 14426 | 27159 |  6644 | 129242 | 0.0511 | 50255019 |
| sailor_moon_comp_1                       |     7 | 23084 | 12800 | 31746 |  8081 | 12994 | 25786 |  6256 | 120754 | 0.0831 | D74045B9 |
| sansara_naga2_comp                       |     9 | 24118 | 13774 | 30929 |  9821 | 14199 | 26412 |  8506 | 127768 | 0.1169 | C70FEA18 |
| sd_gundam_gnext_comp                     |     4 | 25159 | 14182 | 32770 |  9456 | 14808 | 27612 |  6432 | 130423 | 0.0479 | 6D69932E |
| sd_gundam_gx_comp                        |    10 | 21366 | 11622 | 28786 |  7864 | 12286 | 23586 |  5652 | 111172 | 0.1455 | 13F6F07C |
| sd_gundam_x_comp                         |     6 | 21528 | 11843 | 28913 |  8530 | 12736 | 23736 |  7274 | 114566 | 0.1419 | 124A77D8 |
| seiken_densetsu_2_comp                   |     7 | 24525 | 14812 | 30704 | 10624 | 14351 | 26111 |  7846 | 128980 | 0.3821 | 3F9A5966 |
| shadowrun_comp                           |     7 | 22855 | 12763 | 30163 |  8329 | 13436 | 25258 |  5734 | 118545 | 0.3155 | 75D20B18 |
| shima_kousaku_comp                       |    10 | 20494 | 11076 | 28404 |  7379 | 11852 | 23148 |  4786 | 107149 | 0.2196 | B97826C5 |
| shin_megami_tensei2_comp                 |     4 | 23879 | 14110 | 30444 | 10451 | 14658 | 25840 |  7601 | 126987 | 0.0891 | 82F04553 |
| sky_mission_comp                         |     5 | 22375 | 12462 | 31059 |  8550 | 12748 | 25048 |  7632 | 119879 | 0.1222 | A87A5246 |
| slap_stick_comp                          |     5 | 25134 | 14773 | 33207 | 11432 | 14495 | 27671 | 11230 | 137947 | 0.0715 | F4566128 |
| slayers_comp                             |    12 | 21184 | 11594 | 28744 |  7652 | 12230 | 23616 |  5270 | 110302 | 0.1413 | 9DA68A22 |
| sotsugyou_bangai_hen_comp                |     3 | 25132 | 14771 | 33205 | 11430 | 14493 | 27669 | 11228 | 137931 | 0.0717 | 4F745F24 |
| soul_and_sword_comp                      |     6 | 23531 | 13196 | 30625 |  8621 | 13805 | 25796 |  6726 | 122306 | 0.1024 | 4A7448C5 |
| stargate_comp                            |     7 | 22593 | 12599 | 29967 |  8241 | 13283 | 24997 |  5682 | 117369 | 1.9343 | 6174A0CF |
| super_donkey_kong_comp                   |   130 | 23658 | 13533 | 30447 |  9018 | 14401 | 25949 |  6250 | 123386 | 0.6418 | CFB46D6A |
| super_jinsei_game_comp                   |     3 | 24421 | 14167 | 32708 |  9691 | 14422 | 27151 |  6641 | 129204 | 0.0495 | D6574806 |
| super_loopz_comp                         |    19 | 22092 | 12282 | 30215 |  8071 | 12426 | 25070 |  5649 | 115824 | 2.0482 | EA30F4E9 |
| super_mario_rpg_comp                     |    11 | 24435 | 14178 | 32720 |  9699 | 14432 | 27165 |  6650 | 129290 | 0.0499 | 360F25DB |
| super_robot_wars_comp                    |     6 | 23083 | 12800 | 31746 |  8081 | 12993 | 25785 |  6255 | 120749 | 0.0845 | E5939466 |
| syndicate_comp                           |     7 | 23537 | 13174 | 31245 |  8865 | 13669 | 25389 |  5652 | 121538 | 1.2444 | 3A47CBBE |
| tactics_ogre_comp_1                      |     5 | 23726 | 13594 | 30747 |  9014 | 14004 | 26347 |  5571 | 123008 | 0.0750 | 8641BABA |
| tactics_ogre_comp_2                      |     6 | 23083 | 12800 | 31746 |  8081 | 12993 | 25785 |  6255 | 120749 | 0.0838 | B14985F5 |
| tales_of_phantasia_comp                  |    11 | 24321 | 13728 | 32706 |  9022 | 14433 | 27114 |  6527 | 127862 | 0.0949 | 1ED98940 |
| tenchi_souzou_comp                       |     9 | 23086 | 12803 | 31748 |  8084 | 12996 | 25788 |  6258 | 120772 | 0.0839 | C3BD936E |
| tenchi_wo_kurau_comp                     |     5 | 24429 | 14172 | 32714 |  9693 | 14426 | 27159 |  6644 | 129242 | 0.0519 | C5E98D96 |
| time_cop_comp                            |    14 | 20526 | 11114 | 28412 |  7452 | 11936 | 23166 |  5096 | 107716 | 0.2429 | 69BFD84A |
| vortex_comp                              |     7 | 22493 | 12874 | 29956 |  8788 | 12680 | 24850 |  6236 | 117884 | 0.5191 | A27CC256 |
| wizardry5_comp_1                         |     3 | 24178 | 13729 | 32193 |  9775 | 14279 | 26839 |  9429 | 130425 | 0.0654 | 54494D6A |
| wizardry5_comp_2                         |     3 | 25692 | 14401 | 33639 |  9875 | 15170 | 28401 |  8490 | 135671 | 0.0588 | EE6D3436 |
| yatterman_comp                           |    10 | 20492 | 11082 | 28404 |  7380 | 11854 | 23148 |  4786 | 107156 | 0.2163 | E43AD286 |
| zelda_comp_1                             |     4 | 24252 | 14157 | 31052 | 10233 | 14682 | 26820 |  6560 | 127760 | 0.1846 | EF2A502A |
| zelda_comp_2                             |     4 | 24252 | 14157 | 31052 | 10233 | 14682 | 26820 |  6560 | 127760 | 0.1864 | A48D2149 |

## LICENSE

MIT License
