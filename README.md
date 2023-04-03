# SFC Comp

A collection of implementation of compressions mainly used in SFC Games.

Please see [GameList.md](doc/GameList.md) for supported Games.

Compatibility with [LunarCompress](https://fusoya.eludevisibility.org/lc/index.html) is described at [LunarCompress.md](doc/LunarCompress.md).

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
| bahamut_lagoon_comp                      |     9 | 24430 | 14176 | 32718 |  9694 | 14430 | 27163 |  6648 | 129268 | 0.1098 | 3CEAFFE1 |
| bokujou_monogatari_comp                  |     7 | 25004 | 14051 | 33168 |  9453 | 14849 | 27809 |  7225 | 131566 | 0.0819 | 4D976D4D |
| chrono_trigger_comp                      |     9 | 24430 | 14056 | 32718 |  9458 | 14430 | 27163 |  6648 | 128912 | 0.1500 | 0DD1A27C |
| diet_comp                                |    22 | 22245 | 12582 | 30687 |  8203 | 12620 | 24631 |  6489 | 117479 | 0.1616 | 29EF3CCD |
| dokapon_comp                             |     6 | 23983 | 13940 | 31828 |  9434 | 14033 | 26017 |  6689 | 125930 | 0.6151 | 6C1C92FB |
| doom_comp_1                              |     9 | 23696 | 13550 | 30207 |  8950 | 13706 | 25671 |  6797 | 122586 | 0.4377 | E6A4C014 |
| doom_comp_2                              |     7 | 24139 | 14304 | 30205 |  9915 | 14370 | 25669 |  6795 | 125404 | 0.2572 | A0182B5B |
| dq12_comp                                |     5 | 24423 | 14169 | 32710 |  9693 | 14424 | 27153 |  6643 | 129220 | 0.0762 | 9E3F4EF0 |
| dq6_comp                                 |     3 | 25688 | 14401 | 33638 |  9875 | 15169 | 28400 |  8489 | 135663 | 0.0870 | EF916FB5 |
| dragon_knight_4_comp_1                   |     4 | 24678 | 13941 | 32770 |  9433 | 14593 | 27388 |  7217 | 130024 | 0.0826 | 2CA7FBC4 |
| estpolis_biography_comp                  |     5 | 23478 | 13366 | 31256 |  8797 | 13507 | 26323 |  5962 | 122694 | 0.0956 | CC9E43C3 |
| famicom_tantei_club_part_ii_comp         |     5 | 24426 | 14169 | 32713 |  9693 | 14424 | 27158 |  6643 | 129231 | 0.0779 | 23FCFD3F |
| fe3_comp                                 |     4 | 23551 | 13412 | 30818 |  9716 | 14144 | 26259 |  6101 | 124005 | 0.3645 | 2635D057 |
| fe4_comp                                 |     4 | 22270 | 12763 | 29833 |  9019 | 13588 | 24654 |  5436 | 117567 | 0.1778 | 69FDF671 |
| ff5_comp                                 |     5 | 25020 | 14050 | 33168 |  9457 | 14847 | 27812 |  7239 | 131598 | 0.0848 | F94B1331 |
| ff6_comp                                 |     5 | 25002 | 14049 | 33166 |  9451 | 14847 | 27807 |  7223 | 131550 | 0.0810 | 1C4BF2B6 |
| ffusa_comp                               |     6 | 25181 | 15386 | 32282 | 12510 | 14844 | 27388 | 11305 | 138902 | 0.1012 | 555BBE73 |
| front_mission_comp_2                     |     5 | 25900 | 14694 | 34101 |  9511 | 15546 | 28894 |  6641 | 135292 | 0.1191 | FA516D58 |
| hal_comp                                 |     4 | 23589 | 13547 | 30810 | 10181 | 14677 | 26442 |  6556 | 125806 | 0.3861 | 2BA7727B |
| hanjuku_hero_comp                        |    11 | 24418 | 14151 | 32712 |  9692 | 14424 | 27127 |  6649 | 129184 | 0.0859 | 8EC736EB |
| live_a_live_comp_1                       |     5 | 24823 | 15027 | 30716 | 11685 | 14135 | 26503 |  7707 | 130601 | 0.1651 | 3141C171 |
| koei_comp                                |     5 | 21895 | 12134 | 30632 |  8096 | 12369 | 24554 |  6225 | 115910 | 0.2101 | E344CF3E |
| konami_comp_1                            |     5 | 23914 | 14232 | 30474 | 10451 | 14671 | 25860 |  7630 | 127237 | 0.0915 | 3A7D7CB0 |
| konami_comp_2                            |     5 | 23906 | 14139 | 30474 | 10451 | 14671 | 25860 |  7605 | 127111 | 0.0904 | DBF893B8 |
| madara2_comp                             |     4 | 22708 | 13517 | 29768 | 10594 | 14364 | 24851 |  7593 | 123399 | 0.1266 | 3381A7EF |
| marios_super_picross_comp                |     5 | 24423 | 14169 | 32710 |  9693 | 14424 | 27153 |  6643 | 129220 | 0.0741 | 43D4AA0A |
| marvelous_comp                           |     4 | 24252 | 14148 | 31052 | 10233 | 14682 | 26820 |  6560 | 127751 | 1.3319 | 51CFB1CD |
| papuwa_comp                              |     5 | 23761 | 13277 | 30810 |  8715 | 13980 | 26296 |  6503 | 123347 | 0.4522 | 7776B6CB |
| pokemon_gold_comp                        |     4 | 22979 | 13032 | 30593 | 10041 | 14177 | 26125 |  6058 | 123009 | 0.5871 | 07FC9A7F |
| popful_mail_comp                         |     6 | 24033 | 13745 | 30810 |  9053 | 14278 | 26388 |  6501 | 124814 | 0.6177 | 29C3E45A |
| rareware_comp                            |    41 | 22224 | 12737 | 30608 |  8709 | 13202 | 25221 |  5721 | 118463 | 1.9261 | B91B7DA6 |
| rayearth_comp                            |     6 | 25154 | 14655 | 31247 | 10547 | 14686 | 27108 |  8860 | 132263 | 0.1159 | DF56CDD8 |
| rob_northen_comp_1                       |    26 | 21295 | 11784 | 29422 |  7675 | 11970 | 24046 |  5164 | 111382 | 1.3903 | 29D648E8 |
| rob_northen_comp_2                       |    23 | 22773 | 12667 | 30191 |  8110 | 12896 | 25150 |  6517 | 118327 | 0.1467 | F93802D2 |
| rs3_comp_1                               |     5 | 24429 | 14172 | 32714 |  9693 | 14426 | 27159 |  6644 | 129242 | 0.0736 | 50255019 |
| sailor_moon_comp_1                       |     7 | 23084 | 12800 | 31746 |  8081 | 12994 | 25786 |  6256 | 120754 | 0.1107 | D74045B9 |
| sansara_naga2_comp                       |     9 | 24118 | 13774 | 30929 |  9821 | 14199 | 26412 |  8506 | 127768 | 0.1434 | C70FEA18 |
| sd_gundam_gnext_comp                     |     4 | 25159 | 14182 | 32770 |  9456 | 14808 | 27612 |  6432 | 130423 | 0.0786 | 6D69932E |
| sd_gundam_gx_comp                        |    10 | 21366 | 11622 | 28786 |  7864 | 12286 | 23586 |  5652 | 111172 | 0.1912 | 13F6F07C |
| sd_gundam_x_comp                         |     6 | 21528 | 11843 | 28913 |  8530 | 12736 | 23736 |  7274 | 114566 | 0.1833 | 124A77D8 |
| seiken_densetsu_2_comp                   |     7 | 24525 | 14812 | 30704 | 10624 | 14351 | 26111 |  7846 | 128980 | 0.4879 | 3F9A5966 |
| shin_megami_tensei2_comp                 |     4 | 23879 | 14110 | 30444 | 10451 | 14658 | 25840 |  7601 | 126987 | 0.1129 | 82F04553 |
| slap_stick_comp                          |     5 | 25134 | 14773 | 33207 | 11432 | 14495 | 27671 | 11230 | 137947 | 0.1023 | F4566128 |
| slayers_comp                             |    12 | 21184 | 11594 | 28744 |  7652 | 12230 | 23616 |  5270 | 110302 | 0.1889 | 9DA68A22 |
| soul_and_sword_comp                      |     6 | 23531 | 13196 | 30625 |  8621 | 13805 | 25796 |  6726 | 122306 | 0.1240 | 5061A2D4 |
| super_donkey_kong_comp                   |   130 | 23658 | 13533 | 30447 |  9018 | 14401 | 25949 |  6250 | 123386 | 0.7593 | CFB46D6A |
| super_mario_rpg_comp                     |    11 | 24435 | 14178 | 32720 |  9699 | 14432 | 27165 |  6650 | 129290 | 0.0754 | 360F25DB |
| super_robot_wars_comp                    |     6 | 23083 | 12800 | 31746 |  8081 | 12993 | 25785 |  6255 | 120749 | 0.1097 | E5939466 |
| syndicate_comp                           |     7 | 23537 | 13174 | 31245 |  8865 | 13669 | 25389 |  5652 | 121538 | 1.4930 | 3A47CBBE |
| tactics_ogre_comp_1                      |     5 | 23726 | 13594 | 30747 |  9014 | 14004 | 26347 |  5571 | 123008 | 0.1003 | 8641BABA |
| tactics_ogre_comp_2                      |     6 | 23083 | 12800 | 31746 |  8081 | 12993 | 25785 |  6255 | 120749 | 0.1116 | B14985F5 |
| tales_of_phantasia_comp                  |    11 | 24321 | 13728 | 32706 |  9022 | 14433 | 27114 |  6527 | 127862 | 0.1431 | 1ED98940 |
| vortex_comp                              |     7 | 22493 | 12874 | 29956 |  8788 | 12680 | 24850 |  6236 | 117884 | 0.5456 | A27CC256 |
| wizardry5_comp_1                         |     3 | 24178 | 13729 | 32193 |  9775 | 14279 | 26839 |  9429 | 130425 | 0.0962 | 54494D6A |
| wizardry5_comp_2                         |     3 | 25692 | 14401 | 33639 |  9875 | 15170 | 28401 |  8490 | 135671 | 0.0842 | EE6D3436 |
| zelda_comp_1                             |     4 | 24252 | 14157 | 31052 | 10233 | 14682 | 26820 |  6560 | 127760 | 0.2243 | EF2A502A |
| zelda_comp_2                             |     4 | 24252 | 14157 | 31052 | 10233 | 14682 | 26820 |  6560 | 127760 | 0.2243 | A48D2149 |

## LICENSE

MIT License