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

## Tools

### *_comp (e.g. fe4_comp)

This tool compresses the given input files.

#### Usage

```bash
$ ./*_comp [file1] [file2] ...
```

**Note**: This tool is supposed to be used by dragging and dropping files.

#### Sample Output

```text
$ ./fe4_comp fe3_1.4bpp
SFC Compress version 0.2.1 (fe4_comp)

   fe3_1.4bpp
-> fe3_1.4bpp.comp
     elapsed time: 0.053 sec
uncompressed size: 8000h Byte(s)
  compressed size: 56FEh Byte(s)
compression ratio: 67.96%

Press Enter to exit.
```

### bench

This tool performs every compression algorithm for each input file in `<input-dir>`.

#### Usage

```bash
$ ./bench <input-dir>
```

#### Sample Output

| Compression | 2bytes.bin | fe3_1.4bpp | fe3_2.4bpp | ff6.4bpp | ff6.map | lal.event | sample.4bpp | sdk2_1.map | Total Size | Running Time | Hash |
| :--- | :--- | :--- | :--- | :--- | :--- | :--- | :--- | :--- | :--- | :--- | :--- |
| **uncompressed**                         |    2 | 8000 | 8000 | 8000 | 8000 | 441B | 7000 | 7380 |  3279D | ------ | -------- |
| action_pachio_comp                       |    7 | 5F6D | 375D | 7FC8 | 25DF | 385C | 6A14 | 19F6 |  1F8DE | 0.0735 | 19743BD3 |
| addams_family_comp                       |    4 | 5B8D | 33A0 | 7BB2 | 2293 | 34ED | 64AD | 2199 |  1E8A9 | 0.0635 | FDCD9EB3 |
| assault_suits_valken_comp                |    5 | 5F6E | 375D | 727E | 25DE | 385B | 569E | 19F5 |  1D81A | 0.0712 | 4D05E986 |
| bahamut_lagoon_comp                      |    9 | 5F6E | 3760 | 7FCE | 25DE | 385E | 6A1B | 19F8 |  1F8F4 | 0.0504 | 3CEAFFE1 |
| battletech_comp                          |    9 | 5949 | 31DD | 75D5 | 208B | 347E | 62AC | 1668 |  1CF21 | 0.1681 | 9014D4AF |
| bokujou_monogatari_comp                  |    7 | 61AC | 36E3 | 8190 | 24ED | 3A01 | 6CA1 | 1C39 |  201EE | 0.0561 | 4D976D4D |
| bounty_sword_comp                        |   17 | 5D71 | 3B53 | 7200 | 3877 | 3C03 | 5D49 | 2CBE |  2095C | 0.0054 | EEEA45E8 |
| brandish_comp                            |    5 | 5C72 | 32E5 | 727E | 1FEE | 35C9 | 569E | 1863 |  1C692 | 0.1566 | 32710632 |
| burai_comp                               |    3 | 5C82 | 38D2 | 7DB6 | 2EF2 | 33A2 | 657A | 2D08 |  20823 | 0.1221 | FEF1E5C2 |
| cb_chara_wars_comp                       |    6 | 5DE9 | 339B | 7E17 | 21AA | 3746 | 6842 | 186C |  1E93F | 0.0995 | 8947D6F0 |
| chrono_trigger_comp                      |    9 | 5F6E | 36E8 | 7FCE | 24F2 | 385E | 6A1B | 19F8 |  1F790 | 0.0982 | 0DD1A27C |
| danzarb_comp                             |    5 | 5D01 | 35DA | 7B2E | 2458 | 36A4 | 6775 | 167E |  1E6FD | 0.0478 | D5B00688 |
| diet_comp                                |   16 | 56E5 | 3126 | 77DF | 200B | 314C | 6037 | 1959 |  1CAE7 | 0.1058 | 29EF3CCD |
| dokapon_comp                             |    6 | 5DAF | 3674 | 7C54 | 24DA | 36D1 | 65A1 | 1A21 |  1EBEA | 0.4508 | 6C1C92FB |
| doom_comp_1                              |    9 | 5C90 | 34EE | 75FF | 22F6 | 358A | 6447 | 1A8D |  1DEDA | 0.3654 | E6A4C014 |
| doom_comp_2                              |    7 | 5E4B | 37E0 | 75FD | 26BB | 3822 | 6445 | 1A8B |  1E9DC | 0.2043 | A0182B5B |
| doraemon_comp                            |   22 | 503E | 2B66 | 6F18 | 1D49 | 2EBB | 5A7D | 13D8 |  1A537 | 0.2248 | 3EEC5029 |
| dq12_comp                                |    5 | 5F67 | 3759 | 7FC6 | 25DD | 3858 | 6A11 | 19F3 |  1F8C4 | 0.0473 | 9E3F4EF0 |
| dq5_comp_2                               |    5 | 5F67 | 3759 | 7FC6 | 25DD | 3858 | 6A11 | 19F3 |  1F8C4 | 0.0474 | 4F9C5071 |
| dq6_comp                                 |    3 | 6458 | 3841 | 8366 | 2693 | 3B41 | 6EF0 | 2129 |  211EF | 0.0584 | EF916FB5 |
| dragon_knight_4_comp_1                   |    4 | 6066 | 3675 | 8002 | 24D9 | 3901 | 6AFC | 1C31 |  1FBE8 | 0.0513 | 2CA7FBC4 |
| estpolis_biography_comp                  |    5 | 5BB6 | 3436 | 7A18 | 225D | 34C3 | 66D3 | 174A |  1DF46 | 0.0649 | CC9E43C3 |
| famicom_tantei_club_part_ii_comp         |    5 | 5F6A | 3759 | 7FC9 | 25DD | 3858 | 6A16 | 19F3 |  1F8CF | 0.0473 | 23FCFD3F |
| fe3_comp                                 |    4 | 5BFF | 3464 | 7862 | 25F4 | 3740 | 6693 | 17D5 |  1E465 | 0.2505 | DDB2F2D4 |
| fe4_comp                                 |    4 | 56FE | 31DB | 7489 | 233B | 3514 | 604E | 153C |  1CB3F | 0.1393 | 69FDF671 |
| ff5_comp                                 |    5 | 61BC | 36E2 | 8190 | 24F1 | 39FF | 6CA4 | 1C47 |  2020E | 0.0513 | F94B1331 |
| ff6_comp                                 |    5 | 61AA | 36E1 | 818E | 24EB | 39FF | 6C9F | 1C37 |  201DE | 0.0513 | 1C4BF2B6 |
| ffusa_comp                               |    6 | 625D | 3C1A | 7E1A | 30DE | 39FC | 6AFC | 2C29 |  21E96 | 0.0696 | 555BBE73 |
| final_stretch_comp                       |    7 | 5F69 | 375B | 7FC8 | 25DF | 385A | 6A13 | 19F5 |  1F8D4 | 0.0473 | 7D52679A |
| flintstones_comp                         |    8 | 5A7A | 324F | 7DD8 | 2192 | 3267 | 66D0 | 1684 |  1DBF6 | 0.1653 | E2A6853F |
| front_mission_comp_2                     |    5 | 652C | 3966 | 8535 | 2527 | 3CBA | 70DE | 19F1 |  2107C | 0.0926 | FA516D58 |
| gionbana_comp                            |    6 | 5F68 | 375A | 7FC7 | 25DE | 3859 | 6A12 | 19F4 |  1F8CC | 0.0483 | 5C830F84 |
| gokinjo_boukentai_comp                   |   13 | 5D6D | 3B4F | 71FC | 3873 | 3BFF | 5D45 | 2CBA |  2093C | 0.0063 | 7091F696 |
| gun_hazard_comp                          |    5 | 5247 | 2DD9 | 727E | 1FEE | 35C9 | 52ED | 1863 |  1B3AA | 0.2654 | ED5DB7F4 |
| hal_comp                                 |    4 | 5C25 | 34EB | 785A | 27C5 | 3955 | 674A | 199C |  1EB6E | 0.2644 | 61E4F832 |
| hanjuku_hero_comp                        |    B | 5F62 | 3747 | 7FC8 | 25DC | 3858 | 69F7 | 19F9 |  1F8A0 | 0.0542 | 8EC736EB |
| heberekes_popoon_comp                    |   24 | 606A | 37C9 | 79B5 | 2B00 | 38E3 | 68A0 | 27FA |  20689 | 0.6649 | D3B73ED6 |
| ihatovo_monogatari_comp                  |    7 | 5F69 | 375B | 7FC8 | 25DF | 385A | 6A13 | 19F5 |  1F8D4 | 0.0507 | 0391DCB9 |
| jurassic_park_comp                       |    4 | 5B8E | 33A0 | 7BB2 | 2294 | 34EE | 64AE | 219A |  1E8AE | 0.0599 | 34541B52 |
| kamen_rider_sd_comp                      |    5 | 5B0B | 378B | 7D2C | 28A7 | 34B5 | 656E | 23EB |  1F67C | 0.1207 | 01EDB505 |
| kiki_kaikai_comp                         |    5 | 61AF | 36E4 | 8003 | 24ED | 3A01 | 6CA2 | 1C39 |  20064 | 0.0503 | E4A4C89E |
| knights_of_the_round_comp                |    9 | 58DF | 3496 | 74FD | 256C | 342C | 61B4 | 1897 |  1D65E | 0.2115 | 188B914B |
| koei_comp                                |    5 | 5587 | 2F66 | 77A8 | 1FA0 | 3051 | 5FEA | 1851 |  1C4C6 | 0.1134 | E344CF3E |
| konami_comp_1                            |    5 | 5D6A | 3798 | 770A | 28D3 | 394F | 6504 | 1DCE |  1F105 | 0.0634 | 3A7D7CB0 |
| konami_comp_2                            |    5 | 5D62 | 373B | 770A | 28D3 | 394F | 6504 | 1DB5 |  1F087 | 0.0642 | DBF893B8 |
| legend_comp                              |    B | 5A7D | 3298 | 798E | 21BB | 34D2 | 637B | 1916 |  1D9CC | 0.0493 | 96E3F84C |
| lemmings_comp                            |    5 | 5F67 | 3759 | 7FC6 | 25DD | 3858 | 6A11 | 19F3 |  1F8C4 | 0.0491 | 5B1A42EE |
| lennus_2_comp                            |    B | 5D17 | 33D6 | 7E9A | 2216 | 34C6 | 67CD | 18BA |  1E6F5 | 0.1628 | 2EF5A986 |
| live_a_live_comp_1                       |    5 | 60F7 | 3AB3 | 77FC | 2DA5 | 3737 | 6787 | 1E1B |  1FE29 | 0.1375 | BA6526A6 |
| love_quest_comp                          |    6 | 5F6D | 375D | 7FC8 | 25DF | 385C | 6A14 | 19F6 |  1F8DD | 0.0481 | BA9FDAE8 |
| madara2_comp                             |    4 | 58B4 | 34CD | 7448 | 2962 | 381C | 6113 | 1DA9 |  1E207 | 0.0998 | 3381A7EF |
| mahoujin_guru_guru_comp                  |    C | 6104 | 44A0 | 7ACB | 5582 | 45EE | 66F3 | 2E81 |  2515F | 0.0048 | 4AD757B1 |
| maka_maka_comp                           |    5 | 5F6A | 3759 | 7FC9 | 25DD | 3858 | 6A16 | 19F3 |  1F8CF | 0.0489 | 9725EE93 |
| marios_super_picross_comp                |    5 | 5F67 | 3759 | 7FC6 | 25DD | 3858 | 6A11 | 19F3 |  1F8C4 | 0.0467 | 43D4AA0A |
| marvelous_comp                           |    4 | 5EBC | 3744 | 794C | 27F9 | 395A | 68C4 | 19A0 |  1F307 | 0.1474 | 4AC48E73 |
| mujintou_monogatari_comp                 |    5 | 5F33 | 340D | 7FC9 | 21C2 | 3840 | 6A16 | 186B |  1EF91 | 0.0475 | 9AFCB643 |
| odekake_lester_comp                      |    5 | 5F6A | 3759 | 7FC9 | 25DD | 3858 | 6A16 | 19F3 |  1F8CF | 0.0472 | 724C89A4 |
| olivias_mystery_comp                     |    8 | 501E | 2B46 | 6EF8 | 1D29 | 2E9B | 5A5D | 13B8 |  1A43D | 0.2168 | 5FA8A160 |
| oscar_comp                               |    7 | 5954 | 313D | 76EE | 1F58 | 3486 | 64C9 | 14D9 |  1CF06 | 0.1566 | 4C3F30B1 |
| papuwa_comp                              |    5 | 5CD1 | 33DD | 785A | 220B | 369C | 66B8 | 1967 |  1E1D3 | 0.1516 | FC58730D |
| picross_np_comp                          |    5 | 5F6B | 375C | 7FC7 | 25DD | 385A | 6A12 | 19F4 |  1F8D0 | 0.0475 | 06D9CAB2 |
| pokemon_gold_comp                        |    4 | 59C3 | 32E8 | 7781 | 2739 | 3761 | 660D | 17AA |  1E081 | 0.4326 | 6156B892 |
| popful_mail_comp                         |    6 | 5DE1 | 35B1 | 785A | 235D | 37C6 | 6714 | 1965 |  1E78E | 0.0909 | 8E51349B |
| power_piggs_comp                         |    5 | 5E86 | 383C | 7F2E | 271D | 3847 | 69D5 | 18E5 |  1F813 | 0.0443 | 84DF3AE0 |
| rareware_comp                            |   29 | 56EC | 31B3 | 7787 | 2205 | 3391 | 6285 | 1657 |  1CEC1 | 0.5235 | F3FE034F |
| rayearth_comp                            |    6 | 6242 | 393F | 7A0F | 2933 | 395E | 69E4 | 229C |  204A7 | 0.0890 | DF56CDD8 |
| riddick_bowe_boxing_comp                 |    7 | 5769 | 30B0 | 7955 | 2168 | 31CE | 61DA | 1DD2 |  1D457 | 0.0925 | 6C9B46B7 |
| rob_northen_comp_1                       |   1A | 533B | 2E08 | 72C9 | 1DFB | 2EA7 | 5DF0 | 13E8 |  1B2A0 | 0.7413 | A46D379A |
| rob_northen_comp_2                       |   17 | 58F5 | 317B | 75EF | 1FAE | 3260 | 623E | 1975 |  1CE37 | 0.1129 | F93802D2 |
| royal_conquest_comp                      |  1B6 | 538D | 2D9B | 70DF | 1F05 | 30D4 | 5C8F | 169D |  1B6C2 | 0.1883 | 1BAB55C1 |
| rs3_comp_1                               |    5 | 5F6D | 375C | 7FCA | 25DD | 385A | 6A17 | 19F4 |  1F8DA | 0.0461 | 50255019 |
| sailor_moon_comp_1                       |    7 | 5A2C | 3200 | 7C02 | 1F91 | 32C2 | 64BA | 1870 |  1D7B2 | 0.0824 | D74045B9 |
| sansara_naga2_comp                       |    9 | 5E36 | 35CE | 78D1 | 265D | 3777 | 672C | 213A |  1F318 | 0.1131 | C70FEA18 |
| sd_gundam_gnext_comp                     |    4 | 6247 | 3766 | 8002 | 24F0 | 39D8 | 6BDC | 1920 |  1FD77 | 0.0465 | 6D69932E |
| sd_gundam_gx_comp                        |    A | 536E | 2D62 | 7066 | 1EAA | 2FFE | 5C22 | 15F2 |  1B1FC | 0.1648 | 71518908 |
| sd_gundam_x_comp                         |    6 | 5412 | 2E3E | 70DA | 2142 | 31C1 | 5C84 | 1C55 |  1BF0C | 0.1576 | 3E503482 |
| seiken_densetsu_2_comp                   |    7 | 5FCD | 39DC | 77F0 | 2980 | 380F | 65FF | 1EA6 |  1F7D4 | 0.3739 | 3F9A5966 |
| shadowrun_comp                           |    7 | 5947 | 31DB | 75D3 | 2089 | 347C | 62AA | 1666 |  1CF11 | 0.1608 | 6EE99D2E |
| shima_kousaku_comp                       |    A | 500E | 2B2D | 6EF9 | 1CD6 | 2E4C | 5A53 | 12A1 |  1A254 | 0.2138 | BA0850E7 |
| shin_megami_tensei2_comp                 |    4 | 5D47 | 371E | 76EC | 28D3 | 3942 | 64F0 | 1DB1 |  1F00B | 0.0874 | 82F04553 |
| sky_mission_comp                         |    5 | 5767 | 30AE | 7953 | 2166 | 31CC | 61D8 | 1DD0 |  1D447 | 0.0944 | A87A5246 |
| slap_stick_comp                          |    5 | 622E | 39B5 | 81B7 | 2CA8 | 389F | 6C17 | 2BDE |  21ADB | 0.0694 | F4566128 |
| slayers_comp                             |    C | 52BE | 2D4A | 7048 | 1DD2 | 2FC4 | 5C40 | 146E |  1AEA0 | 0.1548 | EA68BBA3 |
| smurfs_comp                              |    7 | 5BA6 | 33AC | 7762 | 2290 | 343D | 645F | 2172 |  1E359 | 0.1100 | 6DFBEF98 |
| soccer_kid_comp                          |    F | 5C96 | 34F3 | 7797 | 22FB | 3590 | 654A | 1AB0 |  1E1B4 | 0.1716 | 8E32C3CA |
| sotsugyou_bangai_hen_comp                |    3 | 622C | 39B3 | 81B5 | 2CA6 | 389D | 6C15 | 2BDC |  21ACB | 0.0705 | 4F745F24 |
| soul_and_sword_comp                      |    6 | 5BEB | 338C | 77A1 | 21AD | 35ED | 64C4 | 1A46 |  1DDC2 | 0.1016 | 4A7448C5 |
| spirou_comp                              |    7 | 58E5 | 316B | 75DF | 1F9E | 3250 | 622E | 1965 |  1CDB7 | 0.1182 | FED6862A |
| stargate_comp                            |    7 | 5841 | 3137 | 750F | 2031 | 33E3 | 61A5 | 1632 |  1CA79 | 0.1577 | 6174A0CF |
| super_4wd_the_baja_comp                  |    8 | 5D74 | 363E | 77BA | 2560 | 3676 | 6568 | 24A8 |  1F15A | 0.1483 | D5E347CD |
| super_bomberman_5_comp                   |    3 | 637E | 3D14 | 7C32 | 38D5 | 4429 | 6AA1 | 3305 |  2376B | 0.0091 | 2862D91F |
| super_donkey_kong_comp                   |   82 | 5C6A | 34DD | 76EF | 233A | 3841 | 655D | 186A |  1E1FA | 0.1377 | CB1196EC |
| super_jinsei_game_comp                   |    3 | 5F65 | 3757 | 7FC4 | 25DB | 3856 | 6A0F | 19F1 |  1F8B4 | 0.0473 | D6574806 |
| super_loopz_comp                         |   13 | 564C | 2FFA | 7607 | 1F87 | 308A | 61EE | 1611 |  1C470 | 0.1227 | 0B362E3A |
| super_mario_rpg_comp                     |    B | 5F73 | 3762 | 7FD0 | 25E3 | 3860 | 6A1D | 19FA |  1F90A | 0.0464 | 360F25DB |
| super_robot_wars_comp                    |    6 | 5A2B | 3200 | 7C02 | 1F91 | 32C1 | 64B9 | 186F |  1D7AD | 0.0822 | E5939466 |
| super_soukoban_comp                      |    3 | 5F66 | 36E0 | 7FC5 | 24EA | 3857 | 6A10 | 19F2 |  1F751 | 0.2214 | 88E348A2 |
| syndicate_comp                           |    7 | 5BF0 | 3375 | 7A0D | 229C | 3565 | 632C | 1613 |  1DAB9 | 0.7518 | 061B0EE9 |
| tactics_ogre_comp_1                      |    5 | 5CAE | 351A | 781B | 2336 | 36B4 | 66EB | 15C3 |  1E080 | 0.0711 | 8641BABA |
| tactics_ogre_comp_2                      |    6 | 5A2B | 3200 | 7C02 | 1F91 | 32C1 | 64B9 | 186F |  1D7AD | 0.0812 | B14985F5 |
| tales_of_phantasia_comp                  |    B | 5F01 | 359D | 7FC2 | 233E | 3861 | 69EA | 197D |  1F371 | 0.0927 | 90B371AB |
| tamolympic_comp                          |    6 | 5639 | 30C2 | 7514 | 200F | 3054 | 5FA2 | 1826 |  1C440 | 0.1663 | 7BFCF590 |
| tenchi_souzou_comp                       |    9 | 5A2E | 3203 | 7C04 | 1F94 | 32C4 | 64BC | 1872 |  1D7C4 | 0.0813 | C3BD936E |
| tenchi_wo_kurau_comp                     |    5 | 5F6D | 375C | 7FCA | 25DD | 385A | 6A17 | 19F4 |  1F8DA | 0.0492 | C5E98D96 |
| time_cop_comp                            |    E | 5024 | 2B4C | 6EFE | 1D30 | 2EA2 | 5A64 | 13BE |  1A470 | 0.2193 | E5F2AF3B |
| vortex_comp                              |    7 | 57DD | 324A | 7504 | 2254 | 3188 | 6112 | 185C |  1CC7C | 0.1493 | A27CC256 |
| wild_guns_comp                           |    5 | 61AF | 36E4 | 7C57 | 24ED | 3A01 | 6C6F | 1C39 |  1FC85 | 0.0531 | DA1D75A4 |
| wizardry5_comp_1                         |    3 | 5E72 | 35A1 | 7DC1 | 262F | 37C7 | 68D7 | 24D5 |  1FD79 | 0.0647 | 54494D6A |
| wizardry5_comp_2                         |    3 | 645C | 3841 | 8367 | 2693 | 3B42 | 6EF1 | 212A |  211F7 | 0.0579 | EE6D3436 |
| wizardry6_comp                           |    5 | 61AD | 36E3 | 818F | 24EB | 3A00 | 6CA0 | 1C38 |  201E7 | 0.0537 | 1DF9E6FE |
| yatterman_comp                           |    A | 500A | 2B36 | 6EF8 | 1CD6 | 2E4E | 5A56 | 129A |  1A256 | 0.2156 | ADD17EC0 |
| zelda_comp_1                             |    4 | 5EBC | 374D | 794C | 27F9 | 395A | 68C4 | 19A0 |  1F310 | 0.1325 | EC940778 |
| zelda_comp_2                             |    4 | 5EBC | 374D | 794C | 27F9 | 395A | 68C4 | 19A0 |  1F310 | 0.1311 | 5633014D |

## LICENSE

MIT License
