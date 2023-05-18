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
| action_pachio_comp                       |    7 | 5F6D | 375D | 7FC8 | 25DF | 385C | 6A14 | 19F6 |  1F8DE | 0.0790 | 4FCA1284 |
| addams_family_comp                       |    4 | 5B8D | 33A0 | 7BB2 | 2293 | 34ED | 64AD | 2199 |  1E8A9 | 0.0656 | 38D09B0C |
| asameshimae_nyanko_comp                  |    8 | 57F9 | 3228 | 7633 | 2654 | 3278 | 61E0 | 19B7 |  1D4BF | 0.7608 | 831A99E2 |
| assault_suits_valken_comp                |    5 | 5F6E | 375D | 727E | 25DE | 385B | 569E | 19F5 |  1D81A | 0.0742 | E614871E |
| bahamut_lagoon_comp                      |    9 | 5E61 | 373B | 7AA8 | 25DE | 379D | 67CB | 19F5 |  1EF88 | 0.1650 | 1BAF87A9 |
| bahamut_lagoon_comp_fast                 |    9 | 5F6E | 3760 | 7FCE | 25DE | 385E | 6A1B | 19F8 |  1F8F4 | 0.0486 | 206A84D5 |
| battle_cross_comp                        |    9 | 558C | 2F69 | 77AC | 1FA4 | 3056 | 5FEE | 1855 |  1C4E7 | 0.1045 | 52A08FAC |
| battletech_comp                          |    9 | 5949 | 31DD | 75D5 | 208B | 347E | 62AC | 1668 |  1CF21 | 0.1683 | 1D369019 |
| bokujou_monogatari_comp                  |    7 | 61AC | 36E3 | 8190 | 24ED | 3A01 | 6CA1 | 1C39 |  201EE | 0.0505 | 4D3AC714 |
| bounty_sword_comp                        |   17 | 5D71 | 3B53 | 7200 | 3877 | 3C03 | 5D49 | 2CBE |  2095C | 0.0056 | EEEA45E8 |
| brandish_comp                            |    5 | 5C72 | 32E5 | 727E | 1FEE | 35C9 | 569E | 1863 |  1C692 | 0.1420 | C5F72E04 |
| cannon_fodder_comp                       |    7 | 5CC5 | 333B | 7E42 | 2291 | 36BF | 67D2 | 1575 |  1E4E0 | 2.1568 | CF337FF8 |
| cb_chara_wars_comp                       |    6 | 5DE9 | 339B | 7E17 | 21AA | 3746 | 6842 | 186C |  1E93F | 0.2079 | 6106D0C8 |
| chrono_trigger_comp                      |    9 | 5EB8 | 36BC | 7D06 | 24F2 | 37E1 | 688F | 19F5 |  1F1DA | 0.3160 | F70281D8 |
| chrono_trigger_comp_fast                 |    9 | 5F6E | 36E8 | 7FCE | 24F2 | 385E | 6A1B | 19F8 |  1F790 | 0.0971 | BFFF850C |
| danzarb_comp                             |    5 | 5D01 | 35DA | 7B2E | 2458 | 36A4 | 6775 | 167E |  1E6FD | 0.0649 | 70DBB17D |
| der_langrisser_comp                      |    3 | 5C82 | 38D2 | 7DB6 | 2EF2 | 33A2 | 657A | 2D08 |  20823 | 0.0948 | D6709F29 |
| derby_stallion_2_comp                    |    6 | 5845 | 322F | 78DA | 21D7 | 326D | 6162 | 1BFA |  1D4F4 | 0.0731 | D2649EB6 |
| diet_comp                                |   16 | 56E4 | 3127 | 77DF | 200C | 314C | 6037 | 1959 |  1CAE8 | 0.1013 | 1180739F |
| dokapon_comp                             |    6 | 5DAF | 3674 | 7C54 | 24DA | 36D1 | 65A1 | 1A21 |  1EBEA | 0.4126 | 336021E2 |
| doom_comp_1                              |    9 | 5C90 | 34EE | 75FF | 22F6 | 358A | 6447 | 1A8D |  1DEDA | 0.3108 | FAF75DA5 |
| doom_comp_2                              |    7 | 5E4B | 37E0 | 75FD | 26BB | 3822 | 6445 | 1A8B |  1E9DC | 0.1924 | 154FEF9E |
| doraemon_comp                            |   22 | 503F | 2B5F | 6F18 | 1D3D | 2EBB | 5A7D | 13F1 |  1A53E | 0.1948 | 2676C957 |
| dq12_comp                                |    5 | 5F67 | 3759 | 7FC6 | 25DD | 3858 | 6A11 | 19F3 |  1F8C4 | 0.0486 | E594BEE6 |
| dq5_comp_2                               |    5 | 5F67 | 3759 | 7FC6 | 25DD | 3858 | 6A11 | 19F3 |  1F8C4 | 0.0485 | 2499B495 |
| dq6_comp                                 |    3 | 6458 | 3841 | 8366 | 2693 | 3B41 | 6EF0 | 2129 |  211EF | 0.0538 | 9D540063 |
| dragon_knight_4_comp                     |    4 | 6066 | 3675 | 8002 | 24D9 | 3901 | 6AFC | 1C31 |  1FBE8 | 0.0561 | 715D1640 |
| estpolis_biography_comp                  |    5 | 5BB6 | 3436 | 7A18 | 225D | 34C3 | 66D3 | 174A |  1DF46 | 0.0644 | 73B10068 |
| famicom_tantei_club_part_ii_comp         |    5 | 5F6A | 3759 | 7FC9 | 25DD | 3858 | 6A16 | 19F3 |  1F8CF | 0.0477 | 02CCA4EB |
| fe3_comp                                 |    4 | 5BFF | 3464 | 7862 | 25F4 | 3740 | 6693 | 17D5 |  1E465 | 0.2030 | 8EFDD535 |
| fe4_comp                                 |    4 | 56FC | 31DE | 7490 | 2336 | 351B | 6051 | 154E |  1CB5E | 0.1569 | 64403DA7 |
| ff5_comp                                 |    5 | 61BC | 36E2 | 8190 | 24F1 | 39FF | 6CA4 | 1C47 |  2020E | 0.0499 | C0B125C1 |
| ff6_comp                                 |    5 | 61AA | 36E1 | 818E | 24EB | 39FF | 6C9F | 1C37 |  201DE | 0.0490 | 79E0D6B6 |
| ffusa_comp                               |    6 | 625D | 3C1A | 7E1A | 30DE | 39FC | 6AFC | 2C29 |  21E96 | 0.0629 | 201C10A1 |
| final_stretch_comp                       |    7 | 5F69 | 375B | 7FC8 | 25DF | 385A | 6A13 | 19F5 |  1F8D4 | 0.0479 | 28996EA6 |
| flintstones_comp                         |    8 | 5A7A | 324F | 7DD8 | 2192 | 3267 | 66D0 | 1684 |  1DBF6 | 0.1284 | 3351C43C |
| front_mission_comp_2                     |    5 | 652C | 3966 | 8535 | 2527 | 3CBA | 70DE | 19F1 |  2107C | 0.0565 | 65B9674A |
| gionbana_comp                            |    6 | 5F68 | 375A | 7FC7 | 25DE | 3859 | 6A12 | 19F4 |  1F8CC | 0.0482 | E5C6DD50 |
| gokinjo_boukentai_comp                   |   13 | 5D6D | 3B4F | 71FC | 3873 | 3BFF | 5D45 | 2CBA |  2093C | 0.0066 | 7091F696 |
| gun_hazard_comp                          |    5 | 5247 | 2DD9 | 727E | 1FEE | 35C9 | 52ED | 1863 |  1B3AA | 0.2454 | DAA12497 |
| hal_comp                                 |    4 | 5C25 | 34EB | 785A | 27C5 | 3955 | 674A | 199C |  1EB6E | 0.2734 | F44CE937 |
| hanjuku_hero_comp                        |    B | 5F62 | 3747 | 7FC8 | 25DC | 3858 | 69F7 | 19F9 |  1F8A0 | 0.0517 | E12BB805 |
| heberekes_popoon_comp                    |   24 | 6051 | 37EB | 79B7 | 2AEE | 38E6 | 6896 | 280F |  20690 | 0.7931 | 362478DF |
| ihatovo_monogatari_comp                  |    7 | 5F69 | 375B | 7FC8 | 25DF | 385A | 6A13 | 19F5 |  1F8D4 | 0.0486 | 605DAFD3 |
| jurassic_park_comp                       |    4 | 5B8E | 33A0 | 7BB2 | 2294 | 34EE | 64AE | 219A |  1E8AE | 0.0641 | A2ABE18E |
| kamen_rider_sd_comp                      |    5 | 5B0B | 378B | 7D2C | 28A6 | 34B5 | 656E | 23EB |  1F67B | 0.0812 | 9997D2D2 |
| kiki_kaikai_comp                         |    5 | 61AF | 36E4 | 8003 | 24ED | 3A01 | 6CA2 | 1C39 |  20064 | 0.0490 | 0C5D6225 |
| knights_of_the_round_comp                |    9 | 58E0 | 3496 | 74FD | 2560 | 3427 | 61B4 | 1897 |  1D64E | 0.2239 | 3DCB2643 |
| koei_comp                                |    5 | 5588 | 2F65 | 77A8 | 1FA0 | 3052 | 5FEA | 1851 |  1C4C7 | 0.1070 | 7F42FC77 |
| konami_comp_1                            |    5 | 5D6A | 3798 | 770A | 28D3 | 394F | 6504 | 1DCE |  1F105 | 0.0844 | D2D0B020 |
| konami_comp_2                            |    5 | 5D62 | 373B | 770A | 28D3 | 394F | 6504 | 1DB5 |  1F087 | 0.0858 | 6B173DBE |
| legend_comp                              |    B | 5A7D | 329B | 798E | 21BB | 34D2 | 637B | 1917 |  1D9D0 | 0.0593 | BCBA0548 |
| lemmings_comp                            |    5 | 5F67 | 3759 | 7FC6 | 25DD | 3858 | 6A11 | 19F3 |  1F8C4 | 0.0487 | 8AA0B85F |
| lennus_2_comp                            |    B | 5D17 | 33D6 | 7E9A | 2216 | 34C6 | 67CD | 18BA |  1E6F5 | 0.1384 | 42E9EE85 |
| live_a_live_comp_1                       |    5 | 60F7 | 3AB3 | 77FC | 2DA5 | 3737 | 6787 | 1E1B |  1FE29 | 0.1664 | 2C4E0E99 |
| love_quest_comp                          |    6 | 5F6D | 375D | 7FC8 | 25DF | 385C | 6A14 | 19F6 |  1F8DD | 0.0485 | EC21F3BF |
| madara2_comp                             |    4 | 58B4 | 34CD | 7448 | 2962 | 381C | 6113 | 1DA9 |  1E207 | 0.0973 | 9E9AC673 |
| mahoujin_guru_guru_comp                  |    C | 6104 | 44A0 | 7ACB | 5582 | 45EE | 66F3 | 2E81 |  2515F | 0.0109 | C08FEEDD |
| maka_maka_comp                           |    5 | 5F6A | 3759 | 7FC9 | 25DD | 3858 | 6A16 | 19F3 |  1F8CF | 0.0482 | 173E0E7A |
| marios_super_picross_comp                |    5 | 5F67 | 3759 | 7FC6 | 25DD | 3858 | 6A11 | 19F3 |  1F8C4 | 0.0487 | 91CB4DB7 |
| marvelous_comp                           |    4 | 5EBC | 3744 | 794C | 27F9 | 395A | 68C4 | 19A0 |  1F307 | 0.1080 | 6D309E04 |
| mujintou_monogatari_comp                 |    6 | 5F34 | 340E | 7FCA | 21C3 | 3841 | 6A17 | 186C |  1EF99 | 0.0529 | 83F39034 |
| nba_jam_comp                             |    8 | 5B00 | 33BF | 7694 | 226D | 31FF | 63A1 | 175D |  1D4C5 | 0.5779 | F2C10868 |
| odekake_lester_comp                      |    5 | 5F6A | 3759 | 7FC9 | 25DD | 3858 | 6A16 | 19F3 |  1F8CF | 0.0479 | DE22BE36 |
| olivias_mystery_comp                     |    8 | 501F | 2B3F | 6EF8 | 1D1D | 2E9B | 5A5D | 13D1 |  1A444 | 0.1963 | 357F9EC4 |
| oscar_comp                               |    7 | 5954 | 313D | 76EE | 1F58 | 3486 | 64C9 | 14D9 |  1CF06 | 0.1469 | C49FE120 |
| pac_in_time_comp                         |    7 | 575B | 30B3 | 7804 | 2022 | 30EF | 6082 | 1553 |  1C6FF | 0.7240 | B8CA1666 |
| papuwa_comp                              |    5 | 5CD1 | 33DD | 785A | 220B | 369C | 66B8 | 1967 |  1E1D3 | 0.1096 | 502D0E93 |
| picross_np_comp                          |    5 | 5F6B | 375C | 7FC7 | 25DD | 385A | 6A12 | 19F4 |  1F8D0 | 0.0477 | E01DA751 |
| pokemon_gold_comp                        |    4 | 59C3 | 32E8 | 7781 | 2739 | 3761 | 660D | 17AA |  1E081 | 0.3263 | 4C92967B |
| popful_mail_comp                         |    6 | 5DE1 | 35B1 | 785A | 235D | 37C6 | 6714 | 1965 |  1E78E | 0.1073 | 25068780 |
| power_piggs_comp                         |    5 | 5E86 | 383C | 7F2E | 271D | 3847 | 69D5 | 18E5 |  1F813 | 0.0495 | 744ACC3A |
| rareware_comp                            |   29 | 56CD | 31E6 | 7789 | 21FD | 3397 | 6289 | 165E |  1CEE0 | 0.3259 | 8FB6BF30 |
| rayearth_comp                            |    6 | 6242 | 393F | 7A0F | 2933 | 395E | 69E4 | 229C |  204A7 | 0.0926 | 7D54FC20 |
| riddick_bowe_boxing_comp                 |    7 | 5769 | 30B0 | 7955 | 2168 | 31CE | 61DA | 1DD2 |  1D457 | 0.0889 | BFA49721 |
| rob_northen_comp_1                       |   1A | 5331 | 2E0E | 72CC | 1E02 | 2EA6 | 5DDA | 13E9 |  1B290 | 0.6597 | F8E57CDA |
| rob_northen_comp_2                       |   17 | 58F5 | 317B | 75EF | 1FAE | 3260 | 623E | 1975 |  1CE37 | 0.1025 | 3B6442BE |
| royal_conquest_comp                      |  1B6 | 538D | 2DA2 | 70DF | 1F32 | 30D4 | 5C8F | 16A6 |  1B6FF | 0.1645 | 5D55ACEE |
| rs3_comp_1                               |    5 | 5F6D | 375C | 7FCA | 25DD | 385A | 6A17 | 19F4 |  1F8DA | 0.0481 | 1288CF07 |
| sailor_moon_comp_1                       |    7 | 5A2B | 3201 | 7C03 | 1F91 | 32C2 | 64BA | 1870 |  1D7B3 | 0.0711 | FFB54D10 |
| sansara_naga2_comp                       |    9 | 5E36 | 35CE | 78D1 | 265D | 3777 | 672C | 213A |  1F318 | 0.0982 | E24D84E8 |
| sd_gundam_gnext_comp                     |    4 | 6247 | 3766 | 8002 | 24F0 | 39D8 | 6BDC | 1920 |  1FD77 | 0.0473 | E5332361 |
| sd_gundam_gx_comp                        |    A | 5310 | 2D72 | 7066 | 1EAE | 3000 | 5C22 | 15EA |  1B1AC | 0.1314 | 893BA799 |
| sd_gundam_x_comp                         |    6 | 540A | 2E29 | 70CB | 213F | 31B6 | 5C7A | 1C45 |  1BEB8 | 0.1286 | 0C515625 |
| seiken_densetsu_2_comp                   |    7 | 5FCD | 39DC | 77F0 | 2980 | 380F | 65FF | 1EA6 |  1F7D4 | 0.4045 | 7AA4D5D8 |
| shadowrun_comp                           |    7 | 5841 | 3137 | 750F | 2031 | 33E3 | 61A5 | 1632 |  1CA79 | 0.1860 | 35318A1A |
| shima_kousaku_comp                       |    A | 500E | 2B2B | 6EF9 | 1CCF | 2E4C | 5A54 | 128A |  1A235 | 0.1813 | 6BA5BB09 |
| shin_megami_tensei2_comp                 |    4 | 5D47 | 371E | 76EC | 28D3 | 3942 | 64F0 | 1DB1 |  1F00B | 0.0845 | 0F8877FD |
| sky_mission_comp                         |    5 | 5767 | 30AE | 7953 | 2166 | 31CC | 61D8 | 1DD0 |  1D447 | 0.0919 | E4E0EE3B |
| slap_stick_comp                          |    5 | 622E | 39B5 | 81B7 | 2CA8 | 389F | 6C17 | 2BDE |  21ADB | 0.0590 | 10F8BC40 |
| slayers_comp                             |    C | 52BE | 2D48 | 703E | 1DCC | 2FDC | 5C38 | 1472 |  1AEA2 | 0.1312 | 6944CB11 |
| smash_tv_comp                            |    7 | 5947 | 31DB | 75D3 | 2089 | 347C | 62AA | 1666 |  1CF11 | 0.1664 | 8C67A01B |
| smurfs_comp                              |    7 | 5BA6 | 33AC | 7762 | 2290 | 343D | 645F | 2172 |  1E359 | 0.0951 | E775B9B5 |
| soccer_kid_comp                          |    F | 5C96 | 34F3 | 7797 | 22FB | 3590 | 654A | 1AB0 |  1E1B4 | 0.1151 | 4732F03B |
| sotsugyou_bangai_hen_comp                |    3 | 622C | 39B3 | 81B5 | 2CA6 | 389D | 6C15 | 2BDC |  21ACB | 0.0592 | AFF33AB6 |
| soul_and_sword_comp                      |    6 | 5BEB | 338C | 77A1 | 21AD | 35ED | 64C4 | 1A46 |  1DDC2 | 0.1277 | 2A481A3B |
| spirou_comp                              |    7 | 58E5 | 316B | 75DF | 1F9E | 3250 | 622E | 1965 |  1CDB7 | 0.1023 | EFDEE25A |
| stargate_comp                            |    7 | 5841 | 3137 | 750F | 2031 | 33E3 | 61A5 | 1632 |  1CA79 | 0.1873 | B9EA80E2 |
| super_4wd_the_baja_comp                  |    8 | 5D74 | 363E | 77BA | 2560 | 3676 | 6568 | 24A8 |  1F15A | 0.1215 | 1C12447C |
| super_bomberman_5_comp                   |    3 | 637E | 3D14 | 7C32 | 38D5 | 4429 | 6AA1 | 3305 |  2376B | 0.0228 | 28B588CB |
| super_donkey_kong_comp                   |   82 | 5C63 | 34C7 | 76F7 | 232D | 3844 | 656D | 1845 |  1E1C6 | 0.2869 | 0CFFFBEB |
| super_jinsei_game_comp                   |    3 | 5F65 | 3757 | 7FC4 | 25DB | 3856 | 6A0F | 19F1 |  1F8B4 | 0.0482 | 013F98F5 |
| super_loopz_comp                         |   13 | 564C | 2FFA | 7607 | 1F87 | 308A | 61EE | 1611 |  1C470 | 0.1235 | F32A8466 |
| super_mario_rpg_comp                     |    B | 5F73 | 3762 | 7FD0 | 25E3 | 3860 | 6A1D | 19FA |  1F90A | 0.0484 | EF5B76EB |
| super_robot_wars_comp                    |    6 | 5A2B | 3200 | 7C02 | 1F91 | 32C1 | 64B9 | 186F |  1D7AD | 0.0716 | 515ADEC4 |
| super_soukoban_comp                      |    3 | 5F66 | 36E0 | 7FC5 | 24EA | 3857 | 6A10 | 19F2 |  1F751 | 0.2037 | 94F58195 |
| syndicate_comp                           |    7 | 5BF0 | 3375 | 7A0D | 229C | 3565 | 632C | 1613 |  1DAB9 | 1.5399 | 24E69F05 |
| tactics_ogre_comp_1                      |    5 | 5CAE | 351A | 781B | 2336 | 36B4 | 66EB | 15C3 |  1E080 | 0.0794 | CDA2F370 |
| tactics_ogre_comp_2                      |    6 | 5A2B | 3200 | 7C02 | 1F91 | 32C1 | 64B9 | 186F |  1D7AD | 0.0709 | FE5BB3B3 |
| tales_of_phantasia_comp                  |    B | 5F01 | 359D | 7FC2 | 233E | 3861 | 69EA | 197D |  1F371 | 0.1028 | 71A09CCF |
| tamolympic_comp                          |    6 | 5638 | 30C1 | 7515 | 200E | 3054 | 5FA2 | 1826 |  1C43E | 0.1459 | AD43BD72 |
| tenchi_souzou_comp                       |    9 | 5A2E | 3203 | 7C04 | 1F94 | 32C4 | 64BC | 1872 |  1D7C4 | 0.0705 | B27FE747 |
| tenchi_wo_kurau_comp                     |    5 | 5F6D | 375C | 7FCA | 25DD | 385A | 6A17 | 19F4 |  1F8DA | 0.0483 | 3D3A46A0 |
| time_cop_comp                            |    E | 502C | 2B4C | 6F05 | 1D29 | 2EA1 | 5A63 | 13D7 |  1A48F | 0.1913 | A3EC0B08 |
| vortex_comp                              |    7 | 57DD | 324A | 7504 | 2254 | 3188 | 6112 | 185C |  1CC7C | 0.1353 | C6BA86B0 |
| wild_guns_comp                           |    5 | 61AF | 36E4 | 7C57 | 24ED | 3A01 | 6C6F | 1C39 |  1FC85 | 0.0631 | 8F738135 |
| wizardry5_comp_1                         |    3 | 5E72 | 35A1 | 7DC1 | 262F | 37C7 | 68D7 | 24D5 |  1FD79 | 0.0520 | C3802A59 |
| wizardry5_comp_2                         |    3 | 645C | 3841 | 8367 | 2693 | 3B42 | 6EF1 | 212A |  211F7 | 0.0529 | AD3FDD58 |
| wizardry6_comp                           |    5 | 61AD | 36E3 | 818F | 24EB | 3A00 | 6CA0 | 1C38 |  201E7 | 0.0526 | A94614DD |
| yatterman_comp                           |    A | 500A | 2B32 | 6EF8 | 1CD2 | 2E4E | 5A56 | 128E |  1A242 | 0.1794 | 9E59E189 |
| zelda_comp_1                             |    4 | 5EBC | 374D | 794C | 27F9 | 395A | 68C4 | 19A0 |  1F310 | 0.1032 | A663D8E6 |
| zelda_comp_2                             |    4 | 5EBC | 374D | 794C | 27F9 | 395A | 68C4 | 19A0 |  1F310 | 0.1012 | CE5E686B |

## LICENSE

MIT License
