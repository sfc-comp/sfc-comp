## Compatibility with Lunar Compress

| Lunar Compress | Corresponding Function           | Note                                                                                                              |
| :------------- | :------------------------------- | :---------------------------------------------------------------------------------------------------------------- |
| LC_LZ1         | zelda_comp_1                     |                                                                                                                   |
| LC_LZ2         | zelda_comp_2                     |                                                                                                                   |
| LC_LZ3         | pokemon_gold_comp                |                                                                                                                   |
| LC_LZ4         | live_a_live_comp_1               | LC_LZ4 does not include the (unused) header value `0x01`.                                                         |
| LC_LZ5         | fe3_comp                         |                                                                                                                   |
| LC_LZ7         | seiken_densetsu_2_comp           |                                                                                                                   |
| LC_LZ8         | super_mario_rpg_comp             |                                                                                                                   |
| LC_LZ9         | estpolis_biography_comp          |                                                                                                                   |
| LC_LZ10        | slap_stick_comp                  |                                                                                                                   |
| LC_LZ11        | bokujou_monogatari_comp          |                                                                                                                   |
| LC_LZ12        | konami_comp_1                    |                                                                                                                   |
| LC_LZ13        | chrono_trigger_comp              | LC_LZ13 uses a tricky technique to improve the compression ratio.                                                 |
| LC_LZ14        | famicom_tantei_club_part_ii_comp |                                                                                                                   |
| LC_LZ15        | vortex_comp                      |                                                                                                                   |
| LC_LZ17        | sailor_moon_comp_1               |                                                                                                                   |
| LC_LZ19        | hal_comp                         |                                                                                                                   |
| LC_LZ20        | marvelous_comp                   | LC_LZ20 (v1.90) seems to have a bug (e.g. It outputs `e4 00 00` when compressing a binary data of `0x00` * 1025). |
