#include <chrono>
#include <filesystem>
#include <functional>

#include "sfc_comp.hpp"

void benchmark(const std::string& path) {
  // Ref
  // - https://stackoverflow.com/questions/664014/what-integer-hash-function-are-good-that-accepts-an-integer-hash-key
  // - boost hash_combine
  const auto hash = [](std::span<const uint8_t> vec) -> size_t {
    uint32_t seed = vec.size();
    for (auto x : vec) {
      x = ((x >> 16) ^ x) * 0x45d9f3b;
      x = ((x >> 16) ^ x) * 0x45d9f3b;
      x = (x >> 16) ^ x;
      seed ^= x + 0x9e3779b9 + (seed << 6) + (seed >> 2);
    }
    return seed;
  };

  using namespace sfc_comp;

  struct CompFunc {
    using func_type = std::function<std::vector<uint8_t>(std::span<const uint8_t>)>;
    CompFunc(const char* name, func_type&& func) : name(name), func(func) {}
    const char* name;
    const func_type& func;
  };

#define C(x) CompFunc(#x, x)
  const auto comps = std::vector<CompFunc>({
    C(action_pachio_comp),
    C(addams_family_comp),
    C(asameshimae_nyanko_comp),
    // C(asameshimae_nyanko_4bpp_comp),
    C(assault_suits_valken_comp),
    C(asuka_kenzan_comp_1),
    C(asuka_kenzan_comp_2),
    C(bahamut_lagoon_comp),
    C(bahamut_lagoon_comp_fast),
    C(battle_cross_comp),
    C(battletech_comp),
    C(bokujou_monogatari_comp),
    C(bounty_sword_comp),
    C(brandish_comp),
    C(cannon_fodder_comp),
    C(cb_chara_wars_comp),
    C(chrono_trigger_comp),
    C(chrono_trigger_comp_fast),
    C(danzarb_comp),
    // C(dekitate_high_school_comp_1),
    // C(dekitate_high_school_comp_2),
    C(der_langrisser_comp),
    C(derby_stallion_2_comp),
    C(diet_comp),
    C(dokapon_comp),
    C(doom_comp_1),
    C(doom_comp_2),
    C(doraemon_comp),
    C(dq12_comp),
    C(dq5_comp_2),
    C(dq6_comp),
    C(dragon_knight_4_comp),
    // C(dragon_knight_4_4bpp_comp),
    C(estpolis_biography_comp),
    C(famicom_tantei_club_part_ii_comp),
    C(fe3_comp),
    C(fe4_comp),
    C(ff5_comp),
    C(ff6_comp),
    C(ffusa_comp),
    C(final_stretch_comp),
    C(flintstones_comp),
    C(front_mission_comp_2),
    C(gionbana_comp),
    C(gokinjo_boukentai_comp),
    C(gun_hazard_comp),
    C(hal_comp),
    C(hanjuku_hero_comp),
    C(heberekes_popoon_comp),
    C(ihatovo_monogatari_comp),
    C(jurassic_park_comp),
    C(kamen_rider_sd_comp),
    // C(keiba_eight_special_2_comp),
    // C(keirin_ou_comp),
    C(kiki_kaikai_comp),
    C(knights_of_the_round_comp),
    C(koei_comp),
    C(konami_comp_1),
    C(konami_comp_2),
    // C(konami_comp_2_r),
    C(legend_comp),
    C(lemmings_comp),
    C(lennus_2_comp),
    C(live_a_live_comp_1),
    C(love_quest_comp),
    C(madara2_comp),
    C(mahoujin_guru_guru_comp),
    C(maka_maka_comp),
    C(marios_super_picross_comp),
    C(marvelous_comp),
    C(mujintou_monogatari_comp),
    C(nba_jam_comp),
    C(odekake_lester_comp),
    C(olivias_mystery_comp),
    C(oscar_comp),
    C(pac_in_time_comp),
    C(papuwa_comp),
    C(picross_np_comp),
    C(pokemon_gold_comp),
    C(popful_mail_comp),
    C(power_piggs_comp),
    C(rareware_comp),
    C(rayearth_comp),
    C(riddick_bowe_boxing_comp),
    C(rob_northen_comp_1),
    C(rob_northen_comp_2),
    C(royal_conquest_comp),
    C(rs3_comp_1),
    C(sailor_moon_comp_1),
    C(sansara_naga2_comp),
    C(sd_gundam_gnext_comp),
    C(sd_gundam_gx_comp),
    C(sd_gundam_x_comp),
    C(seiken_densetsu_2_comp),
    C(shadowrun_comp),
    C(shima_kousaku_comp),
    C(shin_megami_tensei2_comp),
    C(sky_mission_comp),
    C(slap_stick_comp),
    C(slayers_comp),
    C(smash_tv_comp),
    C(smurfs_comp),
    C(soccer_kid_comp),
    C(sotsugyou_bangai_hen_comp),
    C(soul_and_sword_comp),
    C(spirou_comp),
    C(stargate_comp),
    C(super_4wd_the_baja_comp),
    C(super_bomberman_5_comp),
    C(super_donkey_kong_comp),
    // C(super_dunk_star_comp),
    C(super_jinsei_game_comp),
    C(super_loopz_comp),
    C(super_mario_rpg_comp),
    C(super_robot_wars_comp),
    C(super_soukoban_comp),
    C(syndicate_comp),
    C(tactics_ogre_comp_1),
    C(tactics_ogre_comp_2),
    C(tales_of_phantasia_comp),
    C(tamolympic_comp),
    C(tenchi_souzou_comp),
    C(tenchi_wo_kurau_comp),
    C(time_cop_comp),
    C(vortex_comp),
    C(wild_guns_comp),
    C(wizardry5_comp_1),
    C(wizardry5_comp_2),
    C(wizardry6_comp),
    C(yatterman_comp),
    C(zelda_comp_1),
    C(zelda_comp_2)
  });
#undef C

  using namespace std::chrono;

  std::vector<std::vector<uint8_t>> inputs;
  std::vector<std::string> paths;
  std::vector<size_t> orders;

  [&] {
    for (const auto& p : std::filesystem::recursive_directory_iterator(path)) {
      if (p.is_directory()) continue;
      const std::string path = p.path().string();
      orders.push_back(inputs.size());
      inputs.emplace_back(io::load(path));
      paths.emplace_back(p.path().filename().string());
    }
    std::sort(orders.begin(), orders.end(), [&paths](const size_t a, const size_t b) {
      return paths[a] < paths[b];
    });
    printf("| Compression |");
    for (const auto& i : orders) printf(" %s |", paths[i].c_str());
    puts(" Total Size | Running Time | Hash |");
    printf("|");
    for (size_t i = 0; i < inputs.size() + 4; ++i) printf(" :--- |");
    puts("");
    size_t s = 0;
    printf("| %-40s |", "**uncompressed**");
    for (const size_t i : orders) {
      printf(" %4zX |", inputs[i].size());
      s += inputs[i].size();
    }
    printf(" %6zX | ------ | -------- |\n", s);
  }();

  uint32_t total_hash = 0;
  size_t total_size_sum = 0;

  const auto p_comps = comps.begin();
  for (size_t c = 0; c < comps.size(); ++c) {
    const auto& comp = p_comps[c];
    printf("| %-40s | ", comp.name);

    const auto beg = high_resolution_clock::now();
    size_t total_size = 0;

    uint32_t h = 0;
    for (const size_t i : orders) {
      try {
        const auto res = comp.func(inputs[i]);
        const auto sz = res.size();
        h ^= hash(res);
        total_size += sz;
        printf("%4zX | ", sz);
      } catch (const std::exception& e) {
        printf(" ERR | ");
      }
    }
    total_hash ^= h;
    total_size_sum += total_size;
    const auto end = high_resolution_clock::now();
    printf("%6zX | %.4f | %08X |\n", total_size,
            duration_cast<nanoseconds>(end - beg).count() / 1e9, h);
  }
  printf("%zX : %08X\n", total_size_sum, total_hash);
}

int main(int argc, char** argv) {
  using namespace std::chrono;
  if (argc < 2) {
    printf("Usage: %s <input-dir>\n", argv[0]);
    return 1;
  } else {
    const auto beg = high_resolution_clock::now();
    benchmark(argv[1]);
    const auto end = high_resolution_clock::now();
    printf("%.4f seconds.\n", duration_cast<nanoseconds>(end - beg).count() / 1e9);
  }
  return 0;
}
