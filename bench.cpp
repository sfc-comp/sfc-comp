#include <chrono>
#include <filesystem>

#include "sfc_comp.hpp"

#define P(x) std::pair(#x, x)

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

  const auto comps = {
    P(addams_family_comp),
    P(assault_suits_valken_comp),
    P(bahamut_lagoon_comp),
    P(battletech_comp),
    P(bokujou_monogatari_comp),
    P(brandish_comp),
    P(burai_comp),
    P(cb_chara_wars_comp),
    P(chrono_trigger_comp),
    P(danzarb_comp),
    P(diet_comp),
    P(dokapon_comp),
    P(doom_comp_1),
    P(doom_comp_2),
    P(doraemon_comp),
    P(dq12_comp),
    P(dq6_comp),
    P(dragon_knight_4_comp_1),
    P(estpolis_biography_comp),
    P(famicom_tantei_club_part_ii_comp),
    P(fe3_comp),
    P(fe4_comp),
    P(ff5_comp),
    P(ff6_comp),
    P(ffusa_comp),
    P(flintstones_comp),
    P(front_mission_comp_2),
    P(gun_hazard_comp),
    P(hal_comp),
    P(hanjuku_hero_comp),
    P(jurassic_park_comp),
    P(kiki_kaikai_comp),
    P(koei_comp),
    P(konami_comp_1),
    P(konami_comp_2),
    P(legend_comp),
    P(live_a_live_comp_1),
    P(madara2_comp),
    P(marios_super_picross_comp),
    P(marvelous_comp),
    P(olivias_mystery_comp),
    P(papuwa_comp),
    P(pokemon_gold_comp),
    P(popful_mail_comp),
    P(power_piggs_comp),
    P(rareware_comp),
    P(rayearth_comp),
    P(riddick_bowe_boxing_comp),
    P(rob_northen_comp_1),
    P(rob_northen_comp_2),
    P(royal_conquest_comp),
    P(rs3_comp_1),
    P(sailor_moon_comp_1),
    P(sansara_naga2_comp),
    P(sd_gundam_gnext_comp),
    P(sd_gundam_gx_comp),
    P(sd_gundam_x_comp),
    P(seiken_densetsu_2_comp),
    P(shadowrun_comp),
    P(shima_kousaku_comp),
    P(shin_megami_tensei2_comp),
    P(sky_mission_comp),
    P(slap_stick_comp),
    P(slayers_comp),
    P(sotsugyou_bangai_hen_comp),
    P(soul_and_sword_comp),
    P(stargate_comp),
    P(super_donkey_kong_comp),
    P(super_jinsei_game_comp),
    P(super_loopz_comp),
    P(super_mario_rpg_comp),
    P(super_robot_wars_comp),
    P(syndicate_comp),
    P(tactics_ogre_comp_1),
    P(tactics_ogre_comp_2),
    P(tales_of_phantasia_comp),
    P(tenchi_souzou_comp),
    P(tenchi_wo_kurau_comp),
    P(time_cop_comp),
    P(vortex_comp),
    P(wild_guns_comp),
    P(wizardry5_comp_1),
    P(wizardry5_comp_2),
    P(yatterman_comp),
    P(zelda_comp_1),
    P(zelda_comp_2)
  };

  using namespace std::chrono;

  std::vector<std::vector<uint8_t>> inputs;
  std::vector<std::string> paths;
  std::vector<size_t> orders;

  {
    printf("| Compression |");
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
    for (const auto& i : orders) printf(" %s |", paths[i].c_str());
    puts(" Total Size | Running Time | Hash |");

    printf("|");
    for (size_t i = 0; i < inputs.size() + 4; ++i) {
      printf(" :---- |");
    }
    puts("");

    size_t s = 0;
    printf("| %-40s |", "**uncompressed**");
    for (const size_t i : orders) {
      printf(" %5zd |", inputs[i].size());
      s += inputs[i].size();
    }
    printf(" %6zd | ------ | -------- |\n", s);
  }

  uint32_t total_hash = 0;
  size_t total_size_sum = 0;

  const auto p_comps = comps.begin();
  for (size_t c = 0; c < comps.size(); ++c) {
    const auto& comp = p_comps[c];
    printf("| %-40s | ", comp.first);

    const auto beg = high_resolution_clock::now();
    size_t total_size = 0;

    uint32_t h = 0;
    for (const size_t i : orders) {
      try {
        const auto res = (comp.second)(inputs[i]);
        const auto sz = res.size();
        h ^= hash(res);
        total_size += sz;
        printf("%5zd | ", sz);
      } catch (const std::exception& e) {
        printf("ERROR | ");
      }
    }
    total_hash ^= h;
    total_size_sum += total_size;

    const auto end = high_resolution_clock::now();
    printf("%6zd | %.4f | %08X |\n", total_size,
            duration_cast<nanoseconds>(end - beg).count() / 1e9, h);
  }
  printf("%zd : %08X\n", total_size_sum, total_hash);
}

#undef P

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
