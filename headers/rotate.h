#pragma once

#include <string>

extern BYTE                  g_dac_box[256][3];
extern bool                  g_dac_learn;
extern bool                  g_got_real_dac;    // loaddac worked, really got a dac
extern BYTE                  g_old_dac_box[256][3];
extern std::string           g_map_name;
extern bool                  g_map_set;

extern void rotate(int);
extern void save_palette();
extern bool load_palette();
