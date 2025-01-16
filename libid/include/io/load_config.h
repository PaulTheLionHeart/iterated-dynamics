// SPDX-License-Identifier: GPL-3.0-only
//
#pragma once

#include <string>

extern int                   g_cfg_line_nums[];

// Return full path to filename, or empty string if not found.
std::string locate_config_file(const std::string &name);

void load_config();
void load_config(const std::string &cfg_path);
