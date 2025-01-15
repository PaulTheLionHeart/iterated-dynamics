// SPDX-License-Identifier: GPL-3.0-only
//
#pragma once

#include <cstdio>
#include <string>

enum class ItemType
{
    PAR_SET = 0,
    FORMULA = 1,
    L_SYSTEM = 2,
    IFS = 3,
};

bool find_file_item(std::string &filename, char const *item_name, std::FILE **file_ptr, ItemType item_type);

long get_file_entry(ItemType type, char const *title, char const *fn_key_mask, std::string &filename,
    std::string &entry_name);

bool search_for_entry(std::FILE *infile, char const *item_name);