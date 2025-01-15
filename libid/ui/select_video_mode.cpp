// SPDX-License-Identifier: GPL-3.0-only
//
#include "ui/select_video_mode.h"

#include "engine/id_data.h"
#include "helpdefs.h"
#include "io/find_path.h"
#include "io/is_writeable.h"
#include "io/load_config.h"
#include "io/save_file.h"
#include "misc/drivers.h"
#include "misc/ValueSaver.h"
#include "ui/full_screen_choice.h"
#include "ui/get_key_no_help.h"
#include "ui/stop_msg.h"
#include "ui/video_mode.h"

#include <algorithm>
#include <cassert>
#include <cstdio>
#include <cstring>
#include <filesystem>
#include <numeric>
#include <string>
#include <vector>

namespace fs = std::filesystem;

static std::vector<int> s_entry_nums;
static bool s_modes_changed{};

static int check_mode_key(int key, int choice);
static bool ent_less(int lhs, int rhs);
static void update_id_cfg();

static void format_vid_table(int choice, char *buf)
{
    char key_name[5];
    const int idx = s_entry_nums[choice];
    assert(idx < g_video_table_len);
    std::memcpy((char *)&g_video_entry, (char *)&g_video_table[idx],
           sizeof(g_video_entry));
    vid_mode_key_name(g_video_entry.key, key_name);
    std::sprintf(buf, "%-5s %-12s %5d %5d %3d  %.12s %.26s", // 34 chars
        key_name, g_video_entry.driver->get_description().c_str(), g_video_entry.x_dots, g_video_entry.y_dots,
        g_video_entry.colors, g_video_entry.driver->get_name().c_str(), g_video_entry.comment);
}

int select_video_mode(int current_mode)
{
    std::vector<int> attributes;

    attributes.resize(g_video_table_len);
    s_entry_nums.resize(g_video_table_len);
    // init tables
    std::fill(attributes.begin(), attributes.end(), 1);
    std::iota(s_entry_nums.begin(), s_entry_nums.end(), 0);
    std::sort(s_entry_nums.begin(), s_entry_nums.end(), ent_less);

    // pick default mode
    if (current_mode < 0)
    {
        g_video_entry.colors = 256;
    }
    else
    {
        std::memcpy((char *) &g_video_entry, (char *) &g_video_table[current_mode], sizeof(g_video_entry));
    }
    int i;
    for (i = 0; i < g_video_table_len; ++i)  // find default mode
    {
        if (g_video_entry.colors == g_video_table[s_entry_nums[i]].colors
            && (current_mode < 0
                || std::memcmp((char *) &g_video_entry, (char *) &g_video_table[s_entry_nums[i]], sizeof(g_video_entry)) == 0))
        {
            break;
        }
    }
    if (i >= g_video_table_len) // no match, default to first entry
    {
        i = 0;
    }

    {
        ValueSaver saved_tab_mode{g_tab_mode, false};
        ValueSaver saved_help_mode{g_help_mode, HelpLabels::HELP_VIDEO_MODE};
        s_modes_changed = false;
        i = full_screen_choice(ChoiceFlags::HELP, "Select Video Mode",
            "key...name..........xdot..ydot.colr.driver......comment......", nullptr, g_video_table_len,
            nullptr, attributes.data(), 1, 16, 74, i, format_vid_table, nullptr, nullptr, check_mode_key);
    }
    if (i == -1)
    {
        // update id.cfg for new key assignments
        if (s_modes_changed
            && g_bad_config == ConfigStatus::OK
            && stop_msg(StopMsgFlags::CANCEL | StopMsgFlags::NO_BUZZER | StopMsgFlags::INFO_ONLY,
                "Save new function key assignments or cancel changes?") == 0)
        {
            update_id_cfg();
        }
        return -1;
    }
    // picked by function key or ENTER key
    i = (i < 0) ? (-1 - i) : s_entry_nums[i];
    // the selected entry now in g_video_entry
    std::memcpy((char *) &g_video_entry, (char *) &g_video_table[i], sizeof(g_video_entry));

    // copy id.cfg table to resident table, note selected entry
    int k = 0;
    for (i = 0; i < g_video_table_len; ++i)
    {
        if (g_video_table[i].key > 0)
        {
            if (std::memcmp((char *)&g_video_entry, (char *)&g_video_table[i], sizeof(g_video_entry)) == 0)
            {
                k = g_video_table[i].key;
            }
        }
    }
    int ret = k;
    if (k == 0)  // selected entry not a copied (assigned to key) one
    {
        std::memcpy((char *)&g_video_table[MAX_VIDEO_MODES-1],
               (char *)&g_video_entry, sizeof(*g_video_table));
        ret = 1400; // special value for check_vidmode_key
    }

    // update id.cfg for new key assignments
    if (s_modes_changed && g_bad_config == ConfigStatus::OK)
    {
        update_id_cfg();
    }

    return ret;
}

static int check_mode_key(int key, int choice)
{
    int i = check_vid_mode_key(key);
    if (i >= 0)
    {
        return -1-i;
    }
    i = s_entry_nums[choice];
    int ret = 0;
    if ((key == '-' || key == '+')
        && (g_video_table[i].key == 0 || g_video_table[i].key >= 1084))
    {
        if (g_bad_config != ConfigStatus::OK)
        {
            stop_msg("Missing or bad id.cfg file. Can't reassign keys.");
        }
        else
        {
            if (key == '-')
            {
                // deassign key?
                if (g_video_table[i].key >= 1084)
                {
                    g_video_table[i].key = 0;
                    s_modes_changed = true;
                }
            }
            else
            {
                // assign key?
                int j = get_a_key_no_help();
                if (j >= 1084 && j <= 1113)
                {
                    for (int k = 0; k < g_video_table_len; ++k)
                    {
                        if (g_video_table[k].key == j)
                        {
                            g_video_table[k].key = 0;
                            ret = -1; // force redisplay
                        }
                    }
                    g_video_table[i].key = j;
                    s_modes_changed = true;
                }
            }
        }
    }
    return ret;
}

static bool ent_less(const int lhs, const int rhs)
{
    int i = g_video_table[lhs].key;
    if (i == 0)
    {
        i = 9999;
    }
    int j = g_video_table[rhs].key;
    if (j == 0)
    {
        j = 9999;
    }
    return i < j || i == j && lhs < rhs;
}

static void update_id_cfg()
{
    char buf[121];
    const std::string cfg_name = find_path("id.cfg");

    if (!is_writeable(cfg_name))
    {
        std::snprintf(buf, std::size(buf), "Can't write %s", cfg_name.c_str());
        stop_msg(buf);
        return;
    }
    const std::string out_name{(fs::path{cfg_name}.parent_path() / "id.tmp").string()};
    std::FILE *outfile = open_save_file(out_name, "w");
    if (outfile == nullptr)
    {
        std::snprintf(buf, std::size(buf), "Can't create %s", out_name.c_str());
        stop_msg(buf);
        return;
    }
    std::FILE *cfg_file = std::fopen(cfg_name.c_str(), "r");

    int next_mode = 0;
    int line_num = next_mode;
    int next_line_num = g_cfg_line_nums[0];
    while (std::fgets(buf, std::size(buf), cfg_file))
    {
        ++line_num;
        // replace this line?
        if (line_num == next_line_num)
        {
            char key_name[5];
            char colors_buff[10];
            VideoInfo video_entry = g_video_table[next_mode];
            vid_mode_key_name(video_entry.key, key_name);
            std::snprintf(colors_buff, std::size(colors_buff), "%3d", video_entry.colors);
            std::fprintf(outfile, "%-4s,%4d,%5d,%s,%s,%s\n",
                    key_name,
                    video_entry.x_dots,
                    video_entry.y_dots,
                    colors_buff,
                    video_entry.driver->get_name().c_str(),
                    video_entry.comment);
            if (++next_mode >= g_video_table_len)
            {
                next_line_num = 32767;
            }
            else
            {
                next_line_num = g_cfg_line_nums[next_mode];
            }
        }
        else
        {
            std::fputs(buf, outfile);
        }
    }

    std::fclose(cfg_file);
    std::fclose(outfile);
    fs::remove(cfg_name);         // success assumed on these lines
    fs::rename(out_name, cfg_name); // since we checked earlier with access
}

void request_video_mode(int &kbd_char)
{
    driver_stack_screen();
    kbd_char = select_video_mode(g_adapter);
    if (check_vid_mode_key(kbd_char) >= 0) // picked a new mode?
    {
        driver_discard_screen();
    }
    else
    {
        driver_unstack_screen();
    }
}