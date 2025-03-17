// SPDX-License-Identifier: GPL-3.0-only
//
#include <io/loadmap.h>

#include "expected_map.h"
#include "test_data.h"

#include <misc/ValueSaver.h>
#include <ui/rotate.h>

#include <gtest/gtest.h>


TEST(TestValidateLuts, loadMap)
{
    ValueSaver saved_map_name{g_map_name};
    g_map_name = ID_TEST_MAP_DIR "/foo.map";

    const bool result{validate_luts(ID_TEST_MAP_FILE)};

    EXPECT_FALSE(result);
    for (int i = 0; i < 256; ++i)
    {
        EXPECT_EQ(g_expected_map[i][0], g_dac_box[i][0]) << "index " << i << " red";
        EXPECT_EQ(g_expected_map[i][1], g_dac_box[i][1]) << "index " << i << " green";
        EXPECT_EQ(g_expected_map[i][2], g_dac_box[i][2]) << "index " << i << " blue";
    }
}
