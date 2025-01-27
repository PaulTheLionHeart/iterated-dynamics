#include <cmdfiles_test.h>

#include "current_path_saver.h"
#include "test_data.h"

#include <3d.h>
#include <bailout_formula.h>
#include <debug_flags.h>
#include <engine_timer.h>
#include <fractalp.h>
#include <fractype.h>
#include <framain2.h>
#include <history.h>
#include <id.h>
#include <id_data.h>
#include <id_keys.h>
#include <jb.h>
#include <line3d.h>
#include <loadfile.h>
#include <lorenz.h>
#include <make_batch_file.h>
#include <parser.h>
#include <plot3d.h>
#include <rotate.h>
#include <slideshw.h>
#include <soi.h>
#include <sound.h>
#include <stereo.h>
#include <sticky_orbits.h>
#include <stop_msg.h>
#include <trig_fns.h>
#include <value_saver.h>
#include <video_mode.h>

#include <gmock/gmock.h>
#include <gtest/gtest.h>

#include <algorithm>
#include <cctype>
#include <cstring>
#include <string_view>

using namespace testing;

template <typename T>
class ValueUnchanged
{
public:
    ValueUnchanged(const char *label, T &data, T value) :
        m_label(label),
        m_data(data),
        m_saved_value(data),
        m_value(value)
    {
        m_data = value;
    }
    ~ValueUnchanged()
    {
        EXPECT_EQ(m_value, m_data) << m_label;
        m_data = m_saved_value;
    }

private:
    std::string m_label;
    T &m_data;
    T m_saved_value;
    T m_value;
};

#define VALUE_UNCHANGED(var_, value_) ValueUnchanged saved_##var_(#var_, var_, value_)

class TestParameterCommand : public Test
{
protected:
    void SetUp() override;
    void TearDown() override;

    cmdarg_flags exec_cmd_arg(const char *curarg, cmd_file mode);

    cmdarg_flags exec_cmd_arg(const std::string &curarg, cmd_file mode)
    {
        return exec_cmd_arg(curarg.c_str(), mode);
    }
    void expect_stop_msg();

    StrictMock<MockFunction<cmd_arg::StopMsg>> m_stop_msg;
    StrictMock<MockFunction<cmd_arg::Goodbye>> m_goodbye;
    StrictMock<MockFunction<cmd_arg::PrintDoc>> m_print_document;
    cmd_arg::StopMsgFn m_prev_stop_msg;
    cmd_arg::GoodbyeFn m_prev_goodbye;
    cmd_arg::PrintDocFn m_prev_print_document;
    char m_buffer[FILE_MAX_PATH * 2]{};
};

void TestParameterCommand::SetUp()
{
    Test::SetUp();
    m_prev_stop_msg = cmd_arg::get_stop_msg();
    m_prev_goodbye = cmd_arg::get_goodbye();
    m_prev_print_document = cmd_arg::get_print_document();
    cmd_arg::set_stop_msg(m_stop_msg.AsStdFunction());
    cmd_arg::set_goodbye(m_goodbye.AsStdFunction());
    cmd_arg::set_print_document(m_print_document.AsStdFunction());
}

void TestParameterCommand::TearDown()
{
    cmd_arg::set_print_document(m_prev_print_document);
    cmd_arg::set_goodbye(m_prev_goodbye);
    cmd_arg::set_stop_msg(m_prev_stop_msg);
    Test::TearDown();
}

cmdarg_flags TestParameterCommand::exec_cmd_arg(const char *curarg, cmd_file mode)
{
    std::strcpy(m_buffer, curarg);
    return cmdarg(m_buffer, mode);
}

void TestParameterCommand::expect_stop_msg()
{
    EXPECT_CALL(m_stop_msg, Call(stopmsg_flags::NONE, _)).WillOnce(Return(false));
}

TEST_F(TestParameterCommand, parameterTooLong)
{
    expect_stop_msg();

    EXPECT_EQ(
        cmdarg_flags::ERROR, exec_cmd_arg("maximumoftwentycharactersinparametername", cmd_file::SSTOOLS_INI));
}

TEST_F(TestParameterCommand, batchBadArg)
{
    expect_stop_msg();

    EXPECT_EQ(cmdarg_flags::ERROR, exec_cmd_arg("batch=g", cmd_file::SSTOOLS_INI));
}

TEST_F(TestParameterCommand, batchYes)
{
    ValueSaver saved_init_batch{g_init_batch, batch_modes::NONE};

    const cmdarg_flags result = exec_cmd_arg("batch=yes", cmd_file::SSTOOLS_INI);

    EXPECT_EQ(cmdarg_flags::FRACTAL_PARAM | cmdarg_flags::PARAM_3D, result);
    EXPECT_EQ(batch_modes::NORMAL, g_init_batch);
}

TEST_F(TestParameterCommand, batchNo)
{
    ValueSaver saved_init_batch{g_init_batch, batch_modes::NORMAL};

    const cmdarg_flags result = exec_cmd_arg("batch=no", cmd_file::SSTOOLS_INI);

    EXPECT_EQ(cmdarg_flags::FRACTAL_PARAM | cmdarg_flags::PARAM_3D, result);
    EXPECT_EQ(batch_modes::NONE, g_init_batch);
}

TEST_F(TestParameterCommand, batchAfterStartup)
{
    ValueSaver saved_init_batch{g_init_batch, batch_modes::NONE};
    expect_stop_msg();

    EXPECT_EQ(cmdarg_flags::ERROR, exec_cmd_arg("batch=yes", cmd_file::AT_AFTER_STARTUP));
}

TEST_F(TestParameterCommand, maxHistoryNonNumeric)
{
    ValueSaver saved_max_image_history{g_max_image_history, 0};
    expect_stop_msg();

    EXPECT_EQ(cmdarg_flags::ERROR, exec_cmd_arg("maxhistory=yes", cmd_file::SSTOOLS_INI));
}

TEST_F(TestParameterCommand, maxHistoryNegative)
{
    ValueSaver saved_max_image_history{g_max_image_history, 0};
    expect_stop_msg();

    EXPECT_EQ(cmdarg_flags::ERROR, exec_cmd_arg("maxhistory=-10", cmd_file::SSTOOLS_INI));
}

TEST_F(TestParameterCommand, maxHistory)
{
    ValueSaver saved_max_image_history{g_max_image_history, 0};

    const cmdarg_flags result = exec_cmd_arg("maxhistory=10", cmd_file::SSTOOLS_INI);

    EXPECT_EQ(cmdarg_flags::FRACTAL_PARAM | cmdarg_flags::PARAM_3D, result);
    EXPECT_EQ(10, g_max_image_history);
}

TEST_F(TestParameterCommand, maxHistoryAfterStartup)
{
    ValueSaver saved_max_image_history{g_max_image_history, 0};
    expect_stop_msg();

    EXPECT_EQ(cmdarg_flags::ERROR, exec_cmd_arg("maxhistory=10", cmd_file::AT_AFTER_STARTUP));
}

TEST_F(TestParameterCommand, makeDocDefaultFile)
{
    EXPECT_CALL(m_print_document, Call(StrEq("id.txt"), NotNull()));
    EXPECT_CALL(m_goodbye, Call());

    EXPECT_EQ(cmdarg_flags::GOODBYE, exec_cmd_arg("makedoc", cmd_file::SSTOOLS_INI));
}

TEST_F(TestParameterCommand, makeDocCustomFile)
{
    EXPECT_CALL(m_print_document, Call(StrEq("foo.txt"), NotNull()));
    EXPECT_CALL(m_goodbye, Call());

    EXPECT_EQ(cmdarg_flags::GOODBYE, exec_cmd_arg("makedoc=foo.txt", cmd_file::SSTOOLS_INI));
}

TEST_F(TestParameterCommand, makeParTooFewValues)
{
    expect_stop_msg();

    EXPECT_EQ(cmdarg_flags::ERROR, exec_cmd_arg("makepar", cmd_file::SSTOOLS_INI));
}

TEST_F(TestParameterCommand, makeParTooManyValues)
{
    expect_stop_msg();

    EXPECT_EQ(cmdarg_flags::ERROR, exec_cmd_arg("makepar=foo/bar/fmeh", cmd_file::SSTOOLS_INI));
}

// TODO: test makepar with valid arguments

TEST_F(TestParameterCommand, resetBadArg)
{
    ValueSaver saved_escape_exit{g_escape_exit, true};
    expect_stop_msg();

    const cmdarg_flags result = exec_cmd_arg("reset=foo", cmd_file::AT_AFTER_STARTUP);

    EXPECT_EQ(cmdarg_flags::ERROR, result);
    EXPECT_TRUE(g_escape_exit);
}

TEST_F(TestParameterCommand, filenameExtensionTooLong)
{
    ValueSaver saved_gif_filename_mask{g_gif_filename_mask, "*.pot"};
    expect_stop_msg();

    const cmdarg_flags result = exec_cmd_arg("filename=.foobar", cmd_file::AT_AFTER_STARTUP);

    EXPECT_EQ(cmdarg_flags::ERROR, result);
    EXPECT_EQ("*.pot", g_gif_filename_mask);
}

TEST_F(TestParameterCommand, filenameExtension)
{
    ValueSaver saved_gif_filename_mask{g_gif_filename_mask, "*.pot"};

    const cmdarg_flags result = exec_cmd_arg("filename=.gif", cmd_file::AT_AFTER_STARTUP);

    EXPECT_EQ(cmdarg_flags::NONE, result);
    EXPECT_EQ("*.gif", g_gif_filename_mask);
}

TEST_F(TestParameterCommand, filenameValueTooLong)
{
    ValueSaver saved_gif_filename_mask{g_gif_filename_mask, "*.pot"};
    expect_stop_msg();
    const std::string too_long{"filename=" + std::string(FILE_MAX_PATH, 'f') + ".gif"};

    const cmdarg_flags result = exec_cmd_arg(too_long, cmd_file::AT_AFTER_STARTUP);

    EXPECT_EQ(cmdarg_flags::ERROR, result);
    EXPECT_EQ("*.pot", g_gif_filename_mask);
}

TEST_F(TestParameterCommand, mapTooLong)
{
    ValueSaver saved_map_name{g_map_name, "foo.map"};
    expect_stop_msg();
    const std::string too_long{"map=" + std::string(FILE_MAX_PATH, 'f') + ".map"};

    const cmdarg_flags result = exec_cmd_arg(too_long, cmd_file::SSTOOLS_INI);

    EXPECT_EQ(cmdarg_flags::ERROR, result);
    EXPECT_EQ("foo.map", g_map_name);
}

TEST_F(TestParameterCommand, mapSpecifiesSubdir)
{
    ValueSaver saved_map_name{g_map_name, "foo.map"};
    current_path_saver cur_dir(ID_TEST_HOME_DIR);

    const cmdarg_flags result = exec_cmd_arg("map=" ID_TEST_MAP_SUBDIR, cmd_file::SSTOOLS_INI);

    EXPECT_EQ(cmdarg_flags::NONE, result);
    EXPECT_EQ(ID_TEST_MAP_SUBDIR SLASH "foo.map", g_map_name);
}

TEST_F(TestParameterCommand, mapSpecifiesExistingFile)
{
    ValueSaver saved_map_name{g_map_name, "foo.map"};
    current_path_saver cur_dir(ID_TEST_HOME_DIR);

    const cmdarg_flags result =
        exec_cmd_arg("map=" ID_TEST_MAP_SUBDIR SLASH ID_TEST_MAP_FILE, cmd_file::SSTOOLS_INI);

    EXPECT_EQ(cmdarg_flags::NONE, result);
    EXPECT_EQ(ID_TEST_MAP_SUBDIR SLASH ID_TEST_MAP_FILE, g_map_name);
}

TEST_F(TestParameterCommand, adapterDeprecatedValues)
{
    for (const std::string arg : {"egamono", "hgc", "ega", "cga", "mcga", "vga"})
    {
        const cmdarg_flags result = exec_cmd_arg("adapter=" + arg, cmd_file::SSTOOLS_INI);

        EXPECT_EQ(cmdarg_flags::FRACTAL_PARAM | cmdarg_flags::PARAM_3D, result);
    }
}

TEST_F(TestParameterCommand, adapterBadValue)
{
    expect_stop_msg();

    const cmdarg_flags result = exec_cmd_arg("adapter=bad", cmd_file::SSTOOLS_INI);

    EXPECT_EQ(cmdarg_flags::ERROR, result);
}

TEST_F(TestParameterCommand, afiValueIgnored)
{
    const cmdarg_flags result = exec_cmd_arg("afi=anything", cmd_file::SSTOOLS_INI);

    EXPECT_EQ(cmdarg_flags::FRACTAL_PARAM | cmdarg_flags::PARAM_3D, result);
}

TEST_F(TestParameterCommand, textSafeDeprecatedValues)
{
    ValueSaver saved_first_init{g_first_init, true};
    const std::string arg{"textsafe="};

    for (char val : std::string_view{"nybs"})
    {
        const cmdarg_flags result = exec_cmd_arg(arg + val, cmd_file::SSTOOLS_INI);

        EXPECT_EQ(cmdarg_flags::FRACTAL_PARAM | cmdarg_flags::PARAM_3D, result);
    }
}

TEST_F(TestParameterCommand, textSafeInvalidValue)
{
    ValueSaver saved_first_init{g_first_init, true};
    expect_stop_msg();

    const cmdarg_flags result = exec_cmd_arg("textsafe=!", cmd_file::SSTOOLS_INI);

    EXPECT_EQ(cmdarg_flags::ERROR, result);
}

TEST_F(TestParameterCommand, vesaDetectDeprecatedValues)
{
    const std::string arg{"vesadetect="};

    for (char val : std::string_view{"ny"})
    {
        const cmdarg_flags result = exec_cmd_arg(arg + val, cmd_file::SSTOOLS_INI);

        EXPECT_EQ(cmdarg_flags::FRACTAL_PARAM | cmdarg_flags::PARAM_3D, result);
    }
}

TEST_F(TestParameterCommand, vesaDetectInvalidValue)
{
    expect_stop_msg();

    const cmdarg_flags result = exec_cmd_arg("vesadetect=!", cmd_file::SSTOOLS_INI);

    EXPECT_EQ(cmdarg_flags::ERROR, result);
}

TEST_F(TestParameterCommand, biosPaletteDetectDeprecatedValues)
{
    for (char arg : std::string_view{"ny"})
    {
        const cmdarg_flags result = exec_cmd_arg("biospalette=" + std::string{arg}, cmd_file::SSTOOLS_INI);

        EXPECT_EQ(cmdarg_flags::FRACTAL_PARAM | cmdarg_flags::PARAM_3D, result);
    }
}

TEST_F(TestParameterCommand, biosPaletteDetectInvalidValue)
{
    expect_stop_msg();

    const cmdarg_flags result = exec_cmd_arg("biospalette=!", cmd_file::SSTOOLS_INI);

    EXPECT_EQ(cmdarg_flags::ERROR, result);
}

TEST_F(TestParameterCommand, fpuDetectDeprecatedValue)
{
    const cmdarg_flags result = exec_cmd_arg("fpu=387", cmd_file::SSTOOLS_INI);

    EXPECT_EQ(cmdarg_flags::NONE, result);
}

TEST_F(TestParameterCommand, fpuDetectInvalidValue)
{
    expect_stop_msg();

    const cmdarg_flags result = exec_cmd_arg("fpu=487", cmd_file::SSTOOLS_INI);

    EXPECT_EQ(cmdarg_flags::ERROR, result);
}

TEST_F(TestParameterCommand, exitNoAskNo)
{
    ValueSaver saved_escape_exit{g_escape_exit, true};

    const cmdarg_flags result = exec_cmd_arg("exitnoask=n", cmd_file::SSTOOLS_INI);

    EXPECT_FALSE(g_escape_exit);
    EXPECT_EQ(cmdarg_flags::FRACTAL_PARAM | cmdarg_flags::PARAM_3D, result);
}

TEST_F(TestParameterCommand, exitNoAskYes)
{
    ValueSaver saved_escape_exit{g_escape_exit, false};

    const cmdarg_flags result = exec_cmd_arg("exitnoask=y", cmd_file::SSTOOLS_INI);

    EXPECT_TRUE(g_escape_exit);
    EXPECT_EQ(cmdarg_flags::FRACTAL_PARAM | cmdarg_flags::PARAM_3D, result);
}

TEST_F(TestParameterCommand, exitNoAskInvalidValue)
{
    const bool saved_escape_exit{g_escape_exit};
    expect_stop_msg();

    const cmdarg_flags result = exec_cmd_arg("exitnoask=!", cmd_file::SSTOOLS_INI);

    EXPECT_EQ(saved_escape_exit, g_escape_exit);
    EXPECT_EQ(cmdarg_flags::ERROR, result);
}

TEST_F(TestParameterCommand, filenameMask)
{
    ValueSaver saved_gif_filename_mask{g_gif_filename_mask, "*.foo"};

    const cmdarg_flags result = exec_cmd_arg("filename=.pot", cmd_file::AT_CMD_LINE);

    EXPECT_EQ("*.pot", g_gif_filename_mask);
    EXPECT_EQ(cmdarg_flags::NONE, result);
}

// TODO: why does this test cause a crash?
TEST_F(TestParameterCommand, DISABLED_filenameTooLong)
{
    const std::string saved_gif_filename_mask{g_gif_filename_mask};
    const int saved_show_file{g_show_file};
    const std::string saved_browse_name{g_browse_name};
    const std::string too_long{"filename=" + std::string(1024, 'f')};
    expect_stop_msg();

    const cmdarg_flags result = exec_cmd_arg(too_long, cmd_file::AT_CMD_LINE);

    EXPECT_EQ(cmdarg_flags::ERROR, result);
    EXPECT_EQ(saved_gif_filename_mask, g_gif_filename_mask);
    EXPECT_EQ(saved_show_file, g_show_file);
    EXPECT_EQ(saved_browse_name, g_browse_name);
}

TEST_F(TestParameterCommand, videoBadName)
{
    expect_stop_msg();

    const cmdarg_flags result = exec_cmd_arg("video=fmeh", cmd_file::AT_CMD_LINE);

    EXPECT_EQ(cmdarg_flags::ERROR, result);
}

TEST_F(TestParameterCommand, videoNoModes)
{
    expect_stop_msg();
    ValueSaver saved_video_table_len{g_video_table_len, 0};
    ValueSaver saved_init_mode{g_init_mode, 0};

    const cmdarg_flags result = exec_cmd_arg("video=F1", cmd_file::AT_CMD_LINE);

    EXPECT_EQ(cmdarg_flags::ERROR, result);
    EXPECT_EQ(-1, g_init_mode);
}

TEST_F(TestParameterCommand, videoNoMatchingMode)
{
    expect_stop_msg();
    ValueSaver saved_video_table_len{g_video_table_len, 1};
    VIDEOINFO test_mode{};
    test_mode.keynum = ID_KEY_F2;
    ValueSaver saved_video_mode{g_video_table[0], test_mode};
    ValueSaver saved_init_mode{g_init_mode, 0};

    const cmdarg_flags result = exec_cmd_arg("video=F1", cmd_file::AT_CMD_LINE);

    EXPECT_EQ(cmdarg_flags::ERROR, result);
    EXPECT_EQ(-1, g_init_mode);
}

TEST_F(TestParameterCommand, videoMatchingMode)
{
    ValueSaver saved_video_table_len{g_video_table_len, 1};
    VIDEOINFO test_mode{};
    test_mode.keynum = ID_KEY_F1;
    ValueSaver saved_video_mode{g_video_table[0], test_mode};
    ValueSaver saved_init_mode{g_init_mode, 0};

    const cmdarg_flags result = exec_cmd_arg("video=F1", cmd_file::AT_CMD_LINE);

    EXPECT_EQ(cmdarg_flags::FRACTAL_PARAM | cmdarg_flags::PARAM_3D, result);
    EXPECT_EQ(0, g_init_mode);
}

class DACSaver
{
public:
    DACSaver()
    {
        for (int i = 0; i < 256; ++i)
        {
            m_dac_box[i][0] = g_dac_box[i][0];
            m_dac_box[i][1] = g_dac_box[i][1];
            m_dac_box[i][2] = g_dac_box[i][2];
        }
    }
    ~DACSaver()
    {
        for (int i = 0; i < 256; ++i)
        {
            g_dac_box[i][0] = m_dac_box[i][0];
            g_dac_box[i][1] = m_dac_box[i][1];
            g_dac_box[i][2] = m_dac_box[i][2];
        }
    }

private:
    BYTE m_dac_box[256][3];
};

TEST_F(TestParameterCommand, colorsEmptySetsDefaultDAC)
{
    DACSaver saved_dac_box;
    for (int i = 0; i < 256; ++i)
    {
        g_dac_box[i][0] = 0x10;
        g_dac_box[i][1] = 0x20;
        g_dac_box[i][2] = 0x30;
    }

    const cmdarg_flags result = exec_cmd_arg("colors", cmd_file::AT_CMD_LINE);

    EXPECT_EQ(cmdarg_flags::NONE, result);
    EXPECT_EQ(40, g_dac_box[0][0]);
    EXPECT_EQ(40, g_dac_box[0][1]);
    EXPECT_EQ(40, g_dac_box[0][2]);
    EXPECT_EQ(40, g_dac_box[255][0]);
    EXPECT_EQ(40, g_dac_box[255][1]);
    EXPECT_EQ(40, g_dac_box[255][2]);
}

TEST_F(TestParameterCommand, recordColorsInvalidValue)
{
    expect_stop_msg();
    ValueSaver saved_record_colors{g_record_colors, record_colors_mode::none};

    const cmdarg_flags result = exec_cmd_arg("recordcolors=p", cmd_file::AT_CMD_LINE);

    EXPECT_EQ(cmdarg_flags::ERROR, result);
    EXPECT_EQ(record_colors_mode::none, g_record_colors);
}

TEST_F(TestParameterCommand, recordColorsAutomatic)
{
    ValueSaver saved_record_colors{g_record_colors, record_colors_mode::none};

    const cmdarg_flags result = exec_cmd_arg("recordcolors=a", cmd_file::AT_CMD_LINE);

    EXPECT_EQ(cmdarg_flags::NONE, result);
    EXPECT_EQ(record_colors_mode::automatic, g_record_colors);
}

TEST_F(TestParameterCommand, recordColorsComment)
{
    ValueSaver saved_record_colors{g_record_colors, record_colors_mode::none};

    const cmdarg_flags result = exec_cmd_arg("recordcolors=c", cmd_file::AT_CMD_LINE);

    EXPECT_EQ(cmdarg_flags::NONE, result);
    EXPECT_EQ(record_colors_mode::comment, g_record_colors);
}

TEST_F(TestParameterCommand, recordColorsYes)
{
    ValueSaver saved_record_colors{g_record_colors, record_colors_mode::none};

    const cmdarg_flags result = exec_cmd_arg("recordcolors=y", cmd_file::AT_CMD_LINE);

    EXPECT_EQ(cmdarg_flags::NONE, result);
    EXPECT_EQ(record_colors_mode::yes, g_record_colors);
}

TEST_F(TestParameterCommand, maxLineLengthTooSmall)
{
    expect_stop_msg();
    ValueSaver saved_max_line_length{g_max_line_length, 0};
    const std::string arg{"maxlinelength=" + std::to_string(MIN_MAX_LINE_LENGTH - 1)};

    const cmdarg_flags result = exec_cmd_arg(arg, cmd_file::AT_CMD_LINE);

    EXPECT_EQ(cmdarg_flags::ERROR, result);
    EXPECT_EQ(0, g_max_line_length);
}

TEST_F(TestParameterCommand, maxLineLengthTooLarge)
{
    expect_stop_msg();
    ValueSaver saved_max_line_length{g_max_line_length, 0};
    const std::string arg{"maxlinelength=" + std::to_string(MAX_MAX_LINE_LENGTH + 1)};

    const cmdarg_flags result = exec_cmd_arg(arg, cmd_file::AT_CMD_LINE);

    EXPECT_EQ(cmdarg_flags::ERROR, result);
    EXPECT_EQ(0, g_max_line_length);
}

TEST_F(TestParameterCommand, maxLineLength)
{
    ValueSaver saved_max_line_length{g_max_line_length, 0};
    const int value{(MIN_MAX_LINE_LENGTH + MAX_MAX_LINE_LENGTH) / 2};
    const std::string arg{"maxlinelength=" + std::to_string(value)};

    const cmdarg_flags result = exec_cmd_arg(arg, cmd_file::AT_CMD_LINE);

    EXPECT_EQ(cmdarg_flags::NONE, result);
    EXPECT_EQ(value, g_max_line_length);
}

TEST_F(TestParameterCommand, tplusInvalidValue)
{
    expect_stop_msg();

    const cmdarg_flags result = exec_cmd_arg("tplus=!", cmd_file::AT_CMD_LINE);

    EXPECT_EQ(cmdarg_flags::ERROR, result);
}

TEST_F(TestParameterCommand, tplusYes)
{
    const cmdarg_flags result = exec_cmd_arg("tplus=y", cmd_file::AT_CMD_LINE);

    EXPECT_EQ(cmdarg_flags::NONE, result);
}

TEST_F(TestParameterCommand, nonInterlacedInvalidValue)
{
    expect_stop_msg();

    const cmdarg_flags result = exec_cmd_arg("noninterlaced=!", cmd_file::AT_CMD_LINE);

    EXPECT_EQ(cmdarg_flags::ERROR, result);
}

TEST_F(TestParameterCommand, nonInterlacedYes)
{
    const cmdarg_flags result = exec_cmd_arg("noninterlaced=y", cmd_file::AT_CMD_LINE);

    EXPECT_EQ(cmdarg_flags::NONE, result);
}

TEST_F(TestParameterCommand, maxColorResInvalidValue)
{
    expect_stop_msg();

    const cmdarg_flags result = exec_cmd_arg("maxcolorres=3", cmd_file::AT_CMD_LINE);

    EXPECT_EQ(cmdarg_flags::ERROR, result);
}

TEST_F(TestParameterCommand, maxColorResValidValue)
{
    for (int v : {1, 4, 8, 16})
    {
        const std::string arg{"maxcolorres=" + std::to_string(v)};

        const cmdarg_flags result = exec_cmd_arg(arg, cmd_file::AT_CMD_LINE);

        EXPECT_EQ(cmdarg_flags::NONE, result);
    }
}

TEST_F(TestParameterCommand, pixelZoomInvalidValue)
{
    expect_stop_msg();

    const cmdarg_flags result = exec_cmd_arg("pixelzoom=5", cmd_file::AT_CMD_LINE);

    EXPECT_EQ(cmdarg_flags::ERROR, result);
}

TEST_F(TestParameterCommand, pixelZoomValidValue)
{
    const cmdarg_flags result = exec_cmd_arg("pixelzoom=4", cmd_file::AT_CMD_LINE);

    EXPECT_EQ(cmdarg_flags::NONE, result);
}

TEST_F(TestParameterCommand, warnInvalidValue)
{
    expect_stop_msg();
    ValueSaver saved_overwrite_file{g_overwrite_file, false};

    const cmdarg_flags result = exec_cmd_arg("warn=!", cmd_file::AT_CMD_LINE);

    EXPECT_EQ(cmdarg_flags::ERROR, result);
    EXPECT_FALSE(g_overwrite_file);
}

TEST_F(TestParameterCommand, warnNo)
{
    ValueSaver saved_overwrite_file{g_overwrite_file, false};

    const cmdarg_flags result = exec_cmd_arg("warn=n", cmd_file::AT_CMD_LINE);

    EXPECT_EQ(cmdarg_flags::NONE, result);
    EXPECT_TRUE(g_overwrite_file);
}

TEST_F(TestParameterCommand, warnYes)
{
    ValueSaver saved_overwrite_file{g_overwrite_file, true};

    const cmdarg_flags result = exec_cmd_arg("warn=yes", cmd_file::AT_CMD_LINE);

    EXPECT_EQ(cmdarg_flags::NONE, result);
    EXPECT_FALSE(g_overwrite_file);
}

TEST_F(TestParameterCommand, overwriteInvalidValue)
{
    expect_stop_msg();
    ValueSaver saved_overwrite_file{g_overwrite_file, true};

    const cmdarg_flags result = exec_cmd_arg("overwrite=!", cmd_file::AT_CMD_LINE);

    EXPECT_EQ(cmdarg_flags::ERROR, result);
    EXPECT_TRUE(g_overwrite_file);
}

TEST_F(TestParameterCommand, overwriteNo)
{
    ValueSaver saved_overwrite_file{g_overwrite_file, false};

    const cmdarg_flags result = exec_cmd_arg("overwrite=n", cmd_file::AT_CMD_LINE);

    EXPECT_EQ(cmdarg_flags::NONE, result);
    EXPECT_FALSE(g_overwrite_file);
}

TEST_F(TestParameterCommand, overwriteYes)
{
    ValueSaver saved_overwrite_file{g_overwrite_file, true};

    const cmdarg_flags result = exec_cmd_arg("overwrite=yes", cmd_file::AT_CMD_LINE);

    EXPECT_EQ(cmdarg_flags::NONE, result);
    EXPECT_TRUE(g_overwrite_file);
}

TEST_F(TestParameterCommand, gif87aInvalidValue)
{
    expect_stop_msg();

    const cmdarg_flags result = exec_cmd_arg("gif87a=!", cmd_file::AT_CMD_LINE);

    EXPECT_EQ(cmdarg_flags::ERROR, result);
}

TEST_F(TestParameterCommand, gif87aNo)
{
    const cmdarg_flags result = exec_cmd_arg("gif87a=n", cmd_file::AT_CMD_LINE);

    EXPECT_EQ(cmdarg_flags::NONE, result);
}

TEST_F(TestParameterCommand, gif87aYes)
{
    const cmdarg_flags result = exec_cmd_arg("gif87a=yes", cmd_file::AT_CMD_LINE);

    EXPECT_EQ(cmdarg_flags::NONE, result);
}

TEST_F(TestParameterCommand, ditherInvalidValue)
{
    const bool saved_dither_flag{g_dither_flag};
    expect_stop_msg();

    const cmdarg_flags result = exec_cmd_arg("dither=!", cmd_file::AT_CMD_LINE);

    EXPECT_EQ(cmdarg_flags::ERROR, result);
    EXPECT_EQ(saved_dither_flag, g_dither_flag);
}

TEST_F(TestParameterCommand, ditherNo)
{
    ValueSaver saved_dither_flag{g_dither_flag, true};

    const cmdarg_flags result = exec_cmd_arg("dither=n", cmd_file::AT_CMD_LINE);

    EXPECT_EQ(cmdarg_flags::NONE, result);
    EXPECT_FALSE(g_dither_flag);
}

TEST_F(TestParameterCommand, ditherYes)
{
    ValueSaver saved_dither_flag{g_dither_flag, false};

    const cmdarg_flags result = exec_cmd_arg("dither=yes", cmd_file::AT_CMD_LINE);

    EXPECT_EQ(cmdarg_flags::NONE, result);
    EXPECT_TRUE(g_dither_flag);
}

TEST_F(TestParameterCommand, saveTime)
{
    const cmdarg_flags result = exec_cmd_arg("savetime=20", cmd_file::AT_CMD_LINE);

    EXPECT_EQ(cmdarg_flags::NONE, result);
    EXPECT_EQ(20, g_init_save_time);
}

TEST_F(TestParameterCommand, autoKeyInvalidValue)
{
    ValueSaver saved_slides{g_slides, slides_mode::RECORD};
    expect_stop_msg();

    const cmdarg_flags result = exec_cmd_arg("autokey=fmeh", cmd_file::AT_CMD_LINE);

    EXPECT_EQ(cmdarg_flags::ERROR, result);
    EXPECT_EQ(slides_mode::RECORD, g_slides);
}

TEST_F(TestParameterCommand, autoKeyRecord)
{
    ValueSaver saved_slides{g_slides, slides_mode::OFF};

    const cmdarg_flags result = exec_cmd_arg("autokey=record", cmd_file::AT_CMD_LINE);

    EXPECT_EQ(cmdarg_flags::NONE, result);
    EXPECT_EQ(slides_mode::RECORD, g_slides);
}

TEST_F(TestParameterCommand, autoKeyPlay)
{
    ValueSaver saved_slides{g_slides, slides_mode::OFF};

    const cmdarg_flags result = exec_cmd_arg("autokey=play", cmd_file::AT_CMD_LINE);

    EXPECT_EQ(cmdarg_flags::NONE, result);
    EXPECT_EQ(slides_mode::PLAY, g_slides);
}

TEST_F(TestParameterCommand, autoKeyName)
{
    ValueSaver saved_auto_key_name{g_auto_name, "foo.key"};

    const cmdarg_flags result = exec_cmd_arg("autokeyname=baz.key", cmd_file::AT_CMD_LINE);

    EXPECT_EQ(cmdarg_flags::NONE, result);
    EXPECT_EQ("baz.key", g_auto_name);
}

class ParamSaver
{
public:
    ParamSaver(double params[MAX_PARAMS]);
    ~ParamSaver();

private:
    double m_params[MAX_PARAMS];
};

ParamSaver::ParamSaver(double params[MAX_PARAMS])
{
    std::copy(&g_params[0], &g_params[MAX_PARAMS], &m_params[0]);
    std::copy(&params[0], &params[MAX_PARAMS], &g_params[0]);
}

ParamSaver::~ParamSaver()
{
    std::copy(&m_params[0], &m_params[MAX_PARAMS], &g_params[0]);
}

TEST_F(TestParameterCommand, typeSierpinski)
{
    ValueSaver saved_fractal_type{g_fractal_type, fractal_type::LYAPUNOV};
    ValueSaver saved_fractal_specific{g_cur_fractal_specific, nullptr};
    ValueSaver saved_x_min{g_x_min, 111.0};
    ValueSaver saved_x_max{g_x_max, 222.0};
    ValueSaver saved_x_3rd{g_x_3rd, 333.0};
    ValueSaver saved_y_min{g_y_min, 444.0};
    ValueSaver saved_y_max{g_y_max, 555.0};
    ValueSaver saved_y_3rd{g_y_3rd, 666.0};
    double params[MAX_PARAMS]{111.0, 222.0, 333.0, 444.0, 555.0, 666.0, 777.0, 888.0, 999.0, 101010.0};
    ParamSaver saved_params(params);

    const cmdarg_flags result = exec_cmd_arg("type=sierpinski", cmd_file::AT_CMD_LINE);

    EXPECT_EQ(cmdarg_flags::FRACTAL_PARAM, result);
    EXPECT_EQ(fractal_type::SIERPINSKI, g_fractal_type);
    EXPECT_EQ(g_cur_fractal_specific, &g_fractal_specific[+fractal_type::SIERPINSKI]);
    const fractalspecificstuff &fractal{*g_cur_fractal_specific};
    EXPECT_EQ(g_x_min, fractal.xmin);
    EXPECT_EQ(g_x_max, fractal.xmax);
    EXPECT_EQ(g_x_3rd, fractal.xmin);
    EXPECT_EQ(g_y_min, fractal.ymin);
    EXPECT_EQ(g_y_max, fractal.ymax);
    EXPECT_EQ(g_y_3rd, fractal.ymin);
    for (int i = 0; i < 4; ++i)
    {
        EXPECT_EQ(fractal.paramvalue[i], g_params[i]);
    }
}

TEST_F(TestParameterCommand, insideInvalidValue)
{
    expect_stop_msg();
    ValueSaver saved_inside_color{g_inside_color, -9999};

    const cmdarg_flags result = exec_cmd_arg("inside=foo", cmd_file::AT_CMD_LINE);

    EXPECT_EQ(cmdarg_flags::ERROR, result);
    EXPECT_EQ(-9999, g_inside_color);
}

TEST_F(TestParameterCommand, insideZMag)
{
    ValueSaver saved_inside_color{g_inside_color, -9999};

    const cmdarg_flags result = exec_cmd_arg("inside=zmag", cmd_file::AT_CMD_LINE);

    EXPECT_EQ(cmdarg_flags::FRACTAL_PARAM, result);
    EXPECT_EQ(ZMAG, g_inside_color);
}

TEST_F(TestParameterCommand, insideNumber)
{
    ValueSaver saved_inside_color{g_inside_color, -100};

    const cmdarg_flags result = exec_cmd_arg("inside=100", cmd_file::AT_CMD_LINE);

    EXPECT_EQ(cmdarg_flags::FRACTAL_PARAM, result);
    EXPECT_EQ(100, g_inside_color);
}

TEST_F(TestParameterCommand, proximity)
{
    ValueSaver saved_close_proximity{g_close_proximity, -123.0};

    const cmdarg_flags result = exec_cmd_arg("proximity=423", cmd_file::AT_CMD_LINE);

    EXPECT_EQ(cmdarg_flags::FRACTAL_PARAM, result);
    EXPECT_EQ(423.0, g_close_proximity);
}

TEST_F(TestParameterCommand, fillColorNormal)
{
    ValueSaver saved_fill_color{g_fill_color, 999};

    const cmdarg_flags result = exec_cmd_arg("fillcolor=normal", cmd_file::AT_CMD_LINE);

    EXPECT_EQ(cmdarg_flags::FRACTAL_PARAM, result);
    EXPECT_EQ(-1, g_fill_color);
}

TEST_F(TestParameterCommand, fillColorInvalid)
{
    expect_stop_msg();
    ValueSaver saved_fill_color{g_fill_color, 999};

    const cmdarg_flags result = exec_cmd_arg("fillcolor=fmeh", cmd_file::AT_CMD_LINE);

    EXPECT_EQ(cmdarg_flags::ERROR, result);
    EXPECT_EQ(999, g_fill_color);
}

TEST_F(TestParameterCommand, fillColorNumber)
{
    ValueSaver saved_fill_color{g_fill_color, 999};

    const cmdarg_flags result = exec_cmd_arg("fillcolor=100", cmd_file::AT_CMD_LINE);

    EXPECT_EQ(cmdarg_flags::FRACTAL_PARAM, result);
    EXPECT_EQ(100, g_fill_color);
}

TEST_F(TestParameterCommand, finAttractInvalidValue)
{
    expect_stop_msg();
    ValueSaver saved_finite_attractor{g_finite_attractor, true};

    const cmdarg_flags result = exec_cmd_arg("finattract=!", cmd_file::AT_CMD_LINE);

    EXPECT_EQ(cmdarg_flags::ERROR, result);
    EXPECT_TRUE(g_finite_attractor);
}

TEST_F(TestParameterCommand, finAttractYes)
{
    ValueSaver saved_finite_attractor{g_finite_attractor, false};

    const cmdarg_flags result = exec_cmd_arg("finattract=y", cmd_file::AT_CMD_LINE);

    EXPECT_EQ(cmdarg_flags::FRACTAL_PARAM, result);
    EXPECT_TRUE(g_finite_attractor);
}

TEST_F(TestParameterCommand, noBoFInvalidValue)
{
    expect_stop_msg();
    ValueSaver saved_bof_match_book_images{g_bof_match_book_images, true};

    const cmdarg_flags result = exec_cmd_arg("nobof=!", cmd_file::AT_CMD_LINE);

    EXPECT_EQ(cmdarg_flags::ERROR, result);
    EXPECT_TRUE(g_bof_match_book_images);
}

TEST_F(TestParameterCommand, noBoFYes)
{
    ValueSaver saved_bof_match_book_images{g_bof_match_book_images, true};

    const cmdarg_flags result = exec_cmd_arg("nobof=y", cmd_file::AT_CMD_LINE);

    EXPECT_EQ(cmdarg_flags::FRACTAL_PARAM, result);
    EXPECT_FALSE(g_bof_match_book_images);
}

TEST_F(TestParameterCommand, functionSin)
{
    ValueSaver saved_trig_index0{g_trig_index[0], trig_fn::IDENT};
    ValueSaver saved_trig_fns_loaded{g_new_bifurcation_functions_loaded, false};

    const cmdarg_flags result = exec_cmd_arg("function=sin", cmd_file::AT_CMD_LINE);

    EXPECT_EQ(cmdarg_flags::FRACTAL_PARAM, result);
    EXPECT_EQ(trig_fn::SIN, g_trig_index[0]);
    EXPECT_TRUE(g_new_bifurcation_functions_loaded);
}

TEST_F(TestParameterCommand, functionSinCos)
{
    ValueSaver saved_trig_index0{g_trig_index[0], trig_fn::IDENT};
    ValueSaver saved_trig_index1{g_trig_index[1], trig_fn::IDENT};
    ValueSaver saved_trig_fns_loaded{g_new_bifurcation_functions_loaded, false};

    const cmdarg_flags result = exec_cmd_arg("function=sin/cos", cmd_file::AT_CMD_LINE);

    EXPECT_EQ(cmdarg_flags::FRACTAL_PARAM, result);
    EXPECT_EQ(trig_fn::SIN, g_trig_index[0]);
    EXPECT_EQ(trig_fn::COS, g_trig_index[1]);
    EXPECT_TRUE(g_new_bifurcation_functions_loaded);
}

TEST_F(TestParameterCommand, functionSinCosTan)
{
    ValueSaver saved_trig_index0{g_trig_index[0], trig_fn::IDENT};
    ValueSaver saved_trig_index1{g_trig_index[1], trig_fn::IDENT};
    ValueSaver saved_trig_index2{g_trig_index[2], trig_fn::IDENT};
    ValueSaver saved_trig_fns_loaded{g_new_bifurcation_functions_loaded, false};

    const cmdarg_flags result = exec_cmd_arg("function=sin/cos/tan", cmd_file::AT_CMD_LINE);

    EXPECT_EQ(cmdarg_flags::FRACTAL_PARAM, result);
    EXPECT_EQ(trig_fn::SIN, g_trig_index[0]);
    EXPECT_EQ(trig_fn::COS, g_trig_index[1]);
    EXPECT_EQ(trig_fn::TAN, g_trig_index[2]);
    EXPECT_TRUE(g_new_bifurcation_functions_loaded);
}

TEST_F(TestParameterCommand, functionSinCosTanCot)
{
    ValueSaver saved_trig_index0{g_trig_index[0], trig_fn::IDENT};
    ValueSaver saved_trig_index1{g_trig_index[1], trig_fn::IDENT};
    ValueSaver saved_trig_index2{g_trig_index[2], trig_fn::IDENT};
    ValueSaver saved_trig_index3{g_trig_index[3], trig_fn::IDENT};
    ValueSaver saved_trig_fns_loaded{g_new_bifurcation_functions_loaded, false};

    const cmdarg_flags result = exec_cmd_arg("function=sin/cos/tan/cotan", cmd_file::AT_CMD_LINE);

    EXPECT_EQ(cmdarg_flags::FRACTAL_PARAM, result);
    EXPECT_EQ(trig_fn::SIN, g_trig_index[0]);
    EXPECT_EQ(trig_fn::COS, g_trig_index[1]);
    EXPECT_EQ(trig_fn::TAN, g_trig_index[2]);
    EXPECT_EQ(trig_fn::COTAN, g_trig_index[3]);
    EXPECT_TRUE(g_new_bifurcation_functions_loaded);
}

TEST_F(TestParameterCommand, outsideReal)
{
    ValueSaver saved_outside{g_outside_color, -9999};

    const cmdarg_flags result = exec_cmd_arg("outside=real", cmd_file::AT_CMD_LINE);

    EXPECT_EQ(cmdarg_flags::FRACTAL_PARAM, result);
    EXPECT_EQ(REAL, g_outside_color);
}

TEST_F(TestParameterCommand, outsideNumber)
{
    ValueSaver saved_outside{g_outside_color, -9999};

    const cmdarg_flags result = exec_cmd_arg("outside=100", cmd_file::AT_CMD_LINE);

    EXPECT_EQ(cmdarg_flags::FRACTAL_PARAM, result);
    EXPECT_EQ(100, g_outside_color);
}

TEST_F(TestParameterCommand, outsideInvalidName)
{
    expect_stop_msg();
    ValueSaver saved_outside{g_outside_color, -9999};

    const cmdarg_flags result = exec_cmd_arg("outside=zmag", cmd_file::AT_CMD_LINE);

    EXPECT_EQ(cmdarg_flags::ERROR, result);
    EXPECT_EQ(-9999, g_outside_color);
}

TEST_F(TestParameterCommand, bfDigitsNumber)
{
    ValueSaver saved_outside{g_bf_digits, 9999};

    const cmdarg_flags result = exec_cmd_arg("bfdigits=200", cmd_file::AT_CMD_LINE);

    EXPECT_EQ(cmdarg_flags::FRACTAL_PARAM, result);
    EXPECT_EQ(200, g_bf_digits);
}

TEST_F(TestParameterCommand, bfDigitsInvalidValue)
{
    expect_stop_msg();
    ValueSaver saved_outside{g_bf_digits, 9999};

    const cmdarg_flags result = exec_cmd_arg("bfdigits=fmeh", cmd_file::AT_CMD_LINE);

    EXPECT_EQ(cmdarg_flags::ERROR, result);
    EXPECT_EQ(9999, g_bf_digits);
}

TEST_F(TestParameterCommand, bfDigitsTooSmall)
{
    expect_stop_msg();
    ValueSaver saved_outside{g_bf_digits, 9999};

    const cmdarg_flags result = exec_cmd_arg("bfdigits=-1", cmd_file::AT_CMD_LINE);

    EXPECT_EQ(cmdarg_flags::ERROR, result);
    EXPECT_EQ(9999, g_bf_digits);
}

TEST_F(TestParameterCommand, bfDigitsTooLarge)
{
    expect_stop_msg();
    ValueSaver saved_outside{g_bf_digits, 9999};

    const cmdarg_flags result = exec_cmd_arg("bfdigits=2001", cmd_file::AT_CMD_LINE);

    EXPECT_EQ(cmdarg_flags::ERROR, result);
    EXPECT_EQ(9999, g_bf_digits);
}

TEST_F(TestParameterCommand, maxIterNumber)
{
    ValueSaver saved_max_iterations{g_max_iterations, 9999};

    const cmdarg_flags result = exec_cmd_arg("maxiter=20", cmd_file::AT_CMD_LINE);

    EXPECT_EQ(cmdarg_flags::FRACTAL_PARAM, result);
    EXPECT_EQ(20, g_max_iterations);
}

TEST_F(TestParameterCommand, maxIterNumberTooSmall)
{
    expect_stop_msg();
    ValueSaver saved_max_iterations{g_max_iterations, 9999};

    const cmdarg_flags result = exec_cmd_arg("maxiter=1", cmd_file::AT_CMD_LINE);

    EXPECT_EQ(cmdarg_flags::ERROR, result);
    EXPECT_EQ(9999, g_max_iterations);
}

TEST_F(TestParameterCommand, maxIterNotNumber)
{
    expect_stop_msg();
    ValueSaver saved_max_iterations{g_max_iterations, 9999};

    const cmdarg_flags result = exec_cmd_arg("maxiter=fmeh", cmd_file::AT_CMD_LINE);

    EXPECT_EQ(cmdarg_flags::ERROR, result);
    EXPECT_EQ(9999, g_max_iterations);
}

TEST_F(TestParameterCommand, iterIncrIgnored)
{
    const cmdarg_flags result = exec_cmd_arg("iterincr=fmeh", cmd_file::AT_CMD_LINE);

    EXPECT_EQ(cmdarg_flags::NONE, result);
}

TEST_F(TestParameterCommand, passesInvalidValue)
{
    expect_stop_msg();
    ValueSaver saved_user_std_calc_mode{g_user_std_calc_mode, 'Z'};

    const cmdarg_flags result = exec_cmd_arg("passes=!", cmd_file::AT_CMD_LINE);

    EXPECT_EQ(cmdarg_flags::ERROR, result);
    EXPECT_EQ('Z', g_user_std_calc_mode);
}

TEST_F(TestParameterCommand, passesBoundaryTrace)
{
    ValueSaver saved_user_std_calc_mode{g_user_std_calc_mode, 'Z'};

    const cmdarg_flags result = exec_cmd_arg("passes=b", cmd_file::AT_CMD_LINE);

    EXPECT_EQ(cmdarg_flags::FRACTAL_PARAM, result);
    EXPECT_EQ('b', g_user_std_calc_mode);
}

TEST_F(TestParameterCommand, passesSolidGuess3)
{
    ValueSaver saved_user_std_calc_mode{g_user_std_calc_mode, 'Z'};
    ValueSaver saved_stop_pass{g_stop_pass, -1};

    const cmdarg_flags result = exec_cmd_arg("passes=g3", cmd_file::AT_CMD_LINE);

    EXPECT_EQ(cmdarg_flags::FRACTAL_PARAM, result);
    EXPECT_EQ('g', g_user_std_calc_mode);
    EXPECT_EQ(3, g_stop_pass);
}

TEST_F(TestParameterCommand, isMandYes)
{
    ValueSaver saved_is_mandelbrot{g_is_mandelbrot, false};

    const cmdarg_flags result = exec_cmd_arg("ismand=y", cmd_file::AT_CMD_LINE);

    EXPECT_EQ(cmdarg_flags::FRACTAL_PARAM, result);
    EXPECT_TRUE(g_is_mandelbrot);
}

TEST_F(TestParameterCommand, cycleLimitTooLow)
{
    expect_stop_msg();
    ValueSaver saved_cycle_limit{g_init_cycle_limit, 9999};

    const cmdarg_flags result = exec_cmd_arg("cyclelimit=1", cmd_file::AT_CMD_LINE);

    EXPECT_EQ(cmdarg_flags::ERROR, result);
    EXPECT_EQ(9999, g_init_cycle_limit);
}

TEST_F(TestParameterCommand, cycleLimitTooHigh)
{
    expect_stop_msg();
    ValueSaver saved_cycle_limit{g_init_cycle_limit, 9999};

    const cmdarg_flags result = exec_cmd_arg("cyclelimit=257", cmd_file::AT_CMD_LINE);

    EXPECT_EQ(cmdarg_flags::ERROR, result);
    EXPECT_EQ(9999, g_init_cycle_limit);
}

TEST_F(TestParameterCommand, cycleLimitNumber)
{
    ValueSaver saved_cycle_limit{g_init_cycle_limit, 9999};

    const cmdarg_flags result = exec_cmd_arg("cyclelimit=100", cmd_file::AT_CMD_LINE);

    EXPECT_EQ(cmdarg_flags::NONE, result);
    EXPECT_EQ(100, g_init_cycle_limit);
}

TEST_F(TestParameterCommand, cycleRangeNoParams)
{
    ValueSaver saved_cycle_range_lo{g_color_cycle_range_lo, 9999};
    ValueSaver saved_cycle_range_hi{g_color_cycle_range_hi, 9999};

    const cmdarg_flags result = exec_cmd_arg("cyclerange", cmd_file::AT_CMD_LINE);

    EXPECT_EQ(cmdarg_flags::NONE, result);
    EXPECT_EQ(1, g_color_cycle_range_lo);
    EXPECT_EQ(255, g_color_cycle_range_hi);
}

TEST_F(TestParameterCommand, cycleRangeOneParam)
{
    ValueSaver saved_cycle_range_lo{g_color_cycle_range_lo, 9999};
    ValueSaver saved_cycle_range_hi{g_color_cycle_range_hi, 9999};

    const cmdarg_flags result = exec_cmd_arg("cyclerange=10", cmd_file::AT_CMD_LINE);

    EXPECT_EQ(cmdarg_flags::NONE, result);
    EXPECT_EQ(10, g_color_cycle_range_lo);
    EXPECT_EQ(255, g_color_cycle_range_hi);
}

TEST_F(TestParameterCommand, cycleRangeTwoParams)
{
    ValueSaver saved_cycle_range_lo{g_color_cycle_range_lo, 9999};
    ValueSaver saved_cycle_range_hi{g_color_cycle_range_hi, 9999};

    const cmdarg_flags result = exec_cmd_arg("cyclerange=10/20", cmd_file::AT_CMD_LINE);

    EXPECT_EQ(cmdarg_flags::NONE, result);
    EXPECT_EQ(10, g_color_cycle_range_lo);
    EXPECT_EQ(20, g_color_cycle_range_hi);
}

TEST_F(TestParameterCommand, cycleRangeLoTooLow)
{
    expect_stop_msg();
    ValueSaver saved_cycle_range_lo{g_color_cycle_range_lo, 9999};
    ValueSaver saved_cycle_range_hi{g_color_cycle_range_hi, 9999};

    const cmdarg_flags result = exec_cmd_arg("cyclerange=-1/20", cmd_file::AT_CMD_LINE);

    EXPECT_EQ(cmdarg_flags::ERROR, result);
    EXPECT_EQ(9999, g_color_cycle_range_lo);
    EXPECT_EQ(9999, g_color_cycle_range_hi);
}

TEST_F(TestParameterCommand, cycleRangeHiTooHigh)
{
    expect_stop_msg();
    ValueSaver saved_cycle_range_lo{g_color_cycle_range_lo, 9999};
    ValueSaver saved_cycle_range_hi{g_color_cycle_range_hi, 9999};

    const cmdarg_flags result = exec_cmd_arg("cyclerange=10/256", cmd_file::AT_CMD_LINE);

    EXPECT_EQ(cmdarg_flags::ERROR, result);
    EXPECT_EQ(9999, g_color_cycle_range_lo);
    EXPECT_EQ(9999, g_color_cycle_range_hi);
}

TEST_F(TestParameterCommand, rangesInvalidParamCount)
{
    expect_stop_msg();
    ValueSaver saved_log_map_flag(g_log_map_flag, 9999);

    const cmdarg_flags result = exec_cmd_arg("ranges=100/102/foo", cmd_file::AT_CMD_LINE);

    EXPECT_EQ(cmdarg_flags::ERROR, result);
    EXPECT_EQ(9999, g_log_map_flag);
}

TEST_F(TestParameterCommand, rangesNoStriping)
{
    ValueSaver saved_log_map_flag(g_log_map_flag, 9999);

    const cmdarg_flags result = exec_cmd_arg("ranges=10/20/30", cmd_file::AT_CMD_LINE);

    EXPECT_EQ(cmdarg_flags::FRACTAL_PARAM, result);
    EXPECT_EQ(0, g_log_map_flag);
    EXPECT_EQ((std::vector{10, 20, 30}), g_iteration_ranges);
}

TEST_F(TestParameterCommand, rangesOneStripe)
{
    ValueSaver saved_log_map_flag(g_log_map_flag, 9999);

    const cmdarg_flags result = exec_cmd_arg("ranges=10/-20/30", cmd_file::AT_CMD_LINE);

    EXPECT_EQ(cmdarg_flags::FRACTAL_PARAM, result);
    EXPECT_EQ(0, g_log_map_flag);
    EXPECT_EQ((std::vector{10, -1, 20, 30}), g_iteration_ranges);
}

TEST_F(TestParameterCommand, saveNameOnFirstInit)
{
    ValueSaver saved_first_init{g_first_init, true};

    const cmdarg_flags result = exec_cmd_arg("savename=test.gif", cmd_file::AT_CMD_LINE);

    EXPECT_EQ(cmdarg_flags::NONE, result);
    EXPECT_EQ("test.gif", g_save_filename);
}

TEST_F(TestParameterCommand, saveNameAfterStartup)
{
    ValueSaver saved_first_init{g_first_init, false};
    ValueSaver saved_save_filename{g_save_filename, "bar.gif"};

    const cmdarg_flags result = exec_cmd_arg("savename=test.gif", cmd_file::AT_AFTER_STARTUP);

    EXPECT_EQ(cmdarg_flags::NONE, result);
    EXPECT_EQ("test.gif", g_save_filename);
}

TEST_F(TestParameterCommand, saveNameOtherwise)
{
    ValueSaver saved_first_init{g_first_init, false};
    ValueSaver saved_save_filename{g_save_filename, "bar.gif"};

    const cmdarg_flags result = exec_cmd_arg("savename=test.gif", cmd_file::AT_CMD_LINE_SET_NAME);

    EXPECT_EQ(cmdarg_flags::NONE, result);
    EXPECT_EQ("bar.gif", g_save_filename);
}

TEST_F(TestParameterCommand, tweakLZWDeprecated)
{
    const cmdarg_flags result = exec_cmd_arg("tweaklzw=fmeh", cmd_file::AT_AFTER_STARTUP);

    EXPECT_EQ(cmdarg_flags::NONE, result);
}

TEST_F(TestParameterCommand, minStackTooFewValues)
{
    expect_stop_msg();
    ValueSaver saved_soi_min_stack{g_soi_min_stack, 9999};

    const cmdarg_flags result = exec_cmd_arg("minstack=", cmd_file::AT_AFTER_STARTUP);

    EXPECT_EQ(cmdarg_flags::ERROR, result);
    EXPECT_EQ(9999, g_soi_min_stack);
}

TEST_F(TestParameterCommand, minStackTooManyValues)
{
    expect_stop_msg();
    ValueSaver saved_soi_min_stack{g_soi_min_stack, 9999};

    const cmdarg_flags result = exec_cmd_arg("minstack=10/20", cmd_file::AT_AFTER_STARTUP);

    EXPECT_EQ(cmdarg_flags::ERROR, result);
    EXPECT_EQ(9999, g_soi_min_stack);
}

TEST_F(TestParameterCommand, minStackNumber)
{
    ValueSaver saved_soi_min_stack{g_soi_min_stack, 9999};

    const cmdarg_flags result = exec_cmd_arg("minstack=200", cmd_file::AT_AFTER_STARTUP);

    EXPECT_EQ(cmdarg_flags::NONE, result);
    EXPECT_EQ(200, g_soi_min_stack);
}

TEST_F(TestParameterCommand, mathToleranceOneValue)
{
    ValueSaver saved_math_tol0{g_math_tol[0], 9999.0};
    ValueSaver saved_math_tol1{g_math_tol[1], 8888.0};

    const cmdarg_flags result = exec_cmd_arg("mathtolerance=20", cmd_file::AT_AFTER_STARTUP);

    EXPECT_EQ(cmdarg_flags::NONE, result);
    EXPECT_EQ(20.0, g_math_tol[0]);
    EXPECT_EQ(8888.0, g_math_tol[1]);
}

TEST_F(TestParameterCommand, mathToleranceTwoValues)
{
    ValueSaver saved_math_tol0{g_math_tol[0], 9999.0};
    ValueSaver saved_math_tol1{g_math_tol[1], 8888.0};

    const cmdarg_flags result = exec_cmd_arg("mathtolerance=20/30", cmd_file::AT_AFTER_STARTUP);

    EXPECT_EQ(cmdarg_flags::NONE, result);
    EXPECT_EQ(20.0, g_math_tol[0]);
    EXPECT_EQ(30, g_math_tol[1]);
}

TEST_F(TestParameterCommand, mathToleranceSecondValueOnly)
{
    ValueSaver saved_math_tol0{g_math_tol[0], 9999.0};
    ValueSaver saved_math_tol1{g_math_tol[1], 8888.0};

    const cmdarg_flags result = exec_cmd_arg("mathtolerance=/30", cmd_file::AT_AFTER_STARTUP);

    EXPECT_EQ(cmdarg_flags::NONE, result);
    EXPECT_EQ(9999, g_math_tol[0]);
    EXPECT_EQ(30, g_math_tol[1]);
}

static std::string adjust_dir(std::string actual)
{
    if (actual[1] == ':')
    {
        actual[0] = static_cast<char>(std::toupper(static_cast<unsigned char>(actual[0])));
    }
    std::filesystem::path sep;
    sep += std::filesystem::path::preferred_separator;
    if (actual.back() == sep.string().back())
    {
        actual.pop_back();
    }
    return actual;
}

TEST_F(TestParameterCommand, tempDirExisting)
{
    std::filesystem::path start_dir{ID_TEST_DATA_DIR};
    start_dir.make_preferred();
    std::filesystem::path new_dir{ID_TEST_DATA_SUBDIR};
    new_dir.make_preferred();
    ValueSaver saved_temp_dir{g_temp_dir, start_dir.make_preferred().string()};

    const cmdarg_flags result = exec_cmd_arg("tempdir=" + new_dir.string(), cmd_file::AT_AFTER_STARTUP);

    EXPECT_EQ(cmdarg_flags::NONE, result);
    EXPECT_EQ(new_dir.string(), adjust_dir(g_temp_dir));
}

TEST_F(TestParameterCommand, workDirExisting)
{
    std::filesystem::path start_dir{ID_TEST_DATA_DIR};
    start_dir.make_preferred();
    std::filesystem::path new_dir{ID_TEST_DATA_SUBDIR};
    new_dir.make_preferred();
    ValueSaver saved_working_dir{g_working_dir, start_dir.make_preferred().string()};

    const cmdarg_flags result = exec_cmd_arg("workdir=" + new_dir.string(), cmd_file::AT_AFTER_STARTUP);

    EXPECT_EQ(cmdarg_flags::NONE, result);
    EXPECT_EQ(new_dir.string(), adjust_dir(g_working_dir));
}

TEST_F(TestParameterCommand, exitModeDeprecated)
{
    const cmdarg_flags result = exec_cmd_arg("exitmode=fmeh", cmd_file::AT_AFTER_STARTUP);

    EXPECT_EQ(cmdarg_flags::NONE, result);
}

class TextColorSaver
{
public:
    TextColorSaver();
    ~TextColorSaver();

private:
    BYTE m_text_color[31]{};
};

TextColorSaver::TextColorSaver()
{
    std::copy(&g_text_color[0], &g_text_color[31], &m_text_color[0]);
}

TextColorSaver::~TextColorSaver()
{
    std::copy(&m_text_color[0], &m_text_color[31], &g_text_color[0]);
}

TEST_F(TestParameterCommand, textColorsMono)
{
    TextColorSaver saved_text_colors;

    const cmdarg_flags result = exec_cmd_arg("textcolors=mono", cmd_file::AT_AFTER_STARTUP);

    EXPECT_EQ(cmdarg_flags::NONE, result);
    for (int i : {1, 3, 4, 7, 8, 9, 10, 15, 18, 19, 21, 23, 26, 29, 30})
    {
        EXPECT_EQ(BLACK * 16 + WHITE, g_text_color[i]) << i;
    }
    for (int i : {6, 12, 13, 14, 20, 27, 28})
    {
        EXPECT_EQ(WHITE * 16 + BLACK, g_text_color[i]) << i;
    }
    for (int i : {0, 2, 5, 11, 16, 17, 22, 24, 25})
    {
        EXPECT_EQ(BLACK * 16 + L_WHITE, g_text_color[i]) << i;
    }
}

TEST_F(TestParameterCommand, textColorsValues)
{
    TextColorSaver saved_text_colors;

    const cmdarg_flags result = exec_cmd_arg("textcolors=24/28", cmd_file::AT_AFTER_STARTUP);

    EXPECT_EQ(cmdarg_flags::NONE, result);
    EXPECT_EQ(0x24, g_text_color[0]);
    EXPECT_EQ(0x28, g_text_color[1]);
}

TEST_F(TestParameterCommand, textColorsSkippedValues)
{
    TextColorSaver saved_text_colors;
    g_text_color[0] = 0x24;

    const cmdarg_flags result = exec_cmd_arg("textcolors=/28", cmd_file::AT_AFTER_STARTUP);

    EXPECT_EQ(cmdarg_flags::NONE, result);
    EXPECT_EQ(0x24, g_text_color[0]);
    EXPECT_EQ(0x28, g_text_color[1]);
}

TEST_F(TestParameterCommand, potentialOneValue)
{
    ValueSaver saved_potential_params0{g_potential_params[0], 9999.0};
    ValueSaver saved_potential_params1{g_potential_params[1], 8888.0};
    ValueSaver saved_potential_params2{g_potential_params[2], 7777.0};
    ValueSaver saved_potential_16bit{g_potential_16bit, true};

    const cmdarg_flags result = exec_cmd_arg("potential=111", cmd_file::AT_AFTER_STARTUP);

    EXPECT_EQ(cmdarg_flags::FRACTAL_PARAM, result);
    EXPECT_EQ(111.0, g_potential_params[0]);
    EXPECT_EQ(8888.0, g_potential_params[1]);
    EXPECT_EQ(7777.0, g_potential_params[2]);
    EXPECT_FALSE(g_potential_16bit);
}

TEST_F(TestParameterCommand, potentialTwoValues)
{
    ValueSaver saved_potential_params0{g_potential_params[0], 9999.0};
    ValueSaver saved_potential_params1{g_potential_params[1], 8888.0};
    ValueSaver saved_potential_params2{g_potential_params[2], 7777.0};
    ValueSaver saved_potential_16bit{g_potential_16bit, true};

    const cmdarg_flags result = exec_cmd_arg("potential=111/222", cmd_file::AT_AFTER_STARTUP);

    EXPECT_EQ(cmdarg_flags::FRACTAL_PARAM, result);
    EXPECT_EQ(111.0, g_potential_params[0]);
    EXPECT_EQ(222.0, g_potential_params[1]);
    EXPECT_EQ(7777.0, g_potential_params[2]);
    EXPECT_FALSE(g_potential_16bit);
}

TEST_F(TestParameterCommand, potentialThreeValues)
{
    ValueSaver saved_potential_params0{g_potential_params[0], 9999.0};
    ValueSaver saved_potential_params1{g_potential_params[1], 8888.0};
    ValueSaver saved_potential_params2{g_potential_params[2], 7777.0};
    ValueSaver saved_potential_16bit{g_potential_16bit, true};

    const cmdarg_flags result = exec_cmd_arg("potential=111/222/333", cmd_file::AT_AFTER_STARTUP);

    EXPECT_EQ(cmdarg_flags::FRACTAL_PARAM, result);
    EXPECT_EQ(111.0, g_potential_params[0]);
    EXPECT_EQ(222.0, g_potential_params[1]);
    EXPECT_EQ(333.0, g_potential_params[2]);
    EXPECT_FALSE(g_potential_16bit);
}

TEST_F(TestParameterCommand, potential16Bit)
{
    ValueSaver saved_potential_params0{g_potential_params[0], 9999.0};
    ValueSaver saved_potential_params1{g_potential_params[1], 8888.0};
    ValueSaver saved_potential_params2{g_potential_params[2], 7777.0};
    ValueSaver saved_potential_16bit{g_potential_16bit, false};

    const cmdarg_flags result = exec_cmd_arg("potential=111/222/333/16bit", cmd_file::AT_AFTER_STARTUP);

    EXPECT_EQ(cmdarg_flags::FRACTAL_PARAM, result);
    EXPECT_EQ(111.0, g_potential_params[0]);
    EXPECT_EQ(222.0, g_potential_params[1]);
    EXPECT_EQ(333.0, g_potential_params[2]);
    EXPECT_TRUE(g_potential_16bit);
}

TEST_F(TestParameterCommand, paramsOneValue)
{
    ValueSaver saved_params0{g_params[0], 1111.0};
    ValueSaver saved_params1{g_params[1], 2222.0};
    ValueSaver saved_bf_math{g_bf_math, bf_math_type::NONE};

    const cmdarg_flags result = exec_cmd_arg("params=1.0", cmd_file::AT_AFTER_STARTUP);

    EXPECT_EQ(cmdarg_flags::FRACTAL_PARAM, result);
    EXPECT_EQ(1.0, g_params[0]);
    EXPECT_EQ(0.0, g_params[1]);
}

TEST_F(TestParameterCommand, paramsTwoValues)
{
    ValueSaver saved_params0{g_params[0], 1111.0};
    ValueSaver saved_params1{g_params[1], 2222.0};
    ValueSaver saved_params2{g_params[2], 3333.0};
    ValueSaver saved_bf_math{g_bf_math, bf_math_type::NONE};

    const cmdarg_flags result = exec_cmd_arg("params=1.0/2.0", cmd_file::AT_AFTER_STARTUP);

    EXPECT_EQ(cmdarg_flags::FRACTAL_PARAM, result);
    EXPECT_EQ(1.0, g_params[0]);
    EXPECT_EQ(2.0, g_params[1]);
    EXPECT_EQ(0.0, g_params[2]);
}

TEST_F(TestParameterCommand, miimBreadthFirstLeftToRight)
{
    ValueSaver saved_major_method{g_major_method, Major::random_run};
    ValueSaver saved_minor_method{g_inverse_julia_minor_method, Minor::right_first};

    const cmdarg_flags result = exec_cmd_arg("miim=b/l", cmd_file::AT_AFTER_STARTUP);

    EXPECT_EQ(cmdarg_flags::FRACTAL_PARAM, result);
    EXPECT_EQ(Major::breadth_first, g_major_method);
    EXPECT_EQ(Minor::left_first, g_inverse_julia_minor_method);
}

TEST_F(TestParameterCommand, initOrbitPixel)
{
    ValueSaver saved_use_init_orbit{g_use_init_orbit, init_orbit_mode::value};
    ValueSaver saved_init_orbit{g_init_orbit, DComplex{111.0, 222.0}};

    const cmdarg_flags result = exec_cmd_arg("initorbit=pixel", cmd_file::AT_AFTER_STARTUP);

    EXPECT_EQ(cmdarg_flags::FRACTAL_PARAM, result);
    EXPECT_EQ(init_orbit_mode::pixel, g_use_init_orbit);
    EXPECT_EQ(111.0, g_init_orbit.x);
    EXPECT_EQ(222.0, g_init_orbit.y);
}

TEST_F(TestParameterCommand, initOrbitValue)
{
    ValueSaver saved_use_init_orbit{g_use_init_orbit, init_orbit_mode::pixel};
    ValueSaver saved_init_orbit{g_init_orbit, DComplex{111.0, 222.0}};

    const cmdarg_flags result = exec_cmd_arg("initorbit=10/20", cmd_file::AT_AFTER_STARTUP);

    EXPECT_EQ(cmdarg_flags::FRACTAL_PARAM, result);
    EXPECT_EQ(init_orbit_mode::value, g_use_init_orbit);
    EXPECT_EQ(10.0, g_init_orbit.x);
    EXPECT_EQ(20.0, g_init_orbit.y);
}

TEST_F(TestParameterCommand, initOrbitTooFewParameters)
{
    expect_stop_msg();
    ValueSaver saved_use_init_orbit{g_use_init_orbit, init_orbit_mode::pixel};
    ValueSaver saved_init_orbit{g_init_orbit, DComplex{111.0, 222.0}};

    const cmdarg_flags result = exec_cmd_arg("initorbit=10", cmd_file::AT_AFTER_STARTUP);

    EXPECT_EQ(cmdarg_flags::ERROR, result);
    EXPECT_EQ(init_orbit_mode::pixel, g_use_init_orbit);
    EXPECT_EQ(111.0, g_init_orbit.x);
    EXPECT_EQ(222.0, g_init_orbit.y);
}

TEST_F(TestParameterCommand, initOrbitTooFewFloatParameters)
{
    expect_stop_msg();
    ValueSaver saved_use_init_orbit{g_use_init_orbit, init_orbit_mode::pixel};
    ValueSaver saved_init_orbit{g_init_orbit, DComplex{111.0, 222.0}};

    const cmdarg_flags result = exec_cmd_arg("initorbit=10/fmeh", cmd_file::AT_AFTER_STARTUP);

    EXPECT_EQ(cmdarg_flags::ERROR, result);
    EXPECT_EQ(init_orbit_mode::pixel, g_use_init_orbit);
    EXPECT_EQ(111.0, g_init_orbit.x);
    EXPECT_EQ(222.0, g_init_orbit.y);
}

TEST_F(TestParameterCommand, initOrbitTooManyParameters)
{
    expect_stop_msg();
    ValueSaver saved_use_init_orbit{g_use_init_orbit, init_orbit_mode::pixel};
    ValueSaver saved_init_orbit{g_init_orbit, DComplex{111.0, 222.0}};

    const cmdarg_flags result = exec_cmd_arg("initorbit=10/20/30", cmd_file::AT_AFTER_STARTUP);

    EXPECT_EQ(cmdarg_flags::ERROR, result);
    EXPECT_EQ(init_orbit_mode::pixel, g_use_init_orbit);
    EXPECT_EQ(111.0, g_init_orbit.x);
    EXPECT_EQ(222.0, g_init_orbit.y);
}

TEST_F(TestParameterCommand, threeDModeMonocular)
{
    ValueSaver saved_julibrot_3d_mode{g_julibrot_3d_mode, julibrot_3d_mode::LEFT_EYE};

    const cmdarg_flags result = exec_cmd_arg("3dmode=monocular", cmd_file::AT_AFTER_STARTUP);

    EXPECT_EQ(cmdarg_flags::FRACTAL_PARAM, result);
    EXPECT_EQ(julibrot_3d_mode::MONOCULAR, g_julibrot_3d_mode);
}

TEST_F(TestParameterCommand, julibrot3DZDots)
{
    ValueSaver saved_julibrot_z_dots{g_julibrot_z_dots, 9999};
    ValueSaver saved_julibrot_origin_fp{g_julibrot_origin_fp, 111.0f};

    const cmdarg_flags result = exec_cmd_arg("julibrot3d=100", cmd_file::AT_AFTER_STARTUP);

    EXPECT_EQ(cmdarg_flags::FRACTAL_PARAM, result);
    EXPECT_EQ(100, g_julibrot_z_dots);
    EXPECT_EQ(111.0f, g_julibrot_origin_fp);
}

TEST_F(TestParameterCommand, julibrot3DOrigin)
{
    ValueSaver saved_julibrot_z_dots{g_julibrot_z_dots, 9999};
    ValueSaver saved_julibrot_origin{g_julibrot_origin_fp, 111.0f};
    ValueSaver saved_julibrot_depth{g_julibrot_depth_fp, 222.0f};

    const cmdarg_flags result = exec_cmd_arg("julibrot3d=100/10", cmd_file::AT_AFTER_STARTUP);

    EXPECT_EQ(cmdarg_flags::FRACTAL_PARAM, result);
    EXPECT_EQ(100, g_julibrot_z_dots);
    EXPECT_EQ(10.0f, g_julibrot_origin_fp);
    EXPECT_EQ(222.0f, g_julibrot_depth_fp);
}

TEST_F(TestParameterCommand, julibrot3DDepth)
{
    ValueSaver saved_julibrot_z_dots{g_julibrot_z_dots, 9999};
    ValueSaver saved_julibrot_origin{g_julibrot_origin_fp, 111.0f};
    ValueSaver saved_julibrot_depth{g_julibrot_depth_fp, 222.0f};
    ValueSaver saved_julibrot_height{g_julibrot_height_fp, 333.0f};

    const cmdarg_flags result = exec_cmd_arg("julibrot3d=100/10/12", cmd_file::AT_AFTER_STARTUP);

    EXPECT_EQ(cmdarg_flags::FRACTAL_PARAM, result);
    EXPECT_EQ(100, g_julibrot_z_dots);
    EXPECT_EQ(10.0f, g_julibrot_origin_fp);
    EXPECT_EQ(12.0f, g_julibrot_depth_fp);
    EXPECT_EQ(333.0f, g_julibrot_height_fp);
}

TEST_F(TestParameterCommand, julibrot3DHeight)
{
    ValueSaver saved_julibrot_z_dots{g_julibrot_z_dots, 9999};
    ValueSaver saved_julibrot_origin{g_julibrot_origin_fp, 111.0f};
    ValueSaver saved_julibrot_depth{g_julibrot_depth_fp, 222.0f};
    ValueSaver saved_julibrot_height{g_julibrot_height_fp, 333.0f};
    ValueSaver saved_julibrot_width{g_julibrot_width_fp, 444.0f};

    const cmdarg_flags result = exec_cmd_arg("julibrot3d=100/10/12/14", cmd_file::AT_AFTER_STARTUP);

    EXPECT_EQ(cmdarg_flags::FRACTAL_PARAM, result);
    EXPECT_EQ(100, g_julibrot_z_dots);
    EXPECT_EQ(10.0f, g_julibrot_origin_fp);
    EXPECT_EQ(12.0f, g_julibrot_depth_fp);
    EXPECT_EQ(14.0f, g_julibrot_height_fp);
    EXPECT_EQ(444.0f, g_julibrot_width_fp);
}

TEST_F(TestParameterCommand, julibrot3DWidth)
{
    ValueSaver saved_julibrot_z_dots{g_julibrot_z_dots, 9999};
    ValueSaver saved_julibrot_origin{g_julibrot_origin_fp, 111.0f};
    ValueSaver saved_julibrot_depth{g_julibrot_depth_fp, 222.0f};
    ValueSaver saved_julibrot_height{g_julibrot_height_fp, 333.0f};
    ValueSaver saved_julibrot_width{g_julibrot_width_fp, 444.0f};
    ValueSaver saved_julibrot_dist{g_julibrot_dist_fp, 555.0f};

    const cmdarg_flags result = exec_cmd_arg("julibrot3d=100/10/12/14/16", cmd_file::AT_AFTER_STARTUP);

    EXPECT_EQ(cmdarg_flags::FRACTAL_PARAM, result);
    EXPECT_EQ(100, g_julibrot_z_dots);
    EXPECT_EQ(10.0f, g_julibrot_origin_fp);
    EXPECT_EQ(12.0f, g_julibrot_depth_fp);
    EXPECT_EQ(14.0f, g_julibrot_height_fp);
    EXPECT_EQ(16.0f, g_julibrot_width_fp);
    EXPECT_EQ(555.0f, g_julibrot_dist_fp);
}

TEST_F(TestParameterCommand, julibrot3DDistance)
{
    ValueSaver saved_julibrot_z_dots{g_julibrot_z_dots, 9999};
    ValueSaver saved_julibrot_origin{g_julibrot_origin_fp, 111.0f};
    ValueSaver saved_julibrot_depth{g_julibrot_depth_fp, 222.0f};
    ValueSaver saved_julibrot_height{g_julibrot_height_fp, 333.0f};
    ValueSaver saved_julibrot_width{g_julibrot_width_fp, 444.0f};
    ValueSaver saved_julibrot_dist{g_julibrot_dist_fp, 555.0f};

    const cmdarg_flags result = exec_cmd_arg("julibrot3d=100/10/12/14/16/18", cmd_file::AT_AFTER_STARTUP);

    EXPECT_EQ(cmdarg_flags::FRACTAL_PARAM, result);
    EXPECT_EQ(100, g_julibrot_z_dots);
    EXPECT_EQ(10.0f, g_julibrot_origin_fp);
    EXPECT_EQ(12.0f, g_julibrot_depth_fp);
    EXPECT_EQ(14.0f, g_julibrot_height_fp);
    EXPECT_EQ(16.0f, g_julibrot_width_fp);
    EXPECT_EQ(18.0f, g_julibrot_dist_fp);
}

TEST_F(TestParameterCommand, julibrotEyes)
{
    ValueSaver saved_eyes{g_eyes_fp, 111.0f};

    const cmdarg_flags result = exec_cmd_arg("julibroteyes=10", cmd_file::AT_AFTER_STARTUP);

    EXPECT_EQ(cmdarg_flags::FRACTAL_PARAM, result);
    EXPECT_EQ(10.0f, g_eyes_fp);
}

TEST_F(TestParameterCommand, julibrotFromToBadNumberOfValues)
{
    ValueSaver saved_julibrot_x_max{g_julibrot_x_max, 222.0f};
    ValueSaver saved_julibrot_x_min{g_julibrot_x_min, 111.0f};
    ValueSaver saved_julibrot_y_max{g_julibrot_y_max, 444.0f};
    ValueSaver saved_julibrot_y_min{g_julibrot_y_min, 333.0f};

    const cmdarg_flags result = exec_cmd_arg("julibrotfromto=40/30/20/10", cmd_file::AT_AFTER_STARTUP);

    EXPECT_EQ(cmdarg_flags::FRACTAL_PARAM, result);
    EXPECT_EQ(40.0, g_julibrot_x_max);
    EXPECT_EQ(30.0, g_julibrot_x_min);
    EXPECT_EQ(20.0, g_julibrot_y_max);
    EXPECT_EQ(10.0, g_julibrot_y_min);
}

TEST_F(TestParameterCommand, cornersNoValues)
{
    ValueSaver saved_use_center_mag{g_use_center_mag, true};

    const cmdarg_flags result = exec_cmd_arg("corners=", cmd_file::AT_AFTER_STARTUP);

    EXPECT_EQ(cmdarg_flags::NONE, result);
    EXPECT_FALSE(g_use_center_mag);
}

TEST_F(TestParameterCommand, cornersFourValues)
{
    ValueSaver saved_use_center_mag{g_use_center_mag, true};
    ValueSaver saved_x_min{g_x_min, 111.0};
    ValueSaver saved_x_3rd{g_x_3rd, 222.0};
    ValueSaver saved_x_max{g_x_max, 333.0};
    ValueSaver saved_y_min{g_y_min, 444.0};
    ValueSaver saved_y_3rd{g_y_3rd, 555.0};
    ValueSaver saved_y_max{g_y_max, 666.0};

    const cmdarg_flags result = exec_cmd_arg("corners=1/2/3/4", cmd_file::AT_AFTER_STARTUP);

    EXPECT_EQ(cmdarg_flags::FRACTAL_PARAM, result);
    EXPECT_FALSE(g_use_center_mag);
    EXPECT_EQ(1.0, g_x_min);
    EXPECT_EQ(1.0, g_x_3rd);
    EXPECT_EQ(2.0, g_x_max);
    EXPECT_EQ(3.0, g_y_min);
    EXPECT_EQ(3.0, g_y_3rd);
    EXPECT_EQ(4.0, g_y_max);
}

TEST_F(TestParameterCommand, cornersSixValues)
{
    ValueSaver saved_use_center_mag{g_use_center_mag, true};
    ValueSaver saved_x_min{g_x_min, 111.0};
    ValueSaver saved_x_3rd{g_x_3rd, 222.0};
    ValueSaver saved_x_max{g_x_max, 333.0};
    ValueSaver saved_y_min{g_y_min, 444.0};
    ValueSaver saved_y_3rd{g_y_3rd, 555.0};
    ValueSaver saved_y_max{g_y_max, 666.0};

    const cmdarg_flags result = exec_cmd_arg("corners=1/2/3/4/5/6", cmd_file::AT_AFTER_STARTUP);

    EXPECT_EQ(cmdarg_flags::FRACTAL_PARAM, result);
    EXPECT_FALSE(g_use_center_mag);
    EXPECT_EQ(1.0, g_x_min);
    EXPECT_EQ(2.0, g_x_max);
    EXPECT_EQ(3.0, g_y_min);
    EXPECT_EQ(4.0, g_y_max);
    EXPECT_EQ(5.0, g_x_3rd);
    EXPECT_EQ(6.0, g_y_3rd);
}

TEST_F(TestParameterCommand, orbitCornersFourValues)
{
    ValueSaver saved_set_orbit_corners{g_set_orbit_corners, false};
    ValueSaver saved_orbit_corner_min_x{g_orbit_corner_min_x, 999.0};
    ValueSaver saved_orbit_corner_max_x{g_orbit_corner_max_x, 999.0};
    ValueSaver saved_orbit_corner_3_x{g_orbit_corner_3_x, 999.0};
    ValueSaver saved_orbit_corner_min_y{g_orbit_corner_min_y, 999.0};
    ValueSaver saved_orbit_corner_max_y{g_orbit_corner_max_y, 999.0};
    ValueSaver saved_orbit_corner_3_y{g_orbit_corner_3_y, 999.0};

    const cmdarg_flags result = exec_cmd_arg("orbitcorners=1/2/3/4", cmd_file::AT_AFTER_STARTUP);

    EXPECT_EQ(cmdarg_flags::FRACTAL_PARAM, result);
    EXPECT_TRUE(g_set_orbit_corners);
    EXPECT_TRUE(g_keep_screen_coords);
    EXPECT_EQ(1.0, g_orbit_corner_min_x);
    EXPECT_EQ(2.0, g_orbit_corner_max_x);
    EXPECT_EQ(3.0, g_orbit_corner_min_y);
    EXPECT_EQ(4.0, g_orbit_corner_max_y);
    EXPECT_EQ(1.0, g_orbit_corner_3_x);
    EXPECT_EQ(3.0, g_orbit_corner_3_y);
}

TEST_F(TestParameterCommand, orbitCornersSixValues)
{
    ValueSaver saved_set_orbit_corners{g_set_orbit_corners, false};
    ValueSaver saved_orbit_corner_min_x{g_orbit_corner_min_x, 999.0};
    ValueSaver saved_orbit_corner_max_x{g_orbit_corner_max_x, 999.0};
    ValueSaver saved_orbit_corner_3_x{g_orbit_corner_3_x, 999.0};
    ValueSaver saved_orbit_corner_min_y{g_orbit_corner_min_y, 999.0};
    ValueSaver saved_orbit_corner_max_y{g_orbit_corner_max_y, 999.0};
    ValueSaver saved_orbit_corner_3_y{g_orbit_corner_3_y, 999.0};

    const cmdarg_flags result = exec_cmd_arg("orbitcorners=1/2/3/4/5/6", cmd_file::AT_AFTER_STARTUP);

    EXPECT_EQ(cmdarg_flags::FRACTAL_PARAM, result);
    EXPECT_TRUE(g_set_orbit_corners);
    EXPECT_TRUE(g_keep_screen_coords);
    EXPECT_EQ(1.0, g_orbit_corner_min_x);
    EXPECT_EQ(2.0, g_orbit_corner_max_x);
    EXPECT_EQ(3.0, g_orbit_corner_min_y);
    EXPECT_EQ(4.0, g_orbit_corner_max_y);
    EXPECT_EQ(5.0, g_orbit_corner_3_x);
    EXPECT_EQ(6.0, g_orbit_corner_3_y);
}

TEST_F(TestParameterCommand, screenCoordsYes)
{
    ValueSaver saved_keep_screen_coords{g_keep_screen_coords, false};

    const cmdarg_flags result = exec_cmd_arg("screencoords=y", cmd_file::AT_AFTER_STARTUP);

    EXPECT_EQ(cmdarg_flags::FRACTAL_PARAM, result);
    EXPECT_TRUE(g_keep_screen_coords);
}

TEST_F(TestParameterCommand, orbitDrawModeLine)
{
    ValueSaver saved_orbit_draw_mode{g_draw_mode, '!'};

    const cmdarg_flags result = exec_cmd_arg("orbitdrawmode=l", cmd_file::AT_AFTER_STARTUP);

    EXPECT_EQ(cmdarg_flags::FRACTAL_PARAM, result);
    EXPECT_EQ('l', g_draw_mode);
}

TEST_F(TestParameterCommand, orbitDrawModeRectangle)
{
    ValueSaver saved_orbit_draw_mode{g_draw_mode, '!'};

    const cmdarg_flags result = exec_cmd_arg("orbitdrawmode=r", cmd_file::AT_AFTER_STARTUP);

    EXPECT_EQ(cmdarg_flags::FRACTAL_PARAM, result);
    EXPECT_EQ('r', g_draw_mode);
}

TEST_F(TestParameterCommand, orbitDrawModeFunction)
{
    ValueSaver saved_orbit_draw_mode{g_draw_mode, '!'};

    const cmdarg_flags result = exec_cmd_arg("orbitdrawmode=f", cmd_file::AT_AFTER_STARTUP);

    EXPECT_EQ(cmdarg_flags::FRACTAL_PARAM, result);
    EXPECT_EQ('f', g_draw_mode);
}

TEST_F(TestParameterCommand, viewWindowsdDefaults)
{
    ValueSaver saved_view_window{g_view_window, false};
    ValueSaver saved_view_reduction{g_view_reduction, 999.0f};
    ValueSaver saved_screen_aspect{g_screen_aspect, 0.75f};
    ValueSaver saved_final_aspect_ration{g_final_aspect_ratio, 999.0f};
    ValueSaver saved_view_crop{g_view_crop, false};
    ValueSaver saved_view_x_dots{g_view_x_dots, 999};
    ValueSaver saved_view_y_dots{g_view_y_dots, 999};

    const cmdarg_flags result = exec_cmd_arg("viewwindows", cmd_file::AT_AFTER_STARTUP);

    EXPECT_EQ(cmdarg_flags::FRACTAL_PARAM, result);
    EXPECT_TRUE(g_view_window);
    EXPECT_EQ(4.2f, g_view_reduction);
    EXPECT_EQ(0.75f, g_screen_aspect);
    EXPECT_EQ(0.75f, g_final_aspect_ratio);
    EXPECT_TRUE(g_view_crop);
    EXPECT_EQ(0, g_view_x_dots);
    EXPECT_EQ(0, g_view_y_dots);
}

TEST_F(TestParameterCommand, viewWindowsOneValue)
{
    ValueSaver saved_view_window{g_view_window, false};
    ValueSaver saved_view_reduction{g_view_reduction, 999.0f};
    ValueSaver saved_screen_aspect{g_screen_aspect, 0.75f};
    ValueSaver saved_final_aspect_ration{g_final_aspect_ratio, 999.0f};
    ValueSaver saved_view_crop{g_view_crop, false};
    ValueSaver saved_view_x_dots{g_view_x_dots, 999};
    ValueSaver saved_view_y_dots{g_view_y_dots, 999};

    const cmdarg_flags result = exec_cmd_arg("viewwindows=2", cmd_file::AT_AFTER_STARTUP);

    EXPECT_EQ(cmdarg_flags::FRACTAL_PARAM, result);
    EXPECT_EQ(2.0f, g_view_reduction);
}

TEST_F(TestParameterCommand, viewWindowsTwoValues)
{
    ValueSaver saved_view_window{g_view_window, false};
    ValueSaver saved_view_reduction{g_view_reduction, 999.0f};
    ValueSaver saved_screen_aspect{g_screen_aspect, 0.75f};
    ValueSaver saved_final_aspect_ration{g_final_aspect_ratio, 999.0f};
    ValueSaver saved_view_crop{g_view_crop, false};
    ValueSaver saved_view_x_dots{g_view_x_dots, 999};
    ValueSaver saved_view_y_dots{g_view_y_dots, 999};

    const cmdarg_flags result = exec_cmd_arg("viewwindows=2/3", cmd_file::AT_AFTER_STARTUP);

    EXPECT_EQ(cmdarg_flags::FRACTAL_PARAM, result);
    EXPECT_EQ(2.0f, g_view_reduction);
    EXPECT_EQ(3.0f, g_final_aspect_ratio);
}

TEST_F(TestParameterCommand, viewWindowsThreeValues)
{
    ValueSaver saved_view_window{g_view_window, false};
    ValueSaver saved_view_reduction{g_view_reduction, 999.0f};
    ValueSaver saved_screen_aspect{g_screen_aspect, 0.75f};
    ValueSaver saved_final_aspect_ration{g_final_aspect_ratio, 999.0f};
    ValueSaver saved_view_crop{g_view_crop, true};
    ValueSaver saved_view_x_dots{g_view_x_dots, 999};
    ValueSaver saved_view_y_dots{g_view_y_dots, 999};

    const cmdarg_flags result = exec_cmd_arg("viewwindows=2/3/n", cmd_file::AT_AFTER_STARTUP);

    EXPECT_EQ(cmdarg_flags::FRACTAL_PARAM, result);
    EXPECT_EQ(2.0f, g_view_reduction);
    EXPECT_EQ(3.0f, g_final_aspect_ratio);
    EXPECT_FALSE(g_view_crop);
}

TEST_F(TestParameterCommand, viewWindowsFourValues)
{
    ValueSaver saved_view_window{g_view_window, false};
    ValueSaver saved_view_reduction{g_view_reduction, 999.0f};
    ValueSaver saved_screen_aspect{g_screen_aspect, 0.75f};
    ValueSaver saved_final_aspect_ration{g_final_aspect_ratio, 999.0f};
    ValueSaver saved_view_crop{g_view_crop, true};
    ValueSaver saved_view_x_dots{g_view_x_dots, 999};
    ValueSaver saved_view_y_dots{g_view_y_dots, 999};

    const cmdarg_flags result = exec_cmd_arg("viewwindows=2/3/n/800", cmd_file::AT_AFTER_STARTUP);

    EXPECT_EQ(cmdarg_flags::FRACTAL_PARAM, result);
    EXPECT_EQ(2.0f, g_view_reduction);
    EXPECT_EQ(3.0f, g_final_aspect_ratio);
    EXPECT_FALSE(g_view_crop);
    EXPECT_EQ(800, g_view_x_dots);
}

TEST_F(TestParameterCommand, viewWindowsFiveValues)
{
    ValueSaver saved_view_window{g_view_window, false};
    ValueSaver saved_view_reduction{g_view_reduction, 999.0f};
    ValueSaver saved_screen_aspect{g_screen_aspect, 0.75f};
    ValueSaver saved_final_aspect_ration{g_final_aspect_ratio, 999.0f};
    ValueSaver saved_view_crop{g_view_crop, true};
    ValueSaver saved_view_x_dots{g_view_x_dots, 999};
    ValueSaver saved_view_y_dots{g_view_y_dots, 999};

    const cmdarg_flags result = exec_cmd_arg("viewwindows=2/3/n/800/600", cmd_file::AT_AFTER_STARTUP);

    EXPECT_EQ(cmdarg_flags::FRACTAL_PARAM, result);
    EXPECT_EQ(2.0f, g_view_reduction);
    EXPECT_EQ(3.0f, g_final_aspect_ratio);
    EXPECT_FALSE(g_view_crop);
    EXPECT_EQ(800, g_view_x_dots);
    EXPECT_EQ(600, g_view_y_dots);
}

TEST_F(TestParameterCommand, centerMagOn)
{
    ValueSaver saved_use_center_mag{g_use_center_mag, false};
    VALUE_UNCHANGED(g_x_min, 999.0);
    VALUE_UNCHANGED(g_x_max, 999.0);
    VALUE_UNCHANGED(g_x_3rd, 999.0);
    VALUE_UNCHANGED(g_y_min, 999.0);
    VALUE_UNCHANGED(g_y_max, 999.0);
    VALUE_UNCHANGED(g_y_3rd, 999.0);

    const cmdarg_flags result = exec_cmd_arg("center-mag", cmd_file::AT_AFTER_STARTUP);

    EXPECT_EQ(cmdarg_flags::NONE, result);
    EXPECT_TRUE(g_use_center_mag);
}

TEST_F(TestParameterCommand, centerMagThreeValues)
{
    ValueSaver saved_use_center_mag{g_use_center_mag, false};
    ValueSaver saved_x_min{g_x_min, 999.0};
    ValueSaver saved_x_max{g_x_max, 999.0};
    ValueSaver saved_x_3rd{g_x_3rd, 999.0};
    ValueSaver saved_y_min{g_y_min, 999.0};
    ValueSaver saved_y_max{g_y_max, 999.0};
    ValueSaver saved_y_3rd{g_y_3rd, 999.0};

    const cmdarg_flags result = exec_cmd_arg("center-mag=2/4/3", cmd_file::AT_AFTER_STARTUP);

    EXPECT_EQ(cmdarg_flags::FRACTAL_PARAM, result);
    EXPECT_TRUE(g_use_center_mag);
    EXPECT_NEAR(1.555555, g_x_min, 1e-6);
    EXPECT_NEAR(2.444444, g_x_max, 1e-6);
    EXPECT_NEAR(1.555555, g_x_3rd, 1e-6);
    EXPECT_NEAR(3.666666, g_y_min, 1e-6);
    EXPECT_NEAR(4.333333, g_y_max, 1e-6);
    EXPECT_NEAR(3.666666, g_y_3rd, 1e-6);
}

TEST_F(TestParameterCommand, centerMagFourValues)
{
    ValueSaver saved_use_center_mag{g_use_center_mag, false};
    ValueSaver saved_x_min{g_x_min, 999.0};
    ValueSaver saved_x_max{g_x_max, 999.0};
    ValueSaver saved_x_3rd{g_x_3rd, 999.0};
    ValueSaver saved_y_min{g_y_min, 999.0};
    ValueSaver saved_y_max{g_y_max, 999.0};
    ValueSaver saved_y_3rd{g_y_3rd, 999.0};

    const cmdarg_flags result = exec_cmd_arg("center-mag=2/4/3/90", cmd_file::AT_AFTER_STARTUP);

    EXPECT_EQ(cmdarg_flags::FRACTAL_PARAM, result);
    EXPECT_TRUE(g_use_center_mag);
    EXPECT_NEAR(1.995061, g_x_min, 1e-6);
    EXPECT_NEAR(2.004938, g_x_max, 1e-6);
    EXPECT_NEAR(1.995061, g_x_3rd, 1e-6);
    EXPECT_NEAR(3.666666, g_y_min, 1e-6);
    EXPECT_NEAR(4.333333, g_y_max, 1e-6);
    EXPECT_NEAR(3.666666, g_y_3rd, 1e-6);
}

TEST_F(TestParameterCommand, centerMagFiveValues)
{
    ValueSaver saved_use_center_mag{g_use_center_mag, false};
    ValueSaver saved_x_min{g_x_min, 999.0};
    ValueSaver saved_x_max{g_x_max, 999.0};
    ValueSaver saved_x_3rd{g_x_3rd, 999.0};
    ValueSaver saved_y_min{g_y_min, 999.0};
    ValueSaver saved_y_max{g_y_max, 999.0};
    ValueSaver saved_y_3rd{g_y_3rd, 999.0};

    const cmdarg_flags result = exec_cmd_arg("center-mag=2/4/3/90/45", cmd_file::AT_AFTER_STARTUP);

    EXPECT_EQ(cmdarg_flags::FRACTAL_PARAM, result);
    EXPECT_TRUE(g_use_center_mag);
    EXPECT_NEAR(2.232210, g_x_min, 1e-6);
    EXPECT_NEAR(1.767789, g_x_max, 1e-6);
    EXPECT_NEAR(1.760805, g_x_3rd, 1e-6);
    EXPECT_NEAR(3.760805, g_y_min, 1e-6);
    EXPECT_NEAR(4.239194, g_y_max, 1e-6);
    EXPECT_NEAR(3.767789, g_y_3rd, 1e-6);
}

TEST_F(TestParameterCommand, aspectDrift)
{
    ValueSaver saved_aspect_drift{g_aspect_drift, 999.0};

    const cmdarg_flags result = exec_cmd_arg("aspectdrift=12", cmd_file::AT_AFTER_STARTUP);

    EXPECT_EQ(cmdarg_flags::FRACTAL_PARAM, result);
    EXPECT_EQ(12.0f, g_aspect_drift);
}

TEST_F(TestParameterCommand, invertComputedRadius)
{
    ValueSaver saved_inversion0{g_inversion[0], 999.0};
    ValueSaver saved_inversion1{g_inversion[1], 999.0};
    ValueSaver saved_inversion2{g_inversion[2], 999.0};
    ValueSaver saved_invert{g_invert, 999};

    const cmdarg_flags result = exec_cmd_arg("invert=-1", cmd_file::AT_AFTER_STARTUP);

    EXPECT_EQ(cmdarg_flags::FRACTAL_PARAM, result);
    EXPECT_EQ(-1.0, g_inversion[0]);
    EXPECT_EQ(1, g_invert);
}

TEST_F(TestParameterCommand, invertDefaultCenter)
{
    ValueSaver saved_inversion0{g_inversion[0], 999.0};
    ValueSaver saved_inversion1{g_inversion[1], 999.0};
    ValueSaver saved_inversion2{g_inversion[2], 999.0};
    ValueSaver saved_invert{g_invert, 999};

    const cmdarg_flags result = exec_cmd_arg("invert=1", cmd_file::AT_AFTER_STARTUP);

    EXPECT_EQ(cmdarg_flags::FRACTAL_PARAM, result);
    EXPECT_EQ(1.0, g_inversion[0]);
    EXPECT_EQ(1, g_invert);
}

TEST_F(TestParameterCommand, invertRadiusCenter)
{
    ValueSaver saved_inversion0{g_inversion[0], 999.0};
    ValueSaver saved_inversion1{g_inversion[1], 999.0};
    ValueSaver saved_inversion2{g_inversion[2], 999.0};
    ValueSaver saved_invert{g_invert, 999};

    const cmdarg_flags result = exec_cmd_arg("invert=1/100/200", cmd_file::AT_AFTER_STARTUP);

    EXPECT_EQ(cmdarg_flags::FRACTAL_PARAM, result);
    EXPECT_EQ(1.0, g_inversion[0]);
    EXPECT_EQ(100.0, g_inversion[1]);
    EXPECT_EQ(200.0, g_inversion[2]);
    EXPECT_EQ(3, g_invert);
}

TEST_F(TestParameterCommand, oldDemmColorsYes)
{
    ValueSaver saved_old_demm_colors{g_old_demm_colors, false};

    const cmdarg_flags result = exec_cmd_arg("olddemmcolors=y", cmd_file::AT_AFTER_STARTUP);

    EXPECT_EQ(cmdarg_flags::NONE, result);
    EXPECT_TRUE(g_old_demm_colors);
}

TEST_F(TestParameterCommand, askVideoYes)
{
    ValueSaver saved_ask_video{g_ask_video, false};

    const cmdarg_flags result = exec_cmd_arg("askvideo=y", cmd_file::AT_AFTER_STARTUP);

    EXPECT_EQ(cmdarg_flags::NONE, result);
    EXPECT_TRUE(g_ask_video);
}

TEST_F(TestParameterCommand, ramVideoIgnored)
{
    const cmdarg_flags result = exec_cmd_arg("ramvideo=64", cmd_file::AT_AFTER_STARTUP);

    EXPECT_EQ(cmdarg_flags::NONE, result);
}

TEST_F(TestParameterCommand, floatYes)
{
    ValueSaver saved_user_float_flag{g_user_float_flag, false};

    const cmdarg_flags result = exec_cmd_arg("float=y", cmd_file::AT_AFTER_STARTUP);

    EXPECT_EQ(cmdarg_flags::FRACTAL_PARAM | cmdarg_flags::PARAM_3D, result);
    EXPECT_TRUE(g_user_float_flag);
}

TEST_F(TestParameterCommand, fastRestoreYes)
{
    ValueSaver saved_fast_restore{g_fast_restore, false};

    const cmdarg_flags result = exec_cmd_arg("fastrestore=y", cmd_file::AT_AFTER_STARTUP);

    EXPECT_EQ(cmdarg_flags::NONE, result);
    EXPECT_TRUE(g_fast_restore);
}

TEST_F(TestParameterCommand, orgFrmDir)
{
    ValueSaver saved_organize_formulas_search{g_organize_formulas_search, false};
    ValueSaver saved_organize_formulas_dir{g_organize_formulas_dir, "fmeh"};

    const cmdarg_flags result = exec_cmd_arg("orgfrmdir=" ID_TEST_DATA_DIR, cmd_file::AT_AFTER_STARTUP);

    EXPECT_EQ(cmdarg_flags::NONE, result);
    EXPECT_TRUE(g_organize_formulas_search);
    EXPECT_EQ(ID_TEST_DATA_DIR, adjust_dir(g_organize_formulas_dir));
}

TEST_F(TestParameterCommand, biomorph)
{
    ValueSaver saved_user_biomorph_value{g_user_biomorph_value, 999};

    const cmdarg_flags result = exec_cmd_arg("biomorph=52", cmd_file::AT_AFTER_STARTUP);

    EXPECT_EQ(cmdarg_flags::FRACTAL_PARAM, result);
    EXPECT_EQ(52, g_user_biomorph_value);
}

TEST_F(TestParameterCommand, orbitSaveRaw)
{
    ValueSaver saved_orbit_save_flags{g_orbit_save_flags, 0};

    const cmdarg_flags result = exec_cmd_arg("orbitsave=y", cmd_file::AT_AFTER_STARTUP);

    EXPECT_EQ(cmdarg_flags::FRACTAL_PARAM, result);
    EXPECT_EQ(osf_raw, g_orbit_save_flags);
}

TEST_F(TestParameterCommand, orbitSaveMidi)
{
    ValueSaver saved_orbit_save_flags{g_orbit_save_flags, 0};

    const cmdarg_flags result = exec_cmd_arg("orbitsave=s", cmd_file::AT_AFTER_STARTUP);

    EXPECT_EQ(cmdarg_flags::FRACTAL_PARAM, result);
    EXPECT_EQ(osf_midi | osf_raw, g_orbit_save_flags);
}

TEST_F(TestParameterCommand, orbitSaveName)
{
    ValueSaver saved_orbit_save_name{g_orbit_save_name, "foo.txt"};

    const cmdarg_flags result = exec_cmd_arg("orbitsavename=smeagle.txt", cmd_file::AT_AFTER_STARTUP);

    EXPECT_EQ(cmdarg_flags::FRACTAL_PARAM, result);
    EXPECT_EQ("smeagle.txt", g_orbit_save_name);
}

TEST_F(TestParameterCommand, bailOut)
{
    ValueSaver saved_bail_out{g_bail_out, -1};

    const cmdarg_flags result = exec_cmd_arg("bailout=50", cmd_file::AT_AFTER_STARTUP);

    EXPECT_EQ(cmdarg_flags::FRACTAL_PARAM, result);
    EXPECT_EQ(50L, g_bail_out);
}

TEST_F(TestParameterCommand, bailOutTestMod)
{
    ValueSaver saved_bail_out_test{g_bail_out_test, bailouts::Real};

    const cmdarg_flags result = exec_cmd_arg("bailoutest=mod", cmd_file::AT_AFTER_STARTUP);

    EXPECT_EQ(cmdarg_flags::FRACTAL_PARAM, result);
    EXPECT_EQ(bailouts::Mod, g_bail_out_test);
}

TEST_F(TestParameterCommand, bailOutTestBadValue)
{
    expect_stop_msg();
    ValueSaver saved_bail_out_test{g_bail_out_test, bailouts::Real};

    const cmdarg_flags result = exec_cmd_arg("bailoutest=foo", cmd_file::AT_AFTER_STARTUP);

    EXPECT_EQ(cmdarg_flags::ERROR, result);
    EXPECT_EQ(bailouts::Real, g_bail_out_test);
}

TEST_F(TestParameterCommand, symmetryXAxis)
{
    ValueSaver save_force_symmetry{g_force_symmetry, symmetry_type::XY_AXIS};

    const cmdarg_flags result = exec_cmd_arg("symmetry=xaxis", cmd_file::AT_AFTER_STARTUP);

    EXPECT_EQ(cmdarg_flags::FRACTAL_PARAM, result);
    EXPECT_EQ(symmetry_type::X_AXIS, g_force_symmetry);
}

TEST_F(TestParameterCommand, symmetryBadValue)
{
    expect_stop_msg();
    ValueSaver save_force_symmetry{g_force_symmetry, symmetry_type::XY_AXIS};

    const cmdarg_flags result = exec_cmd_arg("symmetry=fmeh", cmd_file::AT_AFTER_STARTUP);

    EXPECT_EQ(cmdarg_flags::ERROR, result);
    EXPECT_EQ(symmetry_type::XY_AXIS, g_force_symmetry);
}

TEST_F(TestParameterCommand, soundYes)
{
    ValueSaver saved_sound_flag{g_sound_flag, 0};

    const cmdarg_flags result = exec_cmd_arg("sound=yes", cmd_file::AT_AFTER_STARTUP);

    EXPECT_EQ(cmdarg_flags::NONE, result);
    EXPECT_EQ(SOUNDFLAG_BEEP | SOUNDFLAG_SPEAKER, g_sound_flag);
}

TEST_F(TestParameterCommand, soundNo)
{
    ValueSaver saved_sound_flag{g_sound_flag, 0};

    const cmdarg_flags result = exec_cmd_arg("sound=no", cmd_file::AT_AFTER_STARTUP);

    EXPECT_EQ(cmdarg_flags::NONE, result);
    EXPECT_EQ(SOUNDFLAG_OFF | SOUNDFLAG_SPEAKER, g_sound_flag);
}

TEST_F(TestParameterCommand, soundBeep)
{
    ValueSaver saved_sound_flag{g_sound_flag, 0};

    const cmdarg_flags result = exec_cmd_arg("sound=b", cmd_file::AT_AFTER_STARTUP);

    EXPECT_EQ(cmdarg_flags::NONE, result);
    EXPECT_EQ(SOUNDFLAG_BEEP | SOUNDFLAG_SPEAKER, g_sound_flag);
}

TEST_F(TestParameterCommand, soundOrbitX)
{
    ValueSaver saved_sound_flag{g_sound_flag, 0};

    const cmdarg_flags result = exec_cmd_arg("sound=x", cmd_file::AT_AFTER_STARTUP);

    EXPECT_EQ(cmdarg_flags::NONE, result);
    EXPECT_EQ(SOUNDFLAG_X | SOUNDFLAG_SPEAKER, g_sound_flag);
}

TEST_F(TestParameterCommand, soundOrbitY)
{
    ValueSaver saved_sound_flag{g_sound_flag, 0};

    const cmdarg_flags result = exec_cmd_arg("sound=y", cmd_file::AT_AFTER_STARTUP);

    EXPECT_EQ(cmdarg_flags::NONE, result);
    EXPECT_EQ(SOUNDFLAG_Y | SOUNDFLAG_SPEAKER, g_sound_flag);
}

TEST_F(TestParameterCommand, soundOrbitZ)
{
    ValueSaver saved_sound_flag{g_sound_flag, 0};

    const cmdarg_flags result = exec_cmd_arg("sound=z", cmd_file::AT_AFTER_STARTUP);

    EXPECT_EQ(cmdarg_flags::NONE, result);
    EXPECT_EQ(SOUNDFLAG_Z | SOUNDFLAG_SPEAKER, g_sound_flag);
}

TEST_F(TestParameterCommand, soundOrbitZMidi)
{
    ValueSaver saved_sound_flag{g_sound_flag, 0};

    const cmdarg_flags result = exec_cmd_arg("sound=z/m", cmd_file::AT_AFTER_STARTUP);

    EXPECT_EQ(cmdarg_flags::NONE, result);
    EXPECT_EQ(SOUNDFLAG_Z | SOUNDFLAG_MIDI, g_sound_flag);
}

TEST_F(TestParameterCommand, soundOrbitZSpeaker)
{
    ValueSaver saved_sound_flag{g_sound_flag, 0};

    const cmdarg_flags result = exec_cmd_arg("sound=z/p", cmd_file::AT_AFTER_STARTUP);

    EXPECT_EQ(cmdarg_flags::NONE, result);
    EXPECT_EQ(SOUNDFLAG_Z | SOUNDFLAG_SPEAKER, g_sound_flag);
}

TEST_F(TestParameterCommand, soundOrbitZQuantized)
{
    ValueSaver saved_sound_flag{g_sound_flag, 0};

    const cmdarg_flags result = exec_cmd_arg("sound=z/q", cmd_file::AT_AFTER_STARTUP);

    EXPECT_EQ(cmdarg_flags::NONE, result);
    EXPECT_EQ(SOUNDFLAG_Z | SOUNDFLAG_QUANTIZED, g_sound_flag);
}

TEST_F(TestParameterCommand, hertzValue)
{
    ValueSaver saved_base_hertz{g_base_hertz, -999};

    const cmdarg_flags result = exec_cmd_arg("hertz=100", cmd_file::AT_AFTER_STARTUP);

    EXPECT_EQ(cmdarg_flags::NONE, result);
    EXPECT_EQ(100, g_base_hertz);
}

TEST_F(TestParameterCommand, volumeValue)
{
    ValueSaver saved_fm_volume{g_fm_volume, -999};

    const cmdarg_flags result = exec_cmd_arg("volume=63", cmd_file::AT_AFTER_STARTUP);

    EXPECT_EQ(cmdarg_flags::NONE, result);
    EXPECT_EQ(63, g_fm_volume);
}

TEST_F(TestParameterCommand, attenuateNo)
{
    ValueSaver saved_hi_attenuation{g_hi_attenuation, -999};

    const cmdarg_flags result = exec_cmd_arg("attenuate=n", cmd_file::AT_AFTER_STARTUP);

    EXPECT_EQ(cmdarg_flags::NONE, result);
    EXPECT_EQ(0, g_hi_attenuation);
}

TEST_F(TestParameterCommand, attenuateLow)
{
    ValueSaver saved_hi_attenuation{g_hi_attenuation, -999};

    const cmdarg_flags result = exec_cmd_arg("attenuate=l", cmd_file::AT_AFTER_STARTUP);

    EXPECT_EQ(cmdarg_flags::NONE, result);
    EXPECT_EQ(1, g_hi_attenuation);
}

TEST_F(TestParameterCommand, attenuateMiddle)
{
    ValueSaver saved_hi_attenuation{g_hi_attenuation, -999};

    const cmdarg_flags result = exec_cmd_arg("attenuate=m", cmd_file::AT_AFTER_STARTUP);

    EXPECT_EQ(cmdarg_flags::NONE, result);
    EXPECT_EQ(2, g_hi_attenuation);
}

TEST_F(TestParameterCommand, attenuateHigh)
{
    ValueSaver saved_hi_attenuation{g_hi_attenuation, -999};

    const cmdarg_flags result = exec_cmd_arg("attenuate=h", cmd_file::AT_AFTER_STARTUP);

    EXPECT_EQ(cmdarg_flags::NONE, result);
    EXPECT_EQ(3, g_hi_attenuation);
}

TEST_F(TestParameterCommand, polyphonyValue)
{
    ValueSaver saved_polyphony{g_polyphony, -999};

    const cmdarg_flags result = exec_cmd_arg("polyphony=1", cmd_file::AT_AFTER_STARTUP);

    EXPECT_EQ(cmdarg_flags::NONE, result);
    EXPECT_EQ(0, g_polyphony);
}

TEST_F(TestParameterCommand, waveTypeValue)
{
    ValueSaver saved_fm_wavetype{g_fm_wavetype, -999};

    const cmdarg_flags result = exec_cmd_arg("wavetype=4", cmd_file::AT_AFTER_STARTUP);

    EXPECT_EQ(cmdarg_flags::NONE, result);
    EXPECT_EQ(4, g_fm_wavetype);
}

TEST_F(TestParameterCommand, attackValue)
{
    ValueSaver saved_fm_attack{g_fm_attack, -999};

    const cmdarg_flags result = exec_cmd_arg("attack=4", cmd_file::AT_AFTER_STARTUP);

    EXPECT_EQ(cmdarg_flags::NONE, result);
    EXPECT_EQ(4, g_fm_attack);
}

TEST_F(TestParameterCommand, decayValue)
{
    ValueSaver saved_fm_decay{g_fm_decay, -999};

    const cmdarg_flags result = exec_cmd_arg("decay=4", cmd_file::AT_AFTER_STARTUP);

    EXPECT_EQ(cmdarg_flags::NONE, result);
    EXPECT_EQ(4, g_fm_decay);
}

TEST_F(TestParameterCommand, sustainValue)
{
    ValueSaver saved_fm_sustain{g_fm_sustain, -999};

    const cmdarg_flags result = exec_cmd_arg("sustain=4", cmd_file::AT_AFTER_STARTUP);

    EXPECT_EQ(cmdarg_flags::NONE, result);
    EXPECT_EQ(4, g_fm_sustain);
}

TEST_F(TestParameterCommand, sReleaseValue)
{
    ValueSaver saved_fm_release{g_fm_release, -999};

    const cmdarg_flags result = exec_cmd_arg("srelease=4", cmd_file::AT_AFTER_STARTUP);

    EXPECT_EQ(cmdarg_flags::NONE, result);
    EXPECT_EQ(4, g_fm_release);
}

TEST_F(TestParameterCommand, scaleMapValues)
{
    ValueSaver saved_scale_map0{g_scale_map[0], -999};
    ValueSaver saved_scale_map1{g_scale_map[1], -999};
    ValueSaver saved_scale_map2{g_scale_map[2], -999};

    const cmdarg_flags result = exec_cmd_arg("scalemap=4/5", cmd_file::AT_AFTER_STARTUP);

    EXPECT_EQ(cmdarg_flags::NONE, result);
    EXPECT_EQ(4, g_scale_map[0]);
    EXPECT_EQ(5, g_scale_map[1]);
    EXPECT_EQ(-999, g_scale_map[2]);
}

TEST_F(TestParameterCommand, periodicityNo)
{
    ValueSaver saved_user_periodicity_value{g_user_periodicity_value, -999};

    const cmdarg_flags result = exec_cmd_arg("periodicity=no", cmd_file::AT_AFTER_STARTUP);

    EXPECT_EQ(cmdarg_flags::FRACTAL_PARAM, result);
    EXPECT_EQ(0, g_user_periodicity_value);
}

TEST_F(TestParameterCommand, periodicityYes)
{
    ValueSaver saved_user_periodicity_value{g_user_periodicity_value, -999};

    const cmdarg_flags result = exec_cmd_arg("periodicity=yes", cmd_file::AT_AFTER_STARTUP);

    EXPECT_EQ(cmdarg_flags::FRACTAL_PARAM, result);
    EXPECT_EQ(1, g_user_periodicity_value);
}

TEST_F(TestParameterCommand, periodicityShow)
{
    ValueSaver saved_user_periodicity_value{g_user_periodicity_value, -999};

    const cmdarg_flags result = exec_cmd_arg("periodicity=show", cmd_file::AT_AFTER_STARTUP);

    EXPECT_EQ(cmdarg_flags::FRACTAL_PARAM, result);
    EXPECT_EQ(-1, g_user_periodicity_value);
}

TEST_F(TestParameterCommand, periodicityValue)
{
    ValueSaver saved_user_periodicity_value{g_user_periodicity_value, -999};

    const cmdarg_flags result = exec_cmd_arg("periodicity=24", cmd_file::AT_AFTER_STARTUP);

    EXPECT_EQ(cmdarg_flags::FRACTAL_PARAM, result);
    EXPECT_EQ(24, g_user_periodicity_value);
}

TEST_F(TestParameterCommand, periodicityShowValue)
{
    ValueSaver saved_user_periodicity_value{g_user_periodicity_value, -999};

    const cmdarg_flags result = exec_cmd_arg("periodicity=-24", cmd_file::AT_AFTER_STARTUP);

    EXPECT_EQ(cmdarg_flags::FRACTAL_PARAM, result);
    EXPECT_EQ(-24, g_user_periodicity_value);
}

TEST_F(TestParameterCommand, periodicityValueClampedHigh)
{
    ValueSaver saved_user_periodicity_value{g_user_periodicity_value, -999};

    const cmdarg_flags result = exec_cmd_arg("periodicity=500", cmd_file::AT_AFTER_STARTUP);

    EXPECT_EQ(cmdarg_flags::FRACTAL_PARAM, result);
    EXPECT_EQ(255, g_user_periodicity_value);
}

TEST_F(TestParameterCommand, periodicityValueClampedLow)
{
    ValueSaver saved_user_periodicity_value{g_user_periodicity_value, -999};

    const cmdarg_flags result = exec_cmd_arg("periodicity=-500", cmd_file::AT_AFTER_STARTUP);

    EXPECT_EQ(cmdarg_flags::FRACTAL_PARAM, result);
    EXPECT_EQ(-255, g_user_periodicity_value);
}

TEST_F(TestParameterCommand, logMapNo)
{
    ValueSaver saved_log_map_auto_calculate{g_log_map_auto_calculate, true};

    const cmdarg_flags result = exec_cmd_arg("logmap=no", cmd_file::AT_AFTER_STARTUP);

    EXPECT_EQ(cmdarg_flags::FRACTAL_PARAM, result);
    EXPECT_FALSE(g_log_map_auto_calculate);
    EXPECT_EQ(0L, g_log_map_flag);
}

TEST_F(TestParameterCommand, logMapYes)
{
    ValueSaver saved_log_map_auto_calculate{g_log_map_auto_calculate, true};

    const cmdarg_flags result = exec_cmd_arg("logmap=yes", cmd_file::AT_AFTER_STARTUP);

    EXPECT_EQ(cmdarg_flags::FRACTAL_PARAM, result);
    EXPECT_FALSE(g_log_map_auto_calculate);
    EXPECT_EQ(1L, g_log_map_flag);
}

TEST_F(TestParameterCommand, logMapOld)
{
    ValueSaver saved_log_map_auto_calculate{g_log_map_auto_calculate, true};

    const cmdarg_flags result = exec_cmd_arg("logmap=old", cmd_file::AT_AFTER_STARTUP);

    EXPECT_EQ(cmdarg_flags::FRACTAL_PARAM, result);
    EXPECT_FALSE(g_log_map_auto_calculate);
    EXPECT_EQ(-1L, g_log_map_flag);
}

TEST_F(TestParameterCommand, logMapValue)
{
    ValueSaver saved_log_map_auto_calculate{g_log_map_auto_calculate, true};

    const cmdarg_flags result = exec_cmd_arg("logmap=15", cmd_file::AT_AFTER_STARTUP);

    EXPECT_EQ(cmdarg_flags::FRACTAL_PARAM, result);
    EXPECT_FALSE(g_log_map_auto_calculate);
    EXPECT_EQ(15L, g_log_map_flag);
}

TEST_F(TestParameterCommand, logModeFly)
{
    ValueSaver saved_log_map_fly_calculate{g_log_map_fly_calculate, -999};
    ValueSaver saved_log_map_auto_calculate{g_log_map_auto_calculate, true};

    const cmdarg_flags result = exec_cmd_arg("logmode=f", cmd_file::AT_AFTER_STARTUP);

    EXPECT_EQ(cmdarg_flags::FRACTAL_PARAM, result);
    EXPECT_EQ(1, g_log_map_fly_calculate);
    EXPECT_FALSE(g_log_map_auto_calculate);
}

TEST_F(TestParameterCommand, logModeTable)
{
    ValueSaver saved_log_map_fly_calculate{g_log_map_fly_calculate, -999};
    ValueSaver saved_log_map_auto_calculate{g_log_map_auto_calculate, true};

    const cmdarg_flags result = exec_cmd_arg("logmode=t", cmd_file::AT_AFTER_STARTUP);

    EXPECT_EQ(cmdarg_flags::FRACTAL_PARAM, result);
    EXPECT_EQ(2, g_log_map_fly_calculate);
    EXPECT_FALSE(g_log_map_auto_calculate);
}

TEST_F(TestParameterCommand, logModeAuto)
{
    ValueSaver saved_log_map_fly_calculate{g_log_map_fly_calculate, -999};
    ValueSaver saved_log_map_auto_calculate{g_log_map_auto_calculate, false};

    const cmdarg_flags result = exec_cmd_arg("logmode=a", cmd_file::AT_AFTER_STARTUP);

    EXPECT_EQ(cmdarg_flags::FRACTAL_PARAM, result);
    EXPECT_EQ(0, g_log_map_fly_calculate);
    EXPECT_TRUE(g_log_map_auto_calculate);
}

TEST_F(TestParameterCommand, debugFlagValueViaDebug)
{
    ValueSaver saved_debug_flag{g_debug_flag, debug_flags::none};
    ValueSaver saved_timer_flag{g_timer_flag, true};

    const cmdarg_flags result = exec_cmd_arg("debug=300", cmd_file::AT_AFTER_STARTUP);

    EXPECT_EQ(cmdarg_flags::NONE, result);
    EXPECT_EQ(debug_flags::prevent_miim, g_debug_flag);
    EXPECT_FALSE(g_timer_flag);
}

TEST_F(TestParameterCommand, debugFlagValue)
{
    ValueSaver saved_debug_flag{g_debug_flag, debug_flags::none};
    ValueSaver saved_timer_flag{g_timer_flag, true};

    const cmdarg_flags result = exec_cmd_arg("debugflag=300", cmd_file::AT_AFTER_STARTUP);

    EXPECT_EQ(cmdarg_flags::NONE, result);
    EXPECT_EQ(debug_flags::prevent_miim, g_debug_flag);
    EXPECT_FALSE(g_timer_flag);
}

TEST_F(TestParameterCommand, debugFlagValueWithTimer)
{
    ValueSaver saved_debug_flag{g_debug_flag, debug_flags::none};
    ValueSaver saved_timer_flag{g_timer_flag, false};

    const cmdarg_flags result = exec_cmd_arg("debugflag=301", cmd_file::AT_AFTER_STARTUP);

    EXPECT_EQ(cmdarg_flags::NONE, result);
    EXPECT_EQ(debug_flags::prevent_miim, g_debug_flag);
    EXPECT_TRUE(g_timer_flag);
}

TEST_F(TestParameterCommand, randomSeedValue)
{
    ValueSaver save_random_seed{g_random_seed, -999};
    ValueSaver save_random_seed_flag{g_random_seed_flag, false};

    const cmdarg_flags result = exec_cmd_arg("rseed=301", cmd_file::AT_AFTER_STARTUP);

    EXPECT_EQ(cmdarg_flags::FRACTAL_PARAM, result);
    EXPECT_EQ(301, g_random_seed);
    EXPECT_TRUE(g_random_seed_flag);
}

TEST_F(TestParameterCommand, orbitDelayValue)
{
    ValueSaver save_orbit_delay{g_orbit_delay, -999};

    const cmdarg_flags result = exec_cmd_arg("orbitdelay=301", cmd_file::AT_AFTER_STARTUP);

    EXPECT_EQ(cmdarg_flags::NONE, result);
    EXPECT_EQ(301, g_orbit_delay);
}

TEST_F(TestParameterCommand, orbitIntervalValue)
{
    ValueSaver save_orbit_interval{g_orbit_interval, -999L};

    const cmdarg_flags result = exec_cmd_arg("orbitinterval=201", cmd_file::AT_AFTER_STARTUP);

    EXPECT_EQ(cmdarg_flags::NONE, result);
    EXPECT_EQ(201, g_orbit_interval);
}

TEST_F(TestParameterCommand, orbitIntervalValueTooSmall)
{
    ValueSaver save_orbit_interval{g_orbit_interval, -999L};

    const cmdarg_flags result = exec_cmd_arg("orbitinterval=0", cmd_file::AT_AFTER_STARTUP);

    EXPECT_EQ(cmdarg_flags::NONE, result);
    EXPECT_EQ(1, g_orbit_interval);
}

TEST_F(TestParameterCommand, orbitIntervalValueTooBig)
{
    ValueSaver save_orbit_interval{g_orbit_interval, -999L};

    const cmdarg_flags result = exec_cmd_arg("orbitinterval=256", cmd_file::AT_AFTER_STARTUP);

    EXPECT_EQ(cmdarg_flags::NONE, result);
    EXPECT_EQ(255, g_orbit_interval);
}

TEST_F(TestParameterCommand, showDotDefault)
{
    ValueSaver save_show_dot{g_show_dot, -999L};
    VALUE_UNCHANGED(g_auto_show_dot, '!');
    VALUE_UNCHANGED(g_size_dot, 10);

    const cmdarg_flags result = exec_cmd_arg("showdot", cmd_file::AT_AFTER_STARTUP);

    EXPECT_EQ(cmdarg_flags::NONE, result);
    EXPECT_EQ(15, g_show_dot);
}

TEST_F(TestParameterCommand, showDotAuto)
{
    ValueSaver save_show_dot{g_show_dot, -999L};
    ValueSaver save_auto_show_dot{g_auto_show_dot, '!'};
    ValueSaver save_size_dot{g_size_dot, -99};

    const cmdarg_flags result = exec_cmd_arg("showdot=a", cmd_file::AT_AFTER_STARTUP);

    EXPECT_EQ(cmdarg_flags::NONE, result);
    EXPECT_EQ(15, g_show_dot);
    EXPECT_EQ('a', g_auto_show_dot);
    EXPECT_EQ(0, g_size_dot);
}

TEST_F(TestParameterCommand, showDotAutoSize)
{
    ValueSaver save_show_dot{g_show_dot, -999L};
    ValueSaver save_auto_show_dot{g_auto_show_dot, '!'};
    ValueSaver save_size_dot{g_size_dot, -1};

    const cmdarg_flags result = exec_cmd_arg("showdot=a/20", cmd_file::AT_AFTER_STARTUP);

    EXPECT_EQ(cmdarg_flags::NONE, result);
    EXPECT_EQ(15, g_show_dot);
    EXPECT_EQ('a', g_auto_show_dot);
    EXPECT_EQ(20, g_size_dot);
}

TEST_F(TestParameterCommand, showOrbitYes)
{
    ValueSaver save_show_orbit{g_start_show_orbit, false};

    const cmdarg_flags result = exec_cmd_arg("showorbit=y", cmd_file::AT_AFTER_STARTUP);

    EXPECT_EQ(cmdarg_flags::NONE, result);
    EXPECT_TRUE(g_start_show_orbit);
}

TEST_F(TestParameterCommand, decompValue)
{
    ValueSaver save_decomp0{g_decomp[0], -99};
    ValueSaver save_decomp1{g_decomp[1], -99};
    VALUE_UNCHANGED(g_bail_out, -99L);

    const cmdarg_flags result = exec_cmd_arg("decomp=16", cmd_file::AT_AFTER_STARTUP);

    EXPECT_EQ(cmdarg_flags::FRACTAL_PARAM, result);
    EXPECT_EQ(16, g_decomp[0]);
    EXPECT_EQ(0, g_decomp[1]);
}

TEST_F(TestParameterCommand, decompTwoValues)
{
    ValueSaver save_decomp0{g_decomp[0], -99};
    ValueSaver save_decomp1{g_decomp[1], -99};
    ValueSaver save_bail_out{g_bail_out, -99L};

    const cmdarg_flags result = exec_cmd_arg("decomp=16/4", cmd_file::AT_AFTER_STARTUP);

    EXPECT_EQ(cmdarg_flags::FRACTAL_PARAM, result);
    EXPECT_EQ(16, g_decomp[0]);
    EXPECT_EQ(4, g_decomp[1]);
    EXPECT_EQ(4, g_bail_out);
}

TEST_F(TestParameterCommand, distEstOneValue)
{
    ValueSaver saved_user_distance_estimator_value{g_user_distance_estimator_value, -99L};
    ValueSaver saved_distance_estimator_width_factor{g_distance_estimator_width_factor, -99};
    ValueSaver saved_distance_estimator_x_dots{g_distance_estimator_x_dots, -99};
    ValueSaver saved_distance_estimator_y_dots{g_distance_estimator_y_dots, -99};

    const cmdarg_flags result = exec_cmd_arg("distest=4", cmd_file::AT_AFTER_STARTUP);

    EXPECT_EQ(cmdarg_flags::FRACTAL_PARAM, result);
    EXPECT_EQ(4L, g_user_distance_estimator_value);
    EXPECT_EQ(71, g_distance_estimator_width_factor);
    EXPECT_EQ(0, g_distance_estimator_x_dots);
    EXPECT_EQ(0, g_distance_estimator_y_dots);
}

TEST_F(TestParameterCommand, distEstTwoValues)
{
    ValueSaver saved_user_distance_estimator_value{g_user_distance_estimator_value, -99L};
    ValueSaver saved_distance_estimator_width_factor{g_distance_estimator_width_factor, -99};
    ValueSaver saved_distance_estimator_x_dots{g_distance_estimator_x_dots, -99};
    ValueSaver saved_distance_estimator_y_dots{g_distance_estimator_y_dots, -99};

    const cmdarg_flags result = exec_cmd_arg("distest=4/5", cmd_file::AT_AFTER_STARTUP);

    EXPECT_EQ(cmdarg_flags::FRACTAL_PARAM, result);
    EXPECT_EQ(4L, g_user_distance_estimator_value);
    EXPECT_EQ(5, g_distance_estimator_width_factor);
    EXPECT_EQ(0, g_distance_estimator_x_dots);
    EXPECT_EQ(0, g_distance_estimator_y_dots);
}

TEST_F(TestParameterCommand, distEstFourValues)
{
    ValueSaver saved_user_distance_estimator_value{g_user_distance_estimator_value, -99L};
    ValueSaver saved_distance_estimator_width_factor{g_distance_estimator_width_factor, -99};
    ValueSaver saved_distance_estimator_x_dots{g_distance_estimator_x_dots, -99};
    ValueSaver saved_distance_estimator_y_dots{g_distance_estimator_y_dots, -99};

    const cmdarg_flags result = exec_cmd_arg("distest=4/5/6/7", cmd_file::AT_AFTER_STARTUP);

    EXPECT_EQ(cmdarg_flags::FRACTAL_PARAM, result);
    EXPECT_EQ(4L, g_user_distance_estimator_value);
    EXPECT_EQ(5, g_distance_estimator_width_factor);
    EXPECT_EQ(6, g_distance_estimator_x_dots);
    EXPECT_EQ(7, g_distance_estimator_y_dots);
}

TEST_F(TestParameterCommand, formulaFileFilename)
{
    ValueSaver saved_formula_filename{g_formula_filename, ""};

    const cmdarg_flags result = exec_cmd_arg("formulafile=foo.frm", cmd_file::AT_AFTER_STARTUP);

    EXPECT_EQ(cmdarg_flags::FRACTAL_PARAM, result);
    EXPECT_EQ("foo.frm", g_formula_filename);
}

TEST_F(TestParameterCommand, formulaName)
{
    ValueSaver saved_formula_name{g_formula_name, ""};

    const cmdarg_flags result = exec_cmd_arg("formulaname=Monongahela", cmd_file::AT_AFTER_STARTUP);

    EXPECT_EQ(cmdarg_flags::FRACTAL_PARAM, result);
    EXPECT_EQ("monongahela", g_formula_name);
}

TEST_F(TestParameterCommand, lFileFilename)
{
    ValueSaver saved_l_system_filename{g_l_system_filename, ""};

    const cmdarg_flags result = exec_cmd_arg("lfile=foo.l", cmd_file::AT_AFTER_STARTUP);

    EXPECT_EQ(cmdarg_flags::FRACTAL_PARAM, result);
    EXPECT_EQ("foo.l", g_l_system_filename);
}

TEST_F(TestParameterCommand, lName)
{
    ValueSaver saved_l_system_name{g_l_system_name, ""};

    const cmdarg_flags result = exec_cmd_arg("lname=Monongahela", cmd_file::AT_AFTER_STARTUP);

    EXPECT_EQ(cmdarg_flags::FRACTAL_PARAM, result);
    EXPECT_EQ("monongahela", g_l_system_name);
}

TEST_F(TestParameterCommand, ifsFileFilename)
{
    ValueSaver saved_ifs_filename{g_ifs_filename, ""};

    const cmdarg_flags result = exec_cmd_arg("ifsfile=foo.ifs", cmd_file::AT_AFTER_STARTUP);

    EXPECT_EQ(cmdarg_flags::FRACTAL_PARAM, result);
    EXPECT_EQ("foo.ifs", g_ifs_filename);
}

TEST_F(TestParameterCommand, ifsName)
{
    ValueSaver saved_ifs_name{g_ifs_name, ""};

    const cmdarg_flags result = exec_cmd_arg("ifs=Monongahela", cmd_file::AT_AFTER_STARTUP);

    EXPECT_EQ(cmdarg_flags::FRACTAL_PARAM, result);
    EXPECT_EQ("monongahela", g_ifs_name);
}

TEST_F(TestParameterCommand, ifs3DName)
{
    ValueSaver saved_ifs_name{g_ifs_name, ""};

    const cmdarg_flags result = exec_cmd_arg("ifs3d=Monongahela", cmd_file::AT_AFTER_STARTUP);

    EXPECT_EQ(cmdarg_flags::FRACTAL_PARAM, result);
    EXPECT_EQ("monongahela", g_ifs_name);
}

TEST_F(TestParameterCommand, parmFile)
{
    ValueSaver saved_command_file{g_command_file, ""};

    const cmdarg_flags result = exec_cmd_arg("parmfile=foo.par", cmd_file::AT_AFTER_STARTUP);

    EXPECT_EQ(cmdarg_flags::FRACTAL_PARAM, result);
    EXPECT_EQ("foo.par", g_command_file);
}

TEST_F(TestParameterCommand, stereoValue)
{
    ValueSaver saved_glasses_type{g_glasses_type, -1};

    const cmdarg_flags result = exec_cmd_arg("stereo=3", cmd_file::AT_AFTER_STARTUP);

    EXPECT_EQ(cmdarg_flags::FRACTAL_PARAM | cmdarg_flags::PARAM_3D, result);
    EXPECT_EQ(3, g_glasses_type);
}

TEST_F(TestParameterCommand, rotation)
{
    ValueSaver saved_x_rot{XROT, -99};
    ValueSaver saved_y_rot{YROT, -99};
    ValueSaver saved_z_rot{ZROT, -99};

    const cmdarg_flags result = exec_cmd_arg("rotation=30/60/90", cmd_file::AT_AFTER_STARTUP);

    EXPECT_EQ(cmdarg_flags::FRACTAL_PARAM | cmdarg_flags::PARAM_3D, result);
    EXPECT_EQ(30, XROT);
    EXPECT_EQ(60, YROT);
    EXPECT_EQ(90, ZROT);
}

TEST_F(TestParameterCommand, perspective)
{
    ValueSaver saved_z_viewer{ZVIEWER, -99};

    const cmdarg_flags result = exec_cmd_arg("perspective=90", cmd_file::AT_AFTER_STARTUP);

    EXPECT_EQ(cmdarg_flags::FRACTAL_PARAM | cmdarg_flags::PARAM_3D, result);
    EXPECT_EQ(90, ZVIEWER);
}

TEST_F(TestParameterCommand, xyShift)
{
    ValueSaver saved_x_shift{XSHIFT, -99};
    ValueSaver saved_y_shift{YSHIFT, -99};

    const cmdarg_flags result = exec_cmd_arg("xyshift=30/90", cmd_file::AT_AFTER_STARTUP);

    EXPECT_EQ(cmdarg_flags::FRACTAL_PARAM | cmdarg_flags::PARAM_3D, result);
    EXPECT_EQ(30, XSHIFT);
    EXPECT_EQ(90, YSHIFT);
}

TEST_F(TestParameterCommand, interocular)
{
    ValueSaver saved_eye_separation{g_eye_separation, -99};

    const cmdarg_flags result = exec_cmd_arg("interocular=90", cmd_file::AT_AFTER_STARTUP);

    EXPECT_EQ(cmdarg_flags::FRACTAL_PARAM | cmdarg_flags::PARAM_3D, result);
    EXPECT_EQ(90, g_eye_separation);
}

TEST_F(TestParameterCommand, converge)
{
    ValueSaver saved_converge_x_adjust{g_converge_x_adjust, -99};

    const cmdarg_flags result = exec_cmd_arg("converge=90", cmd_file::AT_AFTER_STARTUP);

    EXPECT_EQ(cmdarg_flags::FRACTAL_PARAM | cmdarg_flags::PARAM_3D, result);
    EXPECT_EQ(90, g_converge_x_adjust);
}

TEST_F(TestParameterCommand, crop)
{
    ValueSaver saved_red_crop_left{g_red_crop_left, -99};
    ValueSaver saved_red_crop_right{g_red_crop_right, -99};
    ValueSaver saved_blue_crop_left{g_blue_crop_left, -99};
    ValueSaver saved_blue_crop_right{g_blue_crop_right, -99};

    const cmdarg_flags result = exec_cmd_arg("crop=1/2/3/4", cmd_file::AT_AFTER_STARTUP);

    EXPECT_EQ(cmdarg_flags::FRACTAL_PARAM | cmdarg_flags::PARAM_3D, result);
    EXPECT_EQ(1, g_red_crop_left);
    EXPECT_EQ(2, g_red_crop_right);
    EXPECT_EQ(3, g_blue_crop_left);
    EXPECT_EQ(4, g_blue_crop_right);
}

TEST_F(TestParameterCommand, bright)
{
    ValueSaver saved_red_bright{g_red_bright, -99};
    ValueSaver saved_blue_bright{g_blue_bright, -99};

    const cmdarg_flags result = exec_cmd_arg("bright=1/2", cmd_file::AT_AFTER_STARTUP);

    EXPECT_EQ(cmdarg_flags::FRACTAL_PARAM | cmdarg_flags::PARAM_3D, result);
    EXPECT_EQ(1, g_red_bright);
    EXPECT_EQ(2, g_blue_bright);
}

TEST_F(TestParameterCommand, xyAdjust)
{
    ValueSaver saved_adjust_3d_x{g_adjust_3d_x, -99};
    ValueSaver saved_adjust_3d_{g_adjust_3d_y, -99};

    const cmdarg_flags result = exec_cmd_arg("xyadjust=1/2", cmd_file::AT_AFTER_STARTUP);

    EXPECT_EQ(cmdarg_flags::FRACTAL_PARAM | cmdarg_flags::PARAM_3D, result);
    EXPECT_EQ(1, g_adjust_3d_x);
    EXPECT_EQ(2, g_adjust_3d_y);
}

TEST_F(TestParameterCommand, threeDNo)
{
    VALUE_UNCHANGED(g_overlay_3d, true);
    ValueSaver saved_display_3d{g_display_3d, display_3d_modes::YES};

    const cmdarg_flags result = exec_cmd_arg("3d=no", cmd_file::AT_AFTER_STARTUP);

    EXPECT_EQ(cmdarg_flags::PARAM_3D, result);
    EXPECT_EQ(display_3d_modes::NONE, g_display_3d);
}

TEST_F(TestParameterCommand, threeDYes)
{
    ValueSaver saved_overlay_3d{g_overlay_3d, false};
    ValueSaver saved_display_3d{g_display_3d, display_3d_modes::MINUS_ONE};

    const cmdarg_flags result = exec_cmd_arg("3d=yes", cmd_file::AT_AFTER_STARTUP);

    EXPECT_EQ(cmdarg_flags::PARAM_3D | cmdarg_flags::YES_3D, result);
    EXPECT_FALSE(g_overlay_3d);
    EXPECT_EQ(display_3d_modes::YES, g_display_3d);
}

TEST_F(TestParameterCommand, threeDOverlayWithFractal)
{
    ValueSaver saved_overlay_3d{g_overlay_3d, false};
    ValueSaver saved_display_3d{g_display_3d, display_3d_modes::MINUS_ONE};
    ValueSaver saved_calc_status{g_calc_status, calc_status_value::PARAMS_CHANGED};

    const cmdarg_flags result = exec_cmd_arg("3d=overlay", cmd_file::AT_AFTER_STARTUP);

    EXPECT_EQ(cmdarg_flags::PARAM_3D | cmdarg_flags::YES_3D, result);
    EXPECT_TRUE(g_overlay_3d);
    EXPECT_EQ(display_3d_modes::YES, g_display_3d);
}

TEST_F(TestParameterCommand, threeDOverlayNoFractal)
{
    ValueSaver saved_overlay_3d{g_overlay_3d, false};
    ValueSaver saved_display_3d{g_display_3d, display_3d_modes::MINUS_ONE};
    ValueSaver saved_calc_status{g_calc_status, calc_status_value::NO_FRACTAL};

    const cmdarg_flags result = exec_cmd_arg("3d=overlay", cmd_file::AT_AFTER_STARTUP);

    EXPECT_EQ(cmdarg_flags::PARAM_3D | cmdarg_flags::YES_3D, result);
    EXPECT_FALSE(g_overlay_3d);
    EXPECT_EQ(display_3d_modes::YES, g_display_3d);
}

TEST_F(TestParameterCommand, sphereNo)
{
    ValueSaver saved_SPHERE{SPHERE, -99};

    const cmdarg_flags result = exec_cmd_arg("sphere=no", cmd_file::AT_AFTER_STARTUP);

    EXPECT_EQ(cmdarg_flags::PARAM_3D, result);
    EXPECT_EQ(0, SPHERE);
}

TEST_F(TestParameterCommand, sphereYes)
{
    ValueSaver saved_SPHERE{SPHERE, -99};

    const cmdarg_flags result = exec_cmd_arg("sphere=yes", cmd_file::AT_AFTER_STARTUP);

    EXPECT_EQ(cmdarg_flags::PARAM_3D, result);
    EXPECT_EQ(1, SPHERE);
}

TEST_F(TestParameterCommand, scaleXYZTwoValues)
{
    ValueSaver saved_XSCALE{XSCALE, -99};
    ValueSaver saved_YSCALE{YSCALE, -99};
    VALUE_UNCHANGED(ROUGH, -99);

    const cmdarg_flags result = exec_cmd_arg("scalexyz=1/2", cmd_file::AT_AFTER_STARTUP);

    EXPECT_EQ(cmdarg_flags::PARAM_3D, result);
    EXPECT_EQ(1, XSCALE);
    EXPECT_EQ(2, YSCALE);
}

TEST_F(TestParameterCommand, scaleXYZThreeValues)
{
    ValueSaver saved_XSCALE{XSCALE, -99};
    ValueSaver saved_YSCALE{YSCALE, -99};
    ValueSaver saved_ROUGH{ROUGH, -99};

    const cmdarg_flags result = exec_cmd_arg("scalexyz=1/2/3", cmd_file::AT_AFTER_STARTUP);

    EXPECT_EQ(cmdarg_flags::PARAM_3D, result);
    EXPECT_EQ(1, XSCALE);
    EXPECT_EQ(2, YSCALE);
    EXPECT_EQ(3, ROUGH);
}

TEST_F(TestParameterCommand, roughness)
{
    ValueSaver saved_ROUGH{ROUGH, -99};

    const cmdarg_flags result = exec_cmd_arg("roughness=3", cmd_file::AT_AFTER_STARTUP);

    EXPECT_EQ(cmdarg_flags::PARAM_3D, result);
    EXPECT_EQ(3, ROUGH);
}

TEST_F(TestParameterCommand, waterline)
{
    ValueSaver saved_WATERLINE{WATERLINE, -99};

    const cmdarg_flags result = exec_cmd_arg("waterline=3", cmd_file::AT_AFTER_STARTUP);

    EXPECT_EQ(cmdarg_flags::PARAM_3D, result);
    EXPECT_EQ(3, WATERLINE);
}

TEST_F(TestParameterCommand, fillType)
{
    ValueSaver saved_FILLTYPE{FILLTYPE, -99};

    const cmdarg_flags result = exec_cmd_arg("filltype=3", cmd_file::AT_AFTER_STARTUP);

    EXPECT_EQ(cmdarg_flags::PARAM_3D, result);
    EXPECT_EQ(3, FILLTYPE);
}

TEST_F(TestParameterCommand, lightSource)
{
    ValueSaver saved_XLIGHT{XLIGHT, -99};
    ValueSaver saved_YLIGHT{YLIGHT, -99};
    ValueSaver saved_ZLIGHT{ZLIGHT, -99};

    const cmdarg_flags result = exec_cmd_arg("lightsource=1/2/3", cmd_file::AT_AFTER_STARTUP);

    EXPECT_EQ(cmdarg_flags::PARAM_3D, result);
    EXPECT_EQ(1, XLIGHT);
    EXPECT_EQ(2, YLIGHT);
    EXPECT_EQ(3, ZLIGHT);
}

TEST_F(TestParameterCommand, smoothing)
{
    ValueSaver saved_LIGHTAVG{LIGHTAVG, -99};

    const cmdarg_flags result = exec_cmd_arg("smoothing=3", cmd_file::AT_AFTER_STARTUP);

    EXPECT_EQ(cmdarg_flags::PARAM_3D, result);
    EXPECT_EQ(3, LIGHTAVG);
}

TEST_F(TestParameterCommand, latitude)
{
    ValueSaver saved_THETA1{THETA1, -99};
    ValueSaver saved_THETA2{THETA2, -99};

    const cmdarg_flags result = exec_cmd_arg("latitude=1/3", cmd_file::AT_AFTER_STARTUP);

    EXPECT_EQ(cmdarg_flags::PARAM_3D, result);
    EXPECT_EQ(1, THETA1);
    EXPECT_EQ(3, THETA2);
}

TEST_F(TestParameterCommand, longitude)
{
    ValueSaver saved_PHI1{PHI1, -99};
    ValueSaver saved_PHI2{PHI2, -99};

    const cmdarg_flags result = exec_cmd_arg("longitude=1/3", cmd_file::AT_AFTER_STARTUP);

    EXPECT_EQ(cmdarg_flags::PARAM_3D, result);
    EXPECT_EQ(1, PHI1);
    EXPECT_EQ(3, PHI2);
}

TEST_F(TestParameterCommand, radius)
{
    ValueSaver saved_RADIUS{RADIUS, -99};

    const cmdarg_flags result = exec_cmd_arg("radius=1", cmd_file::AT_AFTER_STARTUP);

    EXPECT_EQ(cmdarg_flags::PARAM_3D, result);
    EXPECT_EQ(1, RADIUS);
}

TEST_F(TestParameterCommand, transparentOneValue)
{
    ValueSaver saved_transparent_color_3d0{g_transparent_color_3d[0], -99};
    ValueSaver saved_transparent_color_3d1{g_transparent_color_3d[1], -99};

    const cmdarg_flags result = exec_cmd_arg("transparent=1", cmd_file::AT_AFTER_STARTUP);

    EXPECT_EQ(cmdarg_flags::PARAM_3D, result);
    EXPECT_EQ(1, g_transparent_color_3d[0]);
    EXPECT_EQ(1, g_transparent_color_3d[1]);
}

TEST_F(TestParameterCommand, transparentTwoValues)
{
    ValueSaver saved_transparent_color_3d0{g_transparent_color_3d[0], -99};
    ValueSaver saved_transparent_color_3d1{g_transparent_color_3d[1], -99};

    const cmdarg_flags result = exec_cmd_arg("transparent=1/2", cmd_file::AT_AFTER_STARTUP);

    EXPECT_EQ(cmdarg_flags::PARAM_3D, result);
    EXPECT_EQ(1, g_transparent_color_3d[0]);
    EXPECT_EQ(2, g_transparent_color_3d[1]);
}

TEST_F(TestParameterCommand, previewNo)
{
    ValueSaver saved_preview{g_preview, true};

    const cmdarg_flags result = exec_cmd_arg("preview=no", cmd_file::AT_AFTER_STARTUP);

    EXPECT_EQ(cmdarg_flags::PARAM_3D, result);
    EXPECT_FALSE(g_preview);
}

TEST_F(TestParameterCommand, showBoxNo)
{
    ValueSaver saved_show_box{g_show_box, true};

    const cmdarg_flags result = exec_cmd_arg("showbox=no", cmd_file::AT_AFTER_STARTUP);

    EXPECT_EQ(cmdarg_flags::PARAM_3D, result);
    EXPECT_FALSE(g_show_box);
}

TEST_F(TestParameterCommand, coarse)
{
    ValueSaver saved_preview_factor{g_preview_factor, -99};

    const cmdarg_flags result = exec_cmd_arg("coarse=3", cmd_file::AT_AFTER_STARTUP);

    EXPECT_EQ(cmdarg_flags::PARAM_3D, result);
    EXPECT_EQ(3, g_preview_factor);
}

TEST_F(TestParameterCommand, randomize)
{
    ValueSaver saved_randomize_3d{g_randomize_3d, -99};

    const cmdarg_flags result = exec_cmd_arg("randomize=3", cmd_file::AT_AFTER_STARTUP);

    EXPECT_EQ(cmdarg_flags::PARAM_3D, result);
    EXPECT_EQ(3, g_randomize_3d);
}

TEST_F(TestParameterCommand, ambient)
{
    ValueSaver saved_ambient{g_ambient, -99};

    const cmdarg_flags result = exec_cmd_arg("ambient=3", cmd_file::AT_AFTER_STARTUP);

    EXPECT_EQ(cmdarg_flags::PARAM_3D, result);
    EXPECT_EQ(3, g_ambient);
}

TEST_F(TestParameterCommand, haze)
{
    ValueSaver saved_haze{g_haze, -99};

    const cmdarg_flags result = exec_cmd_arg("haze=3", cmd_file::AT_AFTER_STARTUP);

    EXPECT_EQ(cmdarg_flags::PARAM_3D, result);
    EXPECT_EQ(3, g_haze);
}

TEST_F(TestParameterCommand, fullColorNo)
{
    ValueSaver saved_targa_out{g_targa_out, true};

    const cmdarg_flags result = exec_cmd_arg("fullcolor=no", cmd_file::AT_AFTER_STARTUP);

    EXPECT_EQ(cmdarg_flags::PARAM_3D, result);
    EXPECT_FALSE(g_targa_out);
}

TEST_F(TestParameterCommand, trueColorNo)
{
    ValueSaver saved_truecolor{g_truecolor, true};

    const cmdarg_flags result = exec_cmd_arg("truecolor=no", cmd_file::AT_AFTER_STARTUP);

    EXPECT_EQ(cmdarg_flags::FRACTAL_PARAM | cmdarg_flags::PARAM_3D, result);
    EXPECT_FALSE(g_truecolor);
}

TEST_F(TestParameterCommand, trueModeDefaultColor)
{
    ValueSaver saved_true_mode{g_true_mode, true_color_mode::iterate};

    const cmdarg_flags result = exec_cmd_arg("truemode=d", cmd_file::AT_AFTER_STARTUP);

    EXPECT_EQ(cmdarg_flags::FRACTAL_PARAM | cmdarg_flags::PARAM_3D, result);
    EXPECT_EQ(true_color_mode::default_color, g_true_mode);
}

TEST_F(TestParameterCommand, trueModeIterate)
{
    ValueSaver saved_true_mode{g_true_mode, true_color_mode::default_color};

    const cmdarg_flags result = exec_cmd_arg("truemode=i", cmd_file::AT_AFTER_STARTUP);

    EXPECT_EQ(cmdarg_flags::FRACTAL_PARAM | cmdarg_flags::PARAM_3D, result);
    EXPECT_EQ(true_color_mode::iterate, g_true_mode);
}

TEST_F(TestParameterCommand, trueModeValue)
{
    ValueSaver saved_true_mode{g_true_mode, true_color_mode::default_color};

    const cmdarg_flags result = exec_cmd_arg("truemode=1", cmd_file::AT_AFTER_STARTUP);

    EXPECT_EQ(cmdarg_flags::FRACTAL_PARAM | cmdarg_flags::PARAM_3D, result);
    EXPECT_EQ(true_color_mode::iterate, g_true_mode);
}

TEST_F(TestParameterCommand, useGrayScaleNo)
{
    ValueSaver saved_gray_flag{g_gray_flag, true};

    const cmdarg_flags result = exec_cmd_arg("usegrayscale=n", cmd_file::AT_AFTER_STARTUP);

    EXPECT_EQ(cmdarg_flags::PARAM_3D, result);
    EXPECT_FALSE(g_gray_flag);
}

TEST_F(TestParameterCommand, stereoWidth)
{
    ValueSaver saved_auto_stereo_width{g_auto_stereo_width, true};

    const cmdarg_flags result = exec_cmd_arg("stereowidth=4.5", cmd_file::AT_AFTER_STARTUP);

    EXPECT_EQ(cmdarg_flags::PARAM_3D, result);
    EXPECT_EQ(4.5, g_auto_stereo_width);
}

TEST_F(TestParameterCommand, monitorWidthAliasForStereoWidth)
{
    ValueSaver saved_auto_stereo_width{g_auto_stereo_width, true};

    const cmdarg_flags result = exec_cmd_arg("monitorwidth=4.5", cmd_file::AT_AFTER_STARTUP);

    EXPECT_EQ(cmdarg_flags::PARAM_3D, result);
    EXPECT_EQ(4.5, g_auto_stereo_width);
}

TEST_F(TestParameterCommand, targaOverlayNo)
{
    ValueSaver saved_targa_overlay{g_targa_overlay, true};

    const cmdarg_flags result = exec_cmd_arg("targa_overlay=n", cmd_file::AT_AFTER_STARTUP);

    EXPECT_EQ(cmdarg_flags::PARAM_3D, result);
    EXPECT_FALSE(g_targa_overlay);
}

TEST_F(TestParameterCommand, background)
{
    ValueSaver saved_background_color0{g_background_color[0], 255};
    ValueSaver saved_background_color1{g_background_color[1], 255};
    ValueSaver saved_background_color2{g_background_color[2], 255};

    const cmdarg_flags result = exec_cmd_arg("background=1/2/3", cmd_file::AT_AFTER_STARTUP);

    EXPECT_EQ(cmdarg_flags::PARAM_3D, result);
    EXPECT_EQ(1, g_background_color[0]);
    EXPECT_EQ(2, g_background_color[1]);
    EXPECT_EQ(3, g_background_color[2]);
}

TEST_F(TestParameterCommand, lightNameFirstInit)
{
    ValueSaver saved_first_init{g_first_init, true};
    ValueSaver saved_light_name{g_light_name, "fmeh"};

    const cmdarg_flags result = exec_cmd_arg("lightname=foo", cmd_file::AT_CMD_LINE);

    EXPECT_EQ(cmdarg_flags::NONE, result);
    EXPECT_EQ("foo", g_light_name);
}

TEST_F(TestParameterCommand, lightNameAfterStartup)
{
    ValueSaver saved_first_init{g_first_init, false};
    ValueSaver saved_light_name{g_light_name, "fmeh"};

    const cmdarg_flags result = exec_cmd_arg("lightname=foo", cmd_file::AT_AFTER_STARTUP);

    EXPECT_EQ(cmdarg_flags::NONE, result);
    EXPECT_EQ("foo", g_light_name);
}

TEST_F(TestParameterCommand, lightNameNotSet)
{
    ValueSaver saved_first_init{g_first_init, false};
    VALUE_UNCHANGED(g_light_name, std::string{"fmeh"});

    const cmdarg_flags result = exec_cmd_arg("lightname=foo", cmd_file::AT_CMD_LINE);

    EXPECT_EQ(cmdarg_flags::NONE, result);
}

TEST_F(TestParameterCommand, ray)
{
    ValueSaver saved_raytrace_format{g_raytrace_format, raytrace_formats::none};

    const cmdarg_flags result = exec_cmd_arg("ray=3", cmd_file::AT_CMD_LINE);

    EXPECT_EQ(cmdarg_flags::PARAM_3D, result);
    EXPECT_EQ(raytrace_formats::raw, g_raytrace_format);
}

TEST_F(TestParameterCommand, briefNo)
{
    ValueSaver saved_brief{g_brief, true};

    const cmdarg_flags result = exec_cmd_arg("brief=n", cmd_file::AT_AFTER_STARTUP);

    EXPECT_EQ(cmdarg_flags::PARAM_3D, result);
    EXPECT_FALSE(g_brief);
}

TEST_F(TestParameterCommand, releaseNotAllowed)
{
    expect_stop_msg();

    const cmdarg_flags result = exec_cmd_arg("release=100", cmd_file::AT_AFTER_STARTUP);

    EXPECT_EQ(cmdarg_flags::ERROR, result);
}

TEST_F(TestParameterCommand, curDirNo)
{
    ValueSaver saved_check_cur_dir{g_check_cur_dir, true};

    const cmdarg_flags result = exec_cmd_arg("curdir=n", cmd_file::AT_AFTER_STARTUP);

    EXPECT_EQ(cmdarg_flags::NONE, result);
    EXPECT_FALSE(g_check_cur_dir);
}

TEST_F(TestParameterCommand, virtualNo)
{
    ValueSaver saved_virtual_screens{g_virtual_screens, true};

    const cmdarg_flags result = exec_cmd_arg("virtual=n", cmd_file::AT_AFTER_STARTUP);

    EXPECT_EQ(cmdarg_flags::FRACTAL_PARAM, result);
    EXPECT_FALSE(g_virtual_screens);
}
