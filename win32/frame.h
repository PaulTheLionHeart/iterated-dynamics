#pragma once

#include <array>
#include <string>

#define KEYBUFMAX 80

class Frame
{
public:
    void init(HINSTANCE instance, LPCSTR title);
    void terminate();
    void create_window(int width, int height);
    int get_key_press(bool wait_for_key);
    int pump_messages(bool waitflag);
    void resize(int width, int height);
    void set_keyboard_timeout(int ms);

    void adjust_size(int width, int height);

    HINSTANCE m_instance{};
    HWND m_window{};
    std::string m_title;
    int m_width{};
    int m_height{};
    int m_nc_width{};
    int m_nc_height{};
    HWND m_child{};
    bool m_has_focus{};
    bool m_timed_out{};

    // the keypress buffer
    unsigned int m_key_press_count{};
    unsigned int m_key_press_head{};
    unsigned int m_key_press_tail{};
    std::array<unsigned int, KEYBUFMAX> m_key_press_buffer{};
};

extern Frame g_frame;
