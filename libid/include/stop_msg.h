// SPDX-License-Identifier: GPL-3.0-only
//
#pragma once

#include <string>

// stop_msg() flags
enum class stopmsg_flags
{
    NONE        = 0,
    NO_STACK    = 1,
    CANCEL      = 2,
    NO_BUZZER   = 4,
    FIXED_FONT  = 8,
    INFO_ONLY   = 16
};

inline int operator+(stopmsg_flags value)
{
    return static_cast<int>(value);
}
inline stopmsg_flags operator|(stopmsg_flags lhs, stopmsg_flags rhs)
{
    return static_cast<stopmsg_flags>(+lhs | +rhs);
}
inline bool bit_set(stopmsg_flags flags, stopmsg_flags bit)
{
    return (+flags & +bit) == +bit;
}

bool stop_msg(stopmsg_flags flags, const std::string &msg);

// the most common case
inline bool stop_msg(const std::string &msg)
{
    return stop_msg(stopmsg_flags::NONE, msg);
}
