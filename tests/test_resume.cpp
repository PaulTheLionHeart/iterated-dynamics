// SPDX-License-Identifier: GPL-3.0-only
//
#include <resume.h>

#include <id_data.h>

#include <gtest/gtest.h>

#include <cstdint>

namespace
{

#pragma pack(push, 1)
struct ResumeData
{
    std::int32_t datum1;
    std::int8_t datum2;
    std::uint16_t datum3;
};
#pragma pack(pop)

inline bool operator==(const ResumeData &lhs, const ResumeData &rhs)
{
    return lhs.datum1 == rhs.datum1 && lhs.datum2 == rhs.datum2 && lhs.datum3 == rhs.datum3;
}
inline bool operator!=(const ResumeData &lhs, const ResumeData &rhs)
{
    return !(lhs == rhs);
}

struct TestResume : testing::Test
{
    ~TestResume() override = default;

protected:
    void SetUp() override;
    void TearDown() override;

    ResumeData m_data{};
};

void TestResume::SetUp()
{
    Test::SetUp();
    m_data.datum1 = 123;
    m_data.datum2 = 64;
    m_data.datum3 = 65536;
    alloc_resume(sizeof(ResumeData) + 20, 3);
}

void TestResume::TearDown()
{
    end_resume();
    Test::TearDown();
}

} // namespace

TEST_F(TestResume, allocate)
{
    EXPECT_LE(sizeof(ResumeData) + 20, g_resume_data.size());
    EXPECT_EQ(sizeof(int), g_resume_len);
    EXPECT_EQ(CalcStatus::RESUMABLE, g_calc_status);
}

TEST_F(TestResume, putCountsStoredData)
{
    put_resume(sizeof(ResumeData), &m_data, 0);

    EXPECT_EQ(sizeof(int) + sizeof(ResumeData), g_resume_len);
}

TEST_F(TestResume, startResumeReturnsVersion)
{
    EXPECT_EQ(3, start_resume());
}

TEST_F(TestResume, templatePutSingle)
{
    put_resume_t(m_data);

    EXPECT_EQ(sizeof(int) + sizeof(ResumeData), g_resume_len);
}

TEST_F(TestResume, templatePutMultiple)
{
    put_resume_t(m_data.datum1, m_data.datum2, m_data.datum3);

    EXPECT_EQ(sizeof(int) + sizeof(ResumeData), g_resume_len);
}

TEST_F(TestResume, getReturnsPutData)
{
    put_resume(sizeof(ResumeData), &m_data, 0);
    start_resume();

    ResumeData actual;
    get_resume(sizeof(ResumeData), &actual, 0);

    EXPECT_EQ(m_data, actual);
}

TEST_F(TestResume, templateGetSingle)
{
    put_resume_t(m_data);
    start_resume();

    ResumeData actual;
    get_resume_t(actual);

    EXPECT_EQ(m_data, actual);
}

TEST_F(TestResume, templateGetMultiple)
{
    put_resume_t(m_data);
    start_resume();

    ResumeData actual;
    get_resume_t(actual.datum1, actual.datum2, actual.datum3);

    EXPECT_EQ(m_data, actual);
}
