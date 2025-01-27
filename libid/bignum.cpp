// C routines for bignumbers

/*
Wesley Loewer's Big Numbers.        (C) 1994-95, Wesley B. Loewer

The bignumber format is simply a signed integer of variable length.  The
bytes are stored in reverse order (least significant byte first, most
significant byte last).  The sign bit is the highest bit of the most
significant byte.  Negatives are stored in 2's complement form.  The byte
length of the bignumbers must be a multiple of 4 for 386+ operations, and
a multiple of 2 for 8086/286 and non 80x86 machines.

Some of the arithmetic operations coded here may alter some of the
operands used.  Therefore, take note of the SIDE-EFFECTS listed with each
procedure.  If the value of an operand needs to be retained, just use
copy_bn() first.  This was done for speed sake to avoid unnecessary
copying.  If space is at such a premium that copying it would be
difficult, some of the operations only change the sign of the value.  In
this case, the original could be optained by calling neg_a_bn().

Most of the bignumber routines operate on true integers.  Some of the
procedures, however, are designed specifically for fixed decimal point
operations.  The results and how the results are interpreted depend on
where the implied decimal point is located.  The routines that depend on
where the decimal is located are:  strtobn(), bntostr(), bntoint(), inttobn(),
bntofloat(), floattobn(), inv_bn(), div_bn().  The number of bytes
designated for the integer part may be 1, 2, or 4.


BIGNUMBER FORMAT:
The following is a discription of the bignumber format and associated
variables.  The number is stored in reverse order (Least Significant Byte,
LSB, stored first in memory, Most Significant Byte, MSB, stored last).
Each '_' below represents a block of memory used for arithmetic (1 block =
4 bytes on 386+, 2 bytes on 286-).  All lengths are given in bytes.

LSB                                MSB
  _  _  _  _  _  _  _  _  _  _  _  _
n <---------- g_bn_length --------->
                 g_int_length  ---> <---

  g_bn_length  = the length in bytes of the bignumber
  g_int_length = the number of bytes used to represent the integer part of
              the bignumber.  Possible values are 1, 2, or 4.  This
              determines the largest number that can be represented by
              the bignumber.
                g_int_length = 1, max value = 127.99...
                g_int_length = 2, max value = 32,767.99...
                g_int_length = 4, max value = 2,147,483,647.99...


FULL DOUBLE PRECISION MULTIPLICATION:
( full_mult_bn(), full_square_bn() )

The product of two bignumbers, n1 and n2, will be a result, r, which is
a double wide bignumber.  The integer part will also be twice as wide,
thereby eliminating the possiblity of overflowing the number.

LSB                                                                    MSB
  _  _  _  _  _  _  _  _  _  _  _  _  _  _  _  _  _  _  _  _  _  _  _  _
r <------------------------- 2*g_bn_length ---------------------------->
                                                  2*g_int_length  --->  <---

If this double wide bignumber, r, needs to be converted to a normal,
single width bignumber, this is easily done with pointer arithmetic.  The
converted value starts at r+g_shift_factor (where g_shift_factor =
g_bn_length-g_int_length) and continues for g_bn_length bytes.  The lower order
bytes and the upper integer part of the double wide number can then be
ignored.

LSB                                                                    MSB
  _  _  _  _  _  _  _  _  _  _  _  _  _  _  _  _  _  _  _  _  _  _  _  _
r <------------------------- 2*g_bn_length ---------------------------->
                                                  2*g_int_length  --->  <---
                                 LSB                                  MSB
                r+g_shift_factor   <--------  g_bn_length  ----------->
                                                    g_int_length  ---> <---


PARTIAL PRECISION MULTIPLICATION:
( mult_bn(), square_bn() )

In most cases, full double precision multiplication is not necessary.  The
lower order bytes are usually thrown away anyway.  The non-"full"
multiplication routines only calculate g_r_length bytes in the result.  The
value of g_r_length must be in the range: 2*g_bn_length <= g_r_length < g_bn_length.
The amount by which g_r_length exceeds g_bn_length accounts for the extra bytes
that must be multiplied so that the first g_bn_length bytes are correct.
These extra bytes are referred to in the code as the "padding," that is:
g_r_length=g_bn_length+g_padding.

All three of the values, g_bn_length, g_r_length, and therefore g_padding, must be
multiples of the size of memory blocks being used for arithmetic (2 on
8086/286 and 4 on 386+).  Typically, the g_padding is 2*blocksize.  In the
case where g_bn_length=blocksize, g_padding can only be blocksize to keep
g_r_length from being too big.

The product of two bignumbers, n1 and n2, will then be a result, r, which
is of length g_r_length.  The integer part will be twice as wide, thereby
eliminating the possiblity of overflowing the number.

         LSB                                            MSB
           _  _  _  _  _  _  _  _  _  _  _  _  _  _  _  _
        r  <---- g_r_length = g_bn_length+g_padding ---->
                                   2*g_int_length  --->  <---

If r needs to be converted to a normal, single width bignumber, this is
easily done with pointer arithmetic.  The converted value starts at
r+g_shift_factor (where g_shift_factor = g_padding-g_int_length) and continues for
g_bn_length bytes.  The lower order bytes and the upper integer part of the
double wide number can then be ignored.

         LSB                                            MSB
           _  _  _  _  _  _  _  _  _  _  _  _  _  _  _  _
        r  <---- g_r_length = g_bn_length+g_padding ---->
                                   2*g_int_length  --->  <---
                   LSB                                MSB
   r+g_shift_factor  <--------  g_bn_length  --------->
                                     g_int_length ---> <---
*/

/************************************************************************/
// There are three parts to the bignum library:
//
// 1) bignum.c - initialization, general routines, routines that would
//    not be speeded up much with assembler.
//
// 2) bignuma.asm - hand coded assembler routines.
//
// 3) bignumc.c - portable C versions of routines in bignuma.asm
//
/************************************************************************/
#include "big.h"

#include "port.h"

#include <algorithm>
#include <cfloat>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <string>

/*************************************************************************
* The original bignumber code was written specifically for a Little Endian
* system (80x86).  The following is not particularly efficient, but was
* simple to incorporate.  If speed with a Big Endian machine is critical,
* the bignumber format could be reversed.
**************************************************************************/
#ifdef ACCESS_BY_BYTE
U32 big_access32(BYTE *addr)
{
    return addr[0] | ((U32)addr[1] << 8) | ((U32)addr[2] << 16) | ((U32)addr[3] << 24);
}

U16 big_access16(BYTE *addr)
{
    return (U16)addr[0] | ((U16)addr[1] << 8);
}

S16 big_accessS16(S16 *addr)
{
    return (S16)((BYTE *)addr)[0] | ((S16)((BYTE *)addr)[1] << 8);
}

U32 big_set32(BYTE *addr, U32 val)
{
    addr[0] = (BYTE)(val&0xff);
    addr[1] = (BYTE)((val >> 8)&0xff);
    addr[2] = (BYTE)((val >> 16)&0xff);
    addr[3] = (BYTE)((val >> 24)&0xff);
    return val;
}

U16 big_set16(BYTE *addr, U16 val)
{
    addr[0] = (BYTE)(val&0xff);
    addr[1] = (BYTE)((val >> 8)&0xff);
    return val;
}

S16 big_setS16(S16 *addr, S16 val)
{
    ((BYTE *)addr)[0] = (BYTE)(val&0xff);
    ((BYTE *)addr)[1] = (BYTE)((val >> 8)&0xff);
    return val;
}

#endif

/************************************************************************/
// convert_bn  -- convert bignum numbers from old to new lengths
int convert_bn(bn_t newnum, bn_t old, int newbnlength, int newintlength,
               int oldbnlength, int oldintlength)
{
    int savebnlength;
    int saveintlength;

    // save lengths so not dependent on external environment
    saveintlength = g_int_length;
    savebnlength  = g_bn_length;

    g_int_length     = newintlength;
    g_bn_length      = newbnlength;
    clear_bn(newnum);

    if (newbnlength - newintlength > oldbnlength - oldintlength)
    {

        // This will keep the integer part from overflowing past the array.
        g_bn_length = oldbnlength - oldintlength + std::min(oldintlength, newintlength);

        std::memcpy(newnum+newbnlength-newintlength-oldbnlength+oldintlength,
               old, g_bn_length);
    }
    else
    {
        g_bn_length = newbnlength - newintlength + std::min(oldintlength, newintlength);
        std::memcpy(newnum, old+oldbnlength-oldintlength-newbnlength+newintlength,
               g_bn_length);
    }
    g_int_length = saveintlength;
    g_bn_length  = savebnlength;
    return 0;
}

/********************************************************************/
// bn_hexdump() - for debugging, dumps to stdout

void bn_hexdump(bn_t r)
{
    for (int i = 0; i < g_bn_length; i++)
    {
        std::printf("%02X ", r[i]);
    }
    std::printf("\n");
}

/**********************************************************************/
// strtobn() - converts a string into a bignumer
//   r - pointer to a bignumber
//   s - string in the floating point format [-][digits].[digits]
//   note: the string may not be empty or have extra space and may
//   not use scientific notation (2.3e4).

bn_t strtobn(bn_t r, char *s)
{
    bn_t onesbyte;
    bool signflag = false;
    long longval;

    clear_bn(r);
    onesbyte = r + g_bn_length - g_int_length;

    if (s[0] == '+')    // for + sign
    {
        s++;
    }
    else if (s[0] == '-')    // for neg sign
    {
        signflag = true;
        s++;
    }

    if (std::strchr(s, '.') != nullptr) // is there a decimal point?
    {
        int l = (int) std::strlen(s) - 1;      // start with the last digit
        while (s[l] >= '0' && s[l] <= '9') // while a digit
        {
            *onesbyte = (BYTE)(s[l--] - '0');
            div_a_bn_int(r, 10);
        }

        if (s[l] == '.')
        {
            longval = std::atol(s);
            switch (g_int_length)
            {
                // only 1, 2, or 4 are allowed
            case 1:
                *onesbyte = (BYTE)longval;
                break;
            case 2:
                big_set16(onesbyte, (U16)longval);
                break;
            case 4:
                big_set32(onesbyte, longval);
                break;
            }
        }
    }
    else
    {
        longval = std::atol(s);
        switch (g_int_length)
        {
            // only 1, 2, or 4 are allowed
        case 1:
            *onesbyte = (BYTE)longval;
            break;
        case 2:
            big_set16(onesbyte, (U16)longval);
            break;
        case 4:
            big_set32(onesbyte, longval);
            break;
        }
    }

    if (signflag)
    {
        neg_a_bn(r);
    }

    return r;
}

/********************************************************************/
// std::strlen_needed() - returns string length needed to hold bignumber

int strlen_needed()
{
    int length = 3;

    // first space for integer part
    switch (g_int_length)
    {
    case 1:
        length = 3;  // max 127
        break;
    case 2:
        length = 5;  // max 32767
        break;
    case 4:
        length = 10; // max 2147483647
        break;
    }
    length += g_decimals;  // decimal part
    length += 2;         // decimal point and sign
    length += 4;         // null and a little extra for safety
    return length;
}

/********************************************************************/
// bntostr() - converts a bignumber into a string
//   s - string, must be large enough to hold the number.
//   r - bignumber
//   will covert to a floating point notation
//   SIDE-EFFECT: the bignumber, r, is destroyed.
//                Copy it first if necessary.

char *unsafe_bntostr(char *s, int dec, bn_t r)
{
    int l = 0;
    bn_t onesbyte;
    long longval = 0;

    if (dec == 0)
    {
        dec = g_decimals;
    }
    onesbyte = r + g_bn_length - g_int_length;

    if (is_bn_neg(r))
    {
        neg_a_bn(r);
        *(s++) = '-';
    }
    switch (g_int_length)
    {
        // only 1, 2, or 4 are allowed
    case 1:
        longval = *onesbyte;
        break;
    case 2:
        longval = big_access16(onesbyte);
        break;
    case 4:
        longval = big_access32(onesbyte);
        break;
    }
    strcpy(s, std::to_string(longval).c_str());
    l = (int) std::strlen(s);
    s[l++] = '.';
    for (int d = 0; d < dec; d++)
    {
        *onesbyte = 0;  // clear out highest byte
        mult_a_bn_int(r, 10);
        if (is_bn_zero(r))
        {
            break;
        }
        s[l++] = (BYTE)(*onesbyte + '0');
    }
    s[l] = '\0'; // don't forget nul char

    return s;
}

/*********************************************************************/
//  b = l
//  Converts a long to a bignumber
bn_t inttobn(bn_t r, long longval)
{
    bn_t onesbyte;

    clear_bn(r);
    onesbyte = r + g_bn_length - g_int_length;
    switch (g_int_length)
    {
        // only 1, 2, or 4 are allowed
    case 1:
        *onesbyte = (BYTE)longval;
        break;
    case 2:
        big_set16(onesbyte, (U16)longval);
        break;
    case 4:
        big_set32(onesbyte, longval);
        break;
    }
    return r;
}

/*********************************************************************/
//  l = floor(b), floor rounds down
//  Converts the integer part a bignumber to a long
long bntoint(bn_t n)
{
    bn_t onesbyte;
    long longval = 0;

    onesbyte = n + g_bn_length - g_int_length;
    switch (g_int_length)
    {
        // only 1, 2, or 4 are allowed
    case 1:
        longval = *onesbyte;
        break;
    case 2:
        longval = big_access16(onesbyte);
        break;
    case 4:
        longval = big_access32(onesbyte);
        break;
    }
    return longval;
}

/*********************************************************************/
//  b = f
//  Converts a double to a bignumber
bn_t floattobn(bn_t r, LDBL f)
{
    bn_t onesbyte;
    bool signflag = false;

    clear_bn(r);
    onesbyte = r + g_bn_length - g_int_length;

    if (f < 0)
    {
        signflag = true;
        f = -f;
    }

    switch (g_int_length)
    {
        // only 1, 2, or 4 are allowed
    case 1:
        *onesbyte = (BYTE)f;
        break;
    case 2:
        big_set16(onesbyte, (U16)f);
        break;
    case 4:
        big_set32(onesbyte, (U32)f);
        break;
    }

    f -= (long)f; // keep only the decimal part
    for (int i = g_bn_length-g_int_length-1; i >= 0 && f != 0.0; i--)
    {
        f *= 256;
        r[i] = (BYTE)f;  // keep use the integer part
        f -= (BYTE)f; // now throw away the integer part
    }

    if (signflag)
    {
        neg_a_bn(r);
    }

    return r;
}

/********************************************************************/
// sign(r)
int sign_bn(bn_t n)
{
    return is_bn_neg(n) ? -1 : is_bn_not_zero(n) ? 1 : 0;
}

/********************************************************************/
// r = |n|
bn_t abs_bn(bn_t r, bn_t n)
{
    copy_bn(r, n);
    if (is_bn_neg(r))
    {
        neg_a_bn(r);
    }
    return r;
}

/********************************************************************/
// r = |r|
bn_t abs_a_bn(bn_t r)
{
    if (is_bn_neg(r))
    {
        neg_a_bn(r);
    }
    return r;
}

/********************************************************************/
// r = 1/n
// uses g_bn_tmp1 - g_bn_tmp3 - global temp bignumbers
//  SIDE-EFFECTS:
//      n ends up as |n|    Make copy first if necessary.
bn_t unsafe_inv_bn(bn_t r, bn_t n)
{
    long maxval;
    LDBL f;
    bn_t orig_r;
    bn_t orig_n; // orig_bntmp1 not needed here
    int orig_bnlength;
    int orig_padding;
    int orig_rlength;
    int orig_shiftfactor;

    // use Newton's recursive method for zeroing in on 1/n : r=r(2-rn)

    bool signflag = false;
    if (is_bn_neg(n))
    {
        // will be a lot easier to deal with just positives
        signflag = true;
        neg_a_bn(n);
    }

    f = bntofloat(n);
    if (f == 0) // division by zero
    {
        max_bn(r);
        return r;
    }
    f = 1/f; // approximate inverse
    maxval = (1L << ((g_int_length << 3)-1)) - 1;
    if (f > maxval) // check for overflow
    {
        max_bn(r);
        return r;
    }
    if (f <= -maxval)
    {
        max_bn(r);
        neg_a_bn(r);
        return r;
    }

    // With Newton's Method, there is no need to calculate all the digits
    // every time.  The precision approximately doubles each iteration.
    // Save original values.
    orig_bnlength      = g_bn_length;
    orig_padding       = g_padding;
    orig_rlength       = g_r_length;
    orig_shiftfactor   = g_shift_factor;
    orig_r             = r;
    orig_n             = n;
    // orig_bntmp1        = g_bn_tmp1;

    // calculate new starting values
    g_bn_length = g_int_length + (int)(LDBL_DIG/LOG10_256) + 1; // round up
    if (g_bn_length > orig_bnlength)
    {
        g_bn_length = orig_bnlength;
    }
    calc_lengths();

    // adjust pointers
    r = orig_r + orig_bnlength - g_bn_length;
    // g_bn_tmp1 = orig_bntmp1 + orig_bnlength - g_bn_length;

    floattobn(r, f); // start with approximate inverse
    clear_bn(g_bn_tmp2); // will be used as 1.0 and 2.0

    for (int i = 0; i < 25; i++) // safety net, this shouldn't ever be needed
    {
        // adjust lengths
        g_bn_length <<= 1; // double precision
        if (g_bn_length > orig_bnlength)
        {
            g_bn_length = orig_bnlength;
        }
        calc_lengths();
        r = orig_r + orig_bnlength - g_bn_length;
        n = orig_n + orig_bnlength - g_bn_length;
        // g_bn_tmp1 = orig_bntmp1 + orig_bnlength - g_bn_length;

        unsafe_mult_bn(g_bn_tmp1, r, n); // g_bn_tmp1=rn
        inttobn(g_bn_tmp2, 1);  // g_bn_tmp2 = 1.0
        if (g_bn_length == orig_bnlength && cmp_bn(g_bn_tmp2, g_bn_tmp1+g_shift_factor) == 0)    // if not different
        {
            break;  // they must be the same
        }
        inttobn(g_bn_tmp2, 2); // g_bn_tmp2 = 2.0
        sub_bn(g_bn_tmp3, g_bn_tmp2, g_bn_tmp1+g_shift_factor); // g_bn_tmp3=2-rn
        unsafe_mult_bn(g_bn_tmp1, r, g_bn_tmp3); // g_bn_tmp1=r(2-rn)
        copy_bn(r, g_bn_tmp1+g_shift_factor); // r = g_bn_tmp1
    }

    // restore original values
    g_bn_length   = orig_bnlength;
    g_padding     = orig_padding;
    g_r_length    = orig_rlength;
    g_shift_factor = orig_shiftfactor;
    r             = orig_r;

    if (signflag)
    {
        neg_a_bn(r);
    }
    return r;
}

/********************************************************************/
// r = n1/n2
//      r - result of length g_bn_length
// uses g_bn_tmp1 - g_bn_tmp3 - global temp bignumbers
//  SIDE-EFFECTS:
//      n1, n2 can end up as GARBAGE
//      Make copies first if necessary.
bn_t unsafe_div_bn(bn_t r, bn_t n1, bn_t n2)
{
    int scale1;
    int scale2;
    int i;
    long maxval;
    LDBL a;
    LDBL b;
    LDBL f;

    // first, check for valid data
    a = bntofloat(n1);
    if (a == 0) // division into zero
    {
        clear_bn(r); // return 0
        return r;
    }
    b = bntofloat(n2);
    if (b == 0) // division by zero
    {
        max_bn(r);
        return r;
    }
    f = a/b; // approximate quotient
    maxval = (1L << ((g_int_length << 3)-1)) - 1;
    if (f > maxval) // check for overflow
    {
        max_bn(r);
        return r;
    }
    if (f <= -maxval)
    {
        max_bn(r);
        neg_a_bn(r);
        return r;
    }
    // appears to be ok, do division

    bool sign = false;
    if (is_bn_neg(n1))
    {
        neg_a_bn(n1);
        sign = !sign;
    }
    if (is_bn_neg(n2))
    {
        neg_a_bn(n2);
        sign = !sign;
    }

    // scale n1 and n2 so: |n| >= 1/256
    // scale = (int)(log(1/fabs(a))/LOG_256) = LOG_256(1/|a|)
    i = g_bn_length-1;
    while (i >= 0 && n1[i] == 0)
    {
        i--;
    }
    scale1 = g_bn_length - i - 2;
    if (scale1 < 0)
    {
        scale1 = 0;
    }
    i = g_bn_length-1;
    while (i >= 0 && n2[i] == 0)
    {
        i--;
    }
    scale2 = g_bn_length - i - 2;
    if (scale2 < 0)
    {
        scale2 = 0;
    }

    // shift n1, n2
    // important!, use std::memmove(), not std::memcpy()
    std::memmove(n1+scale1, n1, g_bn_length-scale1); // shift bytes over
    std::memset(n1, 0, scale1);  // zero out the rest
    std::memmove(n2+scale2, n2, g_bn_length-scale2); // shift bytes over
    std::memset(n2, 0, scale2);  // zero out the rest

    unsafe_inv_bn(r, n2);
    unsafe_mult_bn(g_bn_tmp1, n1, r);
    copy_bn(r, g_bn_tmp1+g_shift_factor); // r = g_bn_tmp1

    if (scale1 != scale2)
    {
        // Rescale r back to what it should be.  Overflow has already been checked
        if (scale1 > scale2) // answer is too big, adjust it
        {
            int scale = scale1-scale2;
            std::memmove(r, r+scale, g_bn_length-scale); // shift bytes over
            std::memset(r+g_bn_length-scale, 0, scale);  // zero out the rest
        }
        else if (scale1 < scale2) // answer is too small, adjust it
        {
            int scale = scale2-scale1;
            std::memmove(r+scale, r, g_bn_length-scale); // shift bytes over
            std::memset(r, 0, scale);                 // zero out the rest
        }
        // else scale1 == scale2

    }

    if (sign)
    {
        neg_a_bn(r);
    }

    return r;
}

/********************************************************************/
// sqrt(r)
// uses g_bn_tmp1 - g_bn_tmp6 - global temp bignumbers
//  SIDE-EFFECTS:
//      n ends up as |n|
bn_t sqrt_bn(bn_t r, bn_t n)
{
    int almost_match = 0;
    LDBL f;
    bn_t orig_r;
    bn_t orig_n;
    int orig_bnlength;
    int orig_padding;
    int orig_rlength;
    int orig_shiftfactor;

    // use Newton's recursive method for zeroing in on sqrt(n): r=.5(r+n/r)

    if (is_bn_neg(n))
    {
        // sqrt of a neg, return 0
        clear_bn(r);
        return r;
    }

    f = bntofloat(n);
    if (f == 0) // division by zero will occur
    {
        clear_bn(r); // sqrt(0) = 0
        return r;
    }
    f = sqrtl(f); // approximate square root
    // no need to check overflow

    // With Newton's Method, there is no need to calculate all the digits
    // every time.  The precision approximately doubles each iteration.
    // Save original values.
    orig_bnlength      = g_bn_length;
    orig_padding       = g_padding;
    orig_rlength       = g_r_length;
    orig_shiftfactor   = g_shift_factor;
    orig_r             = r;
    orig_n             = n;

    // calculate new starting values
    g_bn_length = g_int_length + (int)(LDBL_DIG/LOG10_256) + 1; // round up
    if (g_bn_length > orig_bnlength)
    {
        g_bn_length = orig_bnlength;
    }
    calc_lengths();

    // adjust pointers
    r = orig_r + orig_bnlength - g_bn_length;

    floattobn(r, f); // start with approximate sqrt
    copy_bn(g_bn_tmp4, r);

    for (int i = 0; i < 25; i++) // safety net, this shouldn't ever be needed
    {
        // adjust lengths
        g_bn_length <<= 1; // double precision
        if (g_bn_length > orig_bnlength)
        {
            g_bn_length = orig_bnlength;
        }
        calc_lengths();
        r = orig_r + orig_bnlength - g_bn_length;
        n = orig_n + orig_bnlength - g_bn_length;

        copy_bn(g_bn_tmp6, r);
        copy_bn(g_bn_tmp5, n);
        unsafe_div_bn(g_bn_tmp4, g_bn_tmp5, g_bn_tmp6);
        add_a_bn(r, g_bn_tmp4);
        half_a_bn(r);
        if (g_bn_length == orig_bnlength)
        {
            const int comp = std::abs(cmp_bn(r, g_bn_tmp4));
            if (comp < 8)  // if match or almost match
            {
                if (comp < 4  // perfect or near perfect match
                    || almost_match == 1)   // close enough for 2nd time
                {
                    break;
                }
                // this is the first time they almost matched
                almost_match++;
            }
        }
    }

    // restore original values
    g_bn_length   = orig_bnlength;
    g_padding     = orig_padding;
    g_r_length    = orig_rlength;
    g_shift_factor = orig_shiftfactor;
    // cppcheck-suppress uselessAssignmentPtrArg
    r             = orig_r;

    return r;
}

/********************************************************************/
// exp(r)
// uses g_bn_tmp1, g_bn_tmp2, g_bn_tmp3 - global temp bignumbers
bn_t exp_bn(bn_t r, bn_t n)
{
    U16 fact = 1;

    if (is_bn_zero(n))
    {
        inttobn(r, 1);
        return r;
    }

    // use Taylor Series (very slow convergence)
    inttobn(r, 1); // start with r=1.0
    copy_bn(g_bn_tmp2, r);
    while (true)
    {
        // copy n, if n is negative, mult_bn() alters n
        unsafe_mult_bn(g_bn_tmp3, g_bn_tmp2, copy_bn(g_bn_tmp1, n));
        copy_bn(g_bn_tmp2, g_bn_tmp3+g_shift_factor);
        div_a_bn_int(g_bn_tmp2, fact);
        if (!is_bn_not_zero(g_bn_tmp2))
        {
            break; // too small to register
        }
        add_a_bn(r, g_bn_tmp2);
        fact++;
    }
    return r;
}

/********************************************************************/
// ln(r)
// uses g_bn_tmp1 - g_bn_tmp6 - global temp bignumbers
//  SIDE-EFFECTS:
//      n ends up as |n|
bn_t unsafe_ln_bn(bn_t r, bn_t n)
{
    int almost_match = 0;
    long maxval;
    LDBL f;
    bn_t orig_r;
    bn_t orig_n;
    bn_t orig_bntmp5;
    bn_t orig_bntmp4;
    int orig_bnlength;
    int orig_padding;
    int orig_rlength;
    int orig_shiftfactor;

    // use Newton's recursive method for zeroing in on ln(n): r=r+n*exp(-r)-1

    if (is_bn_neg(n) || is_bn_zero(n))
    {
        // error, return largest neg value
        max_bn(r);
        neg_a_bn(r);
        return r;
    }

    f = bntofloat(n);
    f = logl(f); // approximate ln(x)
    maxval = (1L << ((g_int_length << 3)-1)) - 1;
    if (f > maxval) // check for overflow
    {
        max_bn(r);
        return r;
    }
    if (f <= -maxval)
    {
        max_bn(r);
        neg_a_bn(r);
        return r;
    }
    // appears to be ok, do ln

    // With Newton's Method, there is no need to calculate all the digits
    // every time.  The precision approximately doubles each iteration.
    // Save original values.
    orig_bnlength      = g_bn_length;
    orig_padding       = g_padding;
    orig_rlength       = g_r_length;
    orig_shiftfactor   = g_shift_factor;
    orig_r             = r;
    orig_n             = n;
    orig_bntmp5        = g_bn_tmp5;
    orig_bntmp4        = g_bn_tmp4;

    inttobn(g_bn_tmp4, 1); // set before setting new values

    // calculate new starting values
    g_bn_length = g_int_length + (int)(LDBL_DIG/LOG10_256) + 1; // round up
    if (g_bn_length > orig_bnlength)
    {
        g_bn_length = orig_bnlength;
    }
    calc_lengths();

    // adjust pointers
    r = orig_r + orig_bnlength - g_bn_length;
    g_bn_tmp5 = orig_bntmp5 + orig_bnlength - g_bn_length;
    g_bn_tmp4 = orig_bntmp4 + orig_bnlength - g_bn_length;

    floattobn(r, f); // start with approximate ln
    neg_a_bn(r); // -r
    copy_bn(g_bn_tmp5, r); // -r

    for (int i = 0; i < 25; i++) // safety net, this shouldn't ever be needed
    {
        // adjust lengths
        g_bn_length <<= 1; // double precision
        if (g_bn_length > orig_bnlength)
        {
            g_bn_length = orig_bnlength;
        }
        calc_lengths();
        r = orig_r + orig_bnlength - g_bn_length;
        n = orig_n + orig_bnlength - g_bn_length;
        g_bn_tmp5 = orig_bntmp5 + orig_bnlength - g_bn_length;
        g_bn_tmp4 = orig_bntmp4 + orig_bnlength - g_bn_length;
        exp_bn(g_bn_tmp6, r);     // exp(-r)
        unsafe_mult_bn(g_bn_tmp2, g_bn_tmp6, n);  // n*exp(-r)
        sub_a_bn(g_bn_tmp2+g_shift_factor, g_bn_tmp4);   // n*exp(-r) - 1
        sub_a_bn(r, g_bn_tmp2+g_shift_factor);        // -r - (n*exp(-r) - 1)

        if (g_bn_length == orig_bnlength)
        {
            const int comp = std::abs(cmp_bn(r, g_bn_tmp5));
            if (comp < 8)  // if match or almost match
            {
                if (comp < 4  // perfect or near perfect match
                    || almost_match == 1)   // close enough for 2nd time
                {
                    break;
                }
                // this is the first time they almost matched
                almost_match++;
            }
        }
        copy_bn(g_bn_tmp5, r); // -r
    }

    // restore original values
    g_bn_length   = orig_bnlength;
    g_padding     = orig_padding;
    g_r_length    = orig_rlength;
    g_shift_factor = orig_shiftfactor;
    r             = orig_r;
    g_bn_tmp5        = orig_bntmp5;
    g_bn_tmp4        = orig_bntmp4;

    neg_a_bn(r); // -(-r)
    return r;
}

/********************************************************************/
// sincos_bn(r)
// uses g_bn_tmp1 - g_bn_tmp2 - global temp bignumbers
//  SIDE-EFFECTS:
//      n ends up as |n| mod (pi/4)
bn_t unsafe_sincos_bn(bn_t s, bn_t c, bn_t n)
{
    U16 fact = 2;
    bool k = false;

#ifndef CALCULATING_BIG_PI
    // assure range 0 <= x < pi/4

    if (is_bn_zero(n))
    {
        clear_bn(s);    // sin(0) = 0
        inttobn(c, 1);  // cos(0) = 1
        return s;
    }

    bool signcos = false;
    bool signsin = false;
    bool switch_sincos = false;
    if (is_bn_neg(n))
    {
        signsin = !signsin; // sin(-x) = -sin(x), odd; cos(-x) = cos(x), even
        neg_a_bn(n);
    }
    // n >= 0

    double_bn(g_bn_tmp1, g_bn_pi); // 2*pi
    // this could be done with remainders, but it would probably be slower
    while (cmp_bn(n, g_bn_tmp1) >= 0)   // while n >= 2*pi
    {
        sub_a_bn(n, g_bn_tmp1);
    }
    // 0 <= n < 2*pi

    copy_bn(g_bn_tmp1, g_bn_pi); // pi
    if (cmp_bn(n, g_bn_tmp1) >= 0) // if n >= pi
    {
        sub_a_bn(n, g_bn_tmp1);
        signsin = !signsin;
        signcos = !signcos;
    }
    // 0 <= n < pi

    half_bn(g_bn_tmp1, g_bn_pi); // pi/2
    if (cmp_bn(n, g_bn_tmp1) > 0) // if n > pi/2
    {
        sub_bn(n, g_bn_pi, n);   // pi - n
        signcos = !signcos;
    }
    // 0 <= n < pi/2

    half_bn(g_bn_tmp1, g_bn_pi); // pi/2
    half_a_bn(g_bn_tmp1);      // pi/4
    if (cmp_bn(n, g_bn_tmp1) > 0) // if n > pi/4
    {
        half_bn(g_bn_tmp1, g_bn_pi); // pi/2
        sub_bn(n, g_bn_tmp1, n);  // pi/2 - n
        switch_sincos = !switch_sincos;
    }
    // 0 <= n < pi/4

    // this looks redundant, but n could now be zero when it wasn't before
    if (is_bn_zero(n))
    {
        clear_bn(s);    // sin(0) = 0
        inttobn(c, 1);  // cos(0) = 1
        return s;
    }

    // at this point, the double angle trig identities could be used as many
    // times as desired to reduce the range to pi/8, pi/16, etc...  Each time
    // the range is cut in half, the number of iterations required is reduced
    // by "quite a bit."  It's just a matter of testing to see what gives the
    // optimal results.
    // halves = g_bn_length / 10; */ /* this is experimental
    int halves = 1;
    for (int i = 0; i < halves; i++)
    {
        half_a_bn(n);
    }
#endif

    // use Taylor Series (very slow convergence)
    copy_bn(s, n); // start with s=n
    inttobn(c, 1); // start with c=1
    copy_bn(g_bn_tmp1, n); // the current x^n/n!

    while (true)
    {
        // even terms for cosine
        unsafe_mult_bn(g_bn_tmp2, g_bn_tmp1, n);
        copy_bn(g_bn_tmp1, g_bn_tmp2+g_shift_factor);
        div_a_bn_int(g_bn_tmp1, fact++);
        if (!is_bn_not_zero(g_bn_tmp1))
        {
            break; // too small to register
        }
        if (k)   // alternate between adding and subtracting
        {
            add_a_bn(c, g_bn_tmp1);
        }
        else
        {
            sub_a_bn(c, g_bn_tmp1);
        }

        // odd terms for sine
        unsafe_mult_bn(g_bn_tmp2, g_bn_tmp1, n);
        copy_bn(g_bn_tmp1, g_bn_tmp2+g_shift_factor);
        div_a_bn_int(g_bn_tmp1, fact++);
        if (!is_bn_not_zero(g_bn_tmp1))
        {
            break; // too small to register
        }
        if (k)   // alternate between adding and subtracting
        {
            add_a_bn(s, g_bn_tmp1);
        }
        else
        {
            sub_a_bn(s, g_bn_tmp1);
        }
        k = !k; // toggle
#ifdef CALCULATING_BIG_PI
        std::printf("."); // lets you know it's doing something
#endif
    }

#ifndef CALCULATING_BIG_PI
    // now need to undo what was done by cutting angles in half
    inttobn(g_bn_tmp1, 1);
    for (int i = 0; i < halves; i++)
    {
        unsafe_mult_bn(g_bn_tmp2, s, c); // no need for safe mult
        double_bn(s, g_bn_tmp2+g_shift_factor); // sin(2x) = 2*sin(x)*cos(x)
        unsafe_square_bn(g_bn_tmp2, c);
        double_a_bn(g_bn_tmp2+g_shift_factor);
        sub_bn(c, g_bn_tmp2+g_shift_factor, g_bn_tmp1); // cos(2x) = 2*cos(x)*cos(x) - 1
    }

    if (switch_sincos)
    {
        copy_bn(g_bn_tmp1, s);
        copy_bn(s, c);
        copy_bn(c, g_bn_tmp1);
    }
    if (signsin)
    {
        neg_a_bn(s);
    }
    if (signcos)
    {
        neg_a_bn(c);
    }
#endif

    return s; // return sine I guess
}

/********************************************************************/
// atan(r)
// uses g_bn_tmp1 - g_bn_tmp5 - global temp bignumbers
//  SIDE-EFFECTS:
//      n ends up as |n| or 1/|n|
bn_t unsafe_atan_bn(bn_t r, bn_t n)
{
    int almost_match = 0;
    LDBL f;
    bn_t orig_r;
    bn_t orig_n;
    bn_t orig_bn_pi;
    bn_t orig_bntmp3;
    int orig_bnlength;
    int orig_padding;
    int orig_rlength;
    int orig_shiftfactor;

    // use Newton's recursive method for zeroing in on atan(n): r=r-cos(r)(sin(r)-n*cos(r))
    bool signflag = false;
    if (is_bn_neg(n))
    {
        signflag = true;
        neg_a_bn(n);
    }

    // If n is very large, atanl() won't give enough decimal places to be a
    // good enough initial guess for Newton's Method.  If it is larger than
    // say, 1, atan(n) = pi/2 - acot(n) = pi/2 - atan(1/n).

    f = bntofloat(n);
    bool large_arg = f > 1.0;
    if (large_arg)
    {
        unsafe_inv_bn(g_bn_tmp3, n);
        copy_bn(n, g_bn_tmp3);
        f = bntofloat(n);
    }

    clear_bn(g_bn_tmp3); // not really necessary, but makes things more consistent

    // With Newton's Method, there is no need to calculate all the digits
    // every time.  The precision approximately doubles each iteration.
    // Save original values.
    orig_bnlength      = g_bn_length;
    orig_padding       = g_padding;
    orig_rlength       = g_r_length;
    orig_shiftfactor   = g_shift_factor;
    orig_bn_pi         = g_bn_pi;
    orig_r             = r;
    orig_n             = n;
    orig_bntmp3        = g_bn_tmp3;

    // calculate new starting values
    g_bn_length = g_int_length + (int)(LDBL_DIG/LOG10_256) + 1; // round up
    if (g_bn_length > orig_bnlength)
    {
        g_bn_length = orig_bnlength;
    }
    calc_lengths();

    // adjust pointers
    r = orig_r + orig_bnlength - g_bn_length;
    g_bn_pi = orig_bn_pi + orig_bnlength - g_bn_length;
    g_bn_tmp3 = orig_bntmp3 + orig_bnlength - g_bn_length;

    f = atanl(f); // approximate arctangent
    // no need to check overflow

    floattobn(r, f); // start with approximate atan
    copy_bn(g_bn_tmp3, r);

    for (int i = 0; i < 25; i++) // safety net, this shouldn't ever be needed
    {
        // adjust lengths
        g_bn_length <<= 1; // double precision
        if (g_bn_length > orig_bnlength)
        {
            g_bn_length = orig_bnlength;
        }
        calc_lengths();
        r = orig_r + orig_bnlength - g_bn_length;
        n = orig_n + orig_bnlength - g_bn_length;
        g_bn_pi = orig_bn_pi + orig_bnlength - g_bn_length;
        g_bn_tmp3 = orig_bntmp3 + orig_bnlength - g_bn_length;

#ifdef CALCULATING_BIG_PI
        std::printf("\natan() loop #%i, g_bn_length=%i\nsincos() loops\n", i, g_bn_length);
#endif
        unsafe_sincos_bn(g_bn_tmp4, g_bn_tmp5, g_bn_tmp3);   // sin(r), cos(r)
        copy_bn(g_bn_tmp3, r); // restore g_bn_tmp3 from sincos_bn()
        copy_bn(g_bn_tmp1, g_bn_tmp5);
        unsafe_mult_bn(g_bn_tmp2, n, g_bn_tmp1);     // n*cos(r)
        sub_a_bn(g_bn_tmp4, g_bn_tmp2+g_shift_factor); // sin(r) - n*cos(r)
        unsafe_mult_bn(g_bn_tmp1, g_bn_tmp5, g_bn_tmp4); // cos(r) * (sin(r) - n*cos(r))
        sub_a_bn(r, g_bn_tmp1+g_shift_factor); // r - cos(r) * (sin(r) - n*cos(r))

#ifdef CALCULATING_BIG_PI
        putchar('\n');
        bn_hexdump(r);
#endif
        if (g_bn_length == orig_bnlength)
        {
            const int comp = std::abs(cmp_bn(r, g_bn_tmp3));
            if (comp < 8)  // if match or almost match
            {
#ifdef CALCULATING_BIG_PI
                std::printf("atan() loop comp=%i\n", comp);
#endif
                if (comp < 4  // perfect or near perfect match
                    || almost_match == 1)   // close enough for 2nd time
                {
                    break;
                }
                else     // this is the first time they almost matched
                {
                    almost_match++;
                }
            }
#ifdef CALCULATING_BIG_PI
            if (comp >= 8)
            {
                std::printf("atan() loop comp=%i\n", comp);
            }
#endif
        }

        copy_bn(g_bn_tmp3, r); // make a copy for later comparison
    }

    // restore original values
    g_bn_length   = orig_bnlength;
    g_padding     = orig_padding;
    g_r_length    = orig_rlength;
    g_shift_factor = orig_shiftfactor;
    g_bn_pi         = orig_bn_pi;
    r             = orig_r;
    g_bn_tmp3        = orig_bntmp3;

    if (large_arg)
    {
        half_bn(g_bn_tmp3, g_bn_pi);  // pi/2
        sub_a_bn(g_bn_tmp3, r);     // pi/2 - atan(1/n)
        copy_bn(r, g_bn_tmp3);
    }

    if (signflag)
    {
        neg_a_bn(r);
    }
    return r;
}

/********************************************************************/
// atan2(r, ny, nx)
// uses g_bn_tmp1 - g_bn_tmp6 - global temp bigfloats
bn_t unsafe_atan2_bn(bn_t r, bn_t ny, bn_t nx)
{
    int signx;
    int signy;

    signx = sign_bn(nx);
    signy = sign_bn(ny);

    if (signy == 0)
    {
        if (signx < 0)
        {
            copy_bn(r, g_bn_pi); // negative x axis, 180 deg
        }
        else        // signx >= 0    positive x axis, 0
        {
            clear_bn(r);
        }
        return r;
    }
    if (signx == 0)
    {
        copy_bn(r, g_bn_pi); // y axis
        half_a_bn(r);      // +90 deg
        if (signy < 0)
        {
            neg_a_bn(r);    // -90 deg
        }
        return r;
    }

    if (signy < 0)
    {
        neg_a_bn(ny);
    }
    if (signx < 0)
    {
        neg_a_bn(nx);
    }
    unsafe_div_bn(g_bn_tmp6, ny, nx);
    unsafe_atan_bn(r, g_bn_tmp6);
    if (signx < 0)
    {
        sub_bn(r, g_bn_pi, r);
    }
    if (signy < 0)
    {
        neg_a_bn(r);
    }
    return r;
}


/**********************************************************************/
// The rest of the functions are "safe" versions of the routines that
// have side effects which alter the parameters
/**********************************************************************/

/**********************************************************************/
bn_t full_mult_bn(bn_t r, bn_t n1, bn_t n2)
{
    bool sign1 = is_bn_neg(n1);
    bool sign2 = is_bn_neg(n2);
    unsafe_full_mult_bn(r, n1, n2);
    if (sign1)
    {
        neg_a_bn(n1);
    }
    if (sign2)
    {
        neg_a_bn(n2);
    }
    return r;
}

/**********************************************************************/
bn_t mult_bn(bn_t r, bn_t n1, bn_t n2)
{
    bool sign1 = is_bn_neg(n1);
    bool sign2 = is_bn_neg(n2);
    unsafe_mult_bn(r, n1, n2);
    if (sign1)
    {
        neg_a_bn(n1);
    }
    if (sign2)
    {
        neg_a_bn(n2);
    }
    return r;
}

/**********************************************************************/
bn_t full_square_bn(bn_t r, bn_t n)
{
    bool sign = is_bn_neg(n);
    unsafe_full_square_bn(r, n);
    if (sign)
    {
        neg_a_bn(n);
    }
    return r;
}

/**********************************************************************/
bn_t square_bn(bn_t r, bn_t n)
{
    bool sign = is_bn_neg(n);
    unsafe_square_bn(r, n);
    if (sign)
    {
        neg_a_bn(n);
    }
    return r;
}

/**********************************************************************/
bn_t div_bn_int(bn_t r, bn_t n, U16 u)
{
    bool sign = is_bn_neg(n);
    unsafe_div_bn_int(r, n, u);
    if (sign)
    {
        neg_a_bn(n);
    }
    return r;
}

/**********************************************************************/
char *bntostr(char *s, int dec, bn_t r)
{
    return unsafe_bntostr(s, dec, copy_bn(g_bn_tmp_copy2, r));
}

/**********************************************************************/
bn_t inv_bn(bn_t r, bn_t n)
{
    bool sign = is_bn_neg(n);
    unsafe_inv_bn(r, n);
    if (sign)
    {
        neg_a_bn(n);
    }
    return r;
}

/**********************************************************************/
bn_t div_bn(bn_t r, bn_t n1, bn_t n2)
{
    copy_bn(g_bn_tmp_copy1, n1);
    copy_bn(g_bn_tmp_copy2, n2);
    return unsafe_div_bn(r, g_bn_tmp_copy1, g_bn_tmp_copy2);
}

/**********************************************************************/
bn_t ln_bn(bn_t r, bn_t n)
{
    copy_bn(g_bn_tmp_copy1, n); // allows r and n to overlap memory
    unsafe_ln_bn(r, g_bn_tmp_copy1);
    return r;
}

/**********************************************************************/
bn_t sincos_bn(bn_t s, bn_t c, bn_t n)
{
    return unsafe_sincos_bn(s, c, copy_bn(g_bn_tmp_copy1, n));
}

/**********************************************************************/
bn_t atan_bn(bn_t r, bn_t n)
{
    bool sign = is_bn_neg(n);
    unsafe_atan_bn(r, n);
    if (sign)
    {
        neg_a_bn(n);
    }
    return r;
}

/**********************************************************************/
bn_t atan2_bn(bn_t r, bn_t ny, bn_t nx)
{
    copy_bn(g_bn_tmp_copy1, ny);
    copy_bn(g_bn_tmp_copy2, nx);
    unsafe_atan2_bn(r, g_bn_tmp_copy1, g_bn_tmp_copy2);
    return r;
}

bool is_bn_zero(bn_t n)
{
    return !is_bn_not_zero(n);
}
