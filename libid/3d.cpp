//
// A word about the 3D library. Even though this library supports
// three dimensions, the matrices are 4x4 for the following reason.
// With normal 3 dimensional vectors, translation is an ADDITION,
// and rotation is a MULTIPLICATION. A vector {x,y,z} is represented
// as a 4-tuple {x,y,z,1}. It is then possible to define a 4x4
// matrix such that multiplying the vector by the matrix translates
// the vector. This allows combinations of translation and rotation
// to be obtained in a single matrix by multiplying a translation
// matrix and a rotation matrix together. Note that in the code,
// vectors have three components; since the fourth component is
// always 1, that value is not included in the vector variable to
// save space, but the routines make use of the fourth component
// (see vec_mult()). Similarly, the fourth column of EVERY matrix is
// always
//          0
//          0
//          0
//          1
// but currently the C version of a matrix includes this even though
// it could be left out of the data structure and assumed in the
// routines. Vectors are ROW vectors, and are always multiplied with
// matrices FROM THE LEFT (e.g. vector*matrix). Also note the order
// of indices of a matrix is matrix[row][column], and in usual C
// fashion, numbering starts with 0.
//
// TRANSLATION MATRIX =  1     0     0     0
//                       0     1     0     0
//                       0     0     1     0
//                       Tx    Ty    Tz    1
//
// SCALE MATRIX =        Sx    0     0     0
//                       0     Sy    0     0
//                       0     0     Sz    0
//                       0     0     0     1
//
// Rotation about x axis i degrees:
// ROTX(i) =             1     0     0     0
//                       0   cosi  sini    0
//                       0  -sini  cosi    0
//                       0     0     0     1
//
// Rotation about y axis i degrees:
// ROTY(i) =           cosi    0  -sini    0
//                       0     1     0     0
//                     sini    0   cosi    0
//                       0     0     0     1
//
// Rotation about z axis i degrees:
// ROTZ(i) =           cosi  sini    0     0
//                    -sini  cosi    0     0
//                       0     0     1     0
//                       0     0     0     1
//
#include "3d.h"

#include "fixed_pt.h"
#include "fractals.h"
#include "line3d.h"

#include <cfloat>
#include <cmath>
#include <cstring>

int     g_init_3d[20] = { 0 };     // '3d=nn/nn/nn/...' values

// initialize a matrix and set to identity matrix (all 0's, 1's on diagonal)
void identity(MATRIX m)
{
    for (int i = 0 ; i < CMAX; i++)
    {
        for (int j = 0; j < RMAX; j++)
        {
            if (i == j)
            {
                m[j][i] = 1.0;
            }
            else
            {
                m[j][i] = 0.0;
            }
        }
    }
}

// Multiply two matrices
void mat_mul(MATRIX mat1, MATRIX mat2, MATRIX mat3)
{
    // result stored in MATRIX new to avoid problems
    //  in case parameter mat3 == mat2 or mat 1
    MATRIX newmat;
    for (int i = 0; i < 4; i++)
    {
        for (int j = 0; j < 4; j++)
        {
            newmat[j][i] =  mat1[j][0]*mat2[0][i]+
                            mat1[j][1]*mat2[1][i]+
                            mat1[j][2]*mat2[2][i]+
                            mat1[j][3]*mat2[3][i];
        }
    }
    std::memcpy(mat3, newmat, sizeof(newmat));
}

// multiply a matrix by a scalar
void scale(double sx, double sy, double sz, MATRIX m)
{
    MATRIX scale;
    identity(scale);
    scale[0][0] = sx;
    scale[1][1] = sy;
    scale[2][2] = sz;
    mat_mul(m, scale, m);
}

// rotate about X axis
void xrot(double theta, MATRIX m)
{
    MATRIX rot;
    double sintheta;
    double costheta;
    sintheta = std::sin(theta);
    costheta = std::cos(theta);
    identity(rot);
    rot[1][1] = costheta;
    rot[1][2] = -sintheta;
    rot[2][1] = sintheta;
    rot[2][2] = costheta;
    mat_mul(m, rot, m);
}

// rotate about Y axis
void yrot(double theta, MATRIX m)
{
    MATRIX rot;
    double sintheta;
    double costheta;
    sintheta = std::sin(theta);
    costheta = std::cos(theta);
    identity(rot);
    rot[0][0] = costheta;
    rot[0][2] = sintheta;
    rot[2][0] = -sintheta;
    rot[2][2] = costheta;
    mat_mul(m, rot, m);
}

// rotate about Z axis
void zrot(double theta, MATRIX m)
{
    MATRIX rot;
    double sintheta;
    double costheta;
    sintheta = std::sin(theta);
    costheta = std::cos(theta);
    identity(rot);
    rot[0][0] = costheta;
    rot[0][1] = -sintheta;
    rot[1][0] = sintheta;
    rot[1][1] = costheta;
    mat_mul(m, rot, m);
}

// translate
void trans(double tx, double ty, double tz, MATRIX m)
{
    MATRIX trans;
    identity(trans);
    trans[3][0] = tx;
    trans[3][1] = ty;
    trans[3][2] = tz;
    mat_mul(m, trans, m);
}

// cross product  - useful because cross is perpendicular to v and w
int cross_product(VECTOR v, VECTOR w, VECTOR cross)
{
    VECTOR tmp;
    tmp[0] =  v[1]*w[2] - w[1]*v[2];
    tmp[1] =  w[0]*v[2] - v[0]*w[2];
    tmp[2] =  v[0]*w[1] - w[0]*v[1];
    cross[0] = tmp[0];
    cross[1] = tmp[1];
    cross[2] = tmp[2];
    return 0;
}

// normalize a vector to length 1
bool normalize_vector(VECTOR v)
{
    double vlength;
    vlength = dot_product(v, v);

    // bailout if zero vlength
    if (vlength < FLT_MIN || vlength > FLT_MAX)
    {
        return true;
    }
    vlength = std::sqrt(vlength);
    if (vlength < FLT_MIN)
    {
        return true;
    }

    v[0] /= vlength;
    v[1] /= vlength;
    v[2] /= vlength;
    return false;
}

// multiply source vector s by matrix m, result in target t
// used to apply transformations to a vector
int vmult(VECTOR s, MATRIX m, VECTOR t)
{
    VECTOR tmp;
    for (int j = 0; j < CMAX-1; j++)
    {
        tmp[j] = 0.0;
        for (int i = 0; i < RMAX-1; i++)
        {
            tmp[j] += s[i]*m[i][j];
        }
        // vector is really four dimensional with last component always 1
        tmp[j] += m[3][j];
    }
    // set target = tmp. Necessary to use tmp in case source = target
    std::memcpy(t, tmp, sizeof(tmp));
    return 0;
}

// multiply vector s by matrix m, result in s
// use with a function pointer in line3d.c
// must coordinate calling conventions with
// mult_vec in general.asm
void mult_vec(VECTOR s)
{
    VECTOR tmp;
    for (int j = 0; j < CMAX-1; j++)
    {
        tmp[j] = 0.0;
        for (int i = 0; i < RMAX-1; i++)
        {
            tmp[j] += s[i]*g_m[i][j];
        }
        // vector is really four dimensional with last component always 1
        tmp[j] += g_m[3][j];
    }
    // set target = tmp. Necessary to use tmp in case source = target
    std::memcpy(s, tmp, sizeof(tmp));
}

// perspective projection of vector v with respect to viewpoint vector view
int perspective(VECTOR v)
{
    double denom;
    denom = g_view[2] - v[2];

    if (denom >= 0.0)
    {
        v[0] = g_bad_value;   // clipping will catch these values
        v[1] = g_bad_value;   // so they won't plot values BEHIND viewer
        v[2] = g_bad_value;
        return -1;
    }
    v[0] = (v[0]*g_view[2] - g_view[0]*v[2])/denom;
    v[1] = (v[1]*g_view[2] - g_view[1]*v[2])/denom;

    // calculation of z if needed later
    // v[2] =  v[2]/denom;
    return 0;
}

// long version of vmult and perspective combined for speed
int longvmultpersp(
    LVECTOR s, LMATRIX m, LVECTOR t0, LVECTOR t, LVECTOR lview, int bitshift)
{
    // s: source vector
    // m: transformation matrix
    // t0: after transformation, before persp
    // t: target vector
    // lview: perspective viewer coordinates
    // bitshift: fixed point conversion bitshift
    LVECTOR tmp;
    int k;
    g_overflow = false;
    k = CMAX-1;                  // shorten the math if non-perspective and non-illum
    if (lview[2] == 0 && t0[0] == 0)
    {
        k--;
    }

    for (int j = 0; j < k; j++)
    {
        tmp[j] = 0;
        for (int i = 0; i < RMAX-1; i++)
        {
            tmp[j] += multiply(s[i], m[i][j], bitshift);
        }
        // vector is really four dimensional with last component always 1
        tmp[j] += m[3][j];
    }
    if (t0[0]) // first component of  t0 used as flag
    {
        // faster than for loop, if less general
        t0[0] = tmp[0];
        t0[1] = tmp[1];
        t0[2] = tmp[2];
    }
    if (lview[2] != 0)           // perspective 3D
    {

        LVECTOR tmpview;
        long denom;

        denom = lview[2] - tmp[2];
        if (denom >= 0)           // bail out if point is "behind" us
        {
            t[0] = g_bad_value;
            t[0] = t[0] << bitshift;
            t[1] = t[0];
            t[2] = t[0];
            return -1;
        }

        // doing math in this order helps prevent overflow
        tmpview[0] = divide(lview[0], denom, bitshift);
        tmpview[1] = divide(lview[1], denom, bitshift);
        tmpview[2] = divide(lview[2], denom, bitshift);

        tmp[0] = multiply(tmp[0], tmpview[2], bitshift) -
                 multiply(tmpview[0], tmp[2], bitshift);

        tmp[1] = multiply(tmp[1], tmpview[2], bitshift) -
                 multiply(tmpview[1], tmp[2], bitshift);
    }

    // set target = tmp. Necessary to use tmp in case source = target
    // faster than for loop, if less general
    t[0] = tmp[0];
    t[1] = tmp[1];
    t[2] = tmp[2];
    return g_overflow ? 1 : 0;
}

// Long version of perspective. Because of use of fixed point math, there
// is danger of overflow and underflow
int longpersp(LVECTOR lv, LVECTOR lview, int bitshift)
{
    LVECTOR tmpview;
    long denom;
    g_overflow = false;
    denom = lview[2] - lv[2];
    if (denom >= 0)              // bail out if point is "behind" us
    {
        lv[0] = g_bad_value;
        lv[0] = lv[0] << bitshift;
        lv[1] = lv[0];
        lv[2] = lv[0];
        return -1;
    }

    // doing math in this order helps prevent overflow
    tmpview[0] = divide(lview[0], denom, bitshift);
    tmpview[1] = divide(lview[1], denom, bitshift);
    tmpview[2] = divide(lview[2], denom, bitshift);

    lv[0] = multiply(lv[0], tmpview[2], bitshift) -
            multiply(tmpview[0], lv[2], bitshift);

    lv[1] = multiply(lv[1], tmpview[2], bitshift) -
            multiply(tmpview[1], lv[2], bitshift);

    return g_overflow ? 1 : 0;
}

int longvmult(LVECTOR s, LMATRIX m, LVECTOR t, int bitshift)
{
    LVECTOR tmp;
    int k;
    g_overflow = false;
    k = CMAX-1;

    for (int j = 0; j < k; j++)
    {
        tmp[j] = 0;
        for (int i = 0; i < RMAX-1; i++)
        {
            tmp[j] += multiply(s[i], m[i][j], bitshift);
        }
        // vector is really four dimensional with last component always 1
        tmp[j] += m[3][j];
    }

    // set target = tmp. Necessary to use tmp in case source = target
    // faster than for loop, if less general
    t[0] = tmp[0];
    t[1] = tmp[1];
    t[2] = tmp[2];
    return g_overflow ? 1 : 0;
}
