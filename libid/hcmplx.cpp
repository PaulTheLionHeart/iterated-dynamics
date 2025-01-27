// some hyper complex functions
#include "port.h"
#include "prototyp.h"

#include "hcmplx.h"

#include "mpmath.h"

void HComplexMult(DHyperComplex *arg1, DHyperComplex *arg2, DHyperComplex *out)
{
    out->x = arg1->x * arg2->x - arg1->y * arg2->y
             - arg1->z * arg2->z + arg1->t * arg2->t;
    out->y = arg1->y * arg2->x + arg1->x * arg2->y
             - arg1->t * arg2->z - arg1->z * arg2->t;
    out->z = arg1->z * arg2->x - arg1->t * arg2->y
             + arg1->x * arg2->z - arg1->y * arg2->t;
    out->t = arg1->t * arg2->x + arg1->z * arg2->y
             + arg1->y * arg2->z + arg1->x * arg2->t;
}

void HComplexSqr(DHyperComplex *arg, DHyperComplex *out)
{
    out->x = arg->x * arg->x - arg->y * arg->y
             - arg->z * arg->z + arg->t * arg->t;
    out->y = 2 * arg->x * arg->y - 2 * arg->z * arg->t;
    out->z = 2 * arg->z * arg->x - 2 * arg->t * arg->y;
    out->t = 2 * arg->t * arg->x + 2 * arg->z * arg->y;
}

int HComplexInv(DHyperComplex *arg, DHyperComplex *out)
{
    double det;
    double mod;
    double xt_minus_yz;

    det = (sqr(arg->x - arg->t) + sqr(arg->y + arg->z))*
          (sqr(arg->x + arg->t) + sqr(arg->y - arg->z));

    if (det == 0.0)
    {
        return -1;
    }
    mod = sqr(arg->x) + sqr(arg->y) + sqr(arg->z) + sqr(arg->t);
    xt_minus_yz = arg->x * arg->t - arg->y * arg->z;

    out->x = (arg->x * mod - 2 * arg->t * xt_minus_yz)/det;
    out->y = (-arg->y * mod - 2 * arg->z * xt_minus_yz)/det;
    out->z = (-arg->z * mod - 2 * arg->y * xt_minus_yz)/det;
    out->t = (arg->t * mod - 2 * arg->x * xt_minus_yz)/det;
    return 0;
}

void HComplexAdd(DHyperComplex *arg1, DHyperComplex *arg2, DHyperComplex *out)
{
    out->x = arg1->x + arg2->x;
    out->y = arg1->y + arg2->y;
    out->z = arg1->z + arg2->z;
    out->t = arg1->t + arg2->t;
}

void HComplexSub(DHyperComplex *arg1, DHyperComplex *arg2, DHyperComplex *out)
{
    out->x = arg1->x - arg2->x;
    out->y = arg1->y - arg2->y;
    out->z = arg1->z - arg2->z;
    out->t = arg1->t - arg2->t;
}

void HComplexMinus(DHyperComplex *arg1, DHyperComplex *out)
{
    out->x = -arg1->x;
    out->y = -arg1->y;
    out->z = -arg1->z;
    out->t = -arg1->t;
}

// extends the unary function f to *h1
void HComplexTrig0(DHyperComplex *h, DHyperComplex *out)
{
    /* This is the whole beauty of Hypercomplex numbers - *ANY* unary
       complex valued function of a complex variable can easily
       be generalized to hypercomplex numbers */

    DComplex a;
    DComplex b;
    DComplex resulta;
    DComplex resultb;

    // convert to duplex form
    a.x = h->x - h->t;
    a.y = h->y + h->z;
    b.x = h->x + h->t;
    b.y = h->y - h->z;

    // apply function to each part
    CMPLXtrig0(a, resulta);
    CMPLXtrig0(b, resultb);

    // convert back
    out->x = (resulta.x + resultb.x)/2;
    out->y = (resulta.y + resultb.y)/2;
    out->z = (resulta.y - resultb.y)/2;
    out->t = (resultb.x - resulta.x)/2;
}
