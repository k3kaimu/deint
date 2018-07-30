// Written in the D programming language.

/**
This module provides numerical integration by double-exponential (DE) formula.
*/
module deint.deint;

import std.algorithm;
import std.array;
import std.functional;
import std.math;
import std.typecons;
import std.traits;


/**
This template checks that a type `T` is a kind of NumInt.
*/
enum bool isNumericalIntegrate(T)
    = is(typeof((T integ){
        foreach(xw; lockstep(integ.xs, integ.ws)) {}
    }));



/** This structure stores computing points and weights for numerical integration.
 * In this type, a numerical integration is computed as
 * int f(x) dx = sum_i f(xs[i]) ws[i]
 */
struct NumInt(X, W)
{
    /** Constructor
     * Params:
     *     xs = computing points
     *     ws = weights for each computing points
     */
    this(immutable(X)[] xs, immutable(W)[] ws)
    in{
        assert(xs.length == ws.length);
    }
    do {
        _xs = xs;
        _ws = ws;
    }


  const @property
  {
    immutable(X)[] xs() { return _xs; }
    immutable(W)[] ws() { return _ws; }
  }


    /** Execute integration
    */
    auto integrate(Fn)(scope Fn func) const
    {
        alias T = Unqual!(typeof(func(_xs[0]) * _ws[0]));

        T sum = 0;
        foreach(i; 0 .. _xs.length)
            sum += func(_xs[i]) * _ws[i];

        return sum;
    }


  private:
    immutable(X)[] _xs;
    immutable(W)[] _ws;
}


/** This function create a instance of NumInt!(F, F) for the DE formula.
 * It is also known as "Tanh-sinh quadrature".
 * In the DE formula, the integration (1) is converted to (2).
 * (1) int_{xa}^{xb} f(x) dx
 * (2) int_{ta}^{tb} f(g(t)) g'(t) dt
 * 
 * The type of DE formula is automatically decided from the given interval of the integration.
 *
 * Params:
 *      xa = starting value of original integration.
 *      xb = end value of original integration.
 *      isExpDecay = if the integration is formed as int_a^b f(x) exp(-x) dx, this value is Yes. otherwise No.
 *      trapN = division points of trapezoidal quadrature.
 *      ta = starting value of integration transformed by DE-formula.
 *      tb = starting value of integration transformed by DE-formula.
 *
 * Reference:
 * $(HTTP en.wikipedia.org/wiki/Tanh-sinh_quadrature)
*/
NumInt!(F, F) makeDEInt(F)(
    F xa, F xb,
    Flag!"isExpDecay" isExpDecay = No.isExpDecay,
    size_t trapN = 100,
    F ta = -5, F tb = 5)
{
    NumInt!(F, F) _makeParamsImpl(scope F[2] delegate(F) fn)
    {
        immutable(F)[] xs;
        immutable(F)[] ws;
        immutable F h = (tb - ta) / (trapN-1);
        foreach(i; 0 .. trapN) {
            immutable xWt = fn(i * h + ta);
            xs ~= xWt[0];
            ws ~= xWt[1] * h;
        }

        return NumInt!(F, F)(xs, ws);
    }


    if(xa > xb) {
        auto invInt = .makeDEInt!F(xb, xa, isExpDecay, trapN, ta, tb);
        auto newws = invInt.ws.map!"cast(immutable)(-a)".array;
        return NumInt!(F, F)(invInt.xs, newws);
    } else if(xa == -F.infinity && xb != F.infinity) {
        auto tmpInt = .makeDEInt!F(-xb, F.infinity, isExpDecay, trapN, ta, tb);
        auto newxs = tmpInt.xs.map!"cast(immutable)(-a)".array;
        return NumInt!(F, F)(newxs, tmpInt.ws);
    }


    if(xa == -F.infinity && xb == F.infinity){
        assert(!isExpDecay);

        return _makeParamsImpl(delegate(F t){
            immutable F
                sinht = sinh(t),
                x = sinh(PI / 2 * sinht),
                dx = cosh(PI / 2 * sinht) * PI / 2 * cosh(t);

            return cast(F[2])[x, dx];
        });
    }else if(xb == F.infinity) {
        return _makeParamsImpl(delegate(F t){
            if(!isExpDecay){
                real x = exp(PI / 2 * sinh(t)),
                        dx = x * PI / 2 * cosh(t);

                return cast(F[2])[x + xa, dx];
            }else{
                real expmt = exp(-t),
                        x = exp(t - expmt),
                        dx = (1 + expmt) * x;


                return cast(F[2])[x + xa, dx]; 
            }
        });
    }else{
        immutable F diff2 = (xb - xa) / 2,
                    avg2 = (xb + xa) / 2;

        return _makeParamsImpl(delegate(F t){
            immutable F
                cosht = cosh(t),
                sinht = sinh(t),
                x = tanh(PI / 2 * sinht) * diff2 + avg2,
                cosh2 = cosh(PI / 2 * sinht)^^2,
                dx = PI / 2 * cosht / cosh2;

            return cast(F[2])[x, dx * diff2];
        });
    }
}


///
unittest
{
    // integration on [0, 1]
    auto int01 = makeDEInt!real(0, 1);

    // int_0^1 x dx = 0.5
    assert(int01.integrate((real x) => x).approxEqual(0.5));

    // int_0^1 x^^2 dx = 1/3
    assert(int01.integrate((real x) => x^^2).approxEqual(1/3.0));


    // integration on [-inf, inf]
    auto intII = makeDEInt!real(-real.infinity, real.infinity);

    // Gaussian integral
    assert(intII.integrate((real x) => exp(-x^^2)).approxEqual(sqrt(PI)));


    import std.mathspecial;
    // integration on [1, inf] and integrand is expected to decay exponentially
    auto intFI = makeDEInt!real(1, real.infinity, Yes.isExpDecay);

    // compute incomplete gamma function: int_1^inf x * exp(-x) dx = Gamma(2, 1)
    assert(intFI.integrate((real x) => x * exp(-x)).approxEqual(gammaIncompleteCompl(2, 1) * gamma(2)));
}

// Test of \int_a^b exp(x) dx
unittest 
{
    // [a, b]
    real[2][] dataset = [
        [0.0, 1.0],
        [-1.0, 1.0],
        [0.0, -2.0],
        [2.0, 0.0],
        [-10.0, 0.0],
        [0.0, 10.0]
    ];

    foreach(data; dataset) {
        real intDE = makeDEInt!real(data[0], data[1]).integrate((real x) => exp(x));
        real truevalue = exp(data[1]) - exp(data[0]);

        assert(approxEqual(intDE, truevalue));
    }


    // reverse interval
    foreach(data; dataset) {
        //real intDE = intAToB_DE(data[0], data[1], 100, (real x) => exp(x));
        real intDE = makeDEInt!real(data[1], data[0]).integrate((real x) => exp(x));
        real truevalue = -(exp(data[1]) - exp(data[0]));

        assert(approxEqual(intDE, truevalue));
    }
}

// Complementary error function: erfc(a) = 2/sqrt(pi) * int_a^inf e^{-x^2} dx
unittest 
{
    import std.mathspecial;

    real[] dataset = [
        0.0,
        0.01,
        0.1,
        1,
        10,
        100,
    ];

    foreach(a; dataset) {
        auto intDE = makeDEInt!real(a, real.infinity).integrate((real x) => exp(-x^^2));
        intDE *= M_2_SQRTPI;

        auto truevalue = erfc(a);
        assert(approxEqual(intDE, truevalue));
    }
}

// Complementary error function: erfc(a) = -2/sqrt(pi) * int_(-inf)^(-a) e^{-x^2} dx
unittest 
{
    import std.mathspecial;

    real[] dataset = [
        0.0,
        0.01,
        0.1,
        1,
        10,
        100,
    ];

    foreach(a; dataset) {
        auto intDE = makeDEInt!real(-real.infinity, -a).integrate((real x) => exp(-x^^2));
        intDE *= M_2_SQRTPI;

        auto truevalue = erfc(a);
        assert(approxEqual(intDE, truevalue));
    }
}

// int_inf^inf 1/(1+x^2) dx = pi
unittest
{
    auto intDE1 = makeDEInt!real(-real.infinity, real.infinity)
                    .integrate((real x) => 1/(1 + x^^2));
    assert(approxEqual(intDE1, PI));
}

// int exp(-x^2) = sqrt(pi)
unittest
{
    immutable val1 = makeDEInt!real(0, real.infinity, Yes.isExpDecay, 100)
        .integrate((real x) => 1.0/(2*sqrt(x)) * exp(-x) );

    immutable val2 = makeDEInt!real(real.infinity, 0, Yes.isExpDecay, 100)
        .integrate((real x) => -1.0/(2*sqrt(x)) * exp(-x) );

    immutable val3 = makeDEInt!real(-real.infinity, 0, Yes.isExpDecay, 100)
        .integrate((real x) => 1.0/(2*sqrt(-x)) * exp(x) );

    immutable val4 = makeDEInt!real(0, -real.infinity, Yes.isExpDecay, 100)
        .integrate((real x) => -1.0/(2*sqrt(-x)) * exp(x) );

    assert(approxEqual(val1, sqrt(PI)/2));
    assert(approxEqual(val2, sqrt(PI)/2));
    assert(approxEqual(val3, sqrt(PI)/2));
    assert(approxEqual(val4, sqrt(PI)/2));
}

// int exp(-x^2) * 1i = 1i * sqrt(pi)
unittest
{
    import std.complex : Complex;

    immutable val1 = makeDEInt!real(0, real.infinity, Yes.isExpDecay)
        .integrate((real x) => 1.0/(2*sqrt(x)) * exp(-x) * Complex!real(0, 1));

    assert(val1.re.approxEqual(0));
    assert(val1.im.approxEqual(sqrt(PI)/2));
}

// unittest for README.md
unittest
{
    import std.mathspecial;
    import deint.utils;

    // integrate f(x)exp(-x) on (1, inf)
    // Now, we know that the integrand f(x)exp(-x) decay exponentially.
    auto intFI = makeDEInt!real(1, real.infinity, Yes.isExpDecay);

    // incomplete gamma function
    assert(intFI.integrate((real x) => x * exp(-x)).approxEqual(gammaIncompleteCompl(2, 1) * gamma(2)));

    // Also, you can use `withWeight` which pre-computes and stores weight exp(-x).
    // The `withWeight` is useful when integrands have same weights.
    auto intFIW = intFI.withWeight((real x) => exp(-x));

    // incomplete gamma function
    assert(intFIW.integrate((real x) => x).approxEqual(gammaIncompleteCompl(2, 1) * gamma(2)));
    assert(intFIW.integrate((real x) => x^^2).approxEqual(gammaIncompleteCompl(3, 1) * gamma(3)));
    assert(intFIW.integrate((real x) => x^^3).approxEqual(gammaIncompleteCompl(4, 1) * gamma(4)));
}

unittest
{
    // integration on [-inf, inf]
    auto intII = makeDEInt!real(-real.infinity, real.infinity);

    // Gaussian integral
    assert(intII.integrate((real x) => exp(-x^^2)).approxEqual(sqrt(PI)));
}

unittest
{
    // DEInt!real is a struct which computes x_k and w_k in advance.
    auto int01 = makeDEInt!real(0, 1);

    // When f(x) = x, int_0^1 x dx = 0.5
    assert(int01.integrate((real x) => x).approxEqual(0.5));

    // `DEInt!real` is reusable.
    // When f(x) = x^^2, int_0^1 x^^2 dx = 1/3
    assert(int01.integrate((real x) => x^^2).approxEqual(1/3.0));
}