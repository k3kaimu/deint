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
        if (!__ctfe)
        {
            xs.reserve(trapN);
            ws.reserve(trapN);
        }
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
        if(xa == 0 || xb == 0){
            return _makeParamsImpl(delegate(F t){
                immutable F
                    cosht = cosh(t),
                    sinht = sinh(t),
                    epsinht = exp(PI * sinht),
                    x = (xa + xb*epsinht)/(1 + epsinht),
                    cosh2 = cosh(PI / 2 * sinht)^^2,
                    dx = PI / 2 * cosht / cosh2;

                return cast(F[2])[x, dx * (xb - xa)/2];
            });
        }else{
            immutable F diff2 = (xb - xa)/2;
            immutable F avg2 = (xb + xa)/2;

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

    // int_0^1 log(x)/sqrt(x) dx = -4
    assert(int01.integrate((real x) => log(x)/sqrt(x)).approxEqual(-4, 1E-12, 1E-12));

    // int_{-1}^0 log(-x)/sqrt(-x) dx = -4
    assert(makeDEInt!real(-1, 0).integrate((real x) => log(-x)/sqrt(-x)).approxEqual(-4, 1E-12, 1E-12));

    // integration on [-1, 1]
    auto int11 = makeDEInt!real(-1, 1);

    // int_{-1}^1 1/sqrt(1-x^2) dx = pi
    assert(int11.integrate((real x) => 1/(1 + x^^2)).approxEqual(PI/2, 1E-12, 1E-12));

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



/** This function create a instance of NumInt!(F, F) for computing fourtier-type integration by the DE formula.
 * On the DE formula for the fourier-type integration, the integration (1) is converted to (2).
 * (1) int_0^{inf} g(x) sin(omega * x) dx
 * (2) sum_{k=-Nlow}^{Nhigh} g(M*phi(h*k)) sin(M*phi(h*k)/omega) M*phi'(h*k)/omega * h
 * 
 * Params:
 *      isSine = `Yes` if the type of integration is "sine". `No` if the type is "cosine".
 *      omega = angular frequency
 *      stepH = step value of each interval
 *      Nlow = number of negative-valued computing points
 *      Nhigh = number of positive-valued computing points
 *
 * Reference:
 * $(HTTP www.kurims.kyoto-u.ac.jp/~kyodo/kokyuroku/contents/pdf/1040-20.pdf)
*/
NumInt!(F, F) makeDEIntFourier(F)(
    Flag!"isSine" isSine,
    F omega,
    F stepH = 0.1,
    size_t Nlow = 49, size_t Nhigh = 50)
{
    immutable F beta = 1.0L/4;
    immutable F h = stepH;
    immutable F M = PI / h;
    immutable F alpha = beta / sqrt(1 + M*log(1 + M)/(4*PI));

    F[2] phi(F t) {
        F num = t;
        F den = 1 - exp(-2*t - alpha*(1 - exp(-t)) - beta*(exp(t) - 1));
        F dden = (den - 1) * (-2 - alpha*exp(-t) - beta * exp(t));

        F num2 = den - num * dden;
        F den2 = den^^2;

        if(abs(t) < 1E-3) {
            F ddden = (den-1) * (-2 - alpha*exp(-t) - beta * exp(t))^^2
                    + (den-1) * (   + alpha*exp(-t) - beta * exp(t));

            num = 1;
            den = dden;

            num2 = dden - (dden + num * ddden);
            den2 = 2*den*dden;
        }

        return [num/den, num2/den2];
    }

    immutable(F)[] xs;
    immutable(F)[] ws;
    if (!__ctfe)
    {
        xs.reserve(Nlow + Nhigh);
        ws.reserve(Nlow + Nhigh);
    }

    foreach(long i; -cast(long)Nlow .. cast(long)Nhigh) {
        F t = h * i;
        if(isSine) {
            auto p = phi(t);
            xs ~= M * p[0] / omega;
            ws ~= M * p[1] / omega * h;
        } else {
            auto p = phi(t - PI/(2*M));
            xs ~= M * p[0] / omega;
            ws ~= M * p[1] / omega * h;
        }
    }

    return NumInt!(F, F)(xs, ws);
}

unittest {
    auto de = makeDEIntFourier!real(Yes.isSine, 1, 0.2);
    auto dirichlet = de.integrate((real x) => sin(x)/x);

    assert(dirichlet.approxEqual(PI/2));
}
