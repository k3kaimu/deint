// Written in the D programming language.

/**
This module provides numerical integration by double-exponential (DE) formula.
*/
module deint;

import std.algorithm;
import std.array;
import std.math;
import std.typecons;



/** A struct for performing numerical integration by double-exponential (DE) formula.
 * It is also known as "Tanh-sinh quadrature".
 * In DE formula, the integration (1) is converted to (2).
 * (1) int_{xa}^{xb} f(x) dx
 * (2) int_{ta}^{tb} f(g(t)) g'(t) dt
 * 
 * This struct calculates reusable state in advance.
 * The type of DE formula is automatically decided from the given interval of the integration.
 *
 * Reference:
 * $(HTTP en.wikipedia.org/wiki/Tanh-sinh_quadrature)
 */
struct DEInt(F)
{
    /** Initialize an object for computing DE integration.
     * Params:
     * 		xa = starting value of original integration.
     *      xb = end value of original integration.
     *      isExpType = if the integration is formed as int_a^b f(x) exp(-x) dx, this value is Yes. otherwise No.
     *      trapN = division points of trapezoidal quadrature. If this value is 0, 100 is used internally.
     *      ta = starting value of integration transformed by DE-formula. If this value is NaN, -5 is used internally.
     *      tb = starting value of integration transformed by DE-formula. If this value is NaN, +5 is used internally.
     */
    this(
        F xa, F xb,
        Flag!"isExpType" isExpType = No.isExpType,
        size_t trapN = 0,
        F ta = F.nan, F tb = F.nan)
    {
        setParams(xa, xb, isExpType, trapN, ta, tb);
    }


    /// ditto
    void setParams(
        F xa, F xb,
        Flag!"isExpType" isExpType = No.isExpType,
        size_t trapN = 0,
        F ta = F.nan, F tb = F.nan)
    in{
        if(isExpType)
            assert(_xa != -F.infinity);
    }
    do {
        if(trapN == 0) {
            trapN = 100;
        }

        if(ta.isNaN && tb.isNaN) {
            ta = -5;
            tb = 5;
        }

        if(xa > xb) {
            this.setParams(xb, xa, isExpType, trapN, ta, tb);
            _weights = _weights.map!"cast(immutable)(-a)".array;
        } else if(xa == -F.infinity && xb != F.infinity) {
            assert(!isExpType);
            this.setParams(-xb, F.infinity, isExpType, trapN, ta, tb);
            _xs = _xs.map!"cast(immutable)(-a)".array;
        }

        _xa = xa;
        _xb = xb;
        _ta = ta;
        _tb = tb;
        _trapN = trapN;
        _isExpType = cast(bool)isExpType;

        if(_xs !is null)
            return;


        if(xa == -F.infinity && xb == F.infinity){
            assert(!isExpType);

            auto params = _makeParamsImpl(ta, tb, trapN, delegate(F t){
                immutable F
                    sinht = sinh(t),
                    x = sinh(PI / 2 * sinht),
                    dx = cosh(PI / 2 * sinht) * PI / 2 * cosh(t);

                if(!isExpType)
                    return cast(F[2])[x, dx];
                else
                    return cast(F[2])[x, exp(-x) * dx];
            });


            _xs = params[0];
            _weights = params[1];
            _intType = isExpType ? "IIE" : "II";
        }else if(xb == F.infinity) {
            auto params = _makeParamsImpl(ta, tb, trapN, delegate(F t){
                if(!isExpType){
                    real x = exp(PI / 2 * sinh(t)),
                         dx = x * PI / 2 * cosh(t);

                    return cast(F[2])[x + xa, dx];
                }else{
                    real expmt = exp(-t),
                         x = exp(t - expmt),
                         dx = (1 + expmt) * x;


                    return cast(F[2])[x + xa, exp(-(x + xa)) * dx]; 
                }
            });

            _xs = params[0];
            _weights = params[1];
            _intType = isExpType ? "FIE" : "FI";
        }else{
            immutable F diff2 = (xb - xa) / 2,
                        avg2 = (xb + xa) / 2;

            auto params = _makeParamsImpl(ta, tb, trapN, delegate(F t){
                immutable F
                    cosht = cosh(t),
                    sinht = sinh(t),
                    x = tanh(PI / 2 * sinht) * diff2 + avg2,
                    cosh2 = cosh(PI / 2 * sinht)^^2,
                    dx = PI / 2 * cosht / cosh2;

                if(!isExpType)
                    return cast(F[2])[x, dx * diff2];
                else
                    return cast(F[2])[x, exp(-x) * dx * diff2];
            });

            _xs = params[0];
            _weights = params[1];
            _intType = isExpType ? "FFE" : "FF";
        }
    }



    /** Execute integration of func by DE formula.
    */
    F integrate(Fn)(scope Fn func) const
    {
        F sum = 0;
        foreach(i; 0 .. _xs.length)
            sum += func(_xs[i]) * _weights[i];

        return sum;
    }


  @property const
  {
  	/** Return type of integration.
  	 * The First and second charactor of the return value are explain the the starting value or end value of integration is finite value ('F') or infinity 'I'.
  	 * If the return value has the third character and its value is 'E', the integration is formed as int_{xa}^{xb} f(x) exp(-x) dx 
  	 */
    string type() { return _intType; }

    /// Return the set value
    real xa() { return _xa; }

    /// ditto
    real xb() { return _xb; }

    /// ditto
    real ta() { return _ta; }

    /// ditto
    real tb() { return _tb; }

    /// ditto
    size_t trapN() { return _trapN; }

    /// ditto
    bool isExpType() { return _isExpType; }

    /// Return division points (computing points) of trapezoidal quadrature for DE formula.
    immutable(F)[] xs() { return _xs; }

    /// Return weights of each division points.
    immutable(F)[] weights() { return _weights; }
  }


  private:
    string _intType;
    real _xa, _xb, _ta, _tb;
    size_t _trapN;
    bool _isExpType;
    immutable(F)[] _xs;
    immutable(F)[] _weights;

    static
    immutable(F)[][2] _makeParamsImpl(F ta, F tb, size_t trapN, scope F[2] delegate(F) fn)
    {
        immutable(F)[] xs, weights;
        immutable F h = (tb - ta) / (trapN-1);
        foreach(i; 0 .. trapN) {
            immutable xWt = fn(i * h + ta);
            xs ~= xWt[0];
            weights ~= xWt[1] * h;
        }

        return [xs, weights];
    }
}

///
unittest
{
	// integration on [0, 1]
	auto int01 = DEInt!real(0, 1);
	assert(int01.type == "FF");

	// int_0^1 x dx = 0.5
	assert(int01.integrate((real x) => x).approxEqual(0.5));

	// int_0^1 x^^2 dx = 1/3
	assert(int01.integrate((real x) => x^^2).approxEqual(1/3.0));


	// integration on [-inf, inf]
	auto intII = DEInt!real(-real.infinity, real.infinity);
	assert(intII.type == "II");

	// Gaussian integral
	assert(intII.integrate((real x) => exp(-x^^2)).approxEqual(sqrt(PI)));

	import std.mathspecial;
	// integration int_1^inf x * exp(-x) dx = Gamma(2, 1)
	auto intFI = DEInt!real(1, real.infinity, Yes.isExpType);
	assert(intFI.type == "FIE");

	// incomplete gamma function
	assert(intFI.integrate((real x) => x).approxEqual(gammaIncompleteCompl(2, 1) * gamma(2)));
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
        real intDE = DEInt!real(data[0], data[1]).integrate((real x) => exp(x));
        real truevalue = exp(data[1]) - exp(data[0]);

        assert(approxEqual(intDE, truevalue));
    }


    // reverse interval
    foreach(data; dataset) {
        //real intDE = intAToB_DE(data[0], data[1], 100, (real x) => exp(x));
        real intDE = DEInt!real(data[1], data[0]).integrate((real x) => exp(x));
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
        auto intDE = DEInt!real(a, real.infinity).integrate((real x) => exp(-x^^2));
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
        auto intDE = DEInt!real(-real.infinity, -a).integrate((real x) => exp(-x^^2));
        intDE *= M_2_SQRTPI;

        auto truevalue = erfc(a);
        assert(approxEqual(intDE, truevalue));
    }
}

// int_inf^inf 1/(1+x^2) dx = pi
unittest
{
    auto intDE1 = DEInt!real(-real.infinity, real.infinity)
                    .integrate((real x) => 1/(1 + x^^2));
    assert(approxEqual(intDE1, PI));
}

// int exp(-x^2) = sqrt(pi)
unittest
{
    immutable val = DEInt!real(0, real.infinity, Yes.isExpType, 100)
        .integrate((real x) => 1.0/(2*sqrt(x)));

    assert(approxEqual(val, sqrt(PI)/2));
}
