module deint.utils;


import std.math;
import std.traits;
import std.typecons;

import deint.deint;


/** This function create a instance of NumInt from original another and weight function.
 * This function is useful when we know integrands has same weight.
 */
auto withWeight(Int, Fn)(auto ref Int integ, Fn weightFn)
{
    alias X = Unqual!(typeof(integ.xs[0]));
    alias W = Unqual!(typeof( integ.ws[0] * weightFn(X.init) ));

    auto xs = integ.xs;
    auto ws = integ.ws;
    // return IntWithWeight!(X, W)(integ.xs, integ.ws, weightFn);
    immutable(W)[] newws;
    if (!__ctfe)
        newws.reserve(xs.length);
    foreach(i, x; xs)
        newws ~= weightFn(x) * ws[i];

    return NumInt!(X, W)(xs, newws);
}

///
unittest
{
    // integration on [0, inf] by DE formula
    auto deint = makeDEInt!real(0, real.infinity, Yes.isExpDecay);

    // integration on [0, inf] by DE formula with a weight exp(-x).
    auto weighted = deint.withWeight((real x) => exp(-x));

    assert(deint.integrate((real x) => 1.0/(2*sqrt(x)) * exp(-x) )
        .approxEqual(weighted.integrate((real x) => 1.0/(2*sqrt(x)) )));
}


/** This function creates a instance of NumInt for integration by the trapezoidal rule.
 *
 * Reference:
 * $(HTTP en.wikipedia.org/wiki/Trapezoidal_rule)
 */
NumInt!(F, F) makeTrapInt(F)(F xa, F xb, size_t N, Flag!"isPeriodic" isPeriodic = No.isPeriodic)
{
    immutable(F)[] xs;
    immutable(F)[] ws;
    if (!__ctfe)
    {
        xs.reserve(N);
        ws.reserve(N);
    }

    if(!isPeriodic) {
        F h = (xb - xa) / (N-1);
        foreach(i; 0 .. N) {
            xs ~= h*i + xa;
            ws ~= h * (i == 0 || i == N-1 ? 0.5 : 1);
        }
    } else {
        F h = (xb - xa) / N;
        foreach(i; 0 .. N) {
            xs ~= h*i + xa;
            ws ~= h;
        }
    }

    return NumInt!(F, F)(xs, ws);
}

///
unittest
{
    auto trap = makeTrapInt!real(0, 5, 100);
    assert(trap.integrate((real x) => exp(-x^^2))
        .approxEqual(sqrt(PI)/2));
}
