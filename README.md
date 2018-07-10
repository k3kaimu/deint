# Double Exponential Numerical Integration Library

This library provide numerical integration method by Double Exponential (DE) formula.

For example:

```
// integration on [0, 1]
auto int01 = DEInt!real(0, 1);

// int_0^1 x dx = 0.5
assert(int01.integrate((real x) => x).approxEqual(0.5));

// int_0^1 x^^2 dx = 1/3
assert(int01.integrate((real x) => x^^2).approxEqual(1/3.0));


// integration on [-inf, inf]
auto intII = DEInt!real(-real.infinity, real.infinity);

// Gaussian integral
assert(intII.integrate((real x) => exp(-x^^2)).approxEqual(sqrt(PI)));

// integration int_1^inf x * exp(-x) dx = Gamma(2, 1)
auto intFI = DEInt!real(1, real.infinity, Yes.isExpType);

// incomplete gamma function
assert(intFI.integrate((real x) => x).approxEqual(gammaIncompleteCompl(2, 1) * gamma(2)));
```
