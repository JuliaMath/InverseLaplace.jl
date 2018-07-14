# InverseLaplace

*Numerical inverse Laplace transform*

[![](https://img.shields.io/badge/docs-latest-blue.svg)](https://jlapeyre.github.io/InverseLaplace.jl/latest)
[![](https://img.shields.io/badge/docs-stable-blue.svg)](https://jlapeyre.github.io/InverseLaplace.jl/stable)
Linux, OSX: [![Build Status](https://travis-ci.org/jlapeyre/InverseLaplace.jl.svg)](https://travis-ci.org/jlapeyre/InverseLaplace.jl) &nbsp; Windows: [![Build Status](https://ci.appveyor.com/api/projects/status/github/jlapeyre/InverseLaplace.jl?branch=master&svg=true)](https://ci.appveyor.com/project/jlapeyre/inverselaplace-jl) &nbsp; &nbsp; &nbsp;
[![Coverage Status](https://coveralls.io/repos/github/jlapeyre/InverseLaplace.jl/badge.svg?branch=master)](https://coveralls.io/github/jlapeyre/InverseLaplace.jl?branch=master)
[![codecov](https://codecov.io/gh/jlapeyre/InverseLaplace.jl/branch/master/graph/badge.svg)](https://codecov.io/gh/jlapeyre/InverseLaplace.jl)

This implements some numerical methods for computing inverse Laplace transforms in Julia.

InverseLaplace v0.1.0 is the last version that supports Julia v0.6

See the documentation https://jlapeyre.github.io/InverseLaplace.jl/latest .

Note: A small part of `InverseLaplace.jl` depends on `Optim.jl`, which is currently
broken on Julia v0.7. So optimization is disabled for `InverseLaplace.jl` versions
greater than v0.1.0

Note: the last version of this module supporting Julia v0.4 is tagged v0.0.2.
the last version of this module supporting Julia v0.6 is tagged v0.1.0.
