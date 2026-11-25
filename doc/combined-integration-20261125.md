# Combined LIS integration snapshot (2026-11-25)

Base: upstream develop b58dc86ff4d91e2b7d2d3cb7546ceec25187463e.

This isolated draft combines source and regression changes from upstream #259 (includes #256), #257, #258 and fork #170–178. It is not a request to merge those upstream PRs together. Original PR branches are unchanged.

## Build contract

Run `autoreconf -fiv` before configure: the merged test/Makefile.am is authoritative; generated test/Makefile.in intentionally remains at the base version in this experimental snapshot until it can be regenerated and checked with Automake. Do not distribute this snapshot as a ready-to-configure release tarball. No proposed workflow additions from the component PRs are included; existing base CI workflows are untouched. For extended regressions, use the existing manual/develop workflow or run tests directly in a fresh local build.

Suggested smoke build: `autoreconf -fiv && ./configure --enable-test --disable-fortran --disable-f90 && make -j2 && make check`. Repeat with OpenMP, 64-bit LIS_INT and sanitizer/OOM builds as appropriate using the existing project options. The native Krylov API remains int-indexed. No combined build result is asserted by this snapshot commit.

## Included source heads

- 259: 395731244e8db50c40729de18d9b7991ac6d2efb
- 175+257: 9fc727085afe3289fffc1e20383d67bb852469d7
- 258: 963956876799e47b78b77a725e55da0196dbda56
- 170+172+173+174: 388519c63f977b4041ff96e6fce102c643a9eb77
- 178: 2c40c45dddae5871a851aac9df16020cd2676023
- 171: ffed584e6d9b1fba85569fe30bd9aa944187e0f6
- 176: 52824616afb79e56113494a4481e5aa9835119f2
- 177: f28e8604996c9d0a3e0156c6eeff35c07e32af98
