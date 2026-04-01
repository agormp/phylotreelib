# Changelog

---

## [2.2.0] - 2026 April 1

#### Added
- Construction of random trees using Tree.randtree now does proper simulation of either Yule or coalescent trees, so tree structure is more realistic. Default is to create ultrametric tree with exactly the user-specified tree height. Setting the argument rate_sd > 0 will add noise to branch lengths and break ultrametricity. Specifically all branch lengths in the ultrametric tree will be multiplied by a factor X drawn from a lognormal with E[X]=1 and sd=rate_sd.

#### Fixed
- Now consistently uses space, not mixture of space and tabs, to indent elements in NEXUS format

---

## [2.0.0 - 2.2.0]: under construction