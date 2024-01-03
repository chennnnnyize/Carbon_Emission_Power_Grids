# Carbon_Emission_Power_Grids

[![license](https://img.shields.io/github/license/InternLM/lagent.svg)](https://github.com/chennnnnyize/Carbon_Emission_Power_Grids/blob/main/LICENSE.txt)

This is the code repo for the paper "[Contributions of Individual Generators to Nodal Carbon Emissions](https://arxiv.org/abs/2311.03712)"

Authors: Yize Chen, Deepjyoti Deka, Yuanyuan Shi

Accepted at the *15th ACM International Conference on Future and Sustainable Energy Systems (ACM e-Energy 2024)*.

Hong Kong University of Science and Technology, Los Alamos National Lab, University of California San Diego.

**Summary**: An algorithm to calculate each node's exact carbon emission rate via generation mix.

There are growing interests on carbon accounting and clean energy technologies. In the electricity power networks, each
generator can have distinct carbon emission rates. Due to the existence of physical power flows, nodal power consumption is met
by a combination of a set of generators, while such combination is
determined by network topology, generators’ and lines' characteristics, and
power demand. 

This work describes a technique based on physical
power flow model, which can efficiently compute the nodal carbon
emissions contributed by each single generator given the generation and power flow information. 

We also extend the technique to calculate both the nodal average carbon emission and marginal
carbon emission rates. 

<p align="center">
<img width="439" alt="image" src="https://github.com/chennnnnyize/Carbon_Emission_Power_Grids/assets/116547738/7cf3ea3c-2adf-4e1c-a23f-10c77e199fc8">
</p>

## Getting Started
```bash
python marginal_carbon.py
```

**Contact**: yizehchen@lbl.gov
