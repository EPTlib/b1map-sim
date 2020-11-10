---
layout: default
title: Settings
nav_order: 2
description:
permalink: /settings
last_modified_date: 2020-11-06T15:16:00+02:00
---

# Configuration file settings
{: .no_toc .fs-9}

```b1map-sim``` is configured by a .toml file.
{: .fs-6 .fw-300 }

## Table of contents
{: .no_toc .text-delta }

1. TOC
{: toc}

---

In .toml files, the symbol ```#``` marks the rest of the line as a comment, except when inside a string.

## Title and method

```toml
title = "Example"
method = 0
```

- The ```title``` will appear in the log of the application.
- The ```method``` is selected according to the following table.

| Code | Method              |
|------+---------------------|
| 0    | Double angle        |
| 1    | Actual flip-angle   |
| 2    | Bloch--Siegert shift |

## Mesh

```toml
[mesh]
    size = [128,128,8]
    step = [1e-3,1e-3,4e-3] # [m]
```

- ```size``` is the number of voxels of the input data in each direction.
- ```step``` is the size of a voxel in meters.

## Input

```toml
[input]
    body = "body.toml"
    tx-sensitivity = "phantom.h5:/tx_sens" # [T]
    tx-phase = "phantom.h5:/tx_phase" # [rad]
```

- ```body``` is the address of the .toml file containing the parameters of the body. It can be the configuration file itself or another file (useful if the parameters are shared by more cases).
- ```tx-sensitivity``` is the address of the transmit sensitivity (magnitude) in tesla. It must be a dataset in an .h5 file.
- ```tx-phase``` is the address in an .h5 file of the transmit phase in radians. It must be a dataset in an .h5 file.

## Body

```toml
[body]
    materials = "phantom.h5:/materials"
    proton-density = "phantom.h5:/rho"
    longitudinal-relaxation = "phantom.h5:/T1" # [ms]
    transverse-relaxation = "phantom.h5:/T2" # [ms]
```

- ```materials``` is the address of the distribution of materials within the body. It must be a dataset in an .h5 file.
- ```proton-density``` is the address of the list of the proton density of each material in the body. It must be a dataset in an .h5 file.
- ```longitudinal-relaxation``` is the address of the list of the longitudinal relaxations (T1) expressed in millisecond of each material in the body. It must be a dataset in an .h5 file.
- ```transverse-relaxation``` is the address of the list of the transverse relaxations (T2) expressed in millisecond of each material in the body. It must be a dataset in an .h5 file.

## Output

```toml
[output]
    alpha-estimate = "example.h5:/alpha"
    intermediate-images = "example.h5:/imgs"
```

- ```alpha-estimate``` is the address where the estimated flip-angle expressed in radian (or B1+ magnitude expressed in tesla, in the case of Bloch--Siegert shift) will be written. It must be a dataset in an .h5 file.
- ```intermediate-images``` is the root of the address where the intermediate images of the b1-mapping procedure will be written. It must be a dataset in an .h5 file.

Given the above example, the real and imaginary parts of the two intermediate images would be stored in the four datasets: \\
```example.h5:/imgs1/real``` ```example.h5:/imgs1/imag``` \\
```example.h5:/imgs2/real``` ```example.h5:/imgs2/imag```

## Parameters

```toml
[parameter]
    alpha-nominal = 1.04 # [rad]
    TR = 3000.0 # [ms]
    TE = 20.0 # [ms]
    spoiling = 1.0
```

- ```alpha-nominal``` is the nominal flip-angle that would be obtained if the transmit sensitivity was homogeneous.
- ```TR``` is the repetition time of the pulse sequence.
- ```TE``` is the echo time of the pulse sequence.
- ```spoiling``` is a coefficient for transverse magnetization spoiling: 1 is ideal spoiling, whereas 0 is no spoiling.

In addition, the following parameters are specific for the adopted B1-mapping method.

### Double angle (0)

No additional parameters are needed.

### Actual flip-angle (1)

```toml
[parameter]
    TRratio = 5.0
```

- ```TRratio``` is the ratio between the TR of the two intermediate images.

### Bloch--Siegert shift (2)

```toml
[parameter]
    bss-offres = 4.0 # [kHz]
    bss-length = 4.0 # [ms]
```

- ```bss-offres``` is the off-resonance frequency of the Bloch--Siegert pulse expressed in kilohertz.
- ```bss-length``` is the length of the Bloch--Siegert pulse expressed in millisecond.

The implemented Bloch--Siegert shift method assumes that the Bloch--Siegert pulse has only the selected off-resonance frequency. It is a good approximation of the common Fermi pulse.

## Monte Carlo

```toml
[montecarlo]
    samples = 10
    noise = 0.01
```

- ```samples``` is the number of Monte Carlo samples.
- ```noise``` is the inverse of the average SNR of the intermediate images.

This section is optional. If it is present, then a number of noisy output (the Monte Carlo samples) are generated in addition to the noiseless ones.
The value of ```noise``` must be greater than zero, otherwise this section will be ignored.

Given ```output.alpha-estimate = "example.h5:/alpha"```, then the 10 samples are stored in \\
```example.h5:/alpha0``` ```example.h5:/alpha1``` \\
```example.h5:/alpha2``` ```example.h5:/alpha3``` \\
```example.h5:/alpha4``` ```example.h5:/alpha5``` \\
```example.h5:/alpha6``` ```example.h5:/alpha7``` \\
```example.h5:/alpha8``` ```example.h5:/alpha9```
