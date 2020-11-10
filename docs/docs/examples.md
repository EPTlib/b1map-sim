---
layout: default
title: Examples
nav_order: 2
description:
permalink: /examples
last_modified_date: 2020-07-30T10:44:37+02:00
---

# Examples
{: .no_toc .fs-9}

List of examples for each implemented B1-mapping method.
{: .fs-6 .fw-300 }

## Table of contents
{: .no_toc .text-delta }

1. TOC
{:toc}

---

The following examples are based on data simulated of a homogeneous phantom imaged with a birdcage body coil at 64 MHz, that is the Larmor frequency of a 1.5 T scanner.

The .h5 files used in the examples are [homogeneous-phantom-15t.h5](/assets/examples/homogeneous-phantom-15t.h5), containing the geometry of the phantom and the simulated transmit and receive sensitivities, and [homogeneous-phantom-15t-prop.h5](/assets/examples/homogeneous-phantom-15t-prop.h5), containing the physical properties of the phantom.

The description of the imaged body is provided by the .toml file that can be downloaded [here](/assets/examples/phantom.toml).

```toml
{% include_relative phantom.toml %}
```

## Double angle

The configuration file for the double angle example can be downloaded [here](/assets/examples/da.toml).

```toml
{% include_relative da.toml %}
```

The intermediate images and the actual and the estimated flip-angle are reported in the following picture. The results are obtained by applying the double angle method. Before combining the two intermediate images to obtain the flip-angle estimate, noise with a signal to noise ratio (SNR) equal to 200 is added.

![](/b1map-sim/assets/images/b1map-da.png){: style="width:384px"}

## Actual flip-angle

The configuration file for the actual flip-angle example can be downloaded [here](/assets/examples/afi.toml).

```toml
{% include_relative afi.toml %}
```

The intermediate images and the actual and the estimated flip-angle are reported in the following picture. The results are obtained by applying the actual flip-angle method. Before combining the two intermediate images to obtain the flip-angle estimate, noise with a signal to noise ratio (SNR) equal to 200 is added.

![](/b1map-sim/assets/images/b1map-afi.png){: style="width:384px"}

## Bloch--Siegert shift

The configuration file for the Bloch--Siegert shift example can be downloaded [here](/assets/examples/bss.toml).

```toml
{% include_relative bss.toml %}
```

The intermediate images and the actual and the estimated transmit sensitivity magnitude are reported in the following picture. The results are obtained by applying the Bloch--Siegert shift method. Before combining the two intermediate images to obtain the transmit sensitivity magnitude estimate, noise with a signal to noise ratio (SNR) equal to 200 is added.

![](/b1map-sim/assets/images/b1map-bss.png){: style="width:384px"}
