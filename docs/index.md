---
layout: default
title: Home
nav_order: 1
description:
permalink: /
last_modified_date: 2020-11-06T14:56:00+02:00
---

# ```b1map-sim```
{: .fs-9 }

```b1map-sim``` is a generator of virtual B1-mapping images.
{: .fs-6 .fw-300 }

[Get started now](#getting-started){: .btn .btn-primary .fs-5 .mb-4 .mb-md-0 .mr-2 } [View it on GitHub](https://github.com/eptlib/b1map-sim){: .btn .fs-5 .mb-4 .mb-md-0 target="_blank" }

---

## Getting started

The binaries of ```b1map-sim``` for different operative systems can be downloaded [here](https://github.com/EPTlib/b1map-sim/releases).

Alternatively, ```b1map-sim``` can be build directly from its [source](https://github.com/EPTlib/b1map-sim).

```b1map-sim``` is a terminal application that can be run by typing

```bash
b1map-sim <filename>.toml
```

where ```<filename>``` must be replaced with the name of the configuration file for the case of interest.

Details on the set-up of the configuration file are provided in the [Settings](settings) page.

## List of implemented B1-mapping methods

The following methods are currently implemented in ```b1map-sim``` and can be executed from the application:

1. Double angle
1. Actual flip-angle
1. Bloch-Siegert shift

---

## About the project

```b1map-sim``` is &copy; {{ "now" | date: "%Y" }} by [Alessandro Arduino](http://github.com/alessandroarduino), INRiM.

### License

```b1map-sim``` is distributed by an [MIT license](https://github.com/eptlib/b1map-sim/tree/master/LICENSE).
