<p align="center">
  <img src="https://i.imgur.com/Xzvloc7.png" alt="logo" width="500"/>
</p>

> ⚠️ The plugin is in active development and therefore subject to rapid change. Things may change suddenly and without notice, so please be aware that new commits could possibly render a previously working environment unusable. ⚠️

# DS::compare
![License](https://img.shields.io/badge/License-GPLv3-blue.svg)
[![Docker](https://badgen.net/badge/icon/docker?icon=docker&label)](https://https://docker.com/)
[![Generic badge](https://img.shields.io/badge/Uses-C++-<COLOR>.svg)](https://shields.io/)
[![Generic badge](https://img.shields.io/badge/Mitsuba-0.6-<COLOR>.svg)](https://shields.io/)

**DS::compare** is a utility plugin written for the [Mitsuba Renderer](https://github.com/mitsuba-renderer/mitsuba) to assess the effectivity of data structures in a path guiding context. To this end, this plugin samples spherical environment maps, stores the samples in a data structure of choice, and uses this information to approximate a guiding distribution. Using various metrics, the quality of a reconstructed environment map is then analyzed, allowing us to draw conclusions about the respective data structure.
<br>
<br>
Please keep in mind that this plugin only assesses the directional component of the path guiding process.

## Getting Started

### Approach 1: Docker

### Approach 2: Native Environment

## Plugin Structure

### Overview
The plugin system is structured as follows:<br>
`src/utils/dscompare.cpp` - Main plugin file<br>
`include/ds/*` - Plugin headers (Utility, data structures, etc.)<br>
`src/libcore/ds/*` - Data structure impl.<br>

## How to Use

### General Usage

### Adding Environment Maps

### Adding a Data Structure
As the plugin itself uses an intrinsic plugin system for data structures, adding a data structure to the already existing ones is fairly easy.

## License
DS::compare is available under the [GNU GLPv3 license](https://www.gnu.org/licenses/gpl-3.0.html). See [LICENSE](https://github.com/GitThirteen/mitsuba-comparer/blob/main/LICENSE) for the full license text.
