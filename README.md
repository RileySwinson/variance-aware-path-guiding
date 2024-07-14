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
Please keep in mind that this plugin only considers the directional component of the path guiding process.

## Getting Started

### Approach 1: Docker

**Step 1:** Install docker on your host device and create an image from a pre-existing Mitsuba 0.6 dockerfile. The dockerfile used for this project can be found [here](https://github.com/xehoth/mitsuba-docker).
Assuming the command is run from within the folder with the dockerfile:
```
docker build -t <your-image-name> .
```

**Step 2:** Clone this project:
```
git clone https://github.com/GitThirteen/mitsuba-comparer.git
```

**Step 3:** Attach the project on your local device to a fresh docker container via a [*bind mount*](https://docs.docker.com/storage/bind-mounts/):
```
docker run -it --name <your-container-name> -v <host-path>:/home/mitsuba <your-image-name>
```
The `-v <host-path>:<container-path>` flag specifies where Docker should mount the local host directory inside the container. For more information please refer to the Docker documentation. Alternatively, a [*docker volume*](https://docs.docker.com/storage/volumes/) may be used instead. (Untested.)

**Step 4:** At this point, development should be possible by simply running
```
docker start <container-name>
docker attach <container-name>
```
It is highly suggested to develop locally and only use the console attached to the Docker container for running plugin commands, as the syncing between Docker and local environment slows down things like IntelliSense, linters, and other tools immensely. It is also suggested to **not** run this plugin via a bind mount for <ins>production purposes</ins> (e.g., scientific evaluations).

For instructions on how to run the plugin itself, please refer to the [**How to Use**](#how-to-use) section.

### Approach 2: Native Environment

Compiling Mitsuba manually is possible, however quite cumbersome due to outdated libraries and incompatibilities. Please refer to the [Mitsuba documentation](https://www.mitsuba-renderer.org/releases/current/documentation.pdf) to find out how to build Mitsuba on your own OS.

## Plugin Structure

### Overview
The plugin system is structured as follows:<br>
`src/utils/dscompare.cpp` - Main plugin file<br>
`include/ds/*` - Plugin headers (Utility, data structures, etc.)<br>
`src/libcore/ds/*` - Data structure impl.<br>

<a name="how-to-use"></a>
## How to Use

### General Usage

### Adding Environment Maps

### Adding a Data Structure
As the plugin itself uses an intrinsic plugin system for data structures, adding a data structure to the already existing ones is fairly easy.

## License
DS::compare is available under the [GNU GLPv3 license](https://www.gnu.org/licenses/gpl-3.0.html). See [LICENSE](https://github.com/GitThirteen/mitsuba-comparer/blob/main/LICENSE) for the full license text.
