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
`include/mitsuba/ds/*` - Plugin headers (Utility, data structures, etc.)<br>
`src/libcore/ds/*` - Data structure impl.<br>

<a name="how-to-use"></a>
## How to Use

### General Usage

To compile Mitsuba code, simply use
```
scons
```
in the command line.<br>

To run the ds::compare plugin, use
```
mtsutil dscompare
```

Changing parameters via CL arguments is currently unsupported (but planned). Please alter the arguments directly in `dscompare.cpp` in the meantime.

### Adding Environment Maps

As environment maps can be quite huge in terms of file size, they have been excluded from this repository. To test the plugin on your own environment maps, please add them to `data/tests/envmaps`. If the folder does not exist, you are free to create one yourself. If you want to use a custom folder, please make sure to change the path variable in `dscompare.cpp` accordingly.
<br>
<br>
Some of the environment maps used for testing this plugin can also be found here:
- http://benedikt-bitterli.me/
- https://hdri-haven.com/
- https://hdrmaps.com/
- https://www.textures.com/library
- http://dativ.at/
- https://pbrt.org/resources

The environment maps must either be `.hdr` or `.exr` files. Other file formats are not supported.

### Adding a Data Structure
#### Implementation
As the plugin itself uses an intrinsic plugin system for data structures, adding a custom data structure to the already existing ones is fairly easy. Data structures in ds::compare are divided into a `.cpp` and `.h` file, located in `src/libcore/ds` and `include/mitsuba/ds/structures` respectively. Each data structure must extend the `DataStructure` base class and mark itself as visible to the compiler via `MTS_EXPORT_CORE`.
```cpp
struct MTS_EXPORT_CORE MyDataStructure : public DataStructure { ...
```
and implement all functions marked as pure virtual. The plugin will be interacting with the data structure solely through these functions. For more information about each virtual function, take a look at the `ds.h` file located in `include/mitsuba/ds`.

#### Registration
Both implementation (`.cpp`) and header file (`.h`) must be registered in Mitsuba and the ds::compare system for it to be visible.
<br>
- To register the header file, please include the header in `include/mitsuba/ds/include.h`.
- To register the implementation file, please (a) add the file name to the SConscript file located in `src/libcore` and (b) instantiate and attach an instance to the internal cluster in `src/utils/dscompare.cpp` via
```cpp
cluster.attach(new <Class>());
```

## License
DS::compare is available under the [GNU GLPv3 license](https://www.gnu.org/licenses/gpl-3.0.html). See [LICENSE](https://github.com/GitThirteen/mitsuba-comparer/blob/main/LICENSE) for the full license text.
