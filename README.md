# caams
## Purpose
This project is a very simple implementation of concepts described in
*Computer-aided Analysis of Mechanical Systems* [^caams]. It implements
all of the constraints described in the book plus a screw constraint
that was developed by me.

[^caams]: *Computer-aided Analysis of Mechanical Systems*
  Parviz E. Nikravesh
  Prentice Hall, 1988

## Building the project on Ubuntu

### Dependencies
 - boost `sudo apt install libboost-dev`
 - Eigen `sudo apt install libeigen3-dev`
 - OpenGL `sudo apt install libopengl-dev`
 - Freeglut `sudo apt install freeglut3-dev`
 - SDL2 `sudo apt install libsdl2-dev`
 - cmake `sudo apt install cmake`
 - gcc, g++ and make `sudo apt install build-essential`

### Steps to build
  1. Navigate to a directory to store the repository.
  1. `git clone https://github.com/TimKrause2/caams.git`
  2. `cd caams`
  3. `mkdir build`
  4. `cd build`
  5. `cmake ..`
  6. `make`

## Running the example programs
Once the project is built using the steps described above the apps will be
available at `caams/build/apps`. There not meant to be installed. Each program
is in its own subdirectory with the same name. The graphics are simple
shapes rendered with wireframes.

