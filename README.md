> [!WARNING]
> This version of Hermes (version 1, cold ion) is no longer maintained.
> Users are recommended to move to [Hermes-3](https://github.com/bendudson/hermes-3/).

      _   _
     | | | |                              
     | |_| | ___ _ __ _ __ ___   ___  ___ 
     |  _  |/ _ \ '__| '_ ` _ \ / _ \/ __|
     | | | |  __/ |  | | | | | |  __/\__ \
     \_| |_/\___|_|  |_| |_| |_|\___||___/

Hermes plasma edge simulation model. Uses BOUT++ framework, adds finite volume
operators and neutral gas models.

This is Hermes-1, a cold ion drift-reduced model.

License
=======

Full text of the license is in the file LICENSE. If you are using Hermes-1,
please cite the relevant papers.

    Copyright B.Dudson, J.Leddy, University of York, September 2016
              email: benjamin.dudson@york.ac.uk

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <http://www.gnu.org/licenses/>.

Source code
===========

* hermes-1.hxx   -- Header file, main physics model class
* hermes-1.cxx   -- Source code for initialisation and RHS function
* div_ops.cxx    -- Collection of differential operators
* loadmetric.cxx -- Loads metric tensor from grid file
* neutral_rates.cxx -- Calculates sources and sinks
* radiation.hxx  -- Interface to atomic cross-sections
* radiation.cxx  -- Atomic cross-sections calculation

Compiling
=========

To compile, run "make" and specify the location of BOUT++
> $ make BOUT_TOP=/path/to/BOUT/

Test cases
==========

 * driftwave - Drift wave in slab geometry. To run, execute script
   "run-driftwave" in this directory
  > $ ./run-driftwave


Simulations
===========

* linear-device
  linear device simulation (no neutrals)

* circular 
  Large aspect ratio, circular cross-section. Axisymmetric electric
  fields and current evolution

* d3d
  Axisymmetric plasma transport (no electric fields)

* d3d-neutrals
  Axisymmetric plasma transport with neutral gas

* ISTTOK
  This test case is for the large aspect ratio circular cross-section device.
  Grid files are generated by a Python script "circle.py"

