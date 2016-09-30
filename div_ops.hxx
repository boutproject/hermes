/*
  Finite volume discretisations of divergence operators
 
  ***********

    Copyright B.Dudson, J.Leddy, University of York, September 2016
              email: benjamin.dudson@york.ac.uk

    This file is part of Hermes.

    Hermes is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    Hermes is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with Hermes.  If not, see <http://www.gnu.org/licenses/>.
  
 */

#ifndef __DIV_OPS_H__
#define __DIV_OPS_H__

#include <field3d.hxx>
#include <vector3d.hxx>

// Uses flux-conservative form with ZIP interpolation to cell faces
const Field3D Div_n_a_bxGrad_f_B(const Field3D &n, const Field3D &a, const Field3D &f, bool xflux=false, bool yderiv=true);

// This uses a combination of two flux-conservative discretisations to get an operator
// which is almost, but not quite, anti-symmetric
const Field3D Div_n_bxGrad_f_B(const Field3D &n, const Field3D &f);

// Div ( a * Grad_perp(f) )
const Field3D Div_Perp_Lap(const Field3D &a, const Field3D &f);

const Field3D Div_Perp_Lap_x3(const Field3D &a, const Field3D &f, bool xflux=false);

const Field3D Div_Perp_Lap_XYZ(const Field3D &a, const Field3D &f, bool bndryflux=false);

// Divergence of a parallel diffusion k * Grad_par(f)
const Field3D Div_Par_Diffusion(const Field3D &k, const Field3D &f, bool bndry_flux=true);
const Field3D Div_Par_Diffusion_Index(const Field3D &f, bool bndry_flux=true);

// Added Dissipation method. velocity proportional to 3rd derivative of the pressure
const Field3D AddedDissipation(const Field3D &N, const Field3D &P, const Field3D f, bool bndry_flux=true);

const Field3D Div_n_bxGrad_f_B_XPPM(const Field3D &n, const Field3D &f, bool bndry_flux=true, bool poloidal=false, bool positive=false);

const Field3D Div_f_v_XPPM(const Field3D &n, const Vector3D &v, bool bndry_flux=true);

void communicateFluxes(Field3D &f);

const Field3D Div_Perp_Lap_FV(const Field3D &n, const Field3D &f, bool xflux);
const Field3D Div_Perp_Lap_FV_Index(const Field3D &a, const Field3D &f, bool xflux);

const Field3D Div_par_FV(const Field3D &f, const Field3D &v);

// Finite volume parallel divergence of a flow velocity
const Field3D Div_parV_FV(const Field3D &v);

const Field3D D4DY4_FV(const Field3D &d, const Field3D &f, bool bndry_flux=false);

// 4th-order flux conserving term, in index space
const Field3D D4DX4_FV_Index(const Field3D &f, bool bndry_flux=false);

// Div ( k * Grad(f) )
const Field2D Laplace_FV(const Field2D &k, const Field2D &f);

#endif //  __DIV_OPS_H__
