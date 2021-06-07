#include <globals.hxx>
#include <output.hxx>
#include <utils.hxx>
#include "hermes-1.hxx"

#include "loadmetric.hxx"

void LoadMetric(BoutReal Lnorm, BoutReal Bnorm) {
  // Load metric coefficients from the mesh
  Field2D Rxy, Bpxy, Btxy, hthe, sinty;
  auto mesh = bout::globals::mesh;

  GRID_LOAD5(Rxy, Bpxy, Btxy, hthe, sinty); // Load metrics

  Coordinates *coord = mesh->getCoordinates();

  std::string paralleltransform;
  OPTION(Options::getRoot()->getSection("mesh"), paralleltransform, "identity");
   
  if (paralleltransform == "shifted") {
    // Shifted metric
    output << "LoadMetric: Using shifted metric\n";
    sinty = 0.0;
  } else {
    output << "LoadMetric: Using field-aligned metric\n";
  }
  
  // Checking for dpsi and qinty used in BOUT grids
  Field2D dx;
  if(!mesh->get(dx,   "dpsi")) {
    output << "\tUsing dpsi as the x grid spacing\n";
    coord->dx = dx; // Only use dpsi if found
  }else {
    // dx will have been read already from the grid
    output << "\tUsing dx as the x grid spacing\n";
  }
  Field2D qinty;
  if(!mesh->get(qinty, "qinty")) {
    output << "\tUsing qinty as the Z shift\n";
    // mesh.get(zShift = qinty;
  }else {
    // Keep zShift
    output << "\tUsing zShift as the Z shift\n";
  }

  Rxy      /= Lnorm;
  hthe     /= Lnorm;
  sinty    *= SQ(Lnorm)*Bnorm;
  coord->dx /= SQ(Lnorm)*Bnorm;
  
  Bpxy /= Bnorm;
  Btxy /= Bnorm;
  coord->Bxy  /= Bnorm;
  
  // Calculate metric components
  
  BoutReal sbp = 1.0; // Sign of Bp
  if(min(Bpxy, true) < 0.0)
    sbp = -1.0;
  
  coord->g11 = pow((Rxy*Bpxy),2);
  coord->g22 = 1.0 / (pow(hthe,2));
  coord->g33 = (pow(sinty,2))*coord->g11 + (pow(coord->Bxy,2))/coord->g11;
  coord->g12 = 0.0;
  coord->g13 = -sinty*coord->g11;
  coord->g23 = -sbp*Btxy/(hthe*Bpxy*Rxy);
  
  coord->J = hthe / Bpxy;
  
  coord->g_11 = 1.0/coord->g11 + (pow((sinty*Rxy),2));
  coord->g_22 = pow((coord->Bxy*hthe/Bpxy),2);
  coord->g_33 = Rxy*Rxy;
  coord->g_12 = sbp*Btxy*hthe*sinty*Rxy/Bpxy;
  coord->g_13 = sinty*Rxy*Rxy;
  coord->g_23 = sbp*Btxy*hthe*Rxy/Bpxy;
  
  coord->geometry();
}
