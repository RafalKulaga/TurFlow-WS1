#include "VTKStencil.h"

VTKStencil::VTKStencil ( const Parameters & parameters ) : FieldStencil<FlowField> ( parameters ) {
    if(_parameters.geometry.dim == 2){
      pressure = new FLOAT[Nx*Ny];
      velocity = new FLOAT[2*Nx*Ny];
    }
    else{
      pressure = new FLOAT[Nx*Ny*Nz];
      velocity = new FLOAT[3*Nx*Ny*Nz];
    }

    folderName =_parameters.vtk.prefix + std::to_string(_parameters.geometry.dim) + "D_Re" + std::to_string((int)_parameters.flow.Re);
    std::string tmp = "rm -rf "+folderName; // system command to remove directory with its content
    if(std::system(tmp.c_str())); // remove old folder, "if" to avoid -unused flag warning
    mkdir(folderName.c_str(), S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH); //make diretory
}

VTKStencil::~VTKStencil() {
  delete [] pressure; pressure = NULL;
  std::cout << "pressure deleted" << std::endl;
  delete [] velocity; velocity = NULL;
  std::cout << "velocity deleted" << std::endl;
}

void VTKStencil::apply ( FlowField & flowField, int i, int j ){
  int offset = 2;
  int ii = i - offset; int jj = j - offset; // shorter notation including relevant offset
  const int obstacle = flowField.getFlags().getValue(i, j);

  if((obstacle & OBSTACLE_SELF) == 0){ // check if fluid
    flowField.getPressureAndVelocity( pressure[ii + jj*Nx], tempVelocity_2D,  i, j);
    velocity[2*(ii + Nx*jj)] = tempVelocity_2D[0];
    velocity[2*(ii + Nx*jj)+1] = tempVelocity_2D[1];
  }
  else{
    pressure[ii + jj*Nx] = 0.0;
    velocity[2*(ii + Nx*jj)] = 0.0;
    velocity[2*(ii + Nx*jj)+1] = 0.0;
  }
}

void VTKStencil::apply ( FlowField & flowField, int i, int j, int k ){
  int offset = 2;
  int ii = i - offset; int jj = j - offset; int kk = k - offset; // shorter notation including relevant offset
  const int obstacle = flowField.getFlags().getValue(i, j, k);

  if((obstacle & OBSTACLE_SELF) == 0){ // check if fluid
    flowField.getPressureAndVelocity( pressure[ii + jj*Nx + kk*Ny*Nx], tempVelocity_3D,  i, j, k);
    velocity[3*(ii + Nx*jj + Nx*Ny*kk)] = tempVelocity_3D[0];
    velocity[3*(ii + Nx*jj + Nx*Ny*kk)+1] = tempVelocity_3D[1];
    velocity[3*(ii + Nx*jj + Nx*Ny*kk)+2] = tempVelocity_3D[2];
  }
  else{
    pressure[ii + jj*Nx + kk*Ny*Nx] = 0.0;
    velocity[3*(ii + Nx*jj + Nx*Ny*kk)] = 0.0;
    velocity[3*(ii + Nx*jj + Nx*Ny*kk)+1] = 0.0;
    velocity[3*(ii + Nx*jj + Nx*Ny*kk)+2] = 0.0;
  }
}

void VTKStencil::write ( FlowField & flowField, int timeStep){
  // Initializtion
  std::ofstream file;
  std::string fileName =  folderName + "/" +folderName+ "." + std::to_string(timeStep)+".vtk";
  file.open (fileName);

  // Write header
  file << "# vtk DataFile Version 2.0\n";
  file << "Empty\n" << "ASCII\n\n" << "DATASET STRUCTURED_GRID\n";
  file << "DIMENSIONS " << Nx+1 << " " << Ny+1 << " " << Nz+1 << "\n";
  file << "POINTS " << (Nx+1)*(Ny+1)*(Nz+1) << " float\n";

  // Set borders for loops considering the offsets. k borders set in "else"
  int j_start = _parameters.parallel.firstCorner[1] + 2;
  int j_end = j_start + _parameters.parallel.localSize[1] + 1;
  int i_start = _parameters.parallel.firstCorner[0] + 2;
  int i_end = i_start + _parameters.parallel.localSize[0] + 1;

  if( _parameters.geometry.dim == 2){
    for(int j = j_start; j < j_end; ++j){
      for(int i = i_start; i < i_end; ++i){
        file << _parameters.meshsize->getPosX(i,j) << " " << _parameters.meshsize->getPosY(i,j) << " 0\n";
      }
    }

    file << "\n" << "CELL_DATA " << Nx*Ny << "\n";
    file << "SCALARS pressure float 1\n";
    file << "LOOKUP_TABLE default\n";

    for(int i = 0; i < Nx*Ny; ++i){ // single loop as all pressure data is organized in one long array
      file << pressure[i] << "\n";
    }

    file << "\nVECTORS velocity float\n";
      for(int i = 0; i < 2*Nx*Ny; i+=2){ // single loop as all velocity data is organized in one long array
        file << velocity[i] << " "  << velocity[i+1] << " 0\n";
      }
  }
  else{

    int k_start = _parameters.parallel.firstCorner[2] + 2;
    int k_end = k_start + _parameters.parallel.localSize[2] + 1;

    for(int k = k_start; k < k_end; ++k){
      for(int j = j_start; j < j_end; ++j){
        for(int i = i_start; i < i_end; ++i){
          file << _parameters.meshsize->getPosX(i,j,k) << " " << _parameters.meshsize->getPosY(i,j,k) << " " <<  _parameters.meshsize->getPosZ(i,j,k) << "\n";
        }
      }
    }

    file << "\n" << "CELL_DATA " << Nx*Ny*Nz << "\n";
    file << "SCALARS pressure float 1\n";
    file << "LOOKUP_TABLE default\n";

    for(int i = 0; i < Nx*Ny*Nz; ++i){
      file << pressure[i] << "\n";
    }

    file << "\nVECTORS velocity float\n";
      for(int i = 0; i < 3*Nx*Ny*Nz; i+=3){
        file << velocity[i] << " "  << velocity[i+1] << " " << velocity[i+2] << "\n";
      }
  }

  file.close(); // close the file
}
