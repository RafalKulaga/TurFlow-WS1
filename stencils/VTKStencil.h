#ifndef _VTK_STENCIL_H_
#define _VTK_STENCIL_H_

#include "../Definitions.h"
#include "../Parameters.h"
#include "../Stencil.h"
#include "../FlowField.h"
#include <string>
#include <fstream>
#include <sstream>
#include <sys/stat.h>
#include <unistd.h>

/** DONE WS1: Stencil for writting VTK files
 *
 * When iterated with, creates a VTK file.
 */
class VTKStencil : public FieldStencil<FlowField> {

  protected:
        int Nx = _parameters.parallel.localSize[0]; // number of cells declared as sizeX in config file
        int Ny = _parameters.parallel.localSize[1];
        int Nz = _parameters.parallel.localSize[2];


        FLOAT* pressure;  // stores whole pressure field in each timestep
        FLOAT* velocity;  // stores midpoint velocities in all directions
        FLOAT tempVelocity_2D[2]; // temporary storage for midpoint velocites, used to save values from getPressureAndVelocity
        FLOAT tempVelocity_3D[3];

        std::string folderName; // folderName is global as it is used both to create output folder as well as file names

  public:
        /** Constructor
         *
         * @param prefix String with the prefix of the name of the VTK files
         * allocates memory for pressure and velocity fields wrt dimension
         * creates new output folder and if necessary removes existing one with its content
         */
        VTKStencil ( const Parameters & parameters );
        /** Destructor
         *
         * deletes dynamically allocated pressure and velocity fields
         */
        virtual ~VTKStencil();

        /** 2D operation for one position
         *
         * @param flowField State of the flow field
         * @param i Position in the x direction
         * @param j Position in the y direction
         */
        void apply ( FlowField & flowField, int i, int j );

        /** 3D operation for one position
         *
         * @param flowField State of the flow field
         * @param i Position in the x direction
         * @param j Position in the y direction
         * @param k Position in the z direction
         */
        void apply ( FlowField & flowField, int i, int j, int k );

        /** Writes the information to the file
         * @param flowField Flow field to be written
         */
        void write ( FlowField & flowField, int timeStep);

};

#endif
