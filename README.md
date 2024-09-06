The code needs to be compiled against OpenMesh9.0 available at https://www.graphics.rwth-aachen.de/software/openmesh/

The compilation process should generate the Assembly executable.
 A simple dynamical run takes a flag -d and an input json file: ./Assembly -d -i T7.json
 The .json input file contains all the simulation running parameters including subunit design and run parameters.

 Simulations output a data_log.dat text file with several instantaneous data values and a binary file, snapshots.om containing full snapshots of the structure.
 These binary files may be converted to Lammps and VTK formats for visualization and analysis using ./Assembly -c snapshots.om -f folder_for_snapshot_files/

 Simulations are initialized from structures also stored in .om files. A vertices.dat-faces.dat-edges.dat text file triplet may be used to convert to such an .om file using
 ./Assembly -g datfiles_folder/ -o structure.om   

 The folder init_structures/ contains multiple initializer files which may be specified for use in the .json input file, in the "init_file" parameter.
 It is a good idea to use full paths to all files and folders when providing them as command line arguments to ./Assembly 

The .json file specifies the full input and defines the subunit design, assembly conditions, metaparameters and data dump parameters. The following flags are currently available for running:
 
 ./Assembly -d -i input_carlos2.json --> dynamical run
 
 ./Assembly -u -i input_carlos2.json  --> umbrella sampling
 
 ./Assembly -p -i pt_input.json  --> parallel tempering
 
 ./Assembly -w -i input_carlos2.json  --> wham only on existing data
 
 ./Assembly -W -i PT_input.json  --> wham only on existing PT data
 
 ./Assembly -c snapshots.om -f conversions/ --> convert .om 
 
 ./Assembly -c snapshots.om  -b 803215 -f conversions/ --> convert .om, only the specified byte position 
 
 ./Assembly -g datfiles/ -o structure.om --> generate .om from datfiles/vertices.dat, faces.dat, edges.dat
 
 ./Assembly -a -i input_carlos2.json --> analyze; snapshots.om already written in the .json file
 
 ./Assembly -t -i snapshots.om --> analyze topology only based on what's in snapshots.om alone; only simplified mesh can be accessed
 
 ./Assembly -q -i input_carlos2.json --> quench structure to kT=0
 
 ./Assembly -S -i input_carlos2.json --> dynamical run with wall squishing
 
 ./Assembly -T -i input_carlos2.json --> thermodynamic integration
 
 ./Assembly -i omfile.om -b 803215 -A omfile_to_append_to.om --> append snapshot from position b in omfile.om to the file in omfile_to_append_to.om
