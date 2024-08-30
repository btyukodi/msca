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
