#include <string>
#include <stdlib.h>
// -------------------- OpenMesh
#include <OpenMesh/Core/IO/MeshIO.hh>
#include <OpenMesh/Core/Mesh/TriMesh_ArrayKernelT.hh>
typedef OpenMesh::TriMesh_ArrayKernelT<>  MyMesh;

void analyze_all(std::string input_file);

void topology_analyze_all(std::string om_filename);

int count_boundary_edges(MyMesh & mesh);

std::tuple<double, double> compute_asphericity(MyMesh & mesh);

int count_subunits_of_type(MyMesh & mesh, int face_type);

void process_final_structure(MyMesh & mesh, std::string data_folder);

int convert_om_to_dat_snapshots_undeformed_frame(std::string om_filename, std::string dat_base);

int make_capsid_full_specific(std::string om_filename, std::string outfile_base);

int make_capsid_Sigl_specific(std::string om_filename, std::string outfile_base);

int color_faces_Sigl_specific_bak(std::string om_filename, std::string outfile_base);

int print_capsid_data(std::string om_filename, std::string outfile_base);

int print_capsid_data(MyMesh & mesh, std::string outfile_base);

int push_vertices_to_sphere(std::string om_filename, std::string outfile_base);

void defect_analyze_snapshot(MyMesh & mesh, long t, std::ostream& os);

int make_capsid_daichi_specific(std::string om_filename, std::string daichi_design_file, std::string outfile_base);

//int make_capsid_Sigl_specific_merge(std::string om_filename, std::string outfile_base, std::vector<int> & merge_pattern);

int make_capsid_merged_specific(std::string om_filename, std::string outfile_base, std::vector<std::vector<int> > merge_pattern);

int copy_mesh_test(std::string om_filename);

int make_capsid_merged_specific(std::string om_filename, std::string outfile_base, int n_ck_colors, int k);

int color_facets(std::string om_filename, std::string outfile_base);

int create_simple_seeds(std::string om_filename, std::string outfile_base);

void compute_curvature(std::string om_filename, std::string outfile_path);

void count_disclination_spots(std::string input_file);

int make_capsid_botond_specific(std::string om_filename, std::string outfile_base);

double get_max_height(MyMesh & mesh);

int print_mesh_SAT_input(std::string om_filename, std::string outfile_base);

double compute_mean_quadratic_strain(MyMesh & mesh);

int make_capsid_SAT_specific(std::string om_filename, std::string SAT_coloring_file ,std::string outfile_base);