#pragma once
#include <vector>
#include <array>
#include <cmath>
#include <random>
#include <fstream>
#include <Eigen/Dense>
#include <complex>
#include "Parameters.h"

using Link_array = std::vector<double>;
using tensor_index =std::array<int,6>;
using link_index = std::array<int,5>;
using SU3 = Eigen::Matrix<std::complex<double>, 3, 3>;
using SU2 = Eigen::Matrix<std::complex<double>, 2, 2>;
using complex = std::complex<double>;
using SU3_array = std::array<double,18>;
using lattice_index = std::array<int,4>;
using spat_index = std::array<int,3>;
// Main Functions

SU3 get_SU3_at_link(const Link_array&arr, link_index link_index_array); 
void set_link_SU3(Link_array& arr, link_index link_index_array, const SU3& SU3_matrix);
void cold_start_array(Link_array& arr);



// Helpers 
// Tensor index buisness
int flat_index(tensor_index tensor_index_array);
tensor_index tensor_index_array(int index);
void set_array_value(Link_array& arr, tensor_index tensor_index_array,double value); 
double get_array_value(const Link_array& arr, tensor_index tensor_index_array);


// SU3 <--> array

SU3_array SU3_matrix_to_array(const SU3& matrix);
SU3 SU3_array_to_matrix(const SU3_array& array);

// Saving Data
void save_lattice_config(const Link_array& array, const std::string& filename);


// Moving along the lattice
void moveup(lattice_index& lattice_index_array, int link_direction);
void movedown(lattice_index& lattice_index_array, int link_direction);

link_index combine_lattice_index_with_direction(const lattice_index& latice_index_array, int direction);
// Computer Staple Sum 
SU3 compute_staple_sum_at_link(const Link_array& arr, const link_index& link_index_array);



SU2 R_block(const SU3& SU3_matrix);
SU3 R_block_to_SU3(const SU2& SU2_matrix);

SU2 S_block(const SU3& SU3_matrix);
SU3 S_block_to_SU3(const SU2& SU2_matrix);

SU2 T_block(const SU3& SU3_matrix);
SU3 T_block_to_SU3(const SU2& SU2_matrix);

SU2 R_block_to_2by2_unitary(const SU3& SU3_matrix);
SU2 S_block_to_2by2_unitary(const SU3& SU3_matrix);
SU2 T_block_to_2by2_unitary(const SU3& SU3_matrix);


SU3 type_R_heatbath(const Link_array& arr, link_index link_index_array, double beta, SU3 U, SU3 A);
SU3 type_S_heatbath(const Link_array& arr, link_index link_index_array, double beta, SU3 U, SU3 A);
SU3 type_T_heatbath(const Link_array& arr, link_index link_index_array, double beta, SU3 U, SU3 A);
double single_link_heatbath(Link_array& arr, const link_index& link_index_array, double beta);
double heatbath_update(Link_array& arr, double beta);


void check_unitarity(const Link_array& arr);

complex poly_loop_at_spat_coord(const Link_array& arr, int x, int y, int z);
complex compute_correlator(const Link_array& arr, spat_index m, spat_index n);
complex correlator_over_fixed_distance(const Link_array& arr, int r);