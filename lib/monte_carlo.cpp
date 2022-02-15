/*

   Copyright (c) 2006-2010, The Scripps Research Institute

   Licensed under the Apache License, Version 2.0 (the "License");
   you may not use this file except in compliance with the License.
   You may obtain a copy of the License at

       http://www.apache.org/licenses/LICENSE-2.0

   Unless required by applicable law or agreed to in writing, software
   distributed under the License is distributed on an "AS IS" BASIS,
   WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
   See the License for the specific language governing permissions and
   limitations under the License.

   Author: Dr. Oleg Trott <ot14@columbia.edu>, 
           The Olson Lab, 
           The Scripps Research Institute

*/

#include "monte_carlo.h"
#include "coords.h"
#include "mutate.h"
#include "quasi_newton.h"

#include "commonMacros.h"
#include "wrapcl.h"
#include "random.h"
#include <iostream>

#include <fstream>
#include <boost/progress.hpp>
#include <thread>

//#define DISPLAY_ANALYSIS
//#define DATA_DISTRIBUTION_TEST
//#define BUILD_KERNEL_FROM_SOURCE
output_type monte_carlo::operator()(model& m, const precalculate& p, const igrid& ig, const precalculate& p_widened, const igrid& ig_widened, const vec& corner1, const vec& corner2, incrementable* increment_me, rng& generator) const {
	output_container tmp;
	this->operator()(m, tmp, p, ig, p_widened, ig_widened, corner1, corner2, increment_me, generator); // call the version that produces the whole container
	VINA_CHECK(!tmp.empty());
	return tmp.front();
}

bool metropolis_accept(fl old_f, fl new_f, fl temperature, rng& generator) {
	if(new_f < old_f) return true;
	const fl acceptance_probability = std::exp((old_f - new_f) / temperature);
	return random_fl(0, 1, generator) < acceptance_probability;
}

void monte_carlo::single_run(model& m, output_type& out, const precalculate& p, const igrid& ig, rng& generator) const {
	conf_size s = m.get_size();
	change g(s);
	vec authentic_v(1000, 1000, 1000);
	out.e = max_fl;
	output_type current(out);
	quasi_newton quasi_newton_par; quasi_newton_par.max_steps = ssd_par.evals;
	VINA_U_FOR(step, num_steps) {
		output_type candidate(current.c, max_fl);
		mutate_conf(candidate.c, m, mutation_amplitude, generator);
		quasi_newton_par(m, p, ig, candidate, g, hunt_cap);
		if(step == 0 || metropolis_accept(current.e, candidate.e, temperature, generator)) {
			quasi_newton_par(m, p, ig, candidate, g, authentic_v);
			current = candidate;
			if(current.e < out.e)
				out = current;
		}
	}
	quasi_newton_par(m, p, ig, out, g, authentic_v);
}

void monte_carlo::many_runs(model& m, output_container& out, const precalculate& p, const igrid& ig, const vec& corner1, const vec& corner2, sz num_runs, rng& generator) const {
	conf_size s = m.get_size();
	VINA_FOR(run, num_runs) {
		output_type tmp(s, 0);
		tmp.c.randomize(corner1, corner2, generator);
		single_run(m, tmp, p, ig, generator);
		out.push_back(new output_type(tmp));
	}
	out.sort();
}

std::vector<output_type> monte_carlo::cl_to_vina(output_type_cl result_ptr[], int exhaus) const {
	std::vector<output_type> results_vina;
	for (int i = 0; i < exhaus; i++) {
		output_type_cl tmp_cl = result_ptr[i];
		conf tmp_c;
		tmp_c.ligands.resize(1);
		// Position
		for (int j = 0; j < 3; j++)tmp_c.ligands[0].rigid.position[j] = tmp_cl.position[j];
		// Orientation
		qt q(tmp_cl.orientation[0], tmp_cl.orientation[1], tmp_cl.orientation[2], tmp_cl.orientation[3]);
		tmp_c.ligands[0].rigid.orientation = q;
		output_type tmp_vina(tmp_c, tmp_cl.e);
		// torsion
		for (int j = 0; j < tmp_cl.lig_torsion_size; j++)tmp_vina.c.ligands[0].torsions.push_back(tmp_cl.lig_torsion[j]);
		// coords
		for (int j = 0; j < MAX_NUM_OF_ATOMS; j++) {
			vec v_tmp(tmp_cl.coords[j][0], tmp_cl.coords[j][1], tmp_cl.coords[j][2]);
			if (v_tmp[0] * v_tmp[1] * v_tmp[2] != 0) tmp_vina.coords.push_back(v_tmp);
		}
		results_vina.push_back(tmp_vina);
	}
	return results_vina;
}

void monte_carlo::generate_uniform_position(const vec corner1, const vec corner2, std::vector<vec>& uniform_data, int exhaustiveness) const{
	assert(exhaustiveness == 125);
	int counter = 0;
	int n = 5;
	for (int i = 0; i < n; i++) {
		double p0 = (corner1.data[0] - corner2.data[0]) / n * i + corner2.data[0];
		for (int j = 0; j < n; j++) {
			double p1 = (corner1.data[1] - corner2.data[1]) / n * j + corner2.data[1];
			for (int k = 0; k < n; k++) {
				double p2 = (corner1.data[2] - corner2.data[2]) / n * k + corner2.data[2];
				uniform_data[counter] = vec(p0, p1, p2);
				counter++;
			}
		}
	}
	assert(counter == exhaustiveness);
}

#ifndef OPENCL_PART_2
// out is sorted
void monte_carlo::operator()(model& m, output_container& out, const precalculate& p, const igrid& ig, const precalculate& p_widened, const igrid& ig_widened, const vec& corner1, const vec& corner2, incrementable* increment_me, rng& generator) const {
	vec authentic_v(1000, 1000, 1000); // FIXME? this is here to avoid max_fl/max_fl
	conf_size s = m.get_size();
	change g(s);
	output_type tmp(s, 0);
	tmp.c.randomize(corner1, corner2, generator);
	fl best_e = max_fl;
	quasi_newton quasi_newton_par; quasi_newton_par.max_steps = ssd_par.evals;
	output_type origin = tmp;

	VINA_U_FOR(step, num_steps) {
		if (increment_me)
			++(*increment_me);
		output_type candidate = tmp;
		mutate_conf(candidate.c, m, mutation_amplitude, generator);
		quasi_newton_par(m, p, ig, candidate, g, hunt_cap);
		if (step == 0 || metropolis_accept(tmp.e, candidate.e, temperature, generator)) {
			tmp = candidate;

			m.set(tmp.c); // FIXME? useless?

			// FIXME only for very promising ones
			if (tmp.e < best_e || out.size() < num_saved_mins) {
				quasi_newton_par(m, p, ig, tmp, g, authentic_v);
				m.set(tmp.c); // FIXME? useless?
				tmp.coords = m.get_heavy_atom_movable_coords();
				add_to_output_container(out, tmp, min_rmsd, num_saved_mins); // 20 - max size
				if (tmp.e < best_e)
					best_e = tmp.e;
			}
		}
	}
	VINA_CHECK(!out.empty());
	VINA_CHECK(out.front().e <= out.back().e); // make sure the sorting worked in the correct order
}
#else

volatile bool finished = false; //20211119 Glinttsd
void print_process() {
	int count = 0;
	printf("\n");
	do
	{
#ifdef WIN32
		Sleep(100);
#else
		sleep(0.1);
#endif
		printf("\rPerform docking|");
		for (int i = 0; i < count; i++)printf(" ");
		printf("=======");
		for (int i = 0; i < 30 - count; i++)printf(" ");
		printf("|"); fflush(stdout);

		count++;
		count %= 30;
	} while (finished != true);
	printf("\rPerform docking|");
	for (int i = 0; i < 16; i++)printf("=");
	printf("done");
	for (int i = 0; i < 17; i++)printf("=");
	printf("|\n"); fflush(stdout);
}



void monte_carlo::operator()(model& m, output_container& out, const precalculate& p, const igrid& ig, const precalculate& p_widened, const igrid& ig_widened, const vec& corner1, const vec& corner2, incrementable* increment_me, rng& generator) const {
	/**************************************************************************/
	/***************************    OpenCL Init    ****************************/
	/**************************************************************************/


	cl_int err;
	cl_platform_id* platforms;
	cl_device_id* devices;
	cl_context context;
	cl_command_queue queue;
	cl_int gpu_platform_id;
	SetupPlatform(&platforms, &gpu_platform_id);
	SetupDevice(platforms, &devices, gpu_platform_id);
	SetupContext(platforms, devices, &context, 1, gpu_platform_id);
	SetupQueue(&queue, context, devices);
	char* program_file;
	cl_program program_cl;
	cl_program program;
	size_t program_size;


	
	//Read kernel source code
#ifdef BUILD_KERNEL_FROM_SOURCE
	printf("\nBuild kernels from source"); fflush(stdout);
	const std::string default_work_path = ".";
	const std::string include_path = default_work_path + "/OpenCL/inc"; //FIX it
	const std::string addtion = "";

	char* program_file_n[NUM_OF_FILES];
	size_t program_size_n[NUM_OF_FILES];
	std::string file_paths[NUM_OF_FILES] = {	default_work_path + "/OpenCL/src/kernels/code_head.cpp",
												default_work_path + "/OpenCL/src/kernels/mutate_conf.cpp",
												default_work_path + "/OpenCL/src/kernels/matrix.cpp",
												default_work_path + "/OpenCL/src/kernels/quasi_newton.cpp",
												default_work_path + "/OpenCL/src/kernels/kernel2.cl"}; // The order of files is important!

	read_n_file(program_file_n, program_size_n, file_paths, NUM_OF_FILES);
	std::string final_file;
	size_t final_size = NUM_OF_FILES - 1; // count '\n'
	for (int i = 0; i < NUM_OF_FILES; i++) {
		if (i == 0) final_file = program_file_n[0];
		else final_file = final_file + '\n' + (std::string)program_file_n[i];
		final_size += program_size_n[i];
	}
	const char* final_files_char = final_file.data();	

	program_cl = clCreateProgramWithSource(context, 1, (const char**)&final_files_char, &final_size, &err); checkErr(err);
	SetupBuildProgramWithSource(program_cl, NULL, devices, include_path, addtion);
	SaveProgramToBinary(program_cl, "Kernel2_Opt.bin");
#endif

	
	//Display the progress
	printf("\nSearch depth is set to %d",search_depth);
	std::thread console_thread(print_process);
	

	program_cl = SetupBuildProgramWithBinary(context, devices, "Kernel2_Opt.bin");

	err = clUnloadPlatformCompiler(platforms[gpu_platform_id]); checkErr(err);
	//Set kernel arguments
	cl_kernel kernels[1];
	char kernel_name[][50] = { "kernel2" };
	SetupKernel(kernels, program_cl, 1, kernel_name);
	size_t max_wg_size; // max work item within one work group
	size_t max_wi_size[3]; // max work item within each dimension(global)
	err = clGetDeviceInfo(devices[0], CL_DEVICE_MAX_WORK_GROUP_SIZE, sizeof(size_t), &max_wg_size, NULL); checkErr(err);
	err = clGetDeviceInfo(devices[0], CL_DEVICE_MAX_WORK_ITEM_SIZES, 3 * sizeof(size_t), &max_wi_size, NULL); checkErr(err);

	/**************************************************************************/
	/************************    Original Vina code    ************************/
	/**************************************************************************/

	vec authentic_v(1000, 1000, 1000); // FIXME? this is here to avoid max_fl/max_fl
	conf_size s = m.get_size();
	change g(s);
	output_type tmp(s, 0);
	quasi_newton quasi_newton_par; const int quasi_newton_par_max_steps = ssd_par.evals;

	/**************************************************************************/
	/************************    Allocate CPU memory    ***********************/
	/**************************************************************************/

	// Preparing m related data
	m_cl m_cl;
	assert(m.atoms.size() < MAX_NUM_OF_ATOMS);

	for (int i = 0; i < m.atoms.size(); i++) {
		m_cl.atoms[i].types[0] = m.atoms[i].el;// To store 4 atoms types (el, ad, xs, sy)
		m_cl.atoms[i].types[1] = m.atoms[i].ad;
		m_cl.atoms[i].types[2] = m.atoms[i].xs;
		m_cl.atoms[i].types[3] = m.atoms[i].sy;
		for (int j = 0; j < 3; j++) {
			m_cl.atoms[i].coords[j] = m.atoms[i].coords[j];// To store atom coords
		}
	}

	// To store atoms coords
	for (int i = 0; i < m.coords.size(); i++) {
		for (int j = 0; j < 3; j++) {
			m_cl.m_coords.coords[i][j] = m.coords[i].data[j];
		}
	}

	//To store minus forces
	for (int i = 0; i < m.coords.size(); i++) {
		for (int j = 0; j < 3; j++) {
			m_cl.minus_forces.coords[i][j] = m.minus_forces[i].data[j];
		}
	}

	// Preparing ligand data
	assert(m.num_other_pairs() == 0); // m.other_paris is not supported!
	assert(m.ligands.size() == 1); // Only one ligand supported!
	m_cl.ligand.pairs.num_pairs = m.ligands[0].pairs.size();
	for (int i = 0; i < m_cl.ligand.pairs.num_pairs; i++) {
		m_cl.ligand.pairs.type_pair_index[i]	= m.ligands[0].pairs[i].type_pair_index;
		m_cl.ligand.pairs.a[i]					= m.ligands[0].pairs[i].a;
		m_cl.ligand.pairs.b[i]					= m.ligands[0].pairs[i].b;
	}
	m_cl.ligand.begin = m.ligands[0].begin; // 0
	m_cl.ligand.end = m.ligands[0].end; // 29
	ligand m_ligand = m.ligands[0]; // Only support one ligand 
	assert(m_ligand.end < MAX_NUM_OF_ATOMS);

	// Store root node
	m_cl.ligand.rigid.atom_range[0][0] = m_ligand.node.begin;
	m_cl.ligand.rigid.atom_range[0][1] = m_ligand.node.end;
	for (int i = 0; i < 3; i++) m_cl.ligand.rigid.origin[0][i] = m_ligand.node.get_origin()[i];
	for (int i = 0; i < 9; i++) m_cl.ligand.rigid.orientation_m[0][i] = m_ligand.node.get_orientation_m().data[i];
	m_cl.ligand.rigid.orientation_q[0][0] = m_ligand.node.orientation().R_component_1();
	m_cl.ligand.rigid.orientation_q[0][1] = m_ligand.node.orientation().R_component_2();
	m_cl.ligand.rigid.orientation_q[0][2] = m_ligand.node.orientation().R_component_3();
	m_cl.ligand.rigid.orientation_q[0][3] = m_ligand.node.orientation().R_component_4();
	for (int i = 0; i < 3; i++) {m_cl.ligand.rigid.axis[0][i] = 0;m_cl.ligand.rigid.relative_axis[0][i] = 0;m_cl.ligand.rigid.relative_origin[0][i] = 0;}

	// Store children nodes (in depth-first order)
	struct tmp_struct { 
		int start_index = 0;
		int parent_index = 0;
		void store_node(tree<segment>& child_ptr, rigid_cl& rigid) {
			start_index++; // start with index 1, index 0 is root node
			rigid.parent[start_index] = parent_index;
			rigid.atom_range[start_index][0] = child_ptr.node.begin;
			rigid.atom_range[start_index][1] = child_ptr.node.end;
			for (int i = 0; i < 9; i++) rigid.orientation_m[start_index][i] = child_ptr.node.get_orientation_m().data[i];
			rigid.orientation_q[start_index][0] = child_ptr.node.orientation().R_component_1();
			rigid.orientation_q[start_index][1] = child_ptr.node.orientation().R_component_2();
			rigid.orientation_q[start_index][2] = child_ptr.node.orientation().R_component_3();
			rigid.orientation_q[start_index][3] = child_ptr.node.orientation().R_component_4();
			for (int i = 0; i < 3; i++) {
				rigid.origin[start_index][i] = child_ptr.node.get_origin()[i];
				rigid.axis[start_index][i] = child_ptr.node.get_axis()[i];
				rigid.relative_axis[start_index][i] = child_ptr.node.relative_axis[i];
				rigid.relative_origin[start_index][i] = child_ptr.node.relative_origin[i];
			}
			if (child_ptr.children.size() == 0) return;
			else {
				assert(start_index < MAX_NUM_OF_RIGID);
				int parent_index_tmp = start_index;
				for (int i = 0; i < child_ptr.children.size(); i++) {
					this->parent_index = parent_index_tmp; // Update parent index
					this->store_node(child_ptr.children[i], rigid);
				}
			}
		}
	};
	tmp_struct ts;
	for (int i = 0; i < m_ligand.children.size(); i++) {
		ts.parent_index = 0; // Start a new branch, whose parent is 0
		ts.store_node(m_ligand.children[i], m_cl.ligand.rigid);
	}
	m_cl.ligand.rigid.num_children = ts.start_index;

	// set children_map
	for (int i = 0; i < MAX_NUM_OF_RIGID; i++)
		for (int j = 0; j < MAX_NUM_OF_RIGID; j++)
			m_cl.ligand.rigid.children_map[i][j] = false;
	for (int i = 1; i < m_cl.ligand.rigid.num_children + 1; i++) {
		int parent_index = m_cl.ligand.rigid.parent[i];
		m_cl.ligand.rigid.children_map[parent_index][i] = true;
	}
	m_cl.m_num_movable_atoms = m.num_movable_atoms();
	size_t m_cl_size = sizeof(m_cl);

	// Preparing ig related data
	ig_cl* ig_cl_ptr = (ig_cl*)malloc(sizeof(ig_cl));
	ig_cl_ptr->atu = ig.get_atu(); // atu
	ig_cl_ptr->slope = ig.get_slope(); // slope
	std::vector<grid> tmp_grids = ig.get_grids();
	int grid_size = tmp_grids.size();
	assert(GRIDS_SIZE == grid_size); // grid_size has to be 17

	for (int i = 0; i < grid_size; i++) {
		for (int j = 0; j < 3; j++) {
			ig_cl_ptr->grids[i].m_init[j] = tmp_grids[i].m_init[j];
			ig_cl_ptr->grids[i].m_factor[j] = tmp_grids[i].m_factor[j];
			ig_cl_ptr->grids[i].m_dim_fl_minus_1[j] = tmp_grids[i].m_dim_fl_minus_1[j];
			ig_cl_ptr->grids[i].m_factor_inv[j] = tmp_grids[i].m_factor_inv[j];
		}
		if (tmp_grids[i].m_data.dim0() != 0) {
			ig_cl_ptr->grids[i].m_i = tmp_grids[i].m_data.dim0(); assert(MAX_NUM_OF_GRID_MI >= ig_cl_ptr->grids[i].m_i);
			ig_cl_ptr->grids[i].m_j = tmp_grids[i].m_data.dim1(); assert(MAX_NUM_OF_GRID_MJ >= ig_cl_ptr->grids[i].m_j);
			ig_cl_ptr->grids[i].m_k = tmp_grids[i].m_data.dim2(); assert(MAX_NUM_OF_GRID_MK >= ig_cl_ptr->grids[i].m_k);

			//for (int j = 0; j < ig_cl_ptr->grids[i].m_i * ig_cl_ptr->grids[i].m_j * ig_cl_ptr->grids[i].m_k; j++) {
			//	ig_cl_ptr->grids[i].m_data[j] = tmp_grids[i].m_data.m_data[j];
			//}

			int m_i = tmp_grids[i].m_data.dim0();
			int m_j = tmp_grids[i].m_data.dim1();
			int m_k = tmp_grids[i].m_data.dim2();
			for (sz idx0 = 0; idx0 < m_i - 1; idx0++)
				for (sz idx1 = 0; idx1 < m_j - 1; idx1++)
					for (sz idx2 = 0; idx2 < m_k - 1; idx2++) {
						int base = (idx0 + m_i * (idx1 + m_j * idx2)) * 8;
						ig_cl_ptr->grids[i].m_data[base] = tmp_grids[i].m_data(idx0, idx1, idx2); // f000
						ig_cl_ptr->grids[i].m_data[base + 1] = tmp_grids[i].m_data(idx0 + 1, idx1, idx2); // f100
						ig_cl_ptr->grids[i].m_data[base + 2] = tmp_grids[i].m_data(idx0, idx1 + 1, idx2); // f010
						ig_cl_ptr->grids[i].m_data[base + 3] = tmp_grids[i].m_data(idx0 + 1, idx1 + 1, idx2); // f110
						ig_cl_ptr->grids[i].m_data[base + 4] = tmp_grids[i].m_data(idx0, idx1, idx2 + 1); // f001
						ig_cl_ptr->grids[i].m_data[base + 5] = tmp_grids[i].m_data(idx0 + 1, idx1, idx2 + 1); // f101
						ig_cl_ptr->grids[i].m_data[base + 6] = tmp_grids[i].m_data(idx0, idx1 + 1, idx2 + 1); // f011
						ig_cl_ptr->grids[i].m_data[base + 7] = tmp_grids[i].m_data(idx0 + 1, idx1 + 1, idx2 + 1); // f111
					}

		
		}
		else {
			ig_cl_ptr->grids[i].m_i = 0;
			ig_cl_ptr->grids[i].m_j = 0;
			ig_cl_ptr->grids[i].m_k = 0;
		}
	}
	size_t ig_cl_size = sizeof(*ig_cl_ptr);
	
	// Generating random ligand structures
	std::vector<output_type_cl*> rand_molec_struc_vec; rand_molec_struc_vec.resize(thread);
	for (int i = 0; i < thread; i++) {
		rand_molec_struc_vec[i] = (output_type_cl*)malloc(sizeof(output_type_cl));
	}

	//std::vector<output_type_cl> rand_molec_struc_vec; 
	//rand_molec_struc_vec.resize(thread);
	int lig_torsion_size = tmp.c.ligands[0].torsions.size();
	int flex_torsion_size; if (tmp.c.flex.size() != 0) flex_torsion_size = tmp.c.flex[0].torsions.size(); else flex_torsion_size = 0;
	std::vector<vec> uniform_data;
	uniform_data.resize(thread);
	
	for (int i = 0; i < thread; i++) {
		tmp.c.randomize(corner1, corner2, generator); // generate a random structure
		//Generate positions with uniform probability
		//generate_uniform_position(corner1, corner2, uniform_data, exhaustiveness);
		//for (int j = 0; j < 3; j++) rand_molec_struc_vec[i].position[j] = uniform_data[i].data[j];
		for (int j = 0; j < 3; j++) rand_molec_struc_vec[i]->position[j] = tmp.c.ligands[0].rigid.position[j];
		assert(lig_torsion_size < MAX_NUM_OF_LIG_TORSION);
		for (int j = 0; j < lig_torsion_size; j++) rand_molec_struc_vec[i]->lig_torsion[j] = tmp.c.ligands[0].torsions[j];// Only support one ligand
		assert(flex_torsion_size < MAX_NUM_OF_FLEX_TORSION);
		for (int j = 0; j < flex_torsion_size; j++) rand_molec_struc_vec[i]->flex_torsion[j] = tmp.c.flex[0].torsions[j];// Only support one flex

		rand_molec_struc_vec[i]->orientation[0] = (float)tmp.c.ligands[0].rigid.orientation.R_component_1();
		rand_molec_struc_vec[i]->orientation[1] = (float)tmp.c.ligands[0].rigid.orientation.R_component_2();
		rand_molec_struc_vec[i]->orientation[2] = (float)tmp.c.ligands[0].rigid.orientation.R_component_3();
		rand_molec_struc_vec[i]->orientation[3] = (float)tmp.c.ligands[0].rigid.orientation.R_component_4();

		rand_molec_struc_vec[i]->lig_torsion_size = lig_torsion_size;
	}

	// Preaparing p related data
	p_cl p_cl;
	p_cl.m_cutoff_sqr = p.cutoff_sqr();
	p_cl.factor = p.factor;
	p_cl.n = p.n;
	assert(MAX_P_DATA_M_DATA_SIZE > p.data.m_data.size());
	for (int i = 0; i < p.data.m_data.size(); i++) {
		p_cl.m_data[i].factor = p.data.m_data[i].factor;
		assert(FAST_SIZE == p.data.m_data[i].fast.size());
		assert(SMOOTH_SIZE == p.data.m_data[i].smooth.size());
		for (int j = 0; j < FAST_SIZE; j++) {
			p_cl.m_data[i].fast[j] = p.data.m_data[i].fast[j];
		}
		for (int j = 0; j < SMOOTH_SIZE; j++) {
			p_cl.m_data[i].smooth[j][0] = p.data.m_data[i].smooth[j].first;
			p_cl.m_data[i].smooth[j][1] = p.data.m_data[i].smooth[j].second;
		}
	}
	size_t p_cl_size = sizeof(p_cl);

	// Generate random maps
	random_maps* rand_maps = (random_maps*)malloc(sizeof(random_maps));
	for (int i = 0; i < MAX_NUM_OF_RANDOM_MAP; i++) {
		rand_maps->int_map[i] = random_int(0, int(lig_torsion_size), generator);
		rand_maps->pi_map[i] = random_fl(-pi, pi, generator);
	}
	for (int i = 0; i < MAX_NUM_OF_RANDOM_MAP; i++) {
		vec rand_coords = random_inside_sphere(generator);
		for (int j = 0; j < 3 ; j ++) {
			rand_maps->sphere_map[i][j] = rand_coords[j];
		}
	}
	size_t rand_maps_size = sizeof(*rand_maps);


	float hunt_cap_float[3] = {hunt_cap[0], hunt_cap[1], hunt_cap[2]};
	float authentic_v_float[3] = { authentic_v[0],authentic_v[1], authentic_v[2] };
	float mutation_amplitude_float = mutation_amplitude;
	float epsilon_fl_float = epsilon_fl;
	int	total_wi = max_wi_size[0] * max_wi_size[1];
	/**************************************************************************/
	/************************    Allocate GPU memory    ***********************/
	/**************************************************************************/

	cl_mem rand_molec_struc_vec_gpu;
	CreateDeviceBuffer(&rand_molec_struc_vec_gpu, CL_MEM_READ_ONLY, thread* SIZE_OF_MOLEC_STRUC, context);
	for (int i = 0; i < thread; i++) {
		std::vector<float> pos(rand_molec_struc_vec[i]->position, rand_molec_struc_vec[i]->position + 3);
		std::vector<float> ori(rand_molec_struc_vec[i]->orientation, rand_molec_struc_vec[i]->orientation + 4);
		std::vector<float> lig_tor(rand_molec_struc_vec[i]->lig_torsion, rand_molec_struc_vec[i]->lig_torsion + MAX_NUM_OF_LIG_TORSION);
		std::vector<float> flex_tor(rand_molec_struc_vec[i]->flex_torsion, rand_molec_struc_vec[i]->flex_torsion + MAX_NUM_OF_FLEX_TORSION);
		float lig_tor_size = rand_molec_struc_vec[i]->lig_torsion_size;
		err = clEnqueueWriteBuffer(	queue, rand_molec_struc_vec_gpu, false, i * SIZE_OF_MOLEC_STRUC,
									pos.size() * sizeof(float), pos.data(), 0, NULL, NULL); checkErr(err);
		err = clEnqueueWriteBuffer(	queue, rand_molec_struc_vec_gpu, false, i * SIZE_OF_MOLEC_STRUC + pos.size() * sizeof(float),
									ori.size() * sizeof(float), ori.data(), 0, NULL, NULL); checkErr(err);
		err = clEnqueueWriteBuffer(	queue, rand_molec_struc_vec_gpu, false, i * SIZE_OF_MOLEC_STRUC + (pos.size() + ori.size() ) * sizeof(float),
									lig_tor.size() * sizeof(float), lig_tor.data(), 0, NULL, NULL); checkErr(err);
		err = clEnqueueWriteBuffer(queue, rand_molec_struc_vec_gpu, false, i * SIZE_OF_MOLEC_STRUC + (pos.size() + ori.size() + MAX_NUM_OF_LIG_TORSION) * sizeof(float),
									flex_tor.size() * sizeof(float), flex_tor.data(), 0, NULL, NULL); checkErr(err);
		err = clEnqueueWriteBuffer(queue, rand_molec_struc_vec_gpu, false, i * SIZE_OF_MOLEC_STRUC + (pos.size() + ori.size() + MAX_NUM_OF_LIG_TORSION + MAX_NUM_OF_FLEX_TORSION) * sizeof(float),
			sizeof(float), &lig_tor_size, 0, NULL, NULL); checkErr(err);
	}

	cl_mem best_e_gpu;
	CreateDeviceBuffer(&best_e_gpu, CL_MEM_READ_WRITE, thread * sizeof(float), context);
	err = clEnqueueFillBuffer(queue, best_e_gpu, &max_fl, sizeof(float), 0, thread * sizeof(float), 0, NULL, NULL); checkErr(err);

	cl_mem rand_maps_gpu;
	CreateDeviceBuffer(&rand_maps_gpu, CL_MEM_READ_ONLY, rand_maps_size, context);
	err = clEnqueueWriteBuffer(	queue, rand_maps_gpu, false, 0,	rand_maps_size, rand_maps, 0, NULL, NULL); checkErr(err);

	cl_mem hunt_cap_gpu;
	CreateDeviceBuffer(&hunt_cap_gpu, CL_MEM_READ_ONLY, 3 * sizeof(float), context);
	err = clEnqueueWriteBuffer(queue, hunt_cap_gpu, false, 0, 3 * sizeof(float), hunt_cap_float, 0, NULL, NULL); checkErr(err);

	// Preparing m related data
	cl_mem m_cl_gpu;
	CreateDeviceBuffer(&m_cl_gpu, CL_MEM_READ_WRITE, m_cl_size, context);
	err = clEnqueueWriteBuffer(queue, m_cl_gpu, false, 0, m_cl_size, &m_cl, 0, NULL, NULL); checkErr(err);

	// Preparing p related data
	cl_mem p_cl_gpu;
	CreateDeviceBuffer(&p_cl_gpu, CL_MEM_READ_ONLY, p_cl_size, context);
	err = clEnqueueWriteBuffer(queue, p_cl_gpu, false, 0, p_cl_size, &p_cl, 0, NULL, NULL); checkErr(err);

	// Preparing ig related data (cache related data)
	cl_mem ig_cl_gpu;
	CreateDeviceBuffer(&ig_cl_gpu, CL_MEM_READ_ONLY, ig_cl_size, context);
	err = clEnqueueWriteBuffer(queue, ig_cl_gpu, false, 0, ig_cl_size, ig_cl_ptr, 0, NULL, NULL); checkErr(err);

	cl_mem authentic_v_gpu;
	CreateDeviceBuffer(&authentic_v_gpu, CL_MEM_READ_ONLY, 3 * sizeof(float), context);
	err = clEnqueueWriteBuffer(queue, authentic_v_gpu, false, 0, 3 * sizeof(float), authentic_v_float, 0, NULL, NULL); checkErr(err);

	// Preparing result data
	cl_mem results;
	CreateDeviceBuffer(&results, CL_MEM_READ_WRITE | CL_MEM_ALLOC_HOST_PTR, thread * sizeof(output_type_cl), context);
	
	clFinish(queue);

	/**************************************************************************/
	/************************   Set kernel arguments    ***********************/
	/**************************************************************************/
	SetKernelArg(kernels[0], 0, sizeof(cl_mem),		&m_cl_gpu);
	SetKernelArg(kernels[0], 1, sizeof(cl_mem),		&ig_cl_gpu);
	SetKernelArg(kernels[0], 2, sizeof(cl_mem),		&p_cl_gpu);
	SetKernelArg(kernels[0], 3,	sizeof(cl_mem),		&rand_molec_struc_vec_gpu);
	SetKernelArg(kernels[0], 4, sizeof(cl_mem),		&best_e_gpu);
	SetKernelArg(kernels[0], 5, sizeof(int),		&quasi_newton_par_max_steps);
	SetKernelArg(kernels[0], 6, sizeof(unsigned int),	&num_steps);
	SetKernelArg(kernels[0], 7, sizeof(float),		&mutation_amplitude_float);
	SetKernelArg(kernels[0], 8, sizeof(cl_mem),		&rand_maps_gpu); 
	SetKernelArg(kernels[0], 9, sizeof(float),		&epsilon_fl_float);
	SetKernelArg(kernels[0], 10, sizeof(cl_mem),	&hunt_cap_gpu);
	SetKernelArg(kernels[0], 11, sizeof(cl_mem),	&authentic_v_gpu);
	SetKernelArg(kernels[0], 12, sizeof(cl_mem),	&results);
	SetKernelArg(kernels[0], 13, sizeof(int),		&search_depth);
	SetKernelArg(kernels[0], 14, sizeof(int),		&thread); 
	SetKernelArg(kernels[0], 15, sizeof(int),	&total_wi);
	/**************************************************************************/
	/****************************   Start kernel    ***************************/
	/**************************************************************************/
	size_t global_size[2] = {512, 32 };
	size_t local_size[2] = { 16,2 };
	cl_event monte_clarlo_cl;
	err = clEnqueueNDRangeKernel(queue, kernels[0], 2, 0, global_size, local_size, 0, NULL, &monte_clarlo_cl); checkErr(err);

	clWaitForEvents(1, &monte_clarlo_cl);

	finished = true;
	console_thread.join(); // wait the thread finish

	//getchar();

	// Maping result data
 	output_type_cl* result_ptr = (output_type_cl*)clEnqueueMapBuffer(queue, results, CL_TRUE, CL_MAP_READ, 
																	0, thread * sizeof(output_type_cl),
																	0, NULL, NULL, &err); checkErr(err);

	std::vector<output_type> result_vina = cl_to_vina(result_ptr, thread);
	// if empty, something goes wrong in the device part
	if (result_vina.size() == 0) {
		printf("Error in the device part\n"); exit(-1);
	}


	// Unmaping result data
	err = clEnqueueUnmapMemObject(queue, results, result_ptr, 0, NULL, NULL); checkErr(err);

	// Write back to vina
	for (int i = 0; i < thread; i++) {
		add_to_output_container(out, result_vina[i], min_rmsd, num_saved_mins);
	}
	VINA_CHECK(!out.empty());
	VINA_CHECK(out.front().e <= out.back().e);

	/**************************************************************************/
	/*******************************   Finish    ******************************/
	/**************************************************************************/
	// Memory objects release
	err = clReleaseMemObject(m_cl_gpu); checkErr(err);
	err = clReleaseMemObject(ig_cl_gpu); checkErr(err);
	err = clReleaseMemObject(p_cl_gpu); checkErr(err);
	err = clReleaseMemObject(rand_molec_struc_vec_gpu); checkErr(err);
	err = clReleaseMemObject(best_e_gpu); checkErr(err);
	err = clReleaseMemObject(hunt_cap_gpu); checkErr(err);
	err = clReleaseMemObject(authentic_v_gpu); checkErr(err);
	err = clReleaseMemObject(results); checkErr(err);

	free(ig_cl_ptr);
	free(rand_maps);
	for (int i = 0; i < thread; i++)free(rand_molec_struc_vec[i]);



	// Output Analysis
	cl_ulong time_start, time_end;
	double total_time;
	err = clGetEventProfilingInfo(monte_clarlo_cl, CL_PROFILING_COMMAND_START, sizeof(time_start), &time_start, NULL); checkErr(err);
	err = clGetEventProfilingInfo(monte_clarlo_cl, CL_PROFILING_COMMAND_END, sizeof(time_end), &time_end, NULL); checkErr(err);
	total_time = time_end - time_start;
	printf("GPU monte carlo runtime: %0.3f s", (total_time / 1000000000.0));

#ifdef DISPLAY_ANALYSIS
	// Output Analysis
	cl_ulong time_start, time_end;
	double total_time;
	err = clGetEventProfilingInfo(monte_clarlo_cl, CL_PROFILING_COMMAND_START, sizeof(time_start), &time_start, NULL); checkErr(err);
	err = clGetEventProfilingInfo(monte_clarlo_cl, CL_PROFILING_COMMAND_END, sizeof(time_end), &time_end, NULL); checkErr(err);
	total_time = time_end - time_start;
	printf("\nGPU monte carlo runtime = %0.3f s\n", (total_time / 1000000000.0));

	std::ofstream file("gpu_runtime.txt");
	if (file.is_open())
	{
		file << "GPU monte carlo runtime = " << (double)(total_time / 1000000000.0) << " s" << std::endl;
		file.close();
	}

	cl_ulong private_mem_size;
	err = clGetKernelWorkGroupInfo(kernels[0], devices[0], CL_KERNEL_PRIVATE_MEM_SIZE, sizeof(cl_ulong), &private_mem_size, NULL); checkErr(err);
	printf("\nprivate mem used = %f KBytes", (double)private_mem_size / 1024);

	cl_ulong local_mem_size;
	err = clGetKernelWorkGroupInfo(kernels[0], devices[0], CL_KERNEL_LOCAL_MEM_SIZE, sizeof(cl_ulong), &local_mem_size, NULL); checkErr(err);
	printf("\nlocal mem used = %d Bytes", (int)local_mem_size);

	printf("\nconstant mem used = %0.3f MBytes \n", (double)(	p_cl_size + ig_cl_size + thread * SIZE_OF_MOLEC_STRUC) / (1024 * 1024));
#endif
	//getchar();
}
#endif // OPENCL_VERSION
