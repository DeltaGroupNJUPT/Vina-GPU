//#include "kernel2.h"
//void print_output(output_type_cl* a) {
//	for (int i = 0; i < 3 + 4 + MAX_NUM_OF_LIG_TORSION;i++) {
//		if (i < 3)printf("position[%d] = %f\n", i, a->position[i]);
//		else if (i < 7)printf("orientation[%d] = %f\n", i - 3, a->orientation[i - 3]);
//		else if (i < 3 + 4 + MAX_NUM_OF_LIG_TORSION)printf("lig_torsion[%d] = %f\n", i - 7, a->lig_torsion[i - 7]);
//	}
//}
//
//void print_change(change_cl* a) {
//	for (int i = 0; i < 3 + 3 + MAX_NUM_OF_LIG_TORSION; i++) {
//		if (i < 3)printf("position[%d] = %f\n", i, a->position[i]);
//		else if (i < 3 + 3)printf("orientation[%d] = %f\n", i - 3, a->orientation[i - 3]);
//		else if (i < 3 + 3 + MAX_NUM_OF_LIG_TORSION)printf("lig_torsion[%d] = %f\n", i - 6, a->lig_torsion[i - 6]);
//	}
//}

void m_cl_init_with_m_cl(const __global m_cl* m_cl_old, m_cl* m_cl_new) {
	for (int i = 0; i < MAX_NUM_OF_ATOMS; i++)m_cl_new->atoms[i] = m_cl_old->atoms[i];
	m_cl_new->m_coords = m_cl_old->m_coords;
	m_cl_new->minus_forces = m_cl_old->minus_forces;
	m_cl_new->ligand = m_cl_old->ligand;
	m_cl_new->m_num_movable_atoms = m_cl_old->m_num_movable_atoms;
}



void get_heavy_atom_movable_coords( output_type_cl* tmp, const m_cl* m_cl_gpu) {
	int counter = 0;
	for (int i = 0; i < m_cl_gpu->m_num_movable_atoms; i++) {
		if (m_cl_gpu->atoms[i].types[0] != EL_TYPE_H_CL) {
			for (int j = 0; j < 3; j++)tmp->coords[counter][j] = m_cl_gpu->m_coords.coords[i][j];
			counter++;
		}
		else {
			//printf("\n kernel2: removed H atom coords in get_heavy_atom_movable_coords()!");
		}
	}
	//对其余值赋值0
	for (int i = counter; i < MAX_NUM_OF_ATOMS; i++) {
		for (int j = 0; j < 3; j++)tmp->coords[i][j] = 0;
	}
}

// 冒泡排序
//void container_sort(out_container* out) {
//	output_type_cl out_tmp;
//	for (int i = 0; i < out->current_size - 1; i++) {
//		for (int j = 0; j < out->current_size - 1 - i; j++) {
//			if (out->container[j].e > out->container[j + 1].e) {
//				output_type_cl_init_with_output(&out_tmp, &out->container[j]);
//				output_type_cl_init_with_output(&out->container[j], &out->container[j+1]);
//				output_type_cl_init_with_output(&out->container[j + 1], &out_tmp);
//			}
//		}
//	}
//}


void add_to_output_container(out_container* out, const output_type_cl* tmp) {
	if (out->current_size <= MAX_CONTAINER_SIZE_EVERY_WI) {
		out->container[out->current_size - 1] = *tmp;
		out->current_size++;
		container_sort(out);
	}
	else {
		out->container[MAX_CONTAINER_SIZE_EVERY_WI - 1] = *tmp;
		container_sort(out);
	}
}

//生成一个0-1之间的随机数(根据step)
float generate_n(__constant float* pi_map, const int step) {
	return fabs(pi_map[step]) / M_PI;
}

bool metropolis_accept(float old_f, float new_f, float temperature, float n) {
	if (new_f < old_f)return true;
	const float acceptance_probability = exp((old_f - new_f) / temperature);
	bool res = n < acceptance_probability;
	//if (res)printf("\n Accept!"); else printf("\n Refuse!");
	//printf("\n n = %f", n);
	return n < acceptance_probability;
}

void write_back(__global output_type_cl* results, const output_type_cl* best_out) {
	for (int i = 0; i < 3; i++)results->position[i] = best_out->position[i];
	for (int i = 0; i < 4; i++)results->orientation[i] = best_out->orientation[i];
	for (int i = 0; i < MAX_NUM_OF_LIG_TORSION; i++)results->lig_torsion[i] = best_out->lig_torsion[i];
	for (int i = 0; i < MAX_NUM_OF_FLEX_TORSION; i++)results->flex_torsion[i] = best_out->flex_torsion[i];
	results->lig_torsion_size = best_out->lig_torsion_size;
	results->e = best_out->e;
	for (int i = 0; i < MAX_NUM_OF_ATOMS; i++) {
		for (int j = 0; j < 3; j++) {
			results->coords[i][j] = best_out->coords[i][j];
		}
	}
}

__kernel
void kernel2(	__global	m_cl*			m_cl_global,
				__constant	ig_cl*			ig_cl_gpu,
				__constant	p_cl*			p_cl_gpu,
				__constant	float*			rand_molec_struc_vec_gpu,
				__global	float*			best_e_gpu,
							int				bfgs_max_steps,
							unsigned		num_steps,
							float			mutation_amplitude,
				__constant	random_maps*	rand_maps_gpu,
							float			epsilon_fl,
				__global	float*			hunt_cap_gpu,
				__global	float*			authentic_v_gpu,
				__global	output_type_cl	results[],
							int				search_depth,
							int				e,
							int				total_wi
)
{
	int gx = get_global_id(0);
	int gy = get_global_id(1);
	int gs = get_global_size(0);
	int gl = get_global_linear_id();

	float best_e = INFINITY;




	//测试
	//int check_id = 640;
	for (int gll = gl;
			 gll < e;
			 gll += total_wi
		)
	{
		//printf("\n size of double = %d", sizeof(double));
		//printf("\n size of float = %d", sizeof(float));
		//printf("\n mutation_amplitude = %f", mutation_amplitude);


		if (gll % 100 == 0)printf("\nThread %d START", gll);

		m_cl m_cl_gpu;
		m_cl_init_with_m_cl(m_cl_global, &m_cl_gpu);


		__private output_type_cl tmp; // private memory, shared only in work item
		__private change_cl g;
		output_type_cl_init(&tmp, rand_molec_struc_vec_gpu + gll * (SIZE_OF_MOLEC_STRUC / sizeof(float)));
		//printf("\n lig_torsion_size = %f", tmp.lig_torsion_size);
		
		g.lig_torsion_size = tmp.lig_torsion_size;
		//print_ouput_type(&tmp, tmp.lig_torsion_size);//测试
		//测试
		//printf("\n g lig_torsion_size = %f", g.lig_torsion_size);

		// BFGS
		output_type_cl best_out;
		output_type_cl candidate;

		for (int step = 0; step < search_depth; step++) {
			//printf("\n gx = %d, step = %d", gx, step);
			//printf("\n tmp lig_torsion = %f", tmp.lig_torsion_size);
			output_type_cl_init_with_output(&candidate, &tmp);
			//printf("\n candidate lig_torsion = %f", candidate.lig_torsion_size);
			//print_ouput_type(&candidate, tmp.lig_torsion_size);//测试

			int map_index = (step + gll * search_depth) % MAX_NUM_OF_RANDOM_MAP;// 复用random map
			mutate_conf_cl(map_index,
				num_steps,
				&candidate,
				rand_maps_gpu->int_map,
				rand_maps_gpu->sphere_map,
				rand_maps_gpu->pi_map,
				m_cl_gpu.ligand.begin,
				m_cl_gpu.ligand.end,
				m_cl_gpu.atoms,
				&m_cl_gpu.m_coords,
				m_cl_gpu.ligand.rigid.origin[0],
				epsilon_fl,
				mutation_amplitude
			);
			//print_ouput_type(&candidate, candidate.lig_torsion_size);//测试

			float e_before = candidate.e;
			
			bfgs(&candidate,
				&g,
				&m_cl_gpu,
				p_cl_gpu,
				ig_cl_gpu,
				hunt_cap_gpu,
				epsilon_fl,
				bfgs_max_steps
			);
			//printf("\n After Here!");
			//if(e_before < candidate.e)printf("\n bfgs wrong!");

			//printf("\n candidate.e = %f", candidate.e);
			
			
			float n = generate_n(rand_maps_gpu->pi_map, map_index);
			if (n > 1 || n < 0)printf("\n generate_n error!");
			
			if (step == 0 || metropolis_accept(tmp.e, candidate.e, 1.2, n)) {

				output_type_cl_init_with_output(&tmp, &candidate);

				set(&tmp, &m_cl_gpu.ligand.rigid, m_cl_gpu.m_coords.coords,
					m_cl_gpu.atoms, m_cl_gpu.m_num_movable_atoms, epsilon_fl);
				
				if (tmp.e < best_e) {
					//精细搜索BFGS
					bfgs(&tmp,
						&g,
						&m_cl_gpu,
						p_cl_gpu,
						ig_cl_gpu,
						authentic_v_gpu,
						epsilon_fl,
						bfgs_max_steps
					);
					//if (tmp.e > best_e)printf("\n 精细bfgs wrong!");
					// set
					if (tmp.e < best_e) {
						set(&tmp, &m_cl_gpu.ligand.rigid, m_cl_gpu.m_coords.coords,
							m_cl_gpu.atoms, m_cl_gpu.m_num_movable_atoms, epsilon_fl);

						output_type_cl_init_with_output(&best_out, &tmp); // 存储除了coords以外的数据
						get_heavy_atom_movable_coords(&best_out, &m_cl_gpu); // 赋值coords
						best_e = tmp.e;
					}

				}
			}
			
		}

		// 将最优解写回CPU
		write_back(&results[gll], &best_out);
		if (gll % 100 == 0)printf("\nThread %d FINISH", gll);
	}
}