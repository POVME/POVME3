#include <iostream>
#include <iomanip>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include <algorithm>
#include <stdio.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_eigen.h>

int main (int argc, char * argv[]){

	std::istringstream iss(argv[1]);
	int num_frames=0;
	iss>>num_frames;

	if (argc!=(2+num_frames)){

		std::cerr<<"Usage: ./executable num_frames frame1.xyz_arr frame2.xyz_arr ... frameN.xyz_arr \n";
		return 0;
	}

	std::ifstream ifs;

	std::vector <int> master_coord_x_list;
	std::vector <int> master_coord_y_list;
	std::vector <int> master_coord_z_list;

	int max_x=0;
	int max_y=0;
	int max_z=0;
	int min_x=0;
	int min_y=0;
	int min_z=0;

	std::vector <std::string> master_coord_string_list;

	int total_coords=0;
	ifs.open("pocket_coordinates");
	while(!ifs.eof()){
		std::string my_str;
		getline(ifs, my_str);
		if (ifs.eof()) continue;

		std::istringstream iss(my_str);

		double x_coord, y_coord, z_coord;
		int x_coor, y_coor, z_coor;

		iss>>x_coord;
		iss>>y_coord;
		iss>>z_coord;

		x_coor=int(x_coord);
		y_coor=int(y_coord);
		z_coor=int(z_coord);

		master_coord_x_list.push_back(x_coor);
		master_coord_y_list.push_back(y_coor);
		master_coord_z_list.push_back(z_coor);

		if (x_coor>max_x) max_x=x_coor;
		if (y_coor>max_y) max_y=y_coor;
		if (z_coor>max_z) max_z=z_coor;

		++total_coords;
	}
	ifs.close();

	std::sort(master_coord_string_list.begin(), master_coord_string_list.end());

/*	for (std::vector<string>::iterator it=master_coord_string_list.begin(); it!=master_coord_string_list.end(); ++it){
		std::istringstream iss(*it);
		double x, y, z;
		iss>>x;
		iss>>y;
		iss>>z;

		int x1=int(x);
		int y1=int(y);
		int z1=int(z);

		master_coord_x_list.push_back(x1);
                master_coord_y_list.push_back(y1);
                master_coord_z_list.push_back(z1);

          //      if (x_coor>max_x) max_x=x_coor;
            //    if (y_coor>max_y) max_y=y_coor;
              //  if (z_coor>max_z) max_z=z_coor;
		
	}
*/
	min_x=*std::min_element(master_coord_x_list.begin(), master_coord_x_list.end());
	min_y=*std::min_element(master_coord_y_list.begin(), master_coord_y_list.end());
	min_z=*std::min_element(master_coord_z_list.begin(), master_coord_z_list.end());

	std::cerr<<min_x<<'\t'<<min_y<<'\t'<<min_z<<'\n';
	std::cerr<<max_x<<'\t'<<max_y<<'\t'<<max_z<<'\n';

	std::cerr<<num_frames<<'\n';

	int **** total_frame_coords = new int***[num_frames];
	
	for (int i=0; i<num_frames; ++i){
		total_frame_coords[i]=new int**[max_x+1];
		
		for (int j=0; j<=max_x; ++j){
			total_frame_coords[i][j]=new int*[max_y+1];

			for (int k=0; k<=max_y; ++k){
				total_frame_coords[i][j][k]=new int[max_z+1];
			}
		}

	}
	
	for(int i=0; i<num_frames; ++i){
                ifs.open(argv[i+2]);
                while(!ifs.eof()){
                        std::string my_str;
                        getline(ifs, my_str);
                        if (ifs.eof()) continue;
                        
			std::istringstream iss(my_str);

        	        double x_coord, y_coord, z_coord;
	                int x_coor, y_coor, z_coor;

                	iss>>x_coord;
        	        iss>>y_coord;
	                iss>>z_coord;

                	x_coor=int(x_coord);
        	        y_coor=int(y_coord);
	                z_coor=int(z_coord);

			total_frame_coords[i][x_coor][y_coor][z_coor]=1;
		//	std::cerr<<total_frame_coords[i][x_coor][y_coor][z_coor]<<'\t';
		}
                ifs.close();
        }

	std::cerr<<"\nOn to prob computation\n";

	double *** spatial_probabilities=new double**[max_x+1];
	for (int i=min_x; i<=max_x; ++i){
		spatial_probabilities[i]=new double*[max_y+1];
		for (int j=min_y; j<=max_y; ++j){
			spatial_probabilities[i][j]=new double[max_z+1];
			for (int k=min_z; k<=max_z; ++k){
				spatial_probabilities[i][j][k]=0.0;
			}
		}
	}

	std::vector<int>::iterator it_x=master_coord_x_list.begin();
	std::vector<int>::iterator it_y=master_coord_y_list.begin();
	for (std::vector<int>::iterator it_z=master_coord_z_list.begin(); it_z!=master_coord_z_list.end(); ++it_z){
		
		int x_value=*it_x;
		int y_value=*it_y;
		int z_value=*it_z;

		for (int i=0; i<num_frames; ++i){
			if(total_frame_coords[i][x_value][y_value][z_value]==1){
				spatial_probabilities[x_value][y_value][z_value]+=1.0/double(num_frames);
				
			}
		}
		++it_x;
		++it_y;
	}

	double ****** joint_coord_occurence_prob = new double*****[max_x+1];

	int max_nitrogens=0;

	it_x=master_coord_x_list.begin();
	it_y=master_coord_y_list.begin();

	std::cerr<<"\nCalculating Joint Probabilities and Covariance Matrix\n";

	int n_squared=total_coords*total_coords;

	double data_minus_one=num_frames-1;

	double * covariance_matrix=new double[n_squared];
	double * conditional_probabilities=new double[n_squared];

	int two_dim_index=0;

	for (int i=min_x; i<=max_x; ++i){
                joint_coord_occurence_prob[i]=new double****[max_y+1];

                for (int j=min_y; j<=max_y; ++j){
                        joint_coord_occurence_prob[i][j]=new double***[max_z+1];

                        for (int k=min_z; k<=max_z; ++k){
                                joint_coord_occurence_prob[i][j][k]=new double**[max_x+1];

				int temporary_max=0;

				for (int x1=min_x; x1<=max_x; ++x1){


					joint_coord_occurence_prob[i][j][k][x1]=new double*[max_y+1];
	
					for (int y1=min_y; y1<=max_y; ++y1){
	                                        joint_coord_occurence_prob[i][j][k][x1][y1]=new double[max_z+1];
			
						for (int z1=min_z; z1<=max_z; ++z1){
						
							for (int frame_index=0; frame_index<num_frames; ++frame_index){
								
								if (total_frame_coords[frame_index][i][j][k]==total_frame_coords[frame_index][x1][y1][z1]){

									if(total_frame_coords[frame_index][i][j][k]==1){

										joint_coord_occurence_prob[i][j][k][x1][y1][z1]+=1.0/double(num_frames);
			
									}
								}

							}


							double temp_val=joint_coord_occurence_prob[i][j][k][x1][y1][z1];
							if (spatial_probabilities[x1][y1][z1]!=0){
								temp_val/=spatial_probabilities[x1][y1][z1];
								if (i!=x1 || j!=y1 || k!=z1){
									if(temp_val>=0.75) ++temporary_max;
								}
							}	
						}
					}
				if (temporary_max>max_nitrogens) max_nitrogens=temporary_max;
				}
			}	
                }
        }

	std::cerr<<"\nMaximun nitrogen count: "<<max_nitrogens<<'\n';

	//this loop writes an xyz traj file which cylces over N binding pocket coordinates to display
	//binary autocorrelation data given an arbitrary cut-off.
	//At first, arbitrary cut-off if hard-coded to 75%

	it_x=master_coord_x_list.begin();
        it_y=master_coord_y_list.begin();
        for (std::vector<int>::iterator it_z=master_coord_z_list.begin(); it_z!=master_coord_z_list.end(); ++it_z){

                int x1=*it_x;
                int y1=*it_y;
                int z1=*it_z;

		std::vector<int>::iterator it_x2=master_coord_x_list.begin();
                std::vector<int>::iterator it_y2=master_coord_y_list.begin();

		std::cout<<(max_nitrogens+1)<<'\n';
		std::cout<<"probabiliity=0.75%\n";

		std::cout<<"O"<<"\t"<<x1<<'\t'<<y1<<'\t'<<z1<<'\n';

		int nitrogens_this_frame=0;

                for (std::vector<int>::iterator it_z2=master_coord_z_list.begin(); it_z2!=master_coord_z_list.end(); ++it_z2){

                        int x2=*it_x2;
                        int y2=*it_y2;
                        int z2=*it_z2;

			double prob_a_given_b=joint_coord_occurence_prob[x1][y1][z1][x2][y2][z2];

                        if (spatial_probabilities[x2][y2][z2]>0) prob_a_given_b/=spatial_probabilities[x2][y2][z2];

			if(x1==x2 && y1==y2 && z1==z2){
				//do nothing, do not print
			}
			else{
				
				conditional_probabilities[two_dim_index]=prob_a_given_b;

				if (prob_a_given_b>=0.75){
					std::cout<<"N"<<"\t"<<x2<<'\t'<<y2<<'\t'<<z2<<'\n';
					++nitrogens_this_frame;
				}
				else{
					//do not print carbons;
				}
			}

			double covariance_entry=0.0;

			for (int i=0; i<num_frames; ++i){
				double avg_x=spatial_probabilities[x1][y1][z1];
				double avg_y=spatial_probabilities[x2][y2][z2];
				double x_variance=(double(total_frame_coords[i][x1][y1][z1])-avg_x);
				double y_variance=(double(total_frame_coords[i][x2][y2][z2])-avg_y);
				double numerator=y_variance*x_variance;
				covariance_entry+=numerator/data_minus_one;
			}


			covariance_matrix[two_dim_index]=covariance_entry;

			++it_x2;
			++it_y2;

			++two_dim_index;

		}

//		std::cout<<nitrogens_this_frame<<'\n';

		while(nitrogens_this_frame<max_nitrogens){
			if (nitrogens_this_frame==max_nitrogens) break;
			std::cout<<"N"<<'\t'<<"1000"<<'\t'<<"1000"<<'\t'<<"1000\n";
			++nitrogens_this_frame;	
		}
	
		++it_x;
		++it_y;

	}

	std::cerr<<"Solving for the first two principle components...\n";

	gsl_matrix_view cov_mat_view=gsl_matrix_view_array(covariance_matrix, total_coords, total_coords);
	gsl_vector *eval = gsl_vector_alloc(total_coords);
	gsl_matrix *evec= gsl_matrix_alloc(total_coords, total_coords);
	gsl_eigen_symmv_workspace * w= gsl_eigen_symmv_alloc(total_coords);
	gsl_eigen_symmv(&cov_mat_view.matrix, eval, evec, w);
	gsl_eigen_symmv_free(w);
	gsl_eigen_symmv_sort(eval, evec, GSL_EIGEN_SORT_ABS_DESC);

	for (int i=0; i<100; ++i){

		std::string filename="principal_component";
		std::stringstream ss;
		ss<<i;
		filename+=ss.str();
		filename+=".vect";

		double eval_i=gsl_vector_get(eval, i);
		gsl_vector_view evec_i=gsl_matrix_column(evec,i);

		std::cerr<<"eigenvalue "<<i<<" = "<<eval_i<<'\n';
		std::cerr<<filename<<'\n';
		
		FILE * f=fopen(filename.c_str(), "w");

		gsl_vector_fprintf(f, &evec_i.vector, "%g");

		fclose(f);

	}

	gsl_vector_free(eval);
	gsl_matrix_free(evec);

	return 0;
}
