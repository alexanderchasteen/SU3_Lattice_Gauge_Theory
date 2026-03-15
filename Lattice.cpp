    #include "Lattice.h"
    #include "SU3_Sampling.h"
    #include "Parameters.h"
    #include <iostream>



    // Gives the index in an array given the indices in our 7 index tensor for the lattice_data information (at a singular beta)
    // 4 of the indices are the spatiotemporal positions on the lattice, the 5th index the (positive) link direction (of 4) 
    // 6th runs over each of the 18 indices of SU(3) (9 for real, 9 for imaginary). 
    int flat_index(tensor_index tensor_index_array) {
        int i1=tensor_index_array[0];
        int i2=tensor_index_array[1];
        int i3=tensor_index_array[2];
        int i4=tensor_index_array[3];
        int i5=tensor_index_array[4];
        int i6=tensor_index_array[5];
        return (((((i1 * Spatial_Size + i2) * Spatial_Size + i3) * temporal_size + i4) *4  + i5) * 18 + i6);
    }


    // Inverts the flat index back to the position on the tensor
    tensor_index tensor_index_array(int index){
        int i6= index % 18;
        index /= 18;
        int i5 = index % 4;
        index /= 4;
        int i4 = index % temporal_size;
        index /= temporal_size;
        int i3 = index % Spatial_Size;
        index /= Spatial_Size;
        int i2 = index % Spatial_Size;
        index /= Spatial_Size;
        int i1 = index;
        return {i1,i2,i3,i4,i5,i6};
    }


    void set_array_value(Link_array& arr, tensor_index tensor_index_array, double value) {
        int idx = flat_index(tensor_index_array);
        arr[idx] = value;
        return;
    }

    double get_array_value(const Link_array& arr, tensor_index tensor_index_array) {
        int idx = flat_index(tensor_index_array);
        return arr[idx];
    }


    // Now we want to use these functions to get the matrix at a certain link along the lattice. This will take in a 5 index tensor location and we will loop over 
    // the 18 possibilities for the 6th array. Then we will put these in a matrix.


    // First let us come up with functions that both FLATTEN, and MATRIXify information. That is, from SU(3) to array with 18 components and vice versa

    SU3_array SU3_matrix_to_array(const SU3& matrix){
        SU3_array array={};
        int index = 0;
        for (int i=0;i<3; i++){
            for (int j=0;j<3;j++){
                double real_part = matrix(i,j).real();
                array[index] = real_part;
                index += 1;
                double imaginary_part = matrix(i,j).imag();
                array[index] = imaginary_part;
                index += 1;
            }
        }
        return array;
    }

    SU3 SU3_array_to_matrix(const SU3_array& array){
        SU3 SU3_matrix;
        int index=0;
        for (int i=0; i<3; i++){
            for (int j=0; j<3; j++){
                double real_part=array[index];
                index += 1;
                double imaginary_part=array[index];
                index += 1;
                SU3_matrix(i,j)=complex (real_part,imaginary_part);
            }
        }
        return SU3_matrix;
    }


    // Now we want code that gets the SU3 matrix at a given link. The matrix information is stored as an SU3_array and we want to convert it to a matrix. 
    // We use our array --> matrix converter 

    SU3 get_SU3_at_link(const Link_array&arr, link_index link_index_array){
        // extract the real component
        // Initialize a tensor index_array. The algorithm is to append the value 0<=i<=17 int to the link_index 
        // Put the link location in the first 5 elements of the total tensor index array
        tensor_index tensor_index_array = {};
        for (int i=0; i<5; i++){
            tensor_index_array[i] = link_index_array[i];
        }

        // Put the 6th link location in the last element of the total tensor index array (for some general index). 
        // But let us find the real and the imaginary index at the same time and add it Then find this element and add it to an array
        SU3_array SU3array={};
        for (int i=0; i<18;i+=2){
            // Get the real part
            tensor_index_array[5] = i;
            double real_part=get_array_value(arr,tensor_index_array);
            SU3array[i]=real_part;

            // Get the imaginary part
            tensor_index_array[5] = i+1;
            double imaginary_part=get_array_value(arr,tensor_index_array);
            SU3array[i+1]=imaginary_part;
        }
        SU3 SU3_matrix=SU3_array_to_matrix(SU3array);
        return SU3_matrix;
    }

    void set_link_SU3(Link_array& arr, link_index link_index_array, const SU3& SU3_matrix){
        tensor_index tensor_index_array = {};
        for (int i=0; i<5; i++){
            tensor_index_array[i] = link_index_array[i];
        }
        SU3_array SU3array = SU3_matrix_to_array(SU3_matrix);
        for (int i=0; i<18;i++){
            tensor_index_array[5] = i;
            set_array_value(arr, tensor_index_array, SU3array[i]);
        }
        return;
    }

    // Want to start our lattice cold. This means we want to initlaize all the links to I \in SU(3)
    // This is because the action is β/3 Σ_{n ∈ Λ} Σ_{μ < ν} Re(Tr(1-U_{μν}(n))) so at the identity S=0
    // Since in the path integral it is e^-S as the statisticall factor, this correlates to a 0 energy state whichh exists at high temperatures

    void cold_start_array(Link_array& arr){
        for (int i1 = 0; i1 < Spatial_Size; i1++)
        for (int i2 = 0; i2 < Spatial_Size; i2++)
        for (int i3 = 0; i3 < Spatial_Size; i3++)
        for (int i4 = 0; i4 < temporal_size; i4++)
        for (int i5 = 0; i5 < 4; i5++){
            link_index link_index_array = {i1,i2,i3,i4,i5};
            SU3 I;
            I.setIdentity();
            set_link_SU3(arr, link_index_array,I);
        }
        return;
    }


    void moveup(lattice_index& lattice_index_array, int link_direction){
        lattice_index_array[link_direction]+=1;
        if (link_direction==3){
            if (lattice_index_array[link_direction]>=temporal_size) lattice_index_array[link_direction] -= temporal_size;
        }
        else{
            if (lattice_index_array[link_direction]>=Spatial_Size) lattice_index_array[link_direction] -= Spatial_Size;
        }
        return;
    }


    void movedown(lattice_index& lattice_index_array, int link_direction){
        if (lattice_index_array[link_direction] == 0) {
            if (link_direction == 3)
                lattice_index_array[link_direction] = temporal_size - 1;
            else
                lattice_index_array[link_direction] = Spatial_Size - 1;
        } else {
            lattice_index_array[link_direction] -= 1;
        }
        return;
    }



    link_index combine_lattice_index_with_direction(const lattice_index& latice_index_array, int direction){
        link_index total_link_index_array;
        for (int i=0; i<4; i++){
            total_link_index_array[i]=latice_index_array[i];
        }
        total_link_index_array[4]=direction;
        return total_link_index_array;
    }
    

        SU3 compute_staple_sum_at_link(const Link_array& arr, const link_index& link_index_array){
            // The staple sum is given in gattringer. However note that for the bottom part of the loop, this is the staple sum BUT the staple we compute is actually the dagger of the "correct" staple
            // This does not matter since we will take the real part and the trace. So it is equivalent to use this matrix. 
            // SO IMPORTANT NOTE this is physically unmeaningful without includiong the real part of the trace


            // Seperate into the lattice indices (4) and the direction of the link
            lattice_index lattice_index_array = {link_index_array[0], link_index_array[1], link_index_array[2], link_index_array[3]};
            int d=link_index_array[4];
            link_index local_link_index_array;

            // Takes care of the staplesum calculation
            SU3 staple;
            staple.setZero();
            SU3 staple_sum; 
            staple_sum.setZero();
            for (int dperp=0;dperp<4;dperp++){
                if (dperp!=d){
                    // make temperoray lattice_index_array storage because each time we loop we want to start back at the beginning point 
                    lattice_index tmp = lattice_index_array;
                    // Bottom staple
                
                    //Let us grab U_{ν}(n-ν)
                    movedown(tmp,dperp);
                    local_link_index_array = combine_lattice_index_with_direction(tmp,dperp);
                    staple = get_SU3_at_link(arr,local_link_index_array); 
                    
                    //Let us grab U_{-μ}(n+μ-ν)=U_{μ}^dagger(n-ν)
                    local_link_index_array = combine_lattice_index_with_direction(tmp,d);
                    staple =  (get_SU3_at_link(arr,local_link_index_array).adjoint())*staple;

                    //Let us grab U_{-ν}(n+μ)=U_{ν}^dagger}(n+μ-ν)
                    moveup(tmp,d);
                    local_link_index_array = combine_lattice_index_with_direction(tmp,dperp);
                    staple = (get_SU3_at_link(arr,local_link_index_array).adjoint())*staple;
                    
                    staple_sum += staple;

                    // Top staple
                    //Let us grab U_{-ν}(n+ν)=U_{ν}(n)^dagger
                    moveup(tmp,dperp);
                    movedown(tmp,d);
                    
                    assert(tmp == lattice_index_array);
                    local_link_index_array = combine_lattice_index_with_direction(tmp,dperp);
                    staple = get_SU3_at_link(arr,local_link_index_array).adjoint(); 
                    
                    //Let us grab U_{-μ}(n+μ+ν)
                    moveup(tmp,dperp);
                    local_link_index_array = combine_lattice_index_with_direction(tmp,d);
                    staple =  (get_SU3_at_link(arr,local_link_index_array).adjoint())*staple;

                    //Let us grab U_{ν}(n+μ)
                    moveup(tmp,d);
                    movedown(tmp,dperp);
                    local_link_index_array = combine_lattice_index_with_direction(tmp,dperp);
                    staple = get_SU3_at_link(arr,local_link_index_array)*staple;

                    staple_sum += staple;
                    movedown(tmp,d);
                    assert(tmp == lattice_index_array);
                }
            }
            return staple_sum;
        }


    SU2 R_block(const SU3& SU3_matrix){
        SU2 block = SU3_matrix.block<2,2>(0,0);
        return block;
    }

    SU3 R_block_to_SU3(const SU2& SU2_matrix){
        SU3 matrix;
        matrix(0,0) = SU2_matrix(0,0);
        matrix(0,1) = SU2_matrix(0,1);
        matrix(0,2) = 0.0;
        matrix(1,0) = SU2_matrix(1,0);
        matrix(1,1) = SU2_matrix(1,1);
        matrix(1,2) = 0.0;
        matrix(2,0) = 0.0;
        matrix(2,1) = 0.0;
        matrix(2,2) = 1.0; 
        return matrix;    
    }

    SU2 S_block(const SU3& SU3_matrix){
        SU2 block;
        block(0,0) = SU3_matrix(0,0);
        block(0,1) = SU3_matrix(0,2);
        block(1,0) = SU3_matrix(2,0);
        block(1,1) = SU3_matrix(2,2);
        return block;
    }

    SU3 S_block_to_SU3(const SU2& SU2_matrix){
        SU3 matrix;
        matrix(0,0) = SU2_matrix(0,0);
        matrix(0,1) = 0;
        matrix(0,2) = SU2_matrix(0,1);
        matrix(1,0) = 0.0;
        matrix(1,1) = 1.0;
        matrix(1,2) = 0.0;
        matrix(2,0) = SU2_matrix(1,0);
        matrix(2,1) = 0.0;
        matrix(2,2) = SU2_matrix(1,1); 
        return matrix;    
    }

    SU2 T_block(const SU3& SU3_matrix){
        SU2 block = SU3_matrix.block<2,2>(1,1);
        return block;
    }

    SU3 T_block_to_SU3(const SU2& SU2_matrix){
        SU3 matrix;
        matrix(0,0) = 1.0;
        matrix(0,1) = 0;
        matrix(0,2) = 0;
        matrix(1,0) = 0.0;
        matrix(1,1) = SU2_matrix(0,0);
        matrix(1,2) = SU2_matrix(0,1);
        matrix(2,0) = 0.0;
        matrix(2,1) = SU2_matrix(1,0);
        matrix(2,2) = SU2_matrix(1,1); 
        return matrix;    
    }

    SU2 R_block_to_2by2_unitary(const SU3& SU3_matrix){
        // Note that since we want to sample our random matrix with respect to the 2 by 2 block that we have to project to a unitary matrix, we need to compute the determinant of the projection
        // The projection conserves the real trace so it is purely mathematical trick


        // First extract the two by two block
        SU2 block = R_block(SU3_matrix);
        SU2 Hermetian_part = (block+block.adjoint())/2.0;
        SU2 anti_hermetian_part = (block-block.adjoint())/2.0;
        complex b0 = (Hermetian_part.trace())/2.0;
        complex b1 = (anti_hermetian_part*Pauli::sigma_x).trace() / 2.0;
        complex b2 = (anti_hermetian_part*Pauli::sigma_y).trace() / 2.0;
        complex b3 = (anti_hermetian_part*Pauli::sigma_z).trace() / 2.0;
        SU2 matrix = b0*Pauli::I + b1*Pauli::sigma_x + b2*Pauli::sigma_y +b3*Pauli::sigma_z;
        return matrix; 
    }
    // Similarly for S and T type

    SU2 S_block_to_2by2_unitary(const SU3& SU3_matrix){
        SU2 block = S_block(SU3_matrix);
        SU2 Hermetian_part = (block+block.adjoint())/2.0;
        SU2 anti_hermetian_part = (block-block.adjoint())/2.0;
        complex b0 = (Hermetian_part.trace())/2.0;
        complex b1 = (anti_hermetian_part*Pauli::sigma_x).trace() / 2.0;
        complex b2 = (anti_hermetian_part*Pauli::sigma_y).trace() / 2.0;
        complex b3 = (anti_hermetian_part*Pauli::sigma_z).trace() / 2.0;
        SU2 matrix = b0*Pauli::I + b1*Pauli::sigma_x + b2*Pauli::sigma_y +b3*Pauli::sigma_z;
        return matrix; 
    }

    SU2 T_block_to_2by2_unitary(const SU3& SU3_matrix){
        SU2 block = T_block(SU3_matrix);
        SU2 Hermetian_part = (block+block.adjoint())/2.0;
        SU2 anti_hermetian_part = (block-block.adjoint())/2.0;
        complex b0 = (Hermetian_part.trace())/2.0;
        complex b1 = (anti_hermetian_part*Pauli::sigma_x).trace() / 2.0;
        complex b2 = (anti_hermetian_part*Pauli::sigma_y).trace() / 2.0;
        complex b3 = (anti_hermetian_part*Pauli::sigma_z).trace() / 2.0;
        SU2 matrix = b0*Pauli::I + b1*Pauli::sigma_x + b2*Pauli::sigma_y +b3*Pauli::sigma_z;
        return matrix; 
    }



    SU3 type_R_heatbath(const Link_array& arr, link_index link_index_array, double beta, SU3 U, SU3 A){
        constexpr double small = 1e-14;
        // We want to compute the staple sum A once and we will feed into it
        // The link U to update ach time. However after a type R heatbath, we change the link U. So let us have this function
        // Give out the new link U and keep the same format for type S and T so that we can just after each update, get the new link,
        // use the already computed staplesum, and after we do T we will set the link to the last one. 

        // That is 
        // SU3 A = compute_staple_sum_at_link(arr, link_index_array);
        // SU3 U = get_SU3_at_link(arr, link_index_array);
        SU3 UA = U*A;
        SU2 W_block = R_block_to_2by2_unitary(UA);
        double a = std::sqrt(W_block.determinant().real());

        SU2 R_2by2;
        SU3 R_3by3;
        if (a <= small){
            R_2by2 = Random_SU2_generator();
            R_3by3 = R_block_to_SU3(R_2by2);
            U = R_3by3*U; 
        }

        else{
            W_block = W_block/a;
            double adjusted_beta = 2*beta/3;
            SU2 X = SU2_generator(a,adjusted_beta);
            R_2by2 = X*W_block.adjoint();
            R_3by3 = R_block_to_SU3(R_2by2);    
            U = R_3by3 * U; 
        }    
        return U ;
    }

    SU3 type_S_heatbath(const Link_array& arr, link_index link_index_array, double beta, SU3 U, SU3 A){
        constexpr double small = 1e-14;
        SU3 UA = U*A;
        SU2 W_block = S_block_to_2by2_unitary(UA);
        double a = std::sqrt(W_block.determinant().real());
        
        SU2 S_2by2;
        SU3 S_3by3;
        if (a <= small){
            S_2by2 = Random_SU2_generator();
            S_3by3 = S_block_to_SU3(S_2by2);
            U = S_3by3*U; 
        }

        else{
            W_block = W_block/a;
            double adjusted_beta = 2*beta/3;
            SU2 X = SU2_generator(a,adjusted_beta);
            S_2by2 = X*W_block.adjoint();
            S_3by3 = S_block_to_SU3(S_2by2);    
            U = S_3by3 * U; 
        }    
        return U ;
    }

    SU3 type_T_heatbath(const Link_array& arr, link_index link_index_array, double beta, SU3 U, SU3 A){
        constexpr double small = 1e-14;
        SU3 UA = U*A;
        SU2 W_block = T_block_to_2by2_unitary(UA);
        double a = std::sqrt(W_block.determinant().real());

        SU2 T_2by2;
        SU3 T_3by3;
        if (a <= small){
            T_2by2 = Random_SU2_generator();
            T_3by3 = T_block_to_SU3(T_2by2);
            U = T_3by3*U; 
        }

        else{
            W_block = W_block/a;
            double adjusted_beta = 2*beta/3;
            SU2 X = SU2_generator(a,adjusted_beta);
            T_2by2 = X*W_block.adjoint();
            T_3by3 = T_block_to_SU3(T_2by2);    
            U = T_3by3 * U; 
        }    
        return U ;
    }


    double single_link_heatbath(Link_array& arr, const link_index& link_index_array, double beta){
        SU3 U = get_SU3_at_link(arr,link_index_array);
        SU3 A = compute_staple_sum_at_link(arr, link_index_array);
        // Update the link according to R,S, then T
        U = type_R_heatbath(arr, link_index_array, beta, U, A); 
        
        U = type_S_heatbath(arr, link_index_array, beta, U, A);
        U = type_T_heatbath(arr, link_index_array, beta, U, A);
        set_link_SU3(arr, link_index_array, U);
        double action;
        SU3 I; I.setIdentity();

        // PLAQ 
        // Running through the rest of the code, the observable measured is the average plaquette 
        action = ((U*A).trace().real())/3;


        // REAL ACTION
        // Running through the rest of the code, the observable measured is the average action (per plaq.)
        // action = (beta/3)*(18-(U*A).trace().real());
        return action; 
    }



    double heatbath_update(Link_array& arr, double beta){
        double action = 0.0;
        for (int i1 = 0; i1< Spatial_Size; i1++)
        for (int i2 = 0; i2<Spatial_Size; i2++)
        for (int i3 = 0; i3<Spatial_Size; i3++)
        for (int i4 = 0; i4<temporal_size; i4++)
        for (int d = 0; d<4; d++){
            link_index link_index_array = {i1,i2,i3,i4,d};
            action += single_link_heatbath(arr,link_index_array,beta);
        }
        action /= (Nplaq*4);
        // because overcount the same plaquette 4 times when we have these periodic boundary conditions
        return action;  
    }

    // Saving the configurations of the link array 
    void save_lattice_config(const Link_array& array, const std::string& filename) {
        std::ofstream out(filename, std::ios::binary);
        out.write(reinterpret_cast<const char*>(array.data()),
                array.size() * sizeof(double));
    }


    complex poly_loop_at_spat_coord(const Link_array& arr, int x, int y, int z){
        int T = 3;
        SU3 wl_matrix;
        wl_matrix.setIdentity();
        lattice_index lat= {x,y,z,0};
        for (int i = 0; i < temporal_size; i++) {
            link_index idx = combine_lattice_index_with_direction(lat, T);
            wl_matrix = wl_matrix * get_SU3_at_link(arr, idx); // RIGHT multiply
            moveup(lat, T);
            // Periodic boundary
            }
        std::complex<double> tr = wl_matrix.trace();
        return tr;
    }

    complex compute_correlator(const Link_array& arr, spat_index m, spat_index n){
        int x = m[0];
        int y = m[1];
        int z = m[2];
        complex poly1 = poly_loop_at_spat_coord(arr,x,y,z);

        x = n[0];
        y = n[1];
        z = n[2];
        complex poly2 = poly_loop_at_spat_coord(arr,x,y,z);
        complex poly3 = poly1*std::conj(poly2);
        return poly3;
    }

    complex correlator_over_fixed_distance(const Link_array& arr, int r){
        int N = 0;
        complex polyakov_data = 0;

        for (int i=0; i<Spatial_Size; i++)
        for (int j=0; j<Spatial_Size; j++)
        for (int k=0; k<Spatial_Size; k++){
            spat_index m = {i,j,k};
            for (int i1=0; i1<Spatial_Size; i1++)
            for (int j1=0; j1<Spatial_Size; j1++)
            for (int k1=0; k1<Spatial_Size; k1++){
                // minimal image distance
                int dx = (i1-i + Spatial_Size/2) % Spatial_Size - Spatial_Size/2;
                int dy = (j1-j + Spatial_Size/2) % Spatial_Size - Spatial_Size/2;
                int dz = (k1-k + Spatial_Size/2) % Spatial_Size - Spatial_Size/2;

                int d_squared = dx*dx + dy*dy + dz*dz;
                if (r*r == d_squared){
                    spat_index n = {i1,j1,k1};
                    polyakov_data += compute_correlator(arr, m, n);
                    N++;
                }
            }
        }

        if (N != 0) polyakov_data /= static_cast<double>(N);
        else polyakov_data = complex(0.0,0.0);
        return polyakov_data;
}



     
    void check_unitarity(const Link_array& arr) {
        constexpr double eps = 1e-12; // tolerance for numerical errors

        for (int i1 = 0; i1 < Spatial_Size; i1++)
        for (int i2 = 0; i2 < Spatial_Size; i2++)
        for (int i3 = 0; i3 < Spatial_Size; i3++)
        for (int i4 = 0; i4 < temporal_size; i4++)
        for (int d = 0; d < 4; d++) {
            link_index idx = {i1, i2, i3, i4, d};
            SU3 U = get_SU3_at_link(arr, idx);

            // Check determinant
            std::complex<double> det = U.determinant();
            if (std::abs(det - 1.0) > eps) {
                std::cerr << "WARNING: det(U) != 1 at link "
                        << "(" << i1 << "," << i2 << "," << i3 << "," << i4
                        << "), dir=" << d 
                        << " det(U) = " << det << std::endl;
            }
            assert(std::abs(det - 1.0) < eps);

            // Check unitarity: U^dagger * U = I
            SU3 identity_check = U.adjoint() * U;
            for (int i = 0; i < 3; i++)
            for (int j = 0; j < 3; j++) {
                std::complex<double> expected = (i == j) ? 1.0 : 0.0;
                if (std::abs(identity_check(i, j) - expected) > eps) {
                    std::cerr << "WARNING: U^dagger * U != I at link "
                            << "(" << i1 << "," << i2 << "," << i3 << "," << i4
                            << "), dir=" << d 
                            << " element (" << i << "," << j << ") = " 
                            << identity_check(i, j) << std::endl;
                }
                assert(std::abs(identity_check(i, j) - expected) < eps);
            }
        }

        std::cout << "All links are unitary and have determinant 1 within tolerance." << std::endl;
    }

