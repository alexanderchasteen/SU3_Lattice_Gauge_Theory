    #include <iostream>
    #include <fstream>
    #include <iomanip>
    #include "Lattice.h"
    #include "SU3_Sampling.h"
    #include "rng.h"
    #include "Parameters.h"
    #include "statistics.h"

int main(){

    // Make the Raw Data File
    std::string filename1 = "Thermalization_Data.csv";
    std::ofstream Thermfile(filename1);

    if (!Thermfile) {
        std::cerr << "Error creating file!\n";
        return 1;
    }

    Thermfile << std::fixed << std::setprecision(6);

    for (int i=0; i<CONFIG; i++){
        if (i<CONFIG-1) Thermfile << coupling_constants[i] << ",";
        if (i == CONFIG-1) Thermfile << coupling_constants[i];
    }
    Thermfile << "\n";

    std::string filename2 = "RawMC_Data.csv";
    std::ofstream RawMCfile(filename2);

    if (!RawMCfile) {
        std::cerr << "Error creating file!\n";
        return 1;
    }

    RawMCfile << std::fixed << std::setprecision(6);

    for (int i=0; i<CONFIG; i++){
        if (i<CONFIG-1) RawMCfile << coupling_constants[i] << ",";
        if (i == CONFIG-1) RawMCfile << coupling_constants[i];
    }
    RawMCfile << "\n";

    

    std::string filename3 = "Poly_Corr_Real_Data.csv";
    std::ofstream poly(filename3);

    if (!poly) {
        std::cerr << "Error creating file!\n";
        return 1;
    }

    poly << std::fixed << std::setprecision(6);

    for (int i=1; i<maxdistance+1; i++){
        if (i<maxdistance) poly << i << ",";
        if (i == maxdistance) poly << i;
    }
    poly << "\n";

    std::string filename4 = "Poly_Corr_Imag_Data.csv";
    std::ofstream polycomplex(filename4);

    if (!polycomplex) {
        std::cerr << "Error creating file!\n";
        return 1;
    }

    polycomplex << std::fixed << std::setprecision(6);


    std::string filename5 = "Poly_Loop_Real_Data.csv";
    std::ofstream polyloop(filename5);

    if (!polyloop) {
        std::cerr << "Error creating file!\n";
        return 1;
    }

    polyloop << std::fixed << std::setprecision(6);


    
    for (int i=0; i<CONFIG; i++){
        if (i<CONFIG-1) polyloop << coupling_constants[i] << ",";
        if (i == CONFIG-1) polyloop << coupling_constants[i];
    }
    polyloop << "\n";

    std::string filename6 = "Poly_Loop_Complex_Data.csv";
    std::ofstream polyloopcomplex(filename6);

    if (!polyloopcomplex) {
        std::cerr << "Error creating file!\n";
        return 1;
    }

    polyloopcomplex << std::fixed << std::setprecision(6);

     for (int i=0; i<CONFIG; i++){
        if (i<CONFIG-1) polyloopcomplex << coupling_constants[i] << ",";
        if (i == CONFIG-1) polyloopcomplex << coupling_constants[i];
    }
    polyloopcomplex << "\n";

    for (int i=1; i<maxdistance+1; i++){
        if (i<maxdistance) polycomplex << i << ",";
        if (i == maxdistance) polycomplex << i;
    }
    polycomplex << "\n";

    // Make the array of link arrays (that will be helpful for tempering)
    std::array<Link_array, CONFIG> links_at_coupling;
    for (int i=0; i<CONFIG; i++){
        links_at_coupling[i] = Link_array(array_size, 0.0); // initialize vector
        cold_start_array(links_at_coupling[i]);            // cold start
    }
    std::cout << "Done Making" << std::endl;


    double action;

    for (int j = 0; j < thermal_sweeps; j++) {
        for (int i = 0; i < CONFIG; i++) {
            action = heatbath_update(links_at_coupling[i], coupling_constants[i]);
            if (i < CONFIG - 1) Thermfile << action << ",";
            else Thermfile << action;
        }
        Thermfile << "\n";
        if (j % 100 == 0) Thermfile.flush(); 
    }

  
    std::array<double, autocorrelation_sweeps> data;
    std::array<int,CONFIG> integrated_autocorr_times;

    for (int i = 0; i < CONFIG; i++){
         for (int j=0; j<autocorrelation_sweeps; j++) {
            action = heatbath_update(links_at_coupling[i], coupling_constants[i]);
            data[j] = action;
        }
        int tau = tau_int(data);
        std::cout<<tau<<std::endl;
        integrated_autocorr_times[i] = tau;
    }

    for (int j = 0; j < measurment_sweeps; j++) {
        complex PL; 
        for (int i = 0; i < CONFIG; i++) {
            for (int k =0; k <  integrated_autocorr_times[i];k++) action = heatbath_update(links_at_coupling[i], coupling_constants[i]);
            if (i < CONFIG - 1) RawMCfile << action << ",";
            else RawMCfile << action;
              
            // Now we just want to compute the polyakov loop at points m across the lattice and save them to the loop files to plt the scatter plt
            complex z=0;
            double x,y;
            z = poly_loop_at_spat_coord(links_at_coupling[i],0,0,0);
            // for (int i1 =0; i1<Spatial_Size;i1++)
            //     for (int i2 =0; i2<Spatial_Size;i2++)
            //         for (int i3 =0; i3<Spatial_Size;i3++){
            //             z = poly_loop_at_spat_coord(links_at_coupling[i],0,0,0);
            //         }
              
            x=z.real();
            y=z.imag();     
            if (i < CONFIG - 1){
                polyloop << x << ",";
                polyloopcomplex << y <<",";
                } 
            else{
                polyloop <<x;
                polyloopcomplex << y;
                } 
            }
            


        // Now we will fix the configuration to be CONFIG[i] and use that to measure the string tension
        // For lattice spacing a=1: 
        //<W_C>=3(beta/18)^n_A where n_A is the area of the loop measured
        // Similarly ln(<W_C>/3)=-n_tV(n_r)=-n_t (sigma) n_r where sigma is the string tension
        // So we want to keep track of -ln(<W_C>/3)=n_A sigma and show sigma = -ln(beta/18) to first order in beta



        // So for our loops let fix n_t and run through the possible n_r in one of the spatial directions. 
        // Number of points is linear wrt to the size of the lattice
        // for (int r = 1; r <  maxdistance+1; r++){
        //     PL= correlator_over_fixed_distance(links_at_coupling[0],r);
        //     if (r <  maxdistance){
        //         poly << PL.real()<< ",";
        //         polycomplex << PL.imag()<<",";
        //     }
        //     else {
        //         poly << PL.real();
        //         polycomplex << PL.imag();
        //     }
        // }


      
                        

        RawMCfile << "\n";
        poly<<"\n";
        polycomplex<<"\n";
        polyloop<<"\n";
        polyloopcomplex<<"\n";
    }

    Thermfile.close();
    RawMCfile.close();
    poly.close();
    polycomplex.close();
    polyloop.close();        
    polyloopcomplex.close();
    return 0;
}