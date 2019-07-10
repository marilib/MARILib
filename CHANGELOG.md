
Changes from version 1.0.0 are :

1- Function reynolds_number(...) is now based on classical Sutherland viscosity model
   This change impacts slightly friction drag evaluation and all other related variables depending on design process

2- Design variables for Handling Qualities are now conform to the AIAA paper mentionned above : wing_x_root, htp_area, vtp_area
   This change impacts slightly statistical empennage sizing and all other related variables depending on design process

3- Scenarios have been updated to solve a problem of interaction between statistical empennage sizing and HQ based sizing

4- Library atmosphere_grad.py is beeing updated to get rid of non derivable points (in progress)

5- Name of functions "eval_vtp_coupling" and "eval_htp_coupling" changed resp. to "eval_vtp_statistical_sizing" and "eval_htp_statistical_sizing"

6- Add function "eval_aircraft_geom_analysis" in assembly.py

7- Suppress propulsive_architecture as input of Aircraft() (bug)

8- Add sub-folder test_cases in the directory examples

9- Add methods in Aircraft object to print its content, example in "test_cases.turbofan_sequence.py"
    (warning : all figures in standard units : m, kg, s and aggregates)

10- Correction of a bug on battery_price

11- The cost of electricity for an hybrid architecture is now computed separetly and exposed in the attribute elec_cost

12- Functions in math.py have been simplified and modified to allow using Autograd library

13- Propulsive architecture index has been changed from integer to string ("TF", "PTE1")

14- New sub-processes have been added into assembly.py 
