: '
    Parameters for assembly and shear simualtions
'

#!/bin/bash

# Cifras significativas
cs=5;

# Parameters of the model
# Radii of the main and patch particles
r_Parti=0.5;
r_Patch=0.4;

# Main parameters of the simulation
phi=0.6;
CL_concentration=0.1;
N_particles=2000;

# Number of monomers and cross-linkers given concentration an total amount of patchy particles
N_CL=$(echo "scale=0; $CL_concentration * $N_particles" | bc);
N_CL=${N_CL%.*};
N_MO=$(( $N_particles - $N_CL ));

# Compute the size of the box to get he given packing fraction

# Pi
pi=$(echo "scale=$cs; 4*a(1)" | bc -l);
# Ratio of volume 
f_scl=$(echo "scale=$cs; 4/3" | bc);
# Ghost radius to ensure allocation in the box
r_ghost=$(echo "scale=$cs; $r_Parti + $(echo "scale=$cs; $r_Patch / 2" | bc)" | bc);
r_aux=$(echo "scale=$cs; $r_ghost * $r_ghost * $r_ghost" | bc);

# Compute the volume necesary for all the partilces
Vol_ghost=$(echo "scale=$cs; $f_scl * $pi * $r_ghost" | bc);
Vol_CLg=$(echo "scale=$cs; $Vol_ghost * $N_CL" | bc);
Vol_MOg=$(echo "scale=$cs; $Vol_ghost * $N_MO" | bc);
Vol_Totg=$(echo "scale=$cs; $Vol_CLg + $Vol_MOg" | bc);

# Get the total volume needed taking into account the packing fraction
Vol_Tot=$(echo "scale=$cs; $Vol_Totg / $phi" | bc);
# Translate the volume into length of a box
L_real=$(echo "scale=$cs; e( l($Vol_Tot)/3 )" | bc -l );
# Input parameter for LAMMPS script
L=$(echo "scale=$cs; $L_real / 2" | bc);

# Numerical parameters for LAMMPS simulation
steps=1500000;
tstep=0.005;
sstep=10000;
seed1=1234;
seed2=4321;
seed3=3124;

## Variables for shear deformation simulation
tstep_defor=0.001;
sstep_defor=10000;

shear_rate=0.02;
max_strain=12;
Nstep_per_strain=$(echo "scale=$cs; $(echo "scale=$cs; 1 / $shear_rate" | bc) * $(echo "scale=$cs; 1 / $tstep_defor" | bc)" | bc) ;
Nstep_per_strain=${Nstep_per_strain%.*};

shear_it=$(( $max_strain * $Nstep_per_strain));
Nsave=500;
Nave=$(echo "scale=$cs; 1 / $tstep_defor" | bc);
Nave=${Nave%.*};
cycles=4;

relaxTime1=$(( $Nstep_per_strain / 2));
relaxTime2=$(( 2 * $relaxTime1)); 
relaxTime3=$(( 2 * $relaxTime2));
relaxTime4=$(( 2 * $relaxTime3));

## Adapting the variables to be in the directories name
scl=1000;
shear_aux=$(echo "scale=$cs; $scl*$shear_rate" | bc);
shear_aux=${shear_aux%.*};
phi_aux=$(echo "scale=$cs; $scl*$phi" | bc);
phi_aux=${phi_aux%.*};
rt1_aux=$(echo "scale=0; $tstep_defor*$relaxTime1" | bc);
rt1_aux=${rt1_aux%.*};
rt2_aux=$(echo "scale=0; $tstep_defor*$relaxTime2" | bc);
rt2_aux=${rt2_aux%.*};
rt3_aux=$(echo "scale=0; $tstep_defor*$relaxTime3" | bc);
rt3_aux=${rt3_aux%.*};
rt4_aux=$(echo "scale=0; $tstep_defor*$relaxTime4" | bc);
rt4_aux=${rt4_aux%.*};

dir_name=""systemTotalFreePhi"${phi_aux}"CL"${N_CL}"MO"${N_MO}"ShearRate"${shear_aux}"RT1_"${rt1_aux}"RT2_"${rt2_aux}"RT3_"${rt3_aux}"RT4_"${rt4_aux}";

## Create the directory
mkdir -p ${dir_name};

cd ${dir_name};

## Create the README.md file
file_name="README.md";

#echo -e "# Parameters of the simulation\n\n" >> $file_name;
#touch $file_name;
echo -e "# Parameters of the simulation\n\n" >> $file_name;
echo -e "## Model and System parameters\n" >> $file_name;
echo -e ""- Central particle radius: "${r_Parti}" [sigma]"" >> $file_name;
echo -e ""- Patch particle radius: "${r_Patch}" [sigma]"" >> $file_name;
echo -e ""- Cross-Linker Concentration: "${CL_concentration}" >> $file_name;
echo -e ""- Number of particles: "${N_particles}" >> $file_name;
echo -e "\n ## Assembly System values \n" >> $file_name;
echo -e ""- Box Volume: "${Vol_Tot}" [sigma^3]"" >> $file_name;
echo -e ""- Packing fraction: "${phi}" >> $file_name;
echo -e ""- Time step: "${tstep}" [tau]"" >> $file_name;
echo -e ""- Number of time steps: "${steps}" >> $file_name;
echo -e ""- Save every "${sstep}" time step"" >> $file_name;
echo -e ""- Numer of Cross-Linkers: "${N_CL}" >> $file_name;
echo -e ""- Number of Monomers: "${N_MO}" >> $file_name;
echo -e "\n ## Shear Deformation System values \n" >> $file_name;
echo -e ""- Box Volume: "${Vol_Tot}"[sigma^3]"" >> $file_name;
echo -e ""- Time step: "${tstep_defor}" [tau]"" >> $file_name;
echo -e ""- Number of time steps: "${steps}" >> $file_name;
echo -e ""- Save every "${sstep_defor}" time step"" >> $file_name;
echo -e ""- Shear rate: "${shear_rate}" [1/tau]"" >> $file_name;
echo -e ""- Max deformation per cycle: "${max_strain}" >> $file_name;
echo -e ""- Relax steps 1: "${relaxTime1}" >> $file_name;
echo -e ""- Relax steps 2: "${relaxTime2}" >> $file_name;
echo -e ""- Relax steps 3: "${relaxTime3}" >> $file_name;
echo -e ""- Relax steps 4: "${relaxTime4}" >> $file_name;

cd ..; 

## Create the qsub file
cd sim;
file_name="runSim.sge";

rm -f $file_name

# Parameters for the cluster
nodes=4;

echo "#!/bin/bash">> $file_name
echo "#$ -S /bin/bash">> $file_name
echo "#$ -N G_eta${phi}_r${rnn}" >> $file_name
echo "#$ -pe mpich ${nodes}" >> $file_name
echo "#$ -cwd" >> $file_name
echo "#$ -j yes" >> $file_name
echo "" >> $file_name
echo ". /etc/profile.d/modules.sh">> $file_name
echo "module load gcc/8.3.0" >> $file_name
echo "module load openmpi/gcc/64/1.10.1" >> $file_name
#echo "module load lammps/gcc/4may22" >> $file_name
echo "" >> $file_name
echo "mpirun -n ${nodes} /mnt/MD1200B/cferreiro/fbenavides/lammps-2Aug2023/src/lmp_mpi -in in.assembly.lmp -var L ${L} -var NCL ${N_CL} -var NMO ${N_MO} -var seed1 ${seed1} -var seed2 ${seed2} -var seed3 ${seed3} -var steps ${steps} -var tstep ${tstep} -var sstep ${sstep}" >> $file_name


echo "$dir_name"
echo "$N_CL"
echo "$N_MO"
echo "$pi"
echo "$f_scl"
echo "$Vol_ghost"
echo "$Vol_CLg"
echo "$Vol_MOg"
echo "$Vol_Totg"
echo "$Vol_Tot"
echo "$L"
echo "$Nstep_per_strain"
echo "$shear_it"
echo "$Nave"
echo "$relaxTime1"
echo "$relaxTime2"
echo "$relaxTime3"
echo "$relaxTime4"

qsub $file_name
 
