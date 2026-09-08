import os
import csv
import subprocess
import numpy as np
from datetime import datetime
from pathlib import Path
from dataclasses import dataclass, field
from typing import Tuple, List

# =============================================================================
# 1. PARAMETER DATA STRUCTURES (DOCUMENTATION KEPT IN PYTHON)
# =============================================================================

@dataclass
class ContCuConfig:
    """
    Control parameters for microMegas. All documentation is kept here 
    so the generated 'ContCu' file contains strictly values/flags.
    """
    sim_status: int = 0           # Simulation status: 0=New ; 1= restart ; 2= with particles ; 3= with irradiation loops
    mode_deformation: int = 3     # 0=strain rate; 1=make carto; 2=run carto; 3=stress rate; 4=cst stress; 5=fatigue; 6=metalofute; 7=metalofute2; 8=FE mesh load; 9=creep
    echelle: float = 13.5         # Reference scale, i.e. size of elementary screw vectors in BVD.XXX (Burgers vector unit) intially 13.5
    shear: str = 'T'              # True : shear stress resolved on the slip system with highest Schmid factor; False = uniaxial stress
    sigma0: float = 95.0          # Initial stress (MPa)
    deltat0: str = '0.125D-11'    # Elementary time step (s)
    sigma_point: str = '5D10'     # Elementary stress increments (MPa s-1)      (mode_deformation=3,6,7)
    epsilon_point: float = 0.0    # Imposed strain rate (s-1)                   (mode_deformation=0,5,6,7,9)
    raid: float = 1.5             # Apparent Young modulus                      (mode_deformation=0,5,6,7)
    epsmax: str = '3.0D-5'        # Maximum plastic strain per cycles           (mode_deformation=5)
    fsig: str = '1.D0'            # loading sign for the first cycle            (mode_deformation=5,6,7)
    linten: int = 1               # Local line tension definition:  0->Friedel;    1->DeWit;   2->Foreman;    3->Tabulated;  4->Mohles
    gldev: str = 'T'              # Cross slip activation key : True = cross-slip active
    key_nucleation: str = 'F'     # Nucleation activation key : True = nucleation active (control parameters are found the in nucl_init file )
    key_infline: str = 'F'        # PBC induced infinite lines : True = infinite lines are always pinned at least at a point in the volume
    key_crack: str = 'F'          # Add the stress field of a sharp crack in the simulated volume  (T/F=on/off)
    tensile_axis: str = '0 0 1'   # Tensile or compression uniaxial test direction (Miller indices)
    temperature: float = 300.0    # Temperature (Kelvin)
    facteur_depmax: int = 100     # Maximum segment displacement in a simulation step (Burgers vector unit)
    nstat_control: int = 50       # Number of steps accounted for in the simulation control procedure (mode_deformation = 0,5,6,7)
    relax_tl: int = 0             # Number of steps ascribed to initial relaxation when involving line tension only
    relax_int: int = 100          # Number of steps ascribed to initial relaxation under no applied load and with no contact reactions allowed
    relax_reac: int = 100         # Number of steps ascribed to initial relaxation under no applied loading
    nstep: int = 10000            # Number of steps of simulation
    ldis_act: float = 0.25        # Segment length (micron) at which line discretization is systematically tested  (micron).
    ldis_nact: float = 0.25       # Same as above for segments belonging to an inactive slip systems (micron).
    period: int = 10              # Number of steps before force on waiting (quasi-immobile) segments is recalculated
    krc: int = 10                 # Number of waiting steps before long-range contribution to internal stress is recalculated (if Greengard method used)
    l_boite: float = -1.0         # Linear mean size of domains defined used in Greengard's method (micron), it is dynamically defined when negative
    pbci_dim: str = '-1 -1 -1'    # Number of replicas used in Greengard method (in x,y and z directions), symmetric long rang solution is imposed if negative
    ab_length_seuil: float = 5.0  # Maximum segment length neglected in the long range contribution (if Greengard method used) (Echelle unit see line 3)
    dcfj: int = 50                # Minimum Distance at which stress is calculated on the segments connected to a junction (unit: Burgers vector)
    gb: int = 0                   # interface and surface definitions: 0-> inactivates, 1-> domain definition, 2->spherical, 3->regular 3D network (polycristal)
    tauint_limite: float = 100.0  # Critical stress at which segments are considered as in a singular field (MPa)
    kisauve: int = 1000           # Write periodicity of segment configurations and information needed to restart computation
    kstats: int = 200             # Write periodicity of results
    kkim: int = 200               # Write periodicity of the trajectory film
    kpredraw: int = 100           # Periodicity of refreshment of the graphical interface in gmm mode
    shift_rotation: int = 0       # Key for translation and rotation of the simulated volume (see the shift_rotation file)
    iterinfo: int = -6542         # Step of debugging (no debugging if negative)
    sysinfo: int = -2             # Slip system of interest in debugging procedure (needed in simulations with many segments)


@dataclass
class SimulationConfig:
    """Master configuration encompassing geometry and control parameters."""
    target_length_nm: float
    angle_deg: int
    box_dim_x: int
    box_dim_y: int
    box_dim_z: int
    use_graphics: bool = False
    start_coord: Tuple[int, int, int] = (1200, 1200, 1200)
    material_name: str = "Cu"
    cont_cu: ContCuConfig = field(default_factory=ContCuConfig)


# =============================================================================
# 2. SEGMENT GENERATOR
# =============================================================================

class DislocationGenerator:
    def __init__(self):
        """Finalized Generator using relative indexing (-1) and 'ligne' vectors."""
        # Physical lengths (nm)
        self.l_edge = 5.9684    
        self.l_screw = 3.4459   

        # Native BVD Indices
        self.idx_screw = 1   # 0 degrees
        self.idx_60deg = 8   # 60 degrees
        self.idx_edge  = 7   # 90 degrees

        # Static Geometry Vectors ('ligne')
        self.step_edge = {'idx': 7, 'dx': -2, 'dy': 4, 'dz': -2}
        self.step_screw = {'idx': 1, 'dx': -2, 'dy': 0, 'dz': 2}

    def generate_smart_configuration(self, beta_deg, L_target_nm, start_coord, box_dim, output_filepath):
        reference_scale = 1.2183
        cx, cy, cz = start_coord
        os.makedirs(os.path.dirname(output_filepath) or '.', exist_ok=True)

        for i in range(len(box_dim)):
            box_dim[i] = (((box_dim[i]/reference_scale) + 7) // 8) * 8  # For converting nm dimensions to BVD units, ensuring divisibility by 8

        with open(output_filepath, 'w') as f:
            # --- Header ---
            f.write(" 1 1 1 1 1 1 1 1 1 1 1 1\n")
            
            # --- Native Bypass (0, 60, 90) ---
            if beta_deg in [0, 60, 90]:
                if beta_deg == 0:
                    char_idx = self.idx_screw
                    length = int(round(L_target_nm / self.l_screw))
                elif beta_deg == 60:
                    char_idx = self.idx_60deg
                    length = int(round(L_target_nm / self.l_screw)) 
                elif beta_deg == 90:
                    char_idx = self.idx_edge
                    length = int(round(L_target_nm / self.l_edge))

                f.write(f"                   1\n")
                # Unpacked independent axes here
                f.write(f"                {box_dim[0]}                 {box_dim[1]}                 {box_dim[2]}\n")
                
                # Relative indices: 0 means no neighbors (Pinned)
                line = f"    1     {int(cx):>4}      {int(cy):>4}      {int(cz):>4}   {length:>2}       {char_idx}       0     0     0     0     F    0    0\n"
                f.write(line)
                print(f"✅ Generated Native: {beta_deg}°")

            # --- Grouped Staircase Generator (All others) ---
            else:
                beta_rad = np.radians(beta_deg)
                N_edge = int(round((L_target_nm * np.sin(beta_rad)) / self.l_edge))
                N_screw = int(round((L_target_nm * np.cos(beta_rad)) / self.l_screw))
                total_units = N_edge + N_screw

                # Build sequence
                sequence = []
                e_placed, s_placed = 0, 0
                for i in range(total_units):
                    if ((i + 1) * N_edge / total_units - e_placed) > ((i + 1) * N_screw / total_units - s_placed):
                        sequence.append('E')
                        e_placed += 1
                    else:
                        sequence.append('S')
                        s_placed += 1

                # Group blocks
                grouped = []
                curr_type, curr_len = sequence[0], 1
                for s in sequence[1:]:
                    if s == curr_type: curr_len += 1
                    else: 
                        grouped.append((curr_type, curr_len))
                        curr_type, curr_len = s, 1
                grouped.append((curr_type, curr_len))

                f.write(f"                  {len(grouped)}\n")
                f.write(f"                {box_dim[0]}                 {box_dim[1]}                 {box_dim[2]}\n")

                for i, (seg_type, length) in enumerate(grouped):
                    if seg_type == 'E':
                        char_idx, step_dx, step_dy, step_dz = self.step_edge['idx'], self.step_edge['dx'], self.step_edge['dy'], self.step_edge['dz']
                    else:
                        char_idx, step_dx, step_dy, step_dz = self.step_screw['idx'], self.step_screw['dx'], self.step_screw['dy'], self.step_screw['dz']

                    ID = i + 1
                    # Relative indices: 0 for boundaries, -1 for linked connectivity
                    n1 = 0 if ID == 1 else -1
                    n2 = 0 if ID == len(grouped) else -1

                    # Fixed typo in string template ({length:>2} instead of {length:>2])
                    line = f" {ID:>4}     {int(cx):>4}      {int(cy):>4}      {int(cz):>4}   {length:>2}       {int(char_idx)}      {n1:>3}   {n1:>3}    {n2:>2}    {n2:>2}     F    0    0\n"
                    f.write(line)

                    cx += (step_dx * length); cy += (step_dy * length); cz += (step_dz * length)

                print(f"✅ Generated Staircase: {beta_deg}° | {len(grouped)} Groups")


# =============================================================================
# 3. EPHEMERAL EXECUTION WRAPPER
# =============================================================================

class MicroMegasWrapper:
    def __init__(self, base_repo_path: str):
        self.base_repo = Path(base_repo_path).resolve()
        self.generator = DislocationGenerator()

    def run(self, config: SimulationConfig) -> str:
        """Runs the simulation directly in the repository and converts stat.txt to CSV."""
        # Work directly in the base repository
        work_dir = self.base_repo
        
        in_dir = work_dir / "dd" / "in"
        in_dir.mkdir(parents=True, exist_ok=True)
        
        # Generate a specific timestamp for this run
        timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
        
        # 1. Write Configuration Files (directly to the repo) with timestamps
        cont_filename = f"ContCu_{timestamp}"
        self._write_contcu(in_dir / cont_filename, config.cont_cu)
        
        seg_filename = f"Seg_{config.angle_deg}deg_{timestamp}"
        self.generator.generate_smart_configuration(
            beta_deg=config.angle_deg,
            L_target_nm=config.target_length_nm,
            start_coord=config.start_coord,
            box_dim=(config.box_dim_x, config.box_dim_y, config.box_dim_z),
            output_filepath=str(in_dir / seg_filename)
        )
        
        self._write_input_dd(in_dir / "input.dd", config.material_name, cont_filename, seg_filename)

        # 2. Compile and Run
        self._execute_engine(work_dir, config.use_graphics)

        # 3. Parse stat.txt and save as CSV
        return self._export_stat_csv(work_dir, timestamp, config)

    def _write_contcu(self, filepath: Path, c: ContCuConfig):
        """Writes the control file cleanly, strictly values to appease Fortran."""
        lines = [
            f"{c.sim_status}", f"{c.mode_deformation}", f"{c.echelle}", f"{c.shear}",
            f"{c.sigma0}", f"{c.deltat0}", f"{c.sigma_point}", f"{c.epsilon_point}",
            f"{c.raid}", f"{c.epsmax}", f"{c.fsig}", f"{c.linten}", f"{c.gldev}",
            f"{c.key_nucleation}", f"{c.key_infline}", f"{c.key_crack}",
            f"{c.tensile_axis}", f"{c.temperature}", f"{c.facteur_depmax}",
            f"{c.nstat_control}", f"{c.relax_tl}", f"{c.relax_int}", f"{c.relax_reac}",
            f"{c.nstep}", f"{c.ldis_act}", f"{c.ldis_nact}", f"{c.period}",
            f"{c.krc}", f"{c.l_boite}", f"{c.pbci_dim}", f"{c.ab_length_seuil}",
            f"{c.dcfj}", f"{c.gb}", f"{c.tauint_limite}", f"{c.kisauve}",
            f"{c.kstats}", f"{c.kkim}", f"{c.kpredraw}", f"{c.shift_rotation}",
            f"{c.iterinfo}", f"{c.sysinfo}"
        ]
        with open(filepath, 'w') as f:
            f.write("\n".join(lines) + "\n")

    def _write_input_dd(self, filepath: Path, mat: str, cont: str, seg: str):
        """Writes the pointer file with the required header and formatting."""
        input_content = f"""The three files needed to run a microMegas simulation.
These input files must be defined in the directory "dd/in"
In those files are defined all the info needed to run a simulation.

The sequence order must be respected, but the files name can be modified:
>>> The material variables file
>>> The control parameters file
>>> The initial dislocation segment configuration

Your selection of three files must be defined between those two lines.
------------------------------------------------------------------------------
{mat}
{cont}
{seg}
------------------------------------------------------------------------------
"""
        with open(filepath, 'w') as f:
            f.write(input_content)

    def _execute_engine(self, work_dir: Path, use_graphics: bool):
        """Executes the bash script inside dd/exec."""
        exec_dir = work_dir / "dd" / "exec"
        script_name = "./gcomp" if use_graphics else "./comp"
        try:
            subprocess.run(
                ["bash", script_name], 
                cwd=exec_dir, 
                check=True, 
                stdout=subprocess.DEVNULL, # Mutes compiler flooding
                stderr=subprocess.PIPE,
                text=True
            )
        except subprocess.CalledProcessError as e:
            print(f"❌ Execution Failed. Stderr output:\n{e.stderr}")
            raise

    def _export_stat_csv(self, work_dir: Path, timestamp: str, config: SimulationConfig) -> str:
        """Parses the space-separated stat.txt and exports it as a CSV file."""
        stat_file = work_dir / "dd" / "out" / "stat.txt"
        if not stat_file.exists():
            raise FileNotFoundError("Simulation completed but stat.txt was not found.")
            
        csv_filename = f"stat_angle{config.angle_deg}_stress{config.cont_cu.sigma0}_{timestamp}.csv"
        csv_file = work_dir / "dd" / "out" / csv_filename
        
        # Read the raw text file and write out standard comma-separated values
        with open(stat_file, 'r') as infile, open(csv_file, 'w', newline='') as outfile:
            writer = csv.writer(outfile)
            for line in infile:
                # Strip whitespace and split by any number of spaces/tabs
                row = line.strip().split()
                if row:
                    writer.writerow(row)
                    
        return str(csv_file)


# =============================================================================
# 4. EXAMPLE USAGE SCRIPT
# =============================================================================

if __name__ == "__main__":
    # Point this to wherever your local repository clone sits
    REPO_PATH = "../practical_DDD"
    
    wrapper = MicroMegasWrapper(base_repo_path=REPO_PATH)
    
    # We can easily loop over different orientations
    ANGLES = [90] # Pure Edge
    
    for angle in ANGLES:
        # Create a new, unique configuration for this run
        run_config = SimulationConfig(
            target_length_nm=100.0,
            angle_deg=angle, # Unit degrees, can be any float value
            box_dim_x=3000,  # Unit nm, will be adjusted to BVD units in the generator
            box_dim_y=3000,  # Unit nm, will be adjusted to BVD units in the generator
            box_dim_z=3000,  # Unit nm, will be adjusted to BVD units in the generator
        )
        
        # Override specific ContCu variables if needed
        run_config.cont_cu.sigma0 = 53.6
        
        print(f"Starting simulation for angle {angle}°...")
        csv_path = wrapper.run(run_config)
        
        print(f"✅ Saved CSV data to: {csv_path}")
    
    print("\nAll simulations finished directly in the repository.")