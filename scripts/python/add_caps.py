# Modified version by Yanxiang Meng
# 2026-02-19 
# This version is modified to add ACE and NME capping groups to termini of all discontinuous segments,
# not just the first and last residues in a chain


# Original version written by Mohd Ibrahim
# Technical University of Munich
# Email: ibrahim.mohd@tum.de
import numpy as np
import MDAnalysis as mda
import argparse
import warnings
# Suppress specific warnings from MDAnalysis
warnings.filterwarnings("ignore")#, category=UserWarning, module="MDAnalysis.coordinates.PDB")




# ------ Helper functions ------
def create_universe (n_atoms, name, resname, positions, resids, segid):

    u_new = mda.Universe.empty(n_atoms=n_atoms,
                             n_residues=n_atoms,
                             atom_resindex=np.arange (n_atoms),
                             residue_segindex=np.arange (n_atoms),
                             n_segments=n_atoms,
                             trajectory=True) # necessary for adding coordinate


    u_new.add_TopologyAttr('name',   name)
    u_new.add_TopologyAttr('resid', resids)
    u_new.add_TopologyAttr('resname', resname)
    u_new.atoms.positions = positions
    u_new.add_TopologyAttr('segid', n_atoms*[segid])
    u_new.add_TopologyAttr('chainID', n_atoms * [segid])


    return u_new
    
def get_nme_pos (end_residue):

    if "OXT" in end_residue.names:
        index = np.where (end_residue.names == "OXT")[0][0]
        N_position = end_residue.positions [index]
        index_c = np.where (end_residue.names == "C")[0][0]
        carbon_position = end_residue.positions [index_c]
        vector = N_position - carbon_position
        vector /= np.sqrt (sum(vector**2))
        
        C_position = N_position + vector*1.36

        return N_position, C_position

    else:
        # find midpoint of O and CA
        index_o = np.where (end_residue.names == "O")[0][0]
        index_ca = np.where (end_residue.names == "CA")[0][0]

        mid_point = (end_residue.positions [index_o] + end_residue.positions [index_ca] )/2

        # find vector connecting mid_point and C
        index_c = np.where (end_residue.names == "C")[0][0]
        vector  = end_residue.positions [index_c] - mid_point
        vector /= np.sqrt (sum(vector**2))
        N_position = end_residue.positions [index_c] + 1.36* vector
        ##
        C_position = N_position + 1.36*vector
    
    return N_position, C_position 

def get_ace_pos (end_residue):
    
    index_ca = np.where (end_residue.names == "CA")[0][0]
    index_n  = np.where (end_residue.names == "N")[0][0]
    vector   = end_residue.positions [index_n] - end_residue.positions [index_ca] 
    vector  /= np.sqrt (sum(vector**2))

    C1_position = end_residue.positions [index_n] + 1.36*vector

    xa, ya, za =  end_residue.positions [index_ca] 
    xg, yg, zg = C1_position

    # arbritray unit vector
    # create an arbritray orientaiton for the ACE residue
    # does not really matter
    orientation  = np.array([2*np.random.rand () -1, 2*np.random.rand () -1,2*np.random.rand () -1])
    nx, ny, nz =  orientation/np.sqrt (sum(orientation**2))

    ## The carbon and oxygen are placed on the vertices of an equilatrel triangle
    # with another vertex as the Nitrogen atom and the C as the centroid
    # The plane of the triangle is placed in an arbritrary orientation as defined before
    # The orientation does not matter
    ######################################
    x1 = xg - (xa-xg)/2 + np.sqrt (3)*(ny*(za-zg) - nz*(ya-yg))/2
    y1 = yg - (ya-yg)/2 + np.sqrt (3)*(nz*(xa-xg) - nx*(za-zg))/2
    z1 = zg - (za-zg)/2 + np.sqrt (3)*(nx*(ya-yg) - ny*(xa-xg))/2
    
    ## second coordinate
    x2 = xg - (xa-xg)/2 - np.sqrt (3)*(ny*(za-zg) - nz*(ya-yg))/2
    y2 = yg - (ya-yg)/2 - np.sqrt (3)*(nz*(xa-xg) - nx*(za-zg))/2
    z2 = zg - (za-zg)/2 - np.sqrt (3)*(nx*(ya-yg) - ny*(xa-xg))/2

    C2_position = np.array ([x1,y1,z1])
    O_position = np.array ([x2,y2,z2])

    ### rescale distances, the above points may be a bit far apart like 2.1 angstrom but usual bonds are 1.4 or so
    ## Therefore we shrink it
    #  C positinos
    
    vector = C2_position - C1_position
    vector /= np.sqrt (sum (vector**2))
    
    C2_position = C1_position + 1.36*vector

    # O positions
    vector = O_position - C1_position
    vector /= np.sqrt (sum (vector**2))
    
    O_position = C1_position + 1.36*vector
    
    return C1_position, C2_position, O_position


def add_caps(u):
    """Add ACE and NME capping groups to protein termini. For discontinuous segments
    (gaps in residue numbering), caps each continuous run. Returns merged universe."""
    segment_universes = []

    for seg in u.segments:

        chain = u.select_atoms(f"segid {seg.segid}")
        resids = np.sort(np.unique(chain.residues.resids))

        # Split into continuous runs (gap when resid[i] - resid[i-1] > 1)
        runs = []
        run_start = resids[0]
        for i in range(1, len(resids)):
            if resids[i] - resids[i - 1] > 1:
                runs.append((run_start, resids[i - 1]))
                run_start = resids[i]
        runs.append((run_start, resids[-1]))

        run_universes = []
        for first_resid, last_resid in runs:

            # NME before first residue (N-terminal cap)
            first_residue = u.select_atoms(f"segid {seg.segid} and resid {first_resid}")
            ace_positions = get_ace_pos(first_residue)
            ace_names = ["C", "CH3", "O"]
            ace_resid = last_resid + 1  # ACE = last_resid + 1
            ace_universe = create_universe(
                n_atoms=len(ace_positions), name=ace_names,
                resname=len(ace_names) * ["ACE"], positions=ace_positions,
                resids=ace_resid * np.ones(len(ace_names)),
                segid=chain.segids[0]
            )

            # ACE after last residue (C-terminal cap)
            last_residue = u.select_atoms(f"segid {seg.segid} and resid {last_resid}")
            nme_positions = get_nme_pos(last_residue)
            nme_names = ["N", "C"]
            nme_resid = first_resid - 1  # NME = first_resid - 1
            nme_universe = create_universe(
                n_atoms=len(nme_names), name=nme_names,
                resname=len(nme_names) * ["NME"], positions=nme_positions,
                resids=nme_resid * np.ones(len(nme_names)),
                segid=chain.segids[0]
            )

            # Select run atoms (resid first_resid to last_resid)
            run_atoms = u.select_atoms(
                f"segid {seg.segid} and resid {first_resid}:{last_resid}"
            )

            # Remove OXT from last residue if present
            if "OXT" in last_residue.names:
                oxt_idx = np.where(last_residue.names == "OXT")[0][0]
                oxt_atom = last_residue[oxt_idx]
                run_atoms = run_atoms - oxt_atom

            # Merge NME + run + ACE
            u_run = mda.Merge(nme_universe.atoms, run_atoms, ace_universe.atoms)
            run_universes.append(u_run)

        # Merge all runs for this segment
        u_seg = mda.Merge(*(ru.atoms for ru in run_universes))
        segment_universes.append(u_seg)

    ## Join all the universes
    all_uni = mda.Merge(*(seg.atoms for seg in segment_universes))
    return all_uni


# ------ Workflow ------

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Add capping groups ACE and NME to protein termini. Remove the hydrogens from the input pdb file before using this script")
    parser.add_argument('-i', dest='in_file', type=str, default='protein_noh.pdb',help='pdb file')
    parser.add_argument('-o', dest='out_file', type=str, default='protein_noh_cap.pdb',help='output file')

    args      = parser.parse_args()
    in_file   = args.in_file
    out_file  = args.out_file

    u = mda.Universe(in_file)
    all_uni = add_caps(u)
    all_uni.atoms.write(out_file)