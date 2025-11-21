import gemmi
import numpy as np
from . import config
from . import geometry

def detect_interactions_in_structure(structure, pdb_name):
    """
    接收 gemmi.Structure 对象进行核心分析
    """
    results = []
    try:
        if not structure or len(structure) == 0:
            return []

        resolution = structure.resolution if structure.resolution else 0.0
        model = structure[0] # 使用第一个模型
        
        # 核心优化：建立全局 NeighborSearch Grid
        ns = gemmi.NeighborSearch(model, structure.cell, config.DIST_CUTOFF_6A)
        ns.populate(include_h=True)
        
        for chain in model:
            for residue in chain:
                res_name = residue.name
                
                # Main Pi
                if res_name in config.RING_ATOMS:
                    results.extend(_detect_residue(
                        pdb_name, resolution, model, chain, residue, ns, 
                        config.RING_ATOMS[res_name], mode='main'
                    ))
                
                # Trp A
                if res_name in config.TRP_A_ATOMS:
                    results.extend(_detect_residue(
                        pdb_name, resolution, model, chain, residue, ns, 
                        config.TRP_A_ATOMS[res_name], mode='trpA'
                    ))
                    
    except Exception as e:
        raise RuntimeError(f"Core analysis error: {e}")
        
    return results

def _detect_residue(pdb_name, resolution, model, chain, residue, ns, target_atoms, mode):
    hits = []
    # 1. 获取 Pi 原子
    pi_atoms = [atom for atom in residue if atom.name in target_atoms]
    if len(pi_atoms) != len(target_atoms): return []
    
    # 2. 计算几何
    pi_center, pi_center_arr, pi_normal, pi_b_mean = geometry.get_pi_info(pi_atoms)
    alt_pi = pi_atoms[0].altloc
    
    # 3. 搜索 X 原子
    x_candidates = ns.find_atoms(pi_center, alt=alt_pi, radius=config.DIST_CUTOFF_6A)
    
    for x_mark in x_candidates:
        x_cra = x_mark.to_cra(model) # 传入 model
        x_atom = x_cra.atom
        
        if x_atom.element not in config.TARGET_ELEMENTS_X: continue
        x_pos_arr = np.array(x_atom.pos.tolist())
        
        dist_x_pi = geometry.calculate_distance(x_pos_arr, pi_center_arr)
        if dist_x_pi > config.DIST_CUTOFF_6A: continue
        
        xpcn_angle = geometry.calculate_xpcn_angle(x_pos_arr, pi_center_arr, pi_normal)
        
        # 4. 搜索 H 原子
        h_candidates = ns.find_atoms(x_atom.pos, alt=x_atom.altloc, radius=config.DIST_CUTOFF_H)
        
        for h_mark in h_candidates:
            h_cra = h_mark.to_cra(model)
            h_atom = h_cra.atom
            if h_atom.element != config.ELEMENT_H: continue
            
            h_pos_arr = np.array(h_atom.pos.tolist())
            xh_pi_angle = geometry.calculate_xh_picenter_angle(pi_center_arr, x_pos_arr, h_pos_arr)
            theta = geometry.calculate_hudson_theta(pi_center_arr, x_pos_arr, h_pos_arr, pi_normal)
            
            if xh_pi_angle is None or theta is None or xpcn_angle is None: continue
            
            # --- 判定 ---
            is_plevin = (dist_x_pi < 4.5 and xh_pi_angle > 120.0 and xpcn_angle < 25.0)
            
            proj_threshold = None
            if mode == 'trpA': proj_threshold = 1.6
            elif residue.name == 'HIS': proj_threshold = 1.6
            elif residue.name in ['TRP', 'TYR', 'PHE']: proj_threshold = 2.0
            
            is_hudson = False
            proj_dist = None
            
            if proj_threshold and theta <= 40.0 and dist_x_pi <= 4.5:
                proj_dist = geometry.calculate_projection_dist(pi_normal, pi_center_arr, x_pos_arr)
                if proj_dist is not None and proj_dist <= proj_threshold:
                    is_hudson = True
            
            if not is_plevin and not is_hudson: continue
            
            # --- 记录 ---
            hits.append({
                'pdb': pdb_name,
                'resolution': resolution,
                'pi_res': residue.name,
                'pi_id': residue.seqid.num,
                'pi_chain': chain.name,
                'X_atom': x_atom.name,
                'X_res': x_cra.residue.name,
                'X_id': x_cra.residue.seqid.num,
                'X_chain': x_cra.chain.name,
                'H_atom': h_atom.name,
                'dist_X_Pi': round(dist_x_pi, 3),
                'theta': round(theta, 2),
                'angle_XH_Pi': round(xh_pi_angle, 2),
                'angle_XPCN': round(xpcn_angle, 2),
                'proj_dist': round(proj_dist, 3) if proj_dist else None,
                'is_plevin': 1 if is_plevin else 0,
                'is_hudson': 1 if is_hudson else 0,
                'mode': mode
            })
            
    return hits