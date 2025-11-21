import numpy as np
import gemmi

def get_pi_info(atoms):
    """计算 Pi 系统的中心、法向量(SVD拟合)和平均 B-factor"""
    positions = np.array([atom.pos.tolist() for atom in atoms])
    
    # 1. 几何中心
    center_array = np.mean(positions, axis=0)
    pi_center = gemmi.Position(*center_array)
    
    # 2. 平均 B-factor
    b_mean = sum(atom.b_iso for atom in atoms) / len(atoms)
    
    # 3. 法向量 (SVD 拟合平面)
    centered_pos = positions - center_array
    u, s, vh = np.linalg.svd(centered_pos)
    # vh 的最后一行对应最小奇异值，即平面的法向量
    normal_vector = vh[2, :]
    
    return pi_center, center_array, normal_vector, b_mean

def calculate_distance(pos1_array, pos2_array):
    return np.linalg.norm(pos1_array - pos2_array)

def calculate_xpcn_angle(x_pos, pi_center, pi_normal):
    """Plevin: X-PiCenter 向量与法向量的夹角"""
    v_x_pi = pi_center - x_pos
    norm_v = np.linalg.norm(v_x_pi)
    norm_n = np.linalg.norm(pi_normal)
    
    if norm_v == 0 or norm_n == 0: return None

    dot_product = np.dot(v_x_pi, pi_normal)
    cos_theta = np.clip(dot_product / (norm_v * norm_n), -1.0, 1.0)
    angle = np.degrees(np.arccos(cos_theta))
    
    if angle > 90:
        angle = 180 - angle
    return angle

def calculate_xh_picenter_angle(pi_center, x_pos, h_pos):
    """Plevin: XH 向量与 X-PiCenter 向量的夹角"""
    v_hx = h_pos - x_pos
    v_hc = h_pos - pi_center 
    
    norm_hx = np.linalg.norm(v_hx)
    norm_hc = np.linalg.norm(v_hc)
    
    if norm_hx == 0 or norm_hc == 0: return None
    
    cos_theta = np.clip(np.dot(v_hx, v_hc) / (norm_hx * norm_hc), -1.0, 1.0)
    return np.degrees(np.arccos(cos_theta))

def calculate_hudson_theta(pi_center, x_pos, h_pos, normal):
    """Hudson: XH 向量与法向量的夹角 (前提: 投影指向环内)"""
    v_x_pi = pi_center - x_pos
    v_xh = h_pos - x_pos
    
    norm_xpi = np.linalg.norm(v_x_pi)
    if norm_xpi == 0: return None

    # 投影检查
    proj_len = np.dot(v_xh, v_x_pi) / norm_xpi
    
    if proj_len > 0:
        norm_n = np.linalg.norm(normal)
        norm_xh = np.linalg.norm(v_xh)
        if norm_n == 0 or norm_xh == 0: return None
        
        cos_angle = np.clip(np.dot(normal, v_xh) / (norm_n * norm_xh), -1.0, 1.0)
        angle = np.degrees(np.arccos(cos_angle))
        
        if angle >= 90:
            angle = 180 - angle
        return angle
    return None

def calculate_projection_dist(normal, pi_center, x_pos):
    """Hudson: 投影距离"""
    numerator = np.dot(normal, pi_center - x_pos)
    denominator = np.dot(normal, normal)
    if denominator == 0: return None
    
    t = numerator / denominator
    projection_point = x_pos + t * normal
    return np.linalg.norm(projection_point - pi_center)