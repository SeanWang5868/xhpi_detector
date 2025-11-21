import gemmi
import os
import json
from pathlib import Path

# --- 配置文件管理 ---
CONFIG_FILE = Path.home() / ".xhpi_config.json"

def load_saved_mon_lib():
    """加载保存的 Monomer Library 路径"""
    if CONFIG_FILE.exists():
        try:
            with open(CONFIG_FILE, 'r') as f:
                data = json.load(f)
                return data.get("monomer_library_path", None)
        except Exception:
            return None
    return None

def save_mon_lib_path(path):
    """保存 Monomer Library 路径到本地配置"""
    data = {}
    if CONFIG_FILE.exists():
        try:
            with open(CONFIG_FILE, 'r') as f:
                data = json.load(f)
        except Exception:
            data = {}
    
    abs_path = str(Path(path).resolve())
    data["monomer_library_path"] = abs_path
    
    with open(CONFIG_FILE, 'w') as f:
        json.dump(data, f, indent=4)
    
    return abs_path

# --- 路径优先级 ---
_env_path = os.environ.get("GEMMI_MON_LIB_PATH", None)
_saved_path = load_saved_mon_lib()

DEFAULT_MON_LIB_PATH = _env_path if _env_path else _saved_path

# --- 默认参数 ---
# Gemmi H-Change: 0=NoChange, 1=Shift, 2=Remove, 3=ReAdd, 4=ReAddButWater, 5=ReAddKnown
DEFAULT_H_CHANGE = 4 

# --- 核心原子定义 ---
RING_ATOMS = {
    'TRP': {'CD2', 'CE2', 'CE3', 'CZ2', 'CZ3', 'CH2'},
    'TYR': {'CD1', 'CD2', 'CE1', 'CE2', 'CZ', 'CG'},
    'PHE': {'CD1', 'CD2', 'CE1', 'CE2', 'CZ', 'CG'},
    'HIS': {'CE1', 'ND1', 'NE2', 'CG', 'CD2'}
}

TRP_A_ATOMS = {
    'TRP': {'CD1', 'CD2', 'NE1', 'CG', 'CE2'}
}

TARGET_ELEMENTS_X = {gemmi.Element('C'), gemmi.Element('N'), gemmi.Element('O'), gemmi.Element('S')}
ELEMENT_H = gemmi.Element('H')

# 距离阈值
DIST_CUTOFF_6A = 6.0
DIST_CUTOFF_H = 1.3