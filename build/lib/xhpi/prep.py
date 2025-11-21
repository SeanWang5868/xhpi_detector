import gemmi
import os
import logging
from . import config

logger = logging.getLogger("xhpi.prep")

def add_hydrogens_memory(structure, mon_lib_path=None, h_change_val=config.DEFAULT_H_CHANGE):
    """
    使用 Gemmi Topology 在内存中直接加氢。
    :param h_change_val: 0-5 (对应 Gemmi 的操作模式)
    """
    try:
        # 1. 动态统计用到的残基 (加速读取)
        all_codes = set()
        for model in structure:
            for chain in model:
                for residue in chain:
                    if residue.name.strip():
                        all_codes.add(residue.name)
        
        codes_to_load = list(all_codes)

        # 2. 加载单体库
        monlib = gemmi.MonLib()
        
        if mon_lib_path and os.path.exists(mon_lib_path):
            monlib.read_monomer_lib(mon_lib_path, codes_to_load)
        else:
            # 使用 Gemmi 内置库
            pass 

        # 3. 执行加氢 (prepare_topology)
        gemmi.prepare_topology(
            structure, 
            monlib, 
            model_index=0, 
            h_change=h_change_val, 
            reorder=False, 
            ignore_unknown_links=True 
        )
        
        return structure

    except Exception as e:
        logger.error(f"Topology preparation failed: {e}")
        return None