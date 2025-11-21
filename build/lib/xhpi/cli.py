import argparse
import sys
import re
import multiprocessing
import pandas as pd
import logging
import gemmi
from pathlib import Path
from tqdm import tqdm

from . import prep
from . import core
from . import config

# --- 日志配置 ---
logging.basicConfig(
    level=logging.INFO, 
    format='%(asctime)s - %(message)s',
    datefmt='[%H:%M:%S]'
)
logger = logging.getLogger('xhpi')

# --- 辅助函数 ---

def find_files(inputs):
    """
    查找文件逻辑:
    1. 文件: 直接信任并添加。
    2. 目录: 递归查找 ???? + .cif。
    """
    file_list = set()
    pattern = re.compile(r'^[a-zA-Z0-9]{4}\.cif$', re.IGNORECASE)

    for inp in inputs:
        path = Path(inp)
        if path.is_file():
            file_list.add(path.resolve())
        elif path.is_dir():
            for p in path.rglob("*.cif"):
                if pattern.match(p.name):
                    file_list.add(p.resolve())
    return sorted(list(file_list))

def process_one_file(args_packet):
    """Worker: 读取 -> 加氢 -> 检测 -> 保存"""
    filepath, mon_lib_path, file_type, h_mode, out_dir_root = args_packet
    filepath = Path(filepath)
    pdb_name = filepath.stem 
    
    try:
        # 1. 读取
        try:
            structure = gemmi.read_structure(str(filepath))
        except Exception as e:
            return f"Read Error ({pdb_name}): {e}"

        # 2. 加氢
        structure = prep.add_hydrogens_memory(structure, mon_lib_path, h_change_val=h_mode)
        if not structure:
            return f"AddH Failed ({pdb_name})"

        # 3. 检测
        results = core.detect_interactions_in_structure(structure, pdb_name)

        # 4. 保存
        if results:
            df = pd.DataFrame(results)
            
            # 确定输出目录
            if out_dir_root:
                target_dir = Path(out_dir_root)
                target_dir.mkdir(parents=True, exist_ok=True)
            else:
                target_dir = filepath.parent

            # 确定格式
            if file_type == 'csv':
                out_path = target_dir / f"{pdb_name}_xhpi.csv"
                df.to_csv(out_path, index=False)
            else:
                out_path = target_dir / f"{pdb_name}_xhpi.json"
                df.to_json(out_path, orient='records', indent=2)
            
            return None 
        else:
            return None 
            
    except Exception as e:
        return f"Critical Error ({pdb_name}): {e}"

# --- CLI 主入口 ---

def main():
    # 自定义 Usage
    usage_str = """
      %(prog)s --set-mon-lib <DIR>
      %(prog)s <PATH>... [options]"""

    parser = argparse.ArgumentParser(
        prog="xhpi",
        usage=usage_str,
        description="XH-pi interactions detector in protein structures.",
        formatter_class=argparse.RawTextHelpFormatter,
        add_help=False 
    )

    # --- 参数定义 ---

    parser.add_argument(
        '-h', '--help',
        action='help',
        default=argparse.SUPPRESS,
        help="Show this help message and exit."
    )

    parser.add_argument(
        '--set-mon-lib',
        type=str,
        metavar='DIR',
        help="Set the default Monomer Library path permanently and exit."
    )

    # 核心输入 (隐藏 help 以简化显示，已在 usage 体现)
    parser.add_argument(
        'inputs', 
        nargs='*', 
        metavar='PATH',
        help=argparse.SUPPRESS 
    )

    # 详细的 H-mode help
    h_mode_help = (
        "Hydrogen handling mode for Gemmi (Default: 4).\n"
        "  0 = NoChange\n"
        "  1 = Shift\n"
        "  2 = Remove\n"
        "  3 = ReAdd\n"
        "  4 = ReAddButWater\n"
        "  5 = ReAddKnown"
    )
    parser.add_argument(
        '--h-mode',
        type=int,
        default=config.DEFAULT_H_CHANGE,
        choices=[0, 1, 2, 3, 4, 5],
        metavar='NUM',
        help=h_mode_help
    )

    parser.add_argument(
        '--file-type',
        type=str,
        default='json',
        choices=['json', 'csv', 'JSON', 'CSV'], 
        metavar='STR',
        help="File type of Export results (json/csv) (Default: json)."
    )

    parser.add_argument(
        '--out-dir',
        type=str,
        default=None,
        metavar='DIR',
        help="Optional: Directory to save output files (Default: same as input)."
    )

    parser.add_argument(
        '--jobs', 
        type=int, 
        default=1,
        metavar='NUM',
        help="Number of parallel CPU cores to use (Default: 1)."
    )

    parser.add_argument(
        '--mon-lib', 
        type=str, 
        default=config.DEFAULT_MON_LIB_PATH, 
        metavar='DIR',
        help="Override the monomer library path for this run only."
    )

    args = parser.parse_args()

    # ==========================
    #      逻辑控制
    # ==========================

    # --- 配置模式 ---
    if args.set_mon_lib:
        path = Path(args.set_mon_lib)
        if not path.exists():
            logger.error(f"Path not found: {path}")
            sys.exit(1)
        
        saved_path = config.save_mon_lib_path(path)
        print(f"Configuration saved.\nDefault Monomer Library: {saved_path}")
        sys.exit(0)

    # --- 分析模式 ---
    if not args.inputs:
        parser.print_help()
        sys.exit(0)

    file_type_clean = args.file_type.lower()

    logger.info("--- XH-pi Detector ---")
    files = find_files(args.inputs)
    
    if not files:
        logger.error("No valid CIF files found.")
        sys.exit(1)

    # 信息概览
    lib_status = args.mon_lib if args.mon_lib else "Gemmi Default"
    out_status = args.out_dir if args.out_dir else "Source Directory"

    logger.info(f"Targets    : {len(files)} files")
    logger.info(f"Export     : {file_type_clean.upper()} -> {out_status}")
    logger.info(f"H-Mode     : {args.h_mode}")
    logger.info(f"Mon-Lib    : {lib_status}")
    logger.info(f"Jobs       : {args.jobs}")

    # 任务分发
    tasks = [(f, args.mon_lib, file_type_clean, args.h_mode, args.out_dir) for f in files]
    error_logs = []

    try:
        if args.jobs == 1:
            for task in tqdm(tasks, desc="Processing", unit="file"):
                res = process_one_file(task)
                if res: error_logs.append(res)
        else:
            with multiprocessing.Pool(args.jobs) as pool:
                for res in tqdm(pool.imap_unordered(process_one_file, tasks), total=len(tasks), desc="Processing", unit="file"):
                    if res: error_logs.append(res)

    except KeyboardInterrupt:
        print("\n[Aborted by user]")
        sys.exit(1)

    if error_logs:
        logger.warning(f"Finished with {len(error_logs)} errors.")
        for e in error_logs[:3]:
            logger.warning(f"  - {e}")
    else:
        logger.info("Done. All files processed successfully.")

if __name__ == "__main__":
    multiprocessing.freeze_support()
    main()