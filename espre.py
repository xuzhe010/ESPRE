#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
ESPRE: EccDNA Stacking Prediction for Robust Estimation
Version: 1.0.0
Author: User & Gemini
"""

import argparse
import subprocess
import os
import sys
import pandas as pd
import numpy as np
import joblib
import shutil
import warnings
from sklearn.base import BaseEstimator, ClassifierMixin

# 忽略 sklearn 的版本警告
warnings.filterwarnings("ignore")

# ==========================================
# 核心类定义 (必须与训练时一致)
# ==========================================
class StackingWrapper(BaseEstimator, ClassifierMixin):
    """
    ESPRE Stacking Model Wrapper.
    Required for unpickling the trained model.
    """
    def __init__(self, base_models, meta_model, preprocessor):
        self.base_models = base_models
        self.meta_model = meta_model
        self.preprocessor = preprocessor
        self.trained_base_models_ = {}
        self.is_fitted_ = False

    def fit(self, X, y):
        # 仅用于类结构占位，实际预测时使用的是加载的训练好状态
        pass

    def predict_proba(self, X):
        # 实际预测逻辑
        X_proc = self.preprocessor.transform(X)
        level1_preds = pd.DataFrame(index=range(X.shape[0]))
        for name, model in self.trained_base_models_.items():
            level1_preds[name] = model.predict_proba(X_proc)[:, 1]
        return self.meta_model.predict_proba(level1_preds)

    def predict(self, X):
        proba = self.predict_proba(X)[:, 1]
        return (proba >= 0.5).astype(int)

# ==========================================
# 工具函数
# ==========================================
def print_logo():
    logo = r"""
    ███████╗███████╗██████╗ ██████╗ ███████╗
    ██╔════╝██╔════╝██╔══██╗██╔══██╗██╔════╝
    █████╗  ███████╗██████╔╝██████╔╝█████╗  
    ██╔══╝  ╚════██║██╔═══╝ ██╔══██╗██╔══╝  
    ███████╗███████║██║     ██║  ██║███████╗
    ╚══════╝╚══════╝╚═╝     ╚═╝  ╚═╝╚══════╝
    ----------------------------------------
    EccDNA Stacking Prediction for Robust Estimation
    """
    print(logo)

def check_dependencies(model_path, r_script):
    """检查依赖环境"""
    missing = []
    if shutil.which("bedtools") is None: missing.append("bedtools")
    if shutil.which("Rscript") is None: missing.append("Rscript")
    if not os.path.exists(model_path): missing.append(f"Model file ({model_path})")
    if not os.path.exists(r_script): missing.append(f"R script ({r_script})")
    
    if missing:
        sys.exit(f"[Error] Missing dependencies or files: {', '.join(missing)}")

def run_cmd(cmd, step_name):
    """运行Shell命令"""
    print(f"[Running] {step_name}...")
    # print(f"  CMD: {cmd}") # 调试用
    ret = subprocess.call(cmd, shell=True)
    if ret != 0:
        sys.exit(f"[Error] Failed at step: {step_name}")

# ==========================================
# 主流程
# ==========================================
def main():
    print_logo()
    
    # 路径配置
    BASE_DIR = os.path.dirname(os.path.abspath(__file__))
    DB_DIR = os.path.join(BASE_DIR, "database")
    MODEL_DIR = os.path.join(BASE_DIR, "models")
    SCRIPT_DIR = os.path.join(BASE_DIR, "scripts")
    
    # 参数解析
    parser = argparse.ArgumentParser(description="Predict SCZ risk using eccDNA features.")
    parser.add_argument("-i", "--input", required=True, help="Input eccDNA BED file (chr, start, end)")
    parser.add_argument("-g", "--genome", required=True, help="Path to hg38 reference genome (.fa)")
    parser.add_argument("-o", "--output", required=True, help="Output result CSV file")
    parser.add_argument("--tmp_dir", default=None, help="Temporary directory (default: ./espre_tmp_TIMESTAMP)")
    
    args = parser.parse_args()

    # 定义关键文件路径
    WINDOW_BED = os.path.join(DB_DIR, "hg38_1Mb_windows.bed")
    MODEL_FILE = os.path.join(MODEL_DIR, "scz_stacking_model.joblib")
    COLS_FILE = os.path.join(MODEL_DIR, "feature_columns.joblib")
    R_SCRIPT = os.path.join(SCRIPT_DIR, "feature_extract.R")

    # 检查依赖
    check_dependencies(MODEL_FILE, R_SCRIPT)
    
    # 创建临时目录
    if args.tmp_dir is None:
        import time
        args.tmp_dir = f"./espre_tmp_{int(time.time())}"
    if not os.path.exists(args.tmp_dir):
        os.makedirs(args.tmp_dir)

    try:
        # Step 1: Calculate GC Content
        ecc_gc_file = os.path.join(args.tmp_dir, "eccDNA_with_GC.bed")
        cmd_gc = (
            f"bedtools nuc -fi {args.genome} -bed {args.input} "
            f"| sed '1d' | awk -F '\\t' '{{print $1\"\\t\"$2\"\\t\"$3\"\\t\"$6}}' "
            f"> {ecc_gc_file}"
        )
        run_cmd(cmd_gc, "Calculating eccDNA GC content")

        # Step 2: Map to 1Mb Windows
        mapped_file = os.path.join(args.tmp_dir, "mapped_fragments.bed")
        cmd_map = (
            f"bedtools intersect -a {WINDOW_BED} -b {ecc_gc_file} -wa -wb "
            f"| awk '{{print $1\"\\t\"$2\"\\t\"$3\"\\t\"$4\"\\t\"$5\"\\t\"$6\"\\t\"$7}}' "
            f"> {mapped_file}"
        )
        run_cmd(cmd_map, "Mapping to genomic windows")

        # Step 3: Feature Extraction (R)
        feat_csv = os.path.join(args.tmp_dir, "features.csv")
        cmd_r = (
            f"Rscript {R_SCRIPT} "
            f"--input_bed {mapped_file} "
            f"--gc_file {ecc_gc_file} "
            f"--window {WINDOW_BED} "
            f"--output {feat_csv} "
            f"--min_length 100 --max_length 1000"
        )
        run_cmd(cmd_r, "Extracting biological features (Ratio & GC)")

        # Step 4: Prediction
        print("[Predicting] Loading model and computing risk...")
        
        # 加载特征
        if not os.path.exists(feat_csv):
            sys.exit("[Error] Feature file not generated. Check input data format.")
            
        df_feat = pd.read_csv(feat_csv)
        if df_feat.empty:
            sys.exit("[Error] No valid features extracted. Ensure eccDNA mapping to 1Mb windows.")

        # 转换为模型输入格式
        df_ratio = df_feat.set_index('bin_id')[['ratio.corrected']].T
        df_ratio.columns = [f"ratio.corrected_{c}" for c in df_ratio.columns]
        
        df_gc = df_feat.set_index('bin_id')[['GC']].T
        df_gc.columns = [f"GC_{c}" for c in df_gc.columns]
        
        X_input = pd.concat([df_ratio, df_gc], axis=1)

        # 对齐特征列
        train_cols = joblib.load(COLS_FILE)
        X_final = pd.DataFrame(0, index=[0], columns=train_cols)
        common_cols = list(set(X_input.columns) & set(train_cols))
        X_final[common_cols] = X_input[common_cols]

        # 预测
        model = joblib.load(MODEL_FILE)
        prob = model.predict_proba(X_final)[:, 1][0]
        
        # 设定阈值 (建议使用您验证集确定的最佳 cutoff，这里默认 0.5)
        cutoff = 0.5 
        pred_label = "SCZ" if prob >= cutoff else "Normal"
        
        # Step 5: Output
        print("\n" + "="*40)
        print(f"  Sample: {os.path.basename(args.input)}")
        print(f"  Prediction: \033[1;31m{pred_label}\033[0m" if pred_label == "SCZ" else f"  Prediction: \033[1;32m{pred_label}\033[0m")
        print(f"  SCZ Probability: {prob:.4f}")
        print("="*40 + "\n")

        # 保存结果
        res_df = pd.DataFrame({
            'Sample_ID': [os.path.basename(args.input)], 
            'Prediction': [pred_label], 
            'Probability': [prob],
            'Note': ['Flagged High Risk' if prob > 0.8 else '']
        })
        res_df.to_csv(args.output, index=False)
        print(f"[Done] Result saved to: {args.output}")

    except Exception as e:
        import traceback
        traceback.print_exc()
        sys.exit("[Critical Error] Analysis failed.")
    
    finally:
        # 清理
        if os.path.exists(args.tmp_dir):
            shutil.rmtree(args.tmp_dir)

if __name__ == "__main__":
    main()
