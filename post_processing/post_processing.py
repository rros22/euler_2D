import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import os

def visualize_frame(frame: int, data_dir="results/frames"):
    frame_str = f"{frame:05d}"  # e.g., 10 -> '00010'

    def load_box(name):
        path = os.path.join(data_dir, f"step_{frame_str}_{name}.csv")
        return pd.read_csv(path)

    # === Load CSVs for specified frame ===
    df_left   = load_box("box_left")
    df_center = load_box("box_bottom")
    df_right  = load_box("box_right")
    df_top    = load_box("box_top")

    # === Determine dimensions ===
    nx_left   = df_left['i'].max()
    ny_left   = df_left['j'].max()
    nx_center = df_center['i'].max()
    ny_center = df_center['j'].max()
    nx_right  = df_right['i'].max()
    ny_right  = df_right['j'].max()
    nx_top    = df_top['i'].max()
    ny_top    = df_top['j'].max()

    # === Global grid size ===
    nx_total = nx_left + nx_center + nx_right
    ny_total = ny_left + ny_center + ny_top

    # === Offsets ===
    x_offset_left   = 0
    x_offset_center = nx_left
    x_offset_right  = nx_left + nx_center

    y_offset_center = 0
    y_offset_top    = ny_center

    # === Add global indices ===
    df_left['i_global']   = df_left['i'] + x_offset_left
    df_left['j_global']   = df_left['j']

    df_center['i_global'] = df_center['i'] + x_offset_center
    df_center['j_global'] = df_center['j']

    df_right['i_global']  = df_right['i'] + x_offset_right
    df_right['j_global']  = df_right['j']

    df_top['i_global'] = df_top['i'] + x_offset_center + (nx_center - nx_top) // 2
    df_top['j_global'] = df_top['j'] + y_offset_top

    # === Combine all data ===
    df = pd.concat([df_left, df_center, df_right, df_top], ignore_index=True)
    df = df.groupby(['j_global', 'i_global'], as_index=False).mean()

    # === Variables to visualize ===
    variables = ['w_1', 'w_2', 'w_3', 'w_4', 'u_1', 'u_2', 'u_3', 'u_4']

    # === Plot 2D fields ===
    fig, axes = plt.subplots(2, 4, figsize=(22, 10))

    for idx, var in enumerate(variables):
        ax = axes[idx // 4, idx % 4]
        pivot = df.pivot(index='j_global', columns='i_global', values=var)
        im = ax.imshow(pivot, origin='lower', aspect='auto', cmap='gray')
        ax.set_title(var)
        ax.set_xlabel('i')
        ax.set_ylabel('j')
        fig.colorbar(im, ax=ax, shrink=0.8)

    plt.suptitle(f"Frame {frame}", fontsize=16)
    plt.tight_layout()
    plt.show()

visualize_frame(1200)  # will load step_00010_box_*.csv files from "results/"
