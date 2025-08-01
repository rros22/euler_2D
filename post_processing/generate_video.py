import pandas as pd
import matplotlib.pyplot as plt
import subprocess
import os

def save_w1_frames(start, step, end, data_dir="results/frames", out_dir="video_frames"):
    os.makedirs(out_dir, exist_ok=True)

    for frame in range(start, end, step):
        frame_str = f"{frame:05d}"

        def load_box(name):
            path = os.path.join(data_dir, f"step_{frame_str}_{name}.csv")
            return pd.read_csv(path)

        try:
            df_left   = load_box("box_left")
            df_center = load_box("box_bottom")
            df_right  = load_box("box_right")
            df_top    = load_box("box_top")
        except FileNotFoundError:
            print(f"Skipping missing frame: {frame_str}")
            continue

        nx_left   = df_left['i'].max()
        ny_center = df_center['j'].max()
        nx_center = df_center['i'].max()
        nx_right  = df_right['i'].max()
        nx_top    = df_top['i'].max()
        ny_top    = df_top['j'].max()

        # Offsets
        x_offset_left   = 0
        x_offset_center = nx_left
        x_offset_right  = nx_left + nx_center
        y_offset_top    = ny_center

        df_left['i_global']   = df_left['i'] + x_offset_left
        df_left['j_global']   = df_left['j']

        df_center['i_global'] = df_center['i'] + x_offset_center
        df_center['j_global'] = df_center['j']

        df_right['i_global']  = df_right['i'] + x_offset_right
        df_right['j_global']  = df_right['j']

        df_top['i_global']    = df_top['i'] + x_offset_center + (nx_center - nx_top) // 2
        df_top['j_global']    = df_top['j'] + y_offset_top

        df = pd.concat([df_left, df_center, df_right, df_top], ignore_index=True)
        df = df.groupby(['j_global', 'i_global'], as_index=False).mean()

        pivot = df.pivot(index='j_global', columns='i_global', values='w_1')

        plt.figure(figsize=(8, 6))
        plt.imshow(pivot, origin='lower', aspect='auto', cmap='gray')
        plt.colorbar(label='w_1')
        plt.title(f"w_1 at Frame {frame}")
        plt.xlabel("i")
        plt.ylabel("j")
        plt.tight_layout()
        plt.savefig(f"{out_dir}/frame_{frame:05d}.png")
        plt.close()

        print(f"Saved frame {frame}")

save_w1_frames(10,10,150000)

# Call ffmpeg to compile the video
ffmpeg_cmd = [
    "ffmpeg",
    "-framerate", "24",
    "-pattern_type", "glob",
    "-i", "video_frames/frame_*.png",
    "-c:v", "libx264",
    "-pix_fmt", "yuv420p",
    "w1_video.mp4"
]

try:
    subprocess.run(ffmpeg_cmd, check=True)
    print("Video successfully created: w1_video.mp4")
except subprocess.CalledProcessError as e:
    print("Error occurred while creating video:", e)

# video command: ffmpeg -framerate 30 -pattern_type glob -i "video_frames/frame_*.png" -c:v libx264 -pix_fmt yuv420p w1_video.mp4
