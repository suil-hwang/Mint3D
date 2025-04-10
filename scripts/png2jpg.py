import argparse
import os
from PIL import Image

def parse_args():
    parser = argparse.ArgumentParser(description="Convert PNGs in a directory to JPGs")
    parser.add_argument('-i', '--input_dir', type=str, required=True, help="Input directory")
    parser.add_argument('-o', '--output_dir', type=str, help="Output directory (optional)")
    return parser.parse_args()

def convert_png_to_jpg(png_file, jpg_file):
    try:
        with Image.open(png_file) as img:
            if img.mode == 'RGBA' or img.mode == 'LA' or img.mode == 'P': # Check for transparency
                bg = Image.new("RGB", img.size, (255, 255, 255))  # White background
                bg.paste(img, mask=img)  # Use the image itself as a mask
                bg.save(jpg_file, "JPEG")
            elif img.mode == 'RGB':
                img.save(jpg_file, "JPEG") # If no transparency, directly convert
            else:
                print(f"Skipping {png_file}: Unsupported mode {img.mode}") # Handle other modes
                return
            print(f"Converted {png_file} to {jpg_file}")
    except Exception as e:
        print(f"Error converting {png_file}: {e}")

if __name__ == "__main__":
    args = parse_args()
    input_dir = args.input_dir
    output_dir = args.output_dir or input_dir  # Use input_dir if output_dir is not provided

    os.makedirs(output_dir, exist_ok=True)  # Create output dir if it doesn't exist

    for filename in os.listdir(input_dir):
        if filename.lower().endswith(".png"):
            png_file = os.path.join(input_dir, filename)
            jpg_file = os.path.join(output_dir, os.path.splitext(filename)[0] + ".jpg")
            convert_png_to_jpg(png_file, jpg_file)