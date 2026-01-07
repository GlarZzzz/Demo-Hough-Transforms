# filename: create_rgb_bmp.py
from PIL import Image
import sys

def create_rgb_bmp_from_hex(input_hex_path, output_bmp_path, width, height):
    """
    Creates an RGB BMP image from a hex file where each line contains a 24-bit
    RGB value (e.g., 'RRGGBB').
    """
    rgb_pixels = []
    problematic_lines = 0

    try:
        with open(input_hex_path, 'r') as f:
            line_number = 0
            for line_raw in f:
                line_number += 1
                line = line_raw.strip()
                if not line:  # Skip empty lines
                    continue
                try:
                    # Convert hex string (e.g., 'FF8000') to an integer
                    hex_val = int(line, 16)
                    
                    # Extract R, G, B components using bitwise operations
                    # R is in bits 23-16, G in bits 15-8, B in bits 7-0
                    r = (hex_val >> 16) & 0xFF
                    g = (hex_val >> 8) & 0xFF
                    b = hex_val & 0xFF
                    
                    # Append the (R, G, B) tuple to the list
                    rgb_pixels.append((r, g, b))
                    
                except ValueError:
                    print(f"Python Value Error: Failed to parse line {line_number}: '{line}' - Appending black (0,0,0) instead.")
                    # Append a default black pixel if parsing fails
                    rgb_pixels.append((0, 0, 0))
                    problematic_lines += 1
    except FileNotFoundError:
        print(f"Error: Input hex file '{input_hex_path}' not found.")
        sys.exit(1)
    except Exception as e: # Catch other potential IO errors
        print(f"Error reading file '{input_hex_path}': {e}")
        sys.exit(1)

    if problematic_lines > 0:
        print(f"Warning: Encountered {problematic_lines} problematic line(s) that were replaced with black pixels.")

    expected_pixels = width * height
    actual_pixels = len(rgb_pixels)

    if actual_pixels != expected_pixels:
        print(f"Error/Warning: Expected {expected_pixels} pixels, but found {actual_pixels} in hex file.")
        if actual_pixels < expected_pixels:
            padding_count = expected_pixels - actual_pixels
            print(f"Padding with {padding_count} black pixel(s) to meet target size.")
            # Pad with black pixels
            rgb_pixels.extend([(0, 0, 0)] * padding_count)
        elif actual_pixels > expected_pixels:
            print(f"Truncating to {expected_pixels} pixels.")
            rgb_pixels = rgb_pixels[:expected_pixels]

    # Create a new RGB image
    try:
        # Use 'RGB' mode for 24-bit color images
        img = Image.new('RGB', (width, height))
        # putdata expects a sequence of (R, G, B) tuples for RGB mode
        img.putdata(rgb_pixels)
    except Exception as e:
        print(f"Error creating image with PIL: {e}")
        sys.exit(1)

    try:
        img.save(output_bmp_path)
        print(f"Successfully created '{output_bmp_path}' ({width}x{height}) from '{input_hex_path}'")
    except Exception as e:
        print(f"Error saving BMP file '{output_bmp_path}': {e}")
        sys.exit(1)


if __name__ == "__main__":
    if len(sys.argv) != 5:
        print("Usage: python create_rgb_bmp.py <input_hex_file> <output_bmp_file> <width> <height>")
        sys.exit(1)

    input_file_arg = sys.argv[1]
    output_file_arg = sys.argv[2]
    
    try:
        width_arg = int(sys.argv[3])
        height_arg = int(sys.argv[4])
    except ValueError:
        print("Error: Width and height must be integer values.")
        sys.exit(1)

    create_rgb_bmp_from_hex(input_file_arg, output_file_arg, width_arg, height_arg)
