# filename: preprocess_rgb_to_24bit_hex.py (ฉบับแก้ไข)
from PIL import Image
import sys

def preprocess_bmp_rgb_to_24bit_hex(input_bmp_path, output_hex_path, target_width, target_height):
    try:
        img = Image.open(input_bmp_path)
    except FileNotFoundError:
        print(f"Error: Input BMP file '{input_bmp_path}' not found.")
        sys.exit(1)

    # --- START OF CORRECTION ---
    # 1. แปลงภาพเป็นโหมด 'L' (Grayscale 8-bit) เพื่อให้ได้ค่าความสว่าง (luminance) ที่แท้จริง
    #    ไม่ว่าภาพต้นทางจะมาเป็น RGB, RGBA หรือโหมดอะไรก็ตาม
    print(f"Converting image to 8-bit Grayscale ('L' mode) to ensure true grayscale values.")
    img = img.convert('L')
    # --- END OF CORRECTION ---

    # Resize if necessary (ทำหลัง convert เพื่อความแม่นยำ)
    if img.width != target_width or img.height != target_height:
        print(f"Warning: Resizing image from {img.size} to ({target_width}, {target_height})")
        img = img.resize((target_width, target_height), Image.Resampling.LANCZOS)

    # ตอนนี้ pixels จะเป็น list ของค่าความสว่าง 8-bit ค่าเดียว (เช่น 248, 198, ...)
    pixels = list(img.getdata())

    with open(output_hex_path, 'w') as f:
        # gray_value คือค่า 8-bit เช่น 248
        for gray_value in pixels:
            # นำค่า gray_value มาใช้ซ้ำกัน 3 ครั้งเพื่อให้ได้ Grayscale 24-bit (R=G=B)
            r, g, b = gray_value, gray_value, gray_value
            f.write(f"{r:02x}{g:02x}{b:02x}\n")

    print(f"Successfully preprocessed '{input_bmp_path}' ({img.width}x{img.height}) to '{output_hex_path}' with 24-bit TRUE GRAYSCALE hex data.")

if __name__ == "__main__":
    if len(sys.argv) != 5:
        print("Usage: python preprocess_rgb_to_24bit_hex.py <input_bmp_file> <output_hex_file> <width> <height>")
        sys.exit(1)

    input_file = sys.argv[1]
    output_file = sys.argv[2]
    width = int(sys.argv[3])
    height = int(sys.argv[4])
    preprocess_bmp_rgb_to_24bit_hex(input_file, output_file, width, height)