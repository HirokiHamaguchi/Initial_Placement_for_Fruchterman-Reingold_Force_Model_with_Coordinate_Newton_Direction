import cv2
import numpy as np

for idx, it in zip(range(173, 175 + 1), [10, 20, 500]):
    img = cv2.imread(f"Screenshot ({idx}).png")

    # Clip the image
    clip = img[185:1000, 160:800]

    # save the clipped image
    cv2.imwrite(f"{it}.png", clip)
