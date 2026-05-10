import cv2

for idx, it in zip(range(173, 175 + 1), [10, 20, 500]):
    img = cv2.imread(f"Screenshot ({idx}).png")
    if img is None:
        raise FileNotFoundError(f"Screenshot ({idx}).png not found or unreadable")

    # Clip the image
    clip = img[185:1000, 160:800]

    # save the clipped image
    cv2.imwrite(f"{it}.png", clip)
