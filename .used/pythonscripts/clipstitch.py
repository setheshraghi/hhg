import sys
import glob
import cv2

name = sys.argv[1]
images = glob.glob('images/*.png')
images = sorted(images)
fourcc = cv2.VideoWriter_fourcc(*'mp4v')
video = cv2.VideoWriter(
    filename = '../../clips/' + name + '.mp4', fourcc = fourcc, fps = 30.0,
    frameSize = (430, 430)
)

for image in images:
    frame = cv2.imread(image)
    frame = cv2.resize(frame, dsize=(430, 430))
    video.write(frame)

video.release()
