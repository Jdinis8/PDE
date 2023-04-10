from PIL import Image
from moviepy.editor import *
import glob
import os

# Create the frames
frames = []
imgs = glob.glob(
    os.getcwd() + "/graphics/*.png")


#sorting :=)
imgs.sort(key=lambda entry: float(((entry.rsplit('.',1)[0]).rsplit('/',1)[1]).rsplit('_', 1)[1]) )

#making our frames for GIF
for i in imgs:
    new_frame = Image.open(i)
    frames.append(new_frame)
 

##GIF
#frames[0].save(os.getcwd() + '/graphics/simulation.gif', format='GIF',
#               append_images=frames[1:],
#               save_all=True,
#               duration=300, loop=0)

##video
clip = ImageSequenceClip(imgs, fps = 20) 
clip.write_videofile(
    os.getcwd() + "/graphics/simulation.mp4", fps=24)
