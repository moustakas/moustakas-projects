#!/usr/bin/env python

"""
Takes a file name and a legend as input.
"""

from PIL import Image, ImageDraw, ImageFont
import sys

pngfile = sys.argv[1]
galaxy = sys.argv[2]

im = Image.open(pngfile)

draw = ImageDraw.Draw(im)
font = ImageFont.load_default().font
draw.text((15, 15),galaxy,(255,255,255),font=font)

im.save('junk.png')





