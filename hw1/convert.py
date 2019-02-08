

if (False):
    # linearly interpolate the 256 colors given as 8bits
    # RRRGGGBB to 24bit true RGB color
    colormap = []
    for color in range(256):
        red = (color & 0b11100000) >> 5
        green = (color & 0b00011100) >> 2
        blue = color & 0b00000011

        #now scale the colors to 256 bits each
        red = hex(int(red * 255 / 7))
        green = hex(int(green * 255 / 7))
        blue = hex(int(blue * 255 / 3))

        colormap.append((red, green, blue))

    print(colormap)

if (True):
    # convert matlab parula colormap
    # to hex array we can use in colormap_parula.h
    # read in a python array containing the color values as
    # doubles in the range (0,1)
    # example:
    # matlab_parula = [
    #    (0.2422, 0.1504, 0.6603),
    #    (0.2444, 0.1534, 0.6728)]
    from parula import matlab_parula

    matlab_parula[:] = [(hex(int(i[0] * 255)), hex(int(i[1] * 255)), hex(int(i[2] * 255))) for i in matlab_parula]
    print(matlab_parula)
