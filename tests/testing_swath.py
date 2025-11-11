## Obtener perfil topográfico de una línea
from landspy import DEM, SwathProfile, get_line_strings
import numpy as np

#########
input_shp = "C:/Users/geolo/PycharmProjects/landspy/tests/data/in/tunez_swath_profiles.shp"
input_dem = "C:/Users/geolo/PycharmProjects/landspy/tests/data/in/tunez.tif"
input_width = 1500
n_lines = 50
step_size = 30
names_field = "nombre"
output_profiles = "C:/Users/geolo/PycharmProjects/landspy/tests/data/out/perfiles.npy"
########

swath_profiles = []
dem = DEM(input_dem)
lines, names = get_line_strings(input_shp, names_field)
for n, line in enumerate(lines):
    swath_profiles.append(SwathProfile(line, dem, input_width, n_lines, step_size, names[n]))

np.save(output_profiles, np.array(swath_profiles))

import matplotlib.pyplot as plt

figure, ax = plt.subplots()
my_swaths = np.load(output_profiles, allow_pickle=True)
my_swaths[2].draw_swath(ax, True, True, True, True, True, True, 'RAW', True)

plt.show()