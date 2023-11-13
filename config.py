from datetime import timedelta


mpas_freq = timedelta(hours=3)

wrf_freq = {
    1: timedelta(hours=1),
    2: timedelta(hours=1),
    3: timedelta(minutes=30)}

wrf_min_domain = 1
wrf_max_domain = 3

mpas_trim = {
    "west": -150,
    "east": -80,
    "south": 10,
    "north": 50}

wrf_trim = {
    "west": 150,
    "east": 180,
    "south": -50,
    "north": -30}

mpas_slp_smooth_sigma = 2
mpas_height_smooth_sigma = 1

wrf_slp_smooth_sigma = {1: 1, 2: 4, 3: 4}
wrf_height_smooth_sigma = {1: 1, 2: 4, 3: 4}
