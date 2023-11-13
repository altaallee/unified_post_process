import cartopy.crs as ccrs


class RegionDefault:
    def __init__(self, name, projection, central_lon, central_lat, lon_radius,
                 scale, barb_skip, station_priority, hilo_size):
        self.name = name
        self.projection = projection
        self.extent = [central_lon - lon_radius, central_lon + lon_radius,
                       central_lat, central_lat] 
        self.scale = scale
        self.barb_skip = barb_skip
        self.station_priority = station_priority
        self.hilo_size = hilo_size

class StationsDefault():
    def __init__(self, name, lon, lat):
        self.name = name
        self.lon = lon
        self.lat = lat

wrf_regions = {
    1: {
        "nzo" : RegionDefault(
            name="nzo",
            projection=ccrs.LambertConformal(
                central_latitude=-40, central_longitude=170,
                standard_parallels=[-45, -35], cutoff=0),
            central_lon=172.5,
            central_lat=-38,
            lon_radius=30,
            scale=1,
            barb_skip=14,
            station_priority=6,
            hilo_size=50,
        ),
    },
    2: {
        "nzi" : RegionDefault(
            name="nzi",
            projection=ccrs.LambertConformal(
                central_latitude=-40, central_longitude=170,
                standard_parallels=[-45, -35], cutoff=0),
            central_lon=172.5,
            central_lat=-40,
            lon_radius=15,
            scale=2,
            barb_skip=20,
            station_priority=7,
            hilo_size=100,
        ),
        #"ni" : RegionDefault(
        #    name="ni",
        #    projection=ccrs.LambertConformal(
        #        central_latitude=-40, central_longitude=170,
        #        standard_parallels=[-45, -35], cutoff=0),
        #    central_lon=175.25,
        #    central_lat=-37.9,
        #    lon_radius=8,
        #    scale=2,
        #    barb_skip=11,
        #    station_priority=7,
        #    hilo_size=100,
        #),
        "si" : RegionDefault(
            name="si",
            projection=ccrs.LambertConformal(
                central_latitude=-40, central_longitude=170,
                standard_parallels=[-45, -35], cutoff=0),
            central_lon=170.25,
            central_lat=-43.5,
            lon_radius=8,
            scale=2,
            barb_skip=11,
            station_priority=7,
            hilo_size=100,
        ),
        #"ak" : RegionDefault(
        #    name="ak",
        #    projection=ccrs.LambertConformal(
        #        central_latitude=-40, central_longitude=170,
        #        standard_parallels=[-45, -35], cutoff=0),
        #    central_lon=174.7,
        #    central_lat=-36.8,
        #    lon_radius=2,
        #    scale=2,
        #    barb_skip=3,
        #    station_priority=7,
        #    hilo_size=100,
        #),
        #"cc" : RegionDefault(
        #    name="cc",
        #    projection=ccrs.LambertConformal(
        #        central_latitude=-40, central_longitude=170,
        #        standard_parallels=[-45, -35], cutoff=0),
        #    central_lon=172.6,
        #    central_lat=-43.5,
        #    lon_radius=2,
        #    scale=2,
        #    barb_skip=3,
        #    station_priority=7,
        #    hilo_size=100,
        #),
    },
    3: {
        "cci" : RegionDefault(
            name="cci",
            projection=ccrs.LambertConformal(
                central_latitude=-40, central_longitude=170,
                standard_parallels=[-45, -35], cutoff=0),
            central_lon=172.6,
            central_lat=-43.5,
            lon_radius=1,
            scale=3,
            barb_skip=4,
            station_priority=7,
            hilo_size=100,
        ),
    },
}

wrf_sounding_stations = {
    1: {},
    2: {},
    3: {
        "chc" : StationsDefault(
            name="CHC",
            lon=171.53200,
            lat=-43.48940,
        ),
    },
}
