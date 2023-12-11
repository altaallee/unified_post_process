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

mpas_regions = {
    "test" : RegionDefault(
        name="west",
        projection=ccrs.LambertConformal(),
        central_lon=-120,
        central_lat=38,
        lon_radius=20,
        scale="regional",
        barb_skip=10,
        station_priority=2,
        hilo_size=50,
    ),
}

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
    2: {
        "AKL": StationsDefault(name="AKL", lon=174.800000, lat=-37.016667), 
        "BHE": StationsDefault(name="BHE", lon=173.866667, lat=-41.516667),
        "CHT": StationsDefault(name="CHT", lon=176.566667, lat=-43.950000),
        "CHC": StationsDefault(name="CHC", lon=172.516667, lat=-43.483333),
        "DUD": StationsDefault(name="DUD", lon=170.183333, lat=-45.916667),
        "GIS": StationsDefault(name="GIS", lon=177.966667, lat=-38.666667),
        "HKK": StationsDefault(name="HKK", lon=170.966667, lat=-42.716667),
        "IVC": StationsDefault(name="IVC", lon=168.316667, lat=-46.416667),
        "KBZ": StationsDefault(name="KBZ", lon=173.683333, lat=-42.416667),
        "KAT": StationsDefault(name="KAT", lon=173.266667, lat=-35.116667),
        "MFN": StationsDefault(name="MFN", lon=167.916667, lat=-44.666667),
        "NPL": StationsDefault(name="NPL", lon=174.166667, lat=-39.016667),
        "OAM": StationsDefault(name="OAM", lon=171.066667, lat=-44.966667),
        "OHA": StationsDefault(name="OHA", lon=175.383333, lat=-40.200000),
        "PPQ": StationsDefault(name="PPQ", lon=174.966667, lat=-40.883333),
        "ROT": StationsDefault(name="ROT", lon=176.316667, lat=-38.116667),
        "TRG": StationsDefault(name="TRG", lon=176.183333, lat=-37.666667),
        "WLG": StationsDefault(name="WLG", lon=174.800000, lat=-41.316667),
        "WSZ": StationsDefault(name="WSZ", lon=171.583333, lat=-41.750000),
    },
    3: {},
}
