from .geosearch import GeoSearch

_geo = None

def init(hostname, port, username, password):
    global _geo
    _geo = GeoSearch(hostname, port, username, password)

def geo():
    global _geo
    if _geo is None:
        raise RuntimeError('geo.init() must be called before geo.geo()')
    return _geo