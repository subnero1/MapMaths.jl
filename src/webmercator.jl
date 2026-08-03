struct WmX <: ScalarCoordinateSystem end
struct WmY <: ScalarCoordinateSystem end

@simple_vector_coordinate_system(WebMercator, WmX, WmY)
@simple_vector_coordinate_system(WebMercatorAlt, WmX, WmY, GeodeticAltitude)

###

lossless_parent(::WmX) = Lon()
lossless_parent(::WmY) = Lat()
@symmetric lossless_parent(::WmX, ::WmY) = LonLat()
@symmetric lossless_parent(::WmX, ::WmY, ::Alt) = LonLatAlt()

###

cmap_impl(::WmX, (lon,)::Coordinate{Longitude}) = lon / 180
cmap_impl(::Longitude, (wmx,)::Coordinate{WmX}) = 180 * wmx

cmap_impl(::WmY, (lat,)::Coordinate{GeodeticLatitude}) = log(tand((lat + 90) / 2)) / π
cmap_impl(::GeodeticLatitude, (wmy,)::Coordinate{WmY}) = 2 * atand(exp(π * wmy)) - 90
