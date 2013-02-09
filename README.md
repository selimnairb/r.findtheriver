COPYRIGHT
---------
(C) 2013 by the University of North Carolina at Chapel Hill

This program is free software under the GNU General Public License
(>=v2). Read the file COPYING that comes with GRASS for details.


AUTHOR
------
Brian Miles - brian_miles@unc.edu


DESCRIPTION
-----------
r.findtheriver finds the nearest stream pixel to a coordinate pair
using an upstream accumulating area (UAA) raster map.  This is
necessary because the coordinates for streamflow gages are often not
perfectly registered to the topography represented by a digital
elevation model (DEM) map.  This presents a problem when trying to
derive a watershed contributing area using r.water.outlet; if the
streamflow gage does not fall on the stream as represented in the
DEM, r.water.outlet can fail to derive the watershed area.
 
The basic assumption is that the UAA for "stream" pixels will be much
higher than for adjacent "non-stream" pixels.
r.findtheriver attempts to "snap" the coordinates of the
streamflow gage to the "true" stream location by first identifying
stream pixels within a search window, and then selecting the stream
pixel that is closest (cartesian distance) to the input gage
coordinates.  Stream pixels are identified by searching the UAA
raster window for pixels that exceed a threshold.  This is done by
computing the log10 of the UAA value for the pixel corresponding to
the gage coordinates and subtracting from it the log10 of each pixel
in the window; for a given pixel if this difference is greater than
the threshold, the pixel is deemed to be a stream pixel.

r.findtheriver will automatically compute the window and threshold if
they are not supplied by the user.  The window is determined based on
a THRESHOLD_DISTANCE / cell resolution of the UAA map.  The threshold
is determined by computing the log10 of the minimum and maximum
values of the UAA map, and subtracting 1.</p> <p>The closest stream
pixel is printed to standard output.  If no stream pixels were found
nothing is printed.
