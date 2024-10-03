#*************************************************************************
sub interp {
#
#   Purpose:  To linearly interpolate between 2 points to 
#             determine the xval corresponding to the input
#             yval.
#
#   Input:    y1     -- Y-coordinate of first point
#             y2     -- Y-coordinate of second point
#             yval   -- Y-coordinate of point to be determined
#             x1     -- X-coordinate of first point
#             x2     -- X-coordinate of second point
# 
#   Output:  xval    -- X-coordinate calculated by interpolation
#
#*************************************************************************

  my ($y1, $y2, $yval, $x1, $x2) = @_;

#  Determine ratio of yval to total y distance

  $yrange = $y2 - $y1;
  $ydiff = $y2 - $yval;
  $ratio = $ydiff/$yrange;

#  Calculate xval
  $xrange = $x2 - $x1;
  $xdiff = $xrange * $ratio;
  $xval = $x2 - $xdiff;

  return $xval;
}
1;
