
# Function to get a set of points on sine wave
get_sine_points = function(x, vert, amp, per, horiz){
  y = vert + amp * sin((2*pi*x/per) + horiz)
  return(y)
}

# Function to get peak phase (i.e. time at which sine wave is at max) for given amplitude, period and offsets
get_peak_phase = function(vert, amp, per, horiz){
  x = seq(0, 24, length.out = 288)
  sine_points = get_sine_points(x, vert, amp, per, horiz)
  peak_phase = x[which.max(sine_points)]
  return(peak_phase)
}

