function csAverage(y)
    # Calculate a cross sectional average using the trapezoidal method
    yAvg = zeros(Nz) # Initialize average
    for i = 2:Nr                                # Loop through the radial points
        # Add the trapezoids
        yAvg += (r[i]-r[i-1])*(r[i]*y[i:Nr:end] + r[i-1]*y[(i-1):Nr:end])
    end
    # Divide by the area (pi disappears)
    yAvg/Radius^2
end