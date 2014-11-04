function csAverage(y)
    yAvg = zeros(Nz)
    for i = 2:length(r)
        yAvg += (r[i]-r[i-1])*(r[i]*y[i:Nr:end] + r[i-1]*y[(i-1):Nr:end])
    end
    yAvg/Radius^2
end