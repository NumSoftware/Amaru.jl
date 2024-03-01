
"""
    findpeaks(X::AbstractVector, Y::AbstractVector; minprominence=0.0, maxwidth=Inf)

Find peaks in a vector `Y` with respect to `X` and filter according to `minprominence` and `maxwidth` of peaks.

"""
function findpeaks(X::AbstractVector, Y::AbstractVector; minprominence=0.0, maxwidth=Inf)
    
    idxs = Int[] # peak indices
    if minprominence==0
        mi, ma = extrema(Y)
        minprominence = 0.1*(ma-mi)
    end

    n = length(X)
    for i in 2:n-1
        if Y[i]>Y[i-1] && Y[i]>Y[i+1]

            # find prominence and width
            left = right = i
            while left>1 && Y[left]>=Y[left-1]
                left -= 1
            end
            while right<n && Y[right]>=Y[right+1]
                right += 1
            end
            prominence = Y[i] - max(Y[left], Y[right])
            width = X[right]-X[left]

            if prominence>=minprominence && width<=maxwidth
                push!(idxs, i)
            end
        end
    end

    return idxs
end


"""
    findvalleys(X::AbstractVector, Y::AbstractVector; minprominence=0.0, maxwidth=Inf)

Find valleys in a vector `Y` with respect to `X` and filter according to `minprominence` and `maxwidth` of valleys.

"""
function findvalleys(X::AbstractVector, Y::AbstractVector; minprominence=0.0, maxwidth=Inf)
    return findpeaks(X, -Y; minprominence=minprominence, maxwidth=maxwidth)
end