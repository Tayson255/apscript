# This file is part of the LEPErosion package

	using StatsBase
	using Interpolations


	export trapz;
	export Reals;
	export PDFfromData;
	export diam_increment, diam_max, diam_min;
	export intensity_increment, intensity_max, intensity_min;
	export m_increment, m_max, m_min;

	const Reals = Union{Real, AbstractVector{<:Real}}

	diam_increment::Real = 0.01
	diam_max::Real=6.0
	diam_min::Real=0.01

	intensity_increment::Real = 0.01
	intensity_max::Real=10.0
	intensity_min::Real=0.01

	m_increment::Real = 0.01
	m_max::Real = 10.0
	m_min::Real = 0.0

	"""
	trapz(ix::Real, y::Vector) -> Real

	Calculates the integral of `y` using the trapezoidal rule with a constant spacing `ix`.
	"""
	function trapz(ix::Real, y::Vector{<:Real})
		ix * ((y[1] + y[end]) / 2 + sum(y[2:end-1]))
	end

	"""
	trapz(x::Vector, y::Vector) -> Real

	Calculates the integral of `y` with respect to `x` using the trapezoidal rule.
	"""
	function trapz(x::Vector{<:Real}, y::Vector{<:Real})
		@assert length(x) == length(y) "Vectors x and y must have the same length"
		@views h = diff(x)
		@views res = sum((y[1:end-1] + y[2:end]) .* h) / 2
		return res
	end

	"""
    rect(limits::Vector, values::Vector)

	Calculates the integral of a function using the rectangular rule. This function is intended for integrating histograms.
		limits must be a vector representing the boundaries of the intervals over which the function is integrated.
		values must be a vector of function values defined for each subinterval. Its length should be one less than the length of limits.

	# Arguments
	- `limits`: A vector representing the boundaries of the intervals over which the function is integrated.
	- `values`: A vector of function values defined for each subinterval. Its length should be one less than the length of `limits`.

	# Returns
	The approximated integral as the sum of the products of the function values and the corresponding interval lengths.

	# Example
	```julia
	limits = [0.0, 1.0, 2.0, 3.0]
	values = [2.0, 3.0, 4.0]  # One value per interval
	result = rect(limits, values)  # Approximates the integral
	"""
	function rect(limits::Vector{<:Real}, values::Vector{<:Real})
		@views h = diff(limits)
		@views res = sum(values .* h)
		return res
	end

	"""
    PDFfromData(x::Reals, data::Vector{<:Real}; nbins::Int = 20) -> Any

    Calculates the probability density function (PDF) from the provided data.

    # Arguments
    - x: A Real or a vector of Reals at which the PDF is evaluated.
    - data: A vector of real Reals used to build the histogram.
    - nbins: Real of bins for the histogram. Default is 20.

    # Returns
    Returns the interpolated PDF values corresponding to `x`.

    # Details
    The function constructs a histogram from the data using `fit` from StatsBase. It computes the PDF by dividing the bin counts by the product of the total Real of data points and the bin width. Then, it creates an interpolation function (with constant extrapolation to zero outside the histogram range) that is evaluated at `x`.

    # Example
    ```julia
    data = randn(1000)
    x = -3:0.1:3
    pdf_values = PDFfromData(x, data, nbins=30)
    # You can plot the PDF using, for example, Plots.jl:
    # using Plots
    # plot(x, pdf_values)
    ```
    """
	function PDFfromData(x::Reals, data::Vector{<:Real}; nbins::Int = 20)
		h = fit(Histogram, data, nbins=nbins)

		xvals = h.edges[1]
		yvals = [h.weights/(length(data) * (step(xvals))) ; 0.0]

		itp = extrapolate(scale(interpolate(yvals, BSpline(Constant())), xvals),0)

		return itp(x)
	end
