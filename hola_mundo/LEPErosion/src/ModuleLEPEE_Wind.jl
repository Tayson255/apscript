
	using Interpolations

	export Cw_TSR, Cw_Hist, pdf_weibull_wind

	winddata = (mu=17.27, k=2.0)
	function pdf_weibull_wind(v_wind,winddata=winddata)
		k2=winddata.k #shape parameter
		λ=winddata.mu/gamma(1+1/k2) #scale parameter
		return (k2/λ)*((v_wind/λ).^(k2-1)).*exp.(-1*(v_wind/λ).^k2)
	end

	"""
		TSR = Tip Speed Ratio
		Rmax = blade tip radius (m)
		cut_in = cut-in wind speed (m/s)
		cut_out = cut-out wind speed (m/s)
		MinRS =  minimum rotor speed (rpm)
		MaxRS =  maximum rotor speed (rpm)
	"""
	tsrdata = (tsr=3.0,rmax=45,cutin=3,cutout=30,minrs=6,maxrs=16)

	"""
	Cw_TSR(m, tsrdata=tsrdata, pdf_wind(v_wind) = pdf_weibull_wind(v_wind,winddata))
	"""
	function Cw_TSR(m, tsrdata=tsrdata, 
		pdf_wind = v_wind -> pdf_weibull_wind(v_wind, winddata),
		m_increment::Real = m_increment, 
		m_max::Real = m_max, 
		m_min::Real = m_min)
		u1 = (tsrdata.minrs*2*pi/60)*tsrdata.rmax/tsrdata.tsr
		u2 = (tsrdata.maxrs*2*pi/60)*tsrdata.rmax/tsrdata.tsr

		vwind_increment=0.01
		vwinds = collect(tsrdata.cutin:vwind_increment:tsrdata.cutout)

		vwinds1 =vwinds[vwinds .<= u1]
		vwinds2 =vwinds[(vwinds .>= u1) .&& ( vwinds .<= u2)]
		vwinds3 =vwinds[vwinds .>= u2]

		datx = collect(m_min:m_increment:m_max)
		daty =similar(datx)

			for iter in eachindex(datx)
				res1 = trapz(vwind_increment,pdf_wind(vwinds1)) *(u1*tsrdata.tsr/tsrdata.rmax)^(datx[iter]+1)
				res2 = trapz(vwind_increment, ((vwinds2 .^ (datx[iter] +1)) .* pdf_wind(vwinds2)))* (tsrdata.tsr/tsrdata.rmax)^(datx[iter]+1)
				res3 = trapz(vwind_increment,pdf_wind(vwinds3)) *(u2*tsrdata.tsr/tsrdata.rmax)^(datx[iter]+1)
				daty[iter]=res1 + res2+ res3
			end
		itp = LinearInterpolation(datx,daty)
		return itp(m)
	end

	"""
		Cw_Hist(m, limit_w, freq; m_increment, m_max, m_min)

	Calculates the histogram-based wind speed moment for a material characterized by the parameter `m`.

	# Arguments
	- `m`: Material parameter used in the moment calculation. The moment is computed by integrating the wind speed raised to the power of `m+1`.
	- `limit_w`: Array of wind speed interval boundaries.
	- `freq`: Frequency data corresponding to each wind speed interval.
	- `m_increment`: Increment used to generate the range of material parameters (`ms`).
	- `m_max`: Maximum value of the material parameter to be considered.
	- `m_min`: Minimum value of the material parameter to be considered.

	# Details
	1. A parameter array `ms` is generated spanning from `m_min` to `m_max` with a step of `m_increment`.
	2. The frequency data is normalized by dividing `freq` by the result of `rect(limit_w, freq)`, forming a probability density function.
	3. The center of each wind speed interval is computed as the average of consecutive elements in `limit_w`.
	4. For each value in `ms`, the moment is computed by integrating `(w_center)^(m_iter+1)` weighted by the normalized frequency using the function `rect`.
	5. Linear interpolation (using `LinearInterpolation`) is applied to obtain the moment corresponding to the specified parameter `m`.

	# Returns
	A Real representing the interpolated wind speed moment for the material characterized by `m`.
	"""
	function Cw_Hist(m::Reals,
		limit_w::Vector{<:Real},
		freq::Vector{<:Real};
		m_increment::Real = m_increment, 
		m_max::Real = m_max, 
		m_min::Real = m_min
		) 
		ms = collect(m_min:m_increment:m_max)
		freq2 = freq/rect(limit_w,freq) 
		w_center = (limit_w[1:end-1] + limit_w[2:end])/2
		daty = [ rect(limit_w, (w_center.^(miter+1)) .* freq2) for miter in ms ]
		itp = LinearInterpolation(ms,daty)
		return itp(m)
	end
