
	export IfromShapeT, dsdIntensityDiameter, dsdBest, dsdOffShore, dsdMP, Cr_dsd;

	using Distributions

	"""
	beta(d::Reals) -> Reals

	The impingement efficiency is a measure of the fraction of droplets of size d that will impact the rotor blade leading edge material. The impingement efficiency is approximated as follows by Papadakis M, Wong SC, Rachman A, Hung KE, Vu GT, Bidwell CS. Large and small droplet impingement data on airfoils and two simulated ice shapes. Cleveland, OH, USA: Nasa, Glenn Research Center; 2007. NASA/TM-2007-213959: 1-exp(-15*d)

	# Arguments
	- `d::Reals`: Diameter of the droplet in mm.

	# Returns
	- `Reals`: The impingement efficiency varies from 1 to 0 and is applied to the annual droplet counts based on diameter.
	"""
	function beta(d::Reals)
		return 1 .- exp.(-15*d)
	end

	"""
	IfromShapeT(pdf::Function, T::Real, mu::Real) -> Function

	Scales a rain intensity probability density function (pdf) to represent the Real of hours per year for each intensity value.

	# Arguments
	- `pdf::Function`: A function representing the probability density function of rainfall intensities. Note: Although the rain intensity is provided in mm/hr, the pdf returns values in hr/mm.
	- `T::Real`: Total annual rainfall in mm
	- `mu::Real`: Mean intensity value of the pdf.

	# Returns
	An anonymous function that, given a rain intensity `I` (in mm/hr), returns the Real of hours per year associated with that intensity.

	# Explanation
	The function applies a scaling factor (T * mu / (365*24)) to the probability density function `pdf`.
	This transformation ensures that the resulting function yields the Real of hours per year corresponding to each intensity level, distributing the total rainfall across all hours in a year.
	"""
	function IfromShapeT(pdf::Function, T::Real, mu::Real)
		return (I -> pdf(I) * T / mu / (365*24))
	end

	"""
		dsdIntensityDiameter(diameter::Reals, IntensityDistribution::Function, DiameterDistribution::Function; diam_increment::Real, diam_max::Real, diam_min::Real, intensity_increment::Real=intensity_increment, intensity_max::Real=intensity_max, intensity_min::Real=intensity_min) -> Reals

	Computes the integrated intensity-diameter value for a given droplet diameter using the specified intensity
	and diameter distributions. This function creates a range (datx) from diam_min to diam_max using diam_increment,
	computes an integrated value for each diameter via numerical integration (using trapz), and interpolates the result 
	at the specified diameter.

	# Parameters:
	- diameter: The droplet diameter at which the interpolated value is evaluated in mm
	- IntensityDistribution: A function that returns the intensity distribution in hr/(mm/hr).
	- DiameterDistribution: A function that returns the diameter distribution given a diameter, in mm, and intensity in mm/hr.
	- diam_increment: The increment step used for the diameter range in mm.
	- diam_max: The maximum diameter value for the integration in mm.
	- diam_min: The minimum diameter value for the integration in mm.
	- intensity_increment: The increment step used for the intensity integration in mm/hr.
	- intensity_max: The maximum intensity value for the integration in mm/hr.
	- intensity_min: The minimum intensity value for the integration in mm/hr.

	# Returns:
	- Distribution function indicating the Real of droplets per cubic meter per each diameter.

	# Example
	```julia
	IdistributionMiami(I) = IfromShapeT(I -> pdf(LogNormal,I), 1000, mean(LogNormal))
	DDistributionMiami(diameter) = dsdIntensityDiameter(diameter, intensity)
	dsdMiami(diameter) = dsdIntensityDiameter(diameter, Idistribution, DDistributionMiami)
	```
	"""
	function dsdIntensityDiameter(
		diameter::Reals, 
		IntensityDistribution::Function, 
		DiameterDistribution::Function; 
		diam_increment::Real = diam_increment, 
		diam_max::Real = diam_max, 
		diam_min::Real = diam_min, 
		intensity_increment::Real = diam_increment, 
		intensity_max::Real = intensity_max, 
		intensity_min::Real = intensity_min
		)

		intensities = intensity_min:intensity_increment:intensity_max
		datx = diam_min:diam_increment:diam_max
		daty = similar(datx)

		for iter in eachindex(datx)
			daty[iter] = trapz(intensity_increment, (IntensityDistribution(intensities) .* DiameterDistribution(datx[iter], intensities)) )
		end

		itp = LinearInterpolation(datx, daty)
		return itp(diameter)
	end

	"""
	dsdBest(diameter::Reals; T::Real=1000, mu::Real = -0.8, sigma::Real = 1.2) -> Reals
		

	Calculates the Best drop size distribution parametrized by the total annual rain, mean, and standard deviation of the lognormal distribution.

	The default values are given by DNVGL-RP-0573 and are T=1000 mm/year, mu=-0.8, and sigma=1.2.

	# Arguments
	- `diameter::Reals`: Diameter of the droplet in mm.
	- `T::Real`: The total annual rainfall in millimeters (default is 1000).
	- `mu::Real`: The mean of the lognormal distribution (default is -0.8).
	- `sigma::Real`: The standard deviation of the lognormal distribution (default is 1.2).

	# Returns
	- `Real`: A vector representing the Best drop size distribution.

	# Example
	```julia
		dsdBest_Miami(diameter) = dsdBest(diameter,T=1300, mu=-0.1, sigma = 1.6)
	```
	"""
	function dsdBest(
		diameter::Reals; 
		T::Real=1000, 
		mu::Real = -0.8, 
		sigma::Real = 1.2, 
		diam_increment::Real = diam_increment, 
		diam_max::Real = diam_max, 
		diam_min::Real = diam_min, 
		intensity_increment::Real = diam_increment, 
		intensity_max::Real = intensity_max, 
		intensity_min::Real = intensity_min
		)

		pdf_lognormal(I) = pdf(LogNormal(mu,sigma), I) 

		function pdf_weibull_rain(Φ, I)
			k1=2.25 #shape parameter
			λ=1.3*I.^0.232 #scale parameter
			return (k1./λ).*((Φ./λ).^(k1-1)).* exp.(-1 .*(Φ./λ).^k1)
		end

		intensities = intensity_min:intensity_increment:intensity_max
		datx = diam_min:diam_increment:diam_max
		daty = similar(datx)

		for iter in eachindex(datx)
			daty[iter]=trapz(intensity_increment, (intensities.^0.846 .* pdf_lognormal(intensities) .* pdf_weibull_rain(datx[iter],intensities)) )
		end

		daty = daty * 67*6 * exp(-mu-(sigma^2)/2)/(365*24)/pi ./(datx.^3)

		daty = daty * T
		itp = LinearInterpolation(datx,daty)
		return itp(diameter)

	end


	"""
	"""
	function dsdOffShore(
		diameter::Reals; 
		T::Real=1000, 
		mu::Real=-0.8, 
		sigma::Real=1.2, 
		diam_increment::Real = diam_increment, 
		diam_max::Real = diam_max, 
		diam_min::Real = diam_min, 
		intensity_increment::Real = diam_increment, 
		intensity_max::Real = intensity_max, 
		intensity_min::Real = intensity_min
		)

		pdf_lognormal(I) = 3600 * pdf(LogNormal(mu, sigma), I) * T * exp(-mu - (sigma^2) / 2) / (365*24*60*60)
		function pdf_offshore(Φ, I)
			A = 1.0245
			p = 0.1350
			N = 2.8223
			q = -0.0979
			return N .* (I .^ q) .* ((Φ ./ (A .* (I .^ p))) .^ (N .* (I .^ q))) .* exp.(-(Φ ./ (A .* (I .^ p))) .^ (N .* (I .^ q))) ./ Φ
		end

		function W(I)
			67 .* (I .^ 0.846)
		end

		intensities = intensity_min:intensity_increment:intensity_max
		datx = diam_min:diam_increment:diam_max
		daty = similar(datx)

		for iter in eachindex(datx)
			daty[iter] = trapz(intensity_increment, (W(intensities) .* pdf_lognormal(intensities) .* pdf_offshore(datx[iter], intensities)))
		end

		daty = daty .* 6 / pi ./ (datx .^ 3)

		itp = LinearInterpolation(datx, daty)

		return itp(diameter)
	end

	"""
	dsdMP(diameter::Reals; T::Real=1000, mu::Real = -0.8, sigma::Real = 1.2) -> Real
		

	Calculates the Marshal-Palmer drop size distribution parametrized by the total annual rain, mean, and standard deviation of the lognormal distribution.

	The default values are given by DNVGL-RP-0573 and are T=1000 mm/year, mu=-0.8, and sigma=1.2.

	# Arguments
	- `diameter::Reals`: Diameter of the droplet in mm.
	- `T::Real`: The total annual rainfall in millimeters (default is 1000).
	- `mu::Real`: The mean of the lognormal distribution (default is -0.8).
	- `sigma::Real`: The standard deviation of the lognormal distribution (default is 1.2).

	# Returns
	- `Real`: A vector representing the Marshal-Palmer drop size distribution.

	# Example
	```julia
	dsdMP_Miami(diameter) = dsdMP(diameter,T=1300, mu=-0.1, sigma = 1.6)
	```
	"""
	function dsdMP(
		diameter::Reals;
		T::Real=1000, 
		mu::Real = -0.8, 
		sigma::Real = 1.2,
		diam_increment::Real = diam_increment, 
		diam_max::Real = diam_max, 
		diam_min::Real = diam_min, 
		intensity_increment::Real = diam_increment, 
		intensity_max::Real = intensity_max, 
		intensity_min::Real = intensity_min
		)
		pdf_I(I) = 3600 * pdf(LogNormal(mu,sigma), I) * T * exp(-mu-(sigma^2)/2) / (365*24*60*60)
		function pdf_MP(Φ, I)
			C=8000
			A=4.1
			p=-0.21
			a=A * I .^p
			return C .* exp.(-1 .*(a .* Φ))
		end

		intensities = intensity_min:intensity_increment:intensity_max
		datx = diam_min:diam_increment:diam_max
		daty =similar(datx)

		for iter in eachindex(datx)
			daty[iter]=trapz(intensity_increment, (pdf_I(intensities) .* pdf_MP(datx[iter],intensities)) )
		end
			
		itp = LinearInterpolation(datx,daty)
		return itp(diameter)

	end

	"""
	Cr_dsd(dsd::Function, diam_increment::Real=diam_increment, diam_max::Real=diam_max, diam_min::Real=diam_min, beta::Function=beta) -> Real

	Calculates the Cr value for a given drop size distribution (dsd).

	# Arguments
	- `dsd::Function`: Drop size distribution as a function of diameter in drops per cubic meter per mm.
	- `diam_increment::Real`: Increment of diameter (default is `diam_increment`).
	- `diam_max::Real`: Maximum diameter in mm (default is `diam_max`).
	- `diam_min::Real`: Minimum diameter in mm (default is `diam_min`).
	- `beta::Function`: Impingement efficiency as a function of diameter (default is `beta`).

	# Returns
	- `Real`: The Cr value for the given drop size distribution.
	"""
	function Cr_dsd(
		dsd::Function,
		diam_increment::Real = diam_increment,
		diam_max::Real = diam_max,
		diam_min::Real = diam_min,
		beta::Function = beta
		)
		diams = diam_min:diam_increment:diam_max
		trapz(diam_increment, (diams .^ 2) .* beta(diams) .* dsd(diams))
	end
