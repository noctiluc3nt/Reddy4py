#### Constants #####
#' sigma
#'
#' Stefan-Boltzmann constant [W/m^2/K^4]
#' @keywords internal
def sigma():
	return(5.67*10^(-8))


#' cp
#'
#' heat capacity, constant pressure [J/(kg*K)]
#' @keywords internal
def cp():
	return(1005.7)


#' g
#'
#' gravitational accelaration [m/s^2]
#' @keywords internal
def g():
	return(9.8062)


#' clight
#'
#' speed of light [m/s]
#' @keywords internal
def clight():
	return(299792458)


#' csound
#'
#' speed of sound [m/s]
#' @keywords internal
def csound():
	return(343)


#' von Karman constant
#'
#' @keywords internal
def karman():
	return(0.4)


#' R
#'
#' ideal gas constant [J/(mol*K)]
#' @keywords internal
def Runiversal():
	return(8.31451)


#' Rd
#'
#' gas constant for dry air [J/(kg*K)]
#' @keywords internal
def Rd():
	return(287.05)


#' Rv
#'
#' gas constant for water vapor [J/(kg*K)]
#' @keywords internal
def Rv():
	return(461.52)


#' gamma (ratio cp/cv)
#'
#' ratio of specific heat at constant pressure to that at constant volume (i.e. cp/cv = 1004 / 717 = 1.4) 
#' @keywords internal
def cpcv():
	return(1.4)


#' rhoAir
#'
#' density of air [kg/m^3]
#' @keywords internal
def rhoAir():
	return(1.225) #at surface with p0 = 1013.25 hPa and T0 = 288 K 


#' Lv
#'
#' latent heat of vaporization [J/kg]
#' @param temp temperature [K] (optional)
#' @keywords internal
def Lv(temp=None):
	if temp is None:
		return(2264.705*1000)
	else:
		return((3147.5-2.37*temp)*1000) #Lv temperature dependence
		


#' M_H2O
#'
#' Molar mass of water
#' @keywords internal
def M_H2O():
	return(0.01802)


#' M_CO2
#'
#' Molar mass of carbon dioxide
#' @keywords internal
def M_CO2():
	return(0.044)


#' M_CH4
#'
#' Molar mass of methane
#' @keywords internal
def M_CH4():
	return(0.01604)


#' Charnock constant (alpha)
#'
#' Charnock constant (alpha) used to convert friction velocity to surface roughness length
#' @keywords internal
def alpha():
	return(0.016)


#' Earth's radius
#'
#' Earth's radius
#' @keywords internal
def R_earth():
	return(6378137)
