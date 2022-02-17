# Ripley-s-K-funtion
Estimates de value of the linearized version of the Ripley's K-funtion (Ripley, 1977), as proporsed by Bessag (1977), for a group of points in a cartesian plane. Also, will test the significance of the L-value (if the points pattern is aggregated, random or regular), by Monte Carlo simulations.

BesL{base}                							R Documentation

L-funtion

Description:
Estimates de value of the linearized version of the Ripley's K-funtion (Ripley, 1977), as proporsed by Bessag (1977), for a group of points in a cartesian plane. Also, will test the significance of the L-value (if the points pattern is aggregated, random or regular), by Monte Carlo simulations.
Usage:

     BesL(data,focal, area, t, nsim, type='intra', ylim, xlim)

Arguments:
data		the name of the data.frame to be analysed.
focal	the name of the focal specie which will be compared with the others.
area		the total area of the cartesian plane.
t		the radius, which will be centered in the points of the focal specie; points within this radius will be counted.
nsim		the number of Monte Carlo simulations. it is important for the determination of the p-value.
type		the type of analysis to be made. it can be 'intra', for intraspecific analysis or 'inter' for interspecific analysis.
ylim		the y limit of the cartesian plane.
xlim		the x limit of the cartesian plane.

# Details:
The data to be analysed must be organized as a data.frame. It must have a column named 'x' with the x-coordinates of each point; it must have a column named 'y' with the y-coordinate of each point, and must have a column named 'sp' with the name of the specie corresponding to the x y coordinates.
It's important that all species have the same number of xy coordinates. For this, species with less observed xy coordinates must have it's x and y values filled with NA.
Also, the user must have installed the packages 'reshape2' and 'dplyr'.
Ripley's K-function is based on the variance (second-order analysis) of all point-to-point distances in a two dimensional space, and gives a description of the spatial structure of the points. The intraspecific analysis is defined by the expected number of individuals of the same specie as the focal within a radius t, centred in an arbitrary individual. While the interspecific analysis is defined by the expected number of individuals of a different specie of the focal within a radius t, centred in an arbitrary individual.

For the intraspecific K, will be used the following formula (Haase, 1995):
     K(t)= n^-2*area*sum(wij^-1*It(uij))

Where 'n' is the number of individuals in the analysed plot, 'area' is the total area of the plot, 'uij' is the distance between points i and j, and 'It' have a binary value (0,1), it's value is 1 if 'uij' <= t and it's value is 0 if 'uij' > t; 'wij' is a correcting factor for edge effects, it has the value of 1 if the circle centred at i and passing through the point j (with radius of t) is completly inside the plot area, if part of the circle falls outside the plot area the correctos used is, as proposed by Getis & Franklin (1987):
     wij= 1-cos^-1(ei/uij)/pi
     
Where ei is the distance between the point e and the plot boundary.
For the interspecific K, will be used the following formula (Wiegand & Moloney, 2004):
     K12(t)= (area/(n1*n2))*sum(It(uij)/wij)
     
Where 'n1' is the number of individuals of the focal specie in the analysed plot, and 'n2' is the number of the individuals of other specie in the analysed plot. The 'wij' used is the same as the one proposed by Getis & Franklin (1987) explained above.

To determine the L-value, the formula that will be used is and adaptation of the one proposed by Bessag (1977), extracted from Haase(1995):
     L(t)=sqrt(K(t)/pi)-t
     
If L-value is positive, and above the upper limit of the confidence envelope, the clumped ditribution can be assumed, if L-value is negative the patter can be discribed as dispersed or regular. Else the distribution is alleatory.
For the determination of the confidence levels, it will be realized Monte Carlo simulations, sampling de x and y values, for the creation of alleatory distributions, and comparing the L-values with the L-values simulated, extracting a p-value. If p-value is higher than 0.97 the distribution of the points can be considered aggregated, if p-value is lower than 0.03 the distribution of the points can be considered dispersed or regular.

# Value:
A data.frame is returned, with a column 'Focal' with the name of the focal specie chosed by the user; a column 'Specie'  with the name of the specie compared with the focal; a column 'L-value' with the L-value of the observed pattern and a column 'P-value' with the value of p obtained with the comparison with the L-value of the Monte Carlo simulations.

# Warning:
For the correct estimation, none x and y value must fall up the xlim or the ylim. Also, two individuals must not fall on the same point (same x and y values). Also, t must obey the following condition: t < (area/2)1/2 (Dixon, 2006).
   
# Author(s):
André Mouro D'Angioli, e-mail andremourodangioli@gmail.com.

# References:
Besag, J. (1977). Contribution to the discussion of Dr. Ripley’s paper. Journal of the Royal Statistical Society, Series B 39: 193-195.

Dixon, P. M. (2006). Ripley's K function. Encyclopedia of environmetrics.

Getis, A., & Franklin, J. (1987). Second-order neighborhood analysis of mapped point patterns. Ecology, 68(3), 473-477.

Haase, P. (1995). Spatial pattern analysis in ecology based on Ripley's K‐function: Introduction and methods of edge correction. Journal of Vegetation Science, 6(4), 575-582.

Ripley, B. D. (1977). Modelling spatial patterns. Journal of the Royal Statistical Society. Series B (Methodological), 172-212.

Wiegand, T., & A Moloney, K. (2004). Rings, circles, and null‐models for point pattern analysis in ecology. Oikos, 104(2), 209-229.

# Examples:

####creating a perfect uniform distribution####

     x<- rep(1:10, each=10) 
     y<- rep(1:10, times=10)
     sp<-'a'
     plot(x~y)  

####visualizing the distributions#######

     data=data.frame(x,y,sp)

     BesL(data, 'a', 100, 1, 1000, 'intra', ylim=10.1, xlim=10.1)

####creating a perfect uniform distribution with NAs####

     x1<- c(rep(1:10, each=10),NA,NA)
     y1<- c(rep(1:10, times=10),NA,NA)
     sp<-'a'
     plot(x1~y1)
     data1=data.frame(x1,y1,sp)
     BesL(data1, 'a', 100, 1, 1000, 'intra', ylim=10.1, xlim=10.1) ####notice that the L.values are the same of data and data1

########################################

     x<- c(runif(10, 1,4),runif(10, 6,9))
     y<- c(runif(10,1,4), runif(10, 6,9))
     sp<-rep(c('a','b'),each=10)

     data= data.frame(x,y,sp)

     plot(data[data$sp=='a',]$y~data[data$sp=='a',]$x, xlim=c(0,10.1), ylim=c(0,10.1))
     par(new=T)
     plot(data[data$sp=='b',]$y~data[data$sp=='b',]$x, col='red',xlim=c(0,10.1), ylim=c(0,10.1) )

     BesL(data, 'a', 100, 1, 1000, 'intra', ylim=10, xlim=10) ####notice that the specie 'a' is clumped
     BesL(data, 'a', 100, 1, 1000, 'inter', ylim=10, xlim=10) ####notice that the specie 'a' is dispersed in relation to the specie 'b'

######################################
#####Random Distribution##############

     x<- runif(50,1,10)
     y<- runif(50,1,10)
     sp<- rep(c('a','b','c','d','e'), each=10)
     distribution<-data.frame(x,y,sp)

     BesL(distribution, 'a', 100, 1, 100, 'intra', ylim=10.1, xlim=10.1)
     BesL(distribution, 'a', 100, 1, 100, 'inter', ylim=10.1, xlim=10.1)
