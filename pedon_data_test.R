library(soilDB)
library(aqp)
library(sp)

spc_access
str(spc_access)

#subset site data for ID, location, epipedon (change epipedon type to fit your data "umbric, mollic, ochric" if not using MT663 pedons) and create data frame
my_sub <- site(spc_access)[ ,c('pedlabsamnum','latitude_decimal_degrees','longitude_decimal_degrees')]
str(my_sub) 
names(my_sub)
my_sub #note if the x and y columns below are NA, double check to see that the NASIS site table is loaded and that there are coordinates populated then rerun the fetchNASIS() code.