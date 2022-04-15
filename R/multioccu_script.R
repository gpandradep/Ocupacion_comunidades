


# Instalar CamtrapR versión de desarrollo ---------------------------------

library(remotes)
install_github("jniedballa/camtrapR")

# Nota: También requiere tener instalado Rtools (https://cran.r-project.org/bin/windows/Rtools/)

# Librerías ---------------------------------------------------------------

library(camtrapR)
library(tidyverse)
library(nimble)
library(nimbleEcology)
library(tictoc)
library(beepr)


# Datos -------------------------------------------------------------------

# Cargamos la tabla de registros de las especies
registers <-  read.csv("Data/Survey/recordTable_OC.csv")
table(registers$Species)


# Cargamos la tabla de operación de cámaras
CToperation <-  read.csv("Data/Survey/CTtable_OC.csv") 


# Generamos la matríz de operación de las cámaras

camop <- cameraOperation(CTtable= CToperation, # Tabla de operación
                         stationCol= "Station", # Columna que define la estación
                         setupCol= "Setup_date", #Columna fecha de colocación
                         retrievalCol= "Retrieval_date", #Columna fecha de retiro
                         hasProblems= T, # Hubo fallos de cámaras
                         dateFormat= "%Y-%m-%d") # Formato de las fechas

# Generar las historias de detección ---------------------------------------

# Función para generar las historias de detección para todas las especies seleccionadas

DetHist_list <- lapply(unique(registers$Species), FUN = function(x) {
  detectionHistory(
    recordTable         = registers, # Tabla de registros
    camOp                = camop, # Matriz de operación de cámaras
    stationCol           = "Station",
    speciesCol           = "Species",
    recordDateTimeCol    = "DateTimeOriginal",
    recordDateTimeFormat  = "%d/%m/%Y",
    species              = x,     # la función reemplaza x por cada una de las especies
    occasionLength       = 7,  # Colapso de las historias a 7 días
    day1                 = "station", #inicie en la fecha de cada estación
    datesAsOccasionNames = FALSE,
    includeEffort        = TRUE,
    scaleEffort          = TRUE,
    timeZone             = "America/Mexico_City" 
  )}
)

# Se genera una lista con cada historia de detección y el esfuerzo de muestreo, ahora le colocaremos los nombres para saber a cual especie corresponde
names(DetHist_list) <- unique(registers$Species)

# Finalmente creamos una lista nueva donde esten solo las historias de detección
ylist <- lapply(DetHist_list, FUN = function(x) x$detection_history)

# Todas las historias deben tener el mismo número de sitios y de ocasiones de muestreo


# Covariables -------------------------------------------------------------

#Cargamos la base de covariables
covs.data<- read_csv("Data/Covs/selectedcov_nostd180821.csv")

## Seleccionamos las covariables númericas y las estandarizamos
cov.num <- covs.data %>% 
  dplyr::select(where(is.numeric)) %>%
  scale() %>%  # Standardize covariates
  as.data.frame()

# Nota: se debería veríficar correlaciones o multicolinearidad entre las covariables

# Ahora las no numericas
cov.fac <- covs.data %>% 
  dplyr::select(where(is.character))

## Juntamos de nuevo
covs <- data.frame(cov.fac, cov.num)

# Las covariables deben tener el mismo número de sitios (filas) que las historias de detección

identical(nrow(ylist[[1]]), nrow(covs)) 


# Base de datos para los análisis -----------------------------------------

data_list <- list(ylist    = ylist, # Historias de detección
                  siteCovs = covs, # Covariables de sitio
                  obsCovs  = list(effort = DetHist_list[[1]]$effort))  # agregamos el esfuerzo de muestreo como covariable de observación



# Modelo multi-especie en nimble -----------------------------------------

# CamtrapR permite ajustar modelos multi-especie en JAGS y Nimble, nostros vamos a usar Nimble que en teoría puede ser más rápido para correr estos modelos

# Se creará un txt temporal donde estarán las especificaciones del modelo en enfoque Bayesiano
modelfile <- tempfile(fileext = ".txt")

# Generemos el modelo
mod.nimble <- communityModel(data_list, # la lista de datos
                             occuCovs = list(fixed = "Dpop_G"), # La covariables de sitio
                             detCovsObservation = list(fixed = "effort"), #Covariables de observación
                             intercepts = list(det = "ranef", occu = "ranef"),
                             modelFile = modelfile,
                             nimble = TRUE)      # set nimble = TRUE

summary(mod.nimble)


# Corremos el modelo

fit.nimble.comp <- fit(mod.nimble,
                       n.iter = 1000,
                       n.burnin = 500,
                       chains = 3,
                       compile = TRUE, 
);beep(sound = 4)



# Resultados --------------------------------------------------------------

summary(fit.nimble.comp)


plot_effects(mod.nimble,
             fit.nimble.comp,
             submodel = "state")


plot_effects(mod.nimble,
             fit.nimble.comp,
             submodel = "det")


plot_coef(mod.nimble,
          fit.nimble.comp,
          submodel = "state",
          combine = TRUE)


plot_coef(mod.nimble,
          fit.nimble.comp,
          submodel = "det")


plot(fit.nimble.comp)
