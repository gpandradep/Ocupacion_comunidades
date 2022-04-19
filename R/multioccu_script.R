
##%######################################################%##
#                                                          #
####              Estimación de riqueza de              ####
####          especies # Modelos de ocupación           ####
####            de comunidades # CamtrapR y             ####
####        JAGS/NIMBLE # Gabriel Andrade Ponce         ####
#                                                          #
##%######################################################%##



# 1. Cargar librerias  --------------------------------------------------

# Instalar CamtrapR versión de desarrollo ---------------------------------

remotes::install_github("jniedballa/camtrapR")


# Nota: También requiere tener instalado Rtools (https://cran.r-project.org/bin/windows/Rtools/)

## Librerías 

library(camtrapR) # Datos de cámaras y modelos
library(tidyverse) # Manejo de datos y gráficas
library(rjags) # Lenguaje JAGS
library(nimble) # LEnguaje BUGS
library(nimbleEcology) # Nimble enfocado en jerárquicos
library(SpadeR) # Riqueza Chao2
library(tictoc) # Opcional de tiempo
library(beepr) # Opcional de alertas
library(snowfall)


# 2. Cargar datos ------------------------------------------------------------

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
    occasionLength       = 10, # Colapso de las historias a 10 ías
    day1                 = "station", #inicie en la fecha de cada estación
    datesAsOccasionNames = FALSE,
    includeEffort        = TRUE,
    scaleEffort          = TRUE,
    timeZone             = "America/Mexico_City" 
  )}
)

# Se genera una lista con cada historia de detección y el esfuerzo de muestreo, ahora le colocaremos los nombres para saber a cual especie corresponde
names(DetHist_list) <- unique(registers$Species)

# Finalmente creamos una lista nueva donde estén solo las historias de detección
ylist <- lapply(DetHist_list, FUN = function(x) x$detection_history)

# Todas las historias deben tener el mismo número de sitios y de ocasiones de muestreo


## Covariables 

#Cargamos la base de covariables
covars <- read.csv("Data/Covs/stdcovs_OC.csv")

identical(nrow(ylist[[1]]), nrow(covars)) 


# Base de datos para los análisis -----------------------------------------

data_list <- list(ylist    = ylist, # Historias de detección
                  siteCovs = covars, # Covariables de sitio
                  obsCovs  = list(effort = DetHist_list[[1]]$effort))  # agregamos el esfuerzo de muestreo como covariable de observación



# 3. 1 Modelo multi-especie  -----------------------------------------

# CamtrapR permite ajustar modelos multi-especie en JAGS y Nimble, nosotros vamos a usar JAGS ya que la versión de Nimble aun no permite estimar parámetro N de riqueza de especies

# Se creará un txt temporal donde estarán las especificaciones del modelo en enfoque Bayesiano
modelfile <- (fileext = "modoccu.txt")

# Usaremos la función ` communityModel`

# Generemos el modelo
comu_model <- communityModel(data_list, # la lista de datos
                             occuCovs = list(ranef = "Dcrops"), # La covariables de sitio
                             detCovsObservation = list(ranef = "effort"), #Covariables de observación
                             intercepts = list(det = "ranef", occu = "ranef"),
                             augmentation = c(full = 30),# Número aumentado de especies
                             modelFile = modelfile)

summary(comu_model)


# Corremos el modelo

fit.commu <- fit(comu_model,
                 n.iter = 22000, 
                 n.burnin = 2000,
                 thin = 2,
                 chains = 3,
                 cores = 3,
                 quiet = T
);beep(sound = 4)

# Duración 56 min aprox

save(fit.commu, file="results/DR_result.R") # guardamos los resultados para no correr de nuevo

# Resultados --------------------------------------------------------------

results <- summary(fit.commu)[["statistics"]]


plot_effects(comu_model,
             fit.commu,
             submodel = "state")


plot_effects(comu_model,
             fit.commu,
             submodel = "det")


plot_coef(comu_model,
          fit.commu,
          submodel = "state")


plot_coef(comu_model,
          fit.commu,
          submodel = "det")



# Formatear los datos a un vector de frecuencia
abu_Chao <- yaug %>% 
  select(1:nspec) %>%  # seleccionar especies observadas
  t() %>% # trasponer la tabla
  rowSums(. , na.rm = T) %>% # sumar las filas
  as.data.frame()

# Calcular la riqueza con estimadores no paramétricos
chao_sp <- ChaoSpecies(abu_Chao, datatype = "abundance")

NICHao <- chao_sp$Species_table[5,c(1,3,4)] # Extraer valores de IChao
Nocu <- mod_result$summary[862,c(1,4,7)] # Valores del modelo DR

# Unir en un solo dataframe
Nplotdata <- rbind(IChao=NICHao, DR.mod=Nocu) %>% 
  as.data.frame() %>% 
  rownames_to_column(.)

windowsFonts(TNR = windowsFont("Times New Roman")) # Fuentes

# Gráfico para comparar la riqueza estimada
plotN <- ggplot(Nplotdata, aes(x=rowname, y= Estimate, col=rowname))+
  geom_point(aes(shape=rowname),size=3)+
  geom_errorbar(aes(ymin= Nplotdata$`95%Lower`, ymax= Nplotdata$`95%Upper`), width=.3, size=1)+
  labs(x="Estimador de riqueza",y="Número de especies estimado", title = "Diferencia de los estimadores de riqueza")+
  theme_classic()+
  theme(text=element_text(size = 13, family = "TNR"), plot.title = element_text(hjust= 0.5), legend.position = "none")

plotN