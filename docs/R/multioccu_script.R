
##%######################################################%##
#                                                          #
####              Estimación de riqueza de              ####
####          especies # Modelos de ocupación           ####
####            de comunidades # CamtrapR y             ####
####        JAGS/NIMBLE # Gabriel Andrade Ponce         ####
#                                                          #
##%######################################################%##



# 1. Cargar librerías  --------------------------------------------------

# Instalar CamtrapR versión de desarrollo ---------------------------------

remotes::install_github("jniedballa/camtrapR")


# Nota: También requiere tener instalado Rtools (https://cran.r-project.org/bin/windows/Rtools/)

## Librerías 

library(camtrapR) # Datos de cámaras y modelos
library(tidyverse) # Manejo de datos y gráficas
library(rjags) # Lenguaje JAGS
library(nimble) # LEnguaje BUGS
library(nimbleEcology) # Nimble enfocado en jerárquicos
library(bayesplot) # gráficos estimaciones bayesianas
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
                             detCovsObservation = list(fixed = "effort"), #Covariables de observación
                             intercepts = list(det = "ranef", occu = "ranef"),
                             augmentation = c(full = 30),# Número aumentado de especies
                             modelFile = "multmod")

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
load("results/DR_result.R")

# Resultados --------------------------------------------------------------

# Extraemos lo tabla de valores estimados
modresult <- as.data.frame(summary(fit.commu)[["statistics"]])
View(modresult)

# Gráficos de predicción y de coeficientes

# Otra gran ventaja de CamtrapR es que permite gráficar de manera muy sencilla la predicción posterior del modelo. Veamos que pasa con la ocupación de cada especie

plot_effects(comu_model, # El modelo
             fit.commu, # El objeto ajustado
             submodel = "state") # el parámetro de interés

# Ahora con los coeficientes estimados

plot_coef(comu_model,
          fit.commu,
          submodel = "state")

# Realizamos el mismo procedimiento para el submodelo de detección

plot_effects(comu_model,
             fit.commu,
             submodel = "det")


# y sus respectivos coeficientes

plot_coef(comu_model,
          fit.commu,
          submodel = "det")

# Ahora lo que nos interesa. Llevamos todo este viaje para estimar el número de especies en la comunidad

# Valor de Ntotal, es decir del número de especies estimado
(riqueza_est <- modresult["Ntotal",])

# Veamos el gráfico de la distribución posterior
mcmc_areas(fit.commu, # objeto jags
           pars= "Ntotal", # parámetro de interés
           point_est = "mean",
           prob = 0.95) # intervalos de credibilidad



# La estimación no se ve muy bien, hay que verificar los trace plots

mcmc_trace(fit.commu, pars = "Ntotal")

# Debería verse como un cesped, muy probablemente necesitamos muchas mas iteraciones para este modelo

gd <- as.data.frame(gelman.diag(fit.commu,  multivariate = FALSE)[[1]])
gd["Ntotal",]

#La prueba de Gelman-Rubin debe ser ~1 para considerar que hay buena convergencia. Aunque tenemos un valor bueno para Ntotal, hay varios valores de omega con NA, eso puede estar causando los problemas.


# Comparando con métodos clásicos -----------------------------------------


# Formatear los datos a un vector de frecuencia
inci_Chao <- ylist %>%  # historias de captura
  map(~rowSums(.,na.rm = T)) %>% # sumo las detecciones en cada sitio
  reduce(cbind) %>% # unimos las listas
  t() %>% # trasponer la tabla
  as_tibble() %>% #formato tibble
  mutate_if(is.numeric,~(.>=1)*1) %>%  #como es incidencia, formateo a 1 y 0
  rowSums() %>%  # ahora si la suma de las incidencias en cada sitio
  as_tibble() %>% 
 add_row(value= 67, .before = 1) %>%  # el formato requiere que el primer valor sea el número de sitios
  as.matrix() # Requiere formato de matriz



# Calcular la riqueza con estimadores no paramétricos
chao_sp <- ChaoSpecies(inci_Chao, datatype = "incidence_freq")

NIChao <- chao_sp$Species_table[4,c(1,3,4)] # Extraer valores de IChao

Nocu<- mcmc_intervals(fit.commu, pars = "Ntotal", prob = 0.95,prob_outer = 0.99, point_est = "mean")[[1]] %>%  # Extraer valores del bayes plot
  select(m,l,h) %>% # Seleccionar columnas
  rename("Estimate"= m, # Renombrarlas
         "95%Lower"= l,
         "95%Upper"= h)


# Unir en un solo dataframe
Nplotdata <- rbind(IChao=NIChao, DR.mod=Nocu) %>% 
  as.data.frame() %>% 
  rownames_to_column(.)

# Gráfico para comparar la riqueza estimada
plotN <- ggplot(Nplotdata, aes(x=rowname, y= Estimate, col=rowname))+
  geom_point(aes(shape=rowname),size=3)+
  geom_errorbar(aes(ymin= `95%Lower`, ymax= `95%Upper`), width=.3, size=1)+
  labs(x="Estimador de riqueza",y="Número de especies estimado", title = "Diferencia de los estimadores de riqueza")+
  theme_classic()+
  theme(text=element_text(size = 13), plot.title = element_text(hjust= 0.5), legend.position = "none")

plotN
