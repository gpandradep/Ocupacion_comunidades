##%######################################################%##
#                                                          #
####          Estimación de la diversidad por           ####
####  medio de "DiversityOccupancy" # Gabriel Andrade   ####
#                                                          #
##%######################################################%##


# 1. Cargar librerías ----

# Primero instalar la paquetería
# install.packages("DiversityOccupancy") 

library(DiversityOccupancy) 
library(camtrapR)
library(tidyverse)
library(hillR) # Estimar diversidad
library(ggeffects) #gráficas de pred para glm
library(beepr) # Opcional para avisar R termine
library(tictoc) # Opcional para tomar el tiempo de la función

# 2. Cargar los datos ----


# Cargamos la tabla de registros de las especies
registers <-  read.csv("Data/Survey/recordTable_OC.csv")
table(registers$Species)



# Cargamos la tabla de operación de cámaras
CToperation <-  read.csv("Data/Survey/CTtable_OC.csv") 

# Generamos la matríz de operación de las cámaras

camop <- camtrapR::cameraOperation(CTtable= CToperation, # Tabla de operación
                                   stationCol= "Station", # Columna que define la estación
                                   setupCol= "Setup_date", #Columna fecha de colocación
                                   retrievalCol= "Retrieval_date", #Columna fecha de retiro
                                   hasProblems= T, # Hubo fallos de cámaras
                                   dateFormat= "%Y-%m-%d") # Formato de las fechas



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
    occasionLength       = 10,  # Colapso de las historias a 10 días
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

data <- ylist %>% # Lista con los datos
  reduce(cbind) # Unimos las historias de captura en un solo dataframe

# Ahora cargamos las covariables (previamente estandarizadas)
covars <- read.csv("Data/Covs/stdcovs_OC.csv") %>% 
  select(-"...1", -X, -Station, -Cam, -Cluster) #


# 3. Modelos de abundancia ----
# 3.1 Modelos básicos ----
cam_diver <- diversityoccu(pres = data, # La matriz de datos
                           sitecov = covars, # Covariables de sitio
                           obscov = NULL, # no tenemos covariables de observación,
                           spp = 17, # Número de especies
                           form = ~ Effort + Slope ~ Dcrops, # Formula del modelo p(), lambda()
                           dredge = FALSE # En este primer ejemplo no usaremos AIC
)


# Veamos el modelo de la especie 9

cam_diver$models[[9]] # Modelo para la especie 9

# El problema es que no sabemos si todas las especies responden de la misma manera a las covariables que usamos. Debemos hacer una inferencia multi-modelo para escoger el mejor modelo para cada especie

#3.2 Mejores modelos con AIC ----
# Vamos a correr de nuevo la función pero esta vez escogerá el mejor modelo para cada especie
tic()
cam_diver_AIC <- diversityoccu(pres = data, # La matriz de datos
                           sitecov = covars, # Covariables de sitio
                           obscov = NULL, # no tenemos covariables de observación,
                           spp = 17, # Número de especies
                           form = ~ Effort + Slope ~ Dcrops, # Formula del modelo p(), lambda()
                           dredge = TRUE # escoge los mejores modelos con AIC
); beep(sound = 8)
toc() 



# El código puede tardar mucho más dependiendo de la cantidad de datos, por lo que es mejor guardar el objeto en R para no correrlo cada vez que abramos el script.
# save(cam_diver_AIC, file = "results/diver_AIC.R") # guarda el objeto en la carpeta results

# load("results/diver_AIC.R") # Carga los resultados

cam_diver_AIC$models[[9]] # El nuevo modelo para la sp9

# 3.3 Gráfico de predicción ----
# Gráfico de predicción para la especie 3
responseplot.abund( batch = cam_diver_AIC, # objeto creado con diversityoccu
                    spp = 3, # número o nombre de la sp
                    variable= Dcrops # variable  
)

# 4. Modelando la diversodad ---- 
tic()
glm.div <- model.diversity(DivOcc = cam_diver_AIC,
                           method = "h",
                           delta = 2,
                           squared = T
); beep(sound= 8)
toc()

# save(glm.div, file= "results/diver_glm.R")  # guardar los análisis
#load("results/diver_glm.R") cargarlos
AICtab <- glm.div$Table

# 4.1 Respuesta Gráfica de la diversidad a distintas  variables ----
responseplot.diver(glm.div, Dcrops)
responseplot.diver(glm.div, Dpop_G)
responseplot.diver(glm.div, Slope)

# 4.2 Diversidad para cada indice de entropía

# Shanon
cam_diver_sh <- diversityoccu(pres = data, 
                              sitecov = covars, 
                              obscov = NULL, 
                              spp = 17, 
                              form = ~ Effort + Slope ~ Dcrops,
                              dredge = TRUE, 
                              index = "shannon" #<<
)
save(cam_diver_sh, file = "results/shanon.R")


# simpson
cam_diver_sim <- diversityoccu(pres = data, 
                               sitecov = covars, 
                               obscov = NULL, 
                               spp = 17, 
                               form = ~ Effort + Slope ~ Dcrops,                            dredge = TRUE,  
                               index = "simpson" #<<
)
save(cam_diver_sim, file = "results/simpson.R")


# inverso de simpson
cam_diver_inv <- diversityoccu(pres = data, 
                               sitecov = covars, 
                               obscov = NULL, 
                               spp = 17, 
                               form = ~ Effort + Slope ~ Dcrops,                            dredge = TRUE,  
                               index = "invsimpson" #<<
); beep(sound = 8)

save(cam_diver_inv, file = "results/invsimpson.R")

# 5. Calculamos número efectivo de especies ----

# Extraer los datos de abundancia
hill_data <- cam_diver_inv[[4]] %>% 
  select(-h)

# calcular los perfiles de diversidad
q0 <- hill_taxa(hill_data, q=0) 
q1 <- hill_taxa(hill_data, q=1)
q2 <- hill_taxa(hill_data, q=2)

# Unimos en el mismo data.frame
hill_div <- data.frame(q0=q0, q1=q1, q2=q2)

# Juntamos con las covariables
glm_hill <- cbind(hill_div, covars)

# Ajustamos un GLM
glm_q1 <- glm(q1~ Dcrops, family = gaussian, data = glm_hill)

# Gráfica del modelo
plot_q1 <- ggpredict(glm_q1, terms = "Dcrops")
plot(plot_q1)+
  labs(y= "Diversidad q1", x= "Distancia a poblados (estandarizado)")
  theme_classic()
  
plot_q1
