# Estimación de la diversidad por medio de "DiversityOccupancy"
# Gabriel Andrade

# 1. Cargar librerias ----

# Primero instalar la paquetería
# install.packages("DiversityOccupancy") 

library(DiversityOccupancy) 
library(tidyverse)
library(hillR) # Estimar diversidad
library(ggeffects) #gráficas de pred para glm
library(beepr) # Opcional para avisar R termine
library(tictoc) # Opcional para tomar el tiempo de la función

# 2. Cargar los datos ----
files <- list.files("data/", full.names = T, pattern="*.csv") # lista de arhcivos.csv guardados en la carpeta "data
hists <- lapply(files, read.csv) # Leer cada archivo como un csv

# Juntar todas las bases de datos
data <- hists %>% 
  reduce(full_join, by="X") %>% # Aquí los unimos por el nombre de la estación de muestreo
  select(-X) # Elimino la columna de estación
view(data)

# Ahora cargamos las covariables
covs <- read.csv("data/covs/std_covs.csv", sep = ";")

# 3. Modelos de abundancia ----
# 3.1 Modelos básicos ----
cam_diver <- diversityoccu(pres = data, # La matriz de datos
                           sitecov = covs, # Covariables de sitio
                           obscov = NULL, # no tenemos covariables de observación,
                           spp = 16, # Número de especies
                           form = ~ Effort + Slope ~ SATVI, # Formula del modelo p(), lambda()
                           dredge = FALSE # En este primer ejemplo no usaremos AIC
)

# Veamos el modelo de la especie 2

cam_diver$models[[2]] # Modelo para la especie 2

# El problema es que no sabemos si todas las especies responden de la misma manera a las covariables que usamos. Debemos hacer una inferencia multi-modelo para escoger el mejor modelo para cada especie

#3.2 Mejores modelos con AIC ----
# Vamos a correr de nuevo la función pero esta vez escogerá el mejor modelo para cada especie
tic()
cam_diver_AIC <- diversityoccu(pres = data, # La matriz de datos
                           sitecov = covs, # Covariables de sitio
                           obscov = NULL, # no tenemos covariables de observación,
                           spp = 16, # Número de especies
                           form = ~ Effort + Slope ~ SATVI, # Formula del modelo p(), lambda()
                           dredge = TRUE # escoge los mejores modelos con AIC
); beep(sound = 8)
toc() 

# El código puede tardar mucho más dependiendo de la cantidad de datos, por lo que es mejor guardar el objeto en R para no correrlo cada vez que abramos el script.
# save(cam_diver_AIC, file = "results/diver_AIC.R") # guarda el objeto en la carpeta results

# load("results/diver_AIC.R") # Carga los resultados

cam_diver_AIC$models[[2]] # El nuevo modelo para la sp2

# 3.3 Gráfico de predicción ----
# Gráfico de predicción para la especie 11
responseplot.abund( batch = cam_diver_AIC, # objeto creado con diversityoccu
                    spp = 11, # número o nombre de la sp
                    variable= SATVI # variable  
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
responseplot.diver(glm.div, SATVI)
responseplot.diver(glm.div, Dpop_G)
responseplot.diver(glm.div, Slope)

# 4.2 Diversidad para cada indice de entropía

# Shanon
cam_diver_sh <- diversityoccu(pres = data, 
                              sitecov = covs, 
                              obscov = NULL, 
                              spp = 16, 
                              form = ~ Effort + Slope ~ SATVI,
                              dredge = TRUE, 
                              index = "shannon" #<<
)
#save(cam_diver_sh, file = "results/shanon.R")


# simpson
cam_diver_sim <- diversityoccu(pres = data, 
                               sitecov = covs, 
                               obscov = NULL, 
                               spp = 16, 
                               form = ~ Effort + Slope ~ SATVI,                            dredge = TRUE,  
                               index = "simpson" #<<
)
#save(cam_diver_sim, file = "results/simpson.R")


# inverso de simpson
cam_diver_inv <- diversityoccu(pres = data, 
                               sitecov = covs, 
                               obscov = NULL, 
                               spp = 16, 
                               form = ~ Effort + Slope ~ SATVI,                            dredge = TRUE,  
                               index = "invsimpson" #<<
); beep(sound = 8)

#save(cam_diver_inv, file = "results/invsimpson.R")

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
glm_hill <- cbind(hill_div, covs)

# Ajustamos un GLM
glm_q1 <- glm(q1~ Dpop_G, family = gaussian, data = glm_hill)

# Gráfica del modelo
plot_q1 <- ggpredict(glm_q1, terms = "Dpop_G")
plot(plot_q1)+
  labs(y= "Diversidad q1", x= "Distancia a poblados (estandarizado)")
  theme_classic()
