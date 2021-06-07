
## Variables a estudiar en correspondencia con la OMS:
# 1 Volumen; 
# 2 Concentraci�n; 
# 3 Movilidad R�pida
# 4 Movilidad Total (MR + ML); 
# 5 Morfolog�a; 
# 6 Vitalidad
## Extra
# 7 Leucositos 

library(tidyverse)

# todas las variables numéricas
E2 <- read.csv("D:/Mabe/Espermatograma-num-13.01.2020a.csv")

########################################################################################
# Preparar datos
library(haven)
data.esp <- read_sav("D:/Databases/espermograma-SPSS-21-01-2021.sav")
# agregar variable OMS Mov. Total
data.esp$MovTotal <- data.esp$MovRap + data.esp$MovLen
# 1 Volumen; # 2 Concentraci�n; # 3 Movilidad R�pida# 4 Movilidad Total (MR + ML); 
# 5 Morfolog�a; # 6 Vitalidad ## Extra # 7 Leucositos 
data.esp_n <- data.esp %>% select(Volumen, Concentr, MovRap, MovTotal, 
                                  Morfologia, Vitalidad, Leuco2) 
data.esp_n <- drop_na(data.esp_n) # Eliminar filas con valores perdidos
# data.esp_n$Color <- NULL # sin color
# Analizar  variables del dataframe
glimpse(data.esp_n)
resumen_ini <- pastecs::stat.desc(data.esp_n, basic = FALSE, norm = TRUE)
round(resumen_ini, 2)
data.esp_ab <- data.esp_n %>% mutate(ab = MovRap + MovLen)
resumen_ini_ab <- pastecs::stat.desc(data.esp_ab, basic = FALSE, norm = TRUE)
resumen_ini <- psych::describe(resumen_ini)
write_excel_csv(resumen_ini, "D:/Databases/resumen_ini.csv", row.names=T)
write_sav(resumen_ini, "D:/Databases/resumen_ini.sav", compress = FALSE, row.names=TRUE)
write.csv(resumen_ini, "D:/Databases/resumen_ini.csv", row.names = TRUE) 
summary(resumen_ini_ab$ab)
data.esp_n$ab <- data.esp_n$MovRap + data.esp_n$MovLen

psych::describe(data.esp_n)

##### Correlaciones
library(GGally)
ggpairs(data.esp_n)

########################################################################################

boxplot(data.esp$MovRap)
ggplot(data.esp, aes(MovRap, Volumen)) + 
  geom_boxplot() +
  #geom_point() +
  geom_jitter()


###########  Normalizaci?n

# Estandarizar con scale() y convertirlo en dafaframe porque scale devuelve uma matriz
# use the scale() function to z-score standardize a data frame
data.esp_z <- as.data.frame(scale(data.esp_n))

# create normalization function usando m�x y m�n
normalize <- function(x) {
  return ((x - min(x)) / (max(x) - min(x)))
}
# normalize the data con función normalize anterior
E2_n <- as.data.frame(lapply(E2, normalize))

########################################################################################

# Comprebar la tendencia de los datos a la aglomeración


library(factoextra) 
library(ggpubr) 
# realizar un PCA para reducir dimencionalidad
pca_data.esp_z <- prcomp(data.esp_z) 
summary(pca_data.esp_z)
# Gráfico de puntos para ver posibles agrupaciones visualmente
p1 <- fviz_pca_ind(X = pca_data.esp_z,  
                   axes = c(1, 2), # Seleccionar las dimenciones a graficar
                   geom = "point", title = "PCA - datos espermograma", 
                   pallete = "jco") + 
          theme_bw() + 
          theme(legend.position = "bottom")
p1

# Estadístico H (hopkins) para el set de datos espermograma 
library(clustertend) 
set.seed(321) 
hopkins(data = data.esp_z, n = nrow(data.esp_z) - 1) # mientras más cercano a 0 en mejor, cercano a 0,5 más malo
# Visual Assessment of cluster Tendency (VAT)
dist_data.esp_z <- dist(data.esp_z, method = "euclidean")
p3 <- fviz_dist(dist.obj = dist_data.esp_z, show_labels = FALSE) + 
                labs(title = "Datos espermograma") + 
                theme(legend.position = "bottom")
p3


########################################################################################
# Número óptimo de clusters

# Elbow method
fviz_nbclust(x = data.esp_z, FUNcluster = kmeans, method = "wss") + 
  labs(title = "N�mero �ptimo de clusters")
# con NbClust
library(factoextra) 
library(NbClust) 
numero_clusters <- NbClust(data = data.esp_z, distance = "euclidean", 
                           min.nc = 2, max.nc = 10, method = "median", index = "alllong")
fviz_nbclust(numero_clusters)

########################################################################################

# Selección de los mejores algorítmos de cluster

library(clValid) 
library(kohonen) # para clMethods = "som"
library(mclust)  # para clMethods = "model"
comparacion <- clValid(obj = data.esp_z, nClust = 2:6, 
                       clMethods = c("hierarchical", "kmeans", "pam", "model"), 
                       validation = c("stability", "internal")) 
summary(comparacion)

########################################################################################
### Algoritmos de clusters

## K-means

set.seed(101) 
km_clusters <- kmeans(x = data.esp_z, centers = 2, nstart = 50) 
km_clusters
table(km_clusters$cluster)
fviz_cluster(object = km_clusters, data = data.esp_z, show.clust.cent = TRUE, 
             ellipse.type = "convex", star.plot = TRUE, repel = TRUE,
             axes = c(1, 2)) + # para usar las dimensiones de PCA
  theme_bw() + 
  theme(legend.position = "none")
fviz_cluster(km_clusters, data = data.esp_z, ellipse.type = "norm")

# agregar cluster al dataframe
data.esp_z$kc = km_clusters$cluster
# invertir cluster para que sea al m�s grande igual a kmeans
data.esp_z <- data.esp_z %>% mutate(kc = recode(kc, "1" = "2", "2" = "1"))
table(data.esp_z$kc)

## PAM

set.seed(123) 
pam_clusters <- pam(x = data.esp_z, k = 2, metric = "manhattan") 
pam_clusters
fviz_cluster(object = pam_clusters, data = data.esp_z, ellipse.type = "euclid", 
             repel = TRUE) + 
  theme_bw() + 
  theme(legend.position = "none")

## árboles herarquicos

# Matriz de distancias euclídeas 
mat_dist <- dist(x = data.esp_z, method = "euclidean") 
# Dendrogramas con linkage complete y average 
hc_completo <- hclust(d = mat_dist, method = "complete") 
hc_ward.D <- hclust(d = mat_dist, method = "ward.D")  

cor(x = mat_dist, cophenetic(hc_completo))
df <- df %>% mutate(Cl_km = recode(Cl_km, "1" = "2", "2" = "1"))

data.esp_z$hc = cutree(hc_completo, k = 2)

table(cluster)
table(data.esp_z$hc)

hc_completo_clus <- cutree(hc_completo, 2)
hc_average_clus <- cutree(hc_average, 2)
# Crear dendograma y cortar e 2 clusters
set.seed(101) 
fviz_dend(x = hc_completo, k = 2, cex = 0.6) + 
  geom_hline(yintercept = 8.5, linetype = "dashed")
fviz_dend(x = hc_average, k = 2, cex = 0.6) + 
  geom_hline(yintercept = 8.5, linetype = "dashed")

fviz_cluster(object = list(data = data.esp_z[,1:7], cluster = cutree(hc_completo, k = 2)), 
             ellipse.type = "convex", repel = TRUE, show.clust.cent = FALSE) + 
  theme_bw()

predict(hc_completo, data.esp_z[1,])

## hkmeans

library(clustertend)
hkmeans_cluster <- hkmeans(x = data.esp_z[,1:7], hc.metric = "euclidean", 
                           hc.method = "complete", k = 2) 
fviz_cluster(object = hkmeans_cluster, pallete = "jco", repel = TRUE, 
             ellipse.type = "convex", show.clust.cent = TRUE) + 
  theme_bw() + labs(title = "Hierarchical k-means Cluster plot")
summary(hkmeans_cluster)
table(hkmeans_cluster$cluster)
data.esp_z$hkc <- hkmeans_cluster$cluster

### Otros cluster que no uso ----------------#############3
## mclust
model_clustering <- Mclust(data = data.esp_z, G = 2) 
summary(model_clustering)
fviz_mclust(model_clustering, what = "classification", 
            geom = "point", pallete = "jco")
## DBSCAN
library(fpc) 
library(dbscan)
dbscan::kNNdistplot(data.esp_z[,1:7], k = 2)
set.seed(321) 
# DBSCAN con epsilon = 0.15 y minPts = 5 
dbscan_cluster <- fpc::dbscan(data = data.esp_z[,1:7], eps = 2.15, MinPts = 2)
# Resultados de la asignación 
head(dbscan_cluster$cluster)
# Visualización de los clusters 
fviz_cluster(object = dbscan_cluster, data = data.esp_z[,1:7], stand = FALSE, 
             geom = "point", ellipse = FALSE, show.clust.cent = FALSE, 
             pallete = "jco") + theme_bw() + 
  theme(legend.position = "bottom")

########################################################################################

# Plotear 2 variables para ver cómo quedaron algunas relaciones con los colores del cluster
E4<- as.data.frame(data.esp_n)
E4 <- drop_na(E4)
plot(E4$MovRap, E4$Vitalidad, col = hkmeans_cluster$cluster)
plot(E4$MovRap, E4$Vitalidad, col = model_clustering$classification)
par(mfrow=c(1,1))
write_excel_csv(data.esp_n, "D:/Databases/data.esp_n.csv")

########################################################################################
###    Descripción de los Clusters

data.esp_n <- cbind( data.esp_n, data.esp_z[,8:10])
psych::describeBy(data.esp_n[,c(1:7,10)], group = data.esp_n$hkc)

#le pedimos un peque resumen 

resumen_post <- data.esp_n %>%
  mutate(Cl_km = km_clusters$cluster) %>%
  group_by(Cl_km) %>%
  summarise_all("mean")

write_excel_csv(resumen_post, "D:/Databases/resumen_post.csv")

# Correlacionar resultados

cor(as.numeric(data.esp_n$hkc),as.numeric(data.esp_n$kc))

# ARI Cl_km con Cl_h - Computes the adjusted Rand index comparing two classifications.
library(mclust)
adjustedRandIndex(as.numeric(data.esp_n$hkc),as.numeric(data.esp_n$kc))

# Reorganizar los niveles para que en df_h el 1 sea el de referencia
#df_h$Cl_h <- relevel(df_h$Cl_h, ref = "1") 

#niveles de los factores
levels(as.factor(data.esp_z$hkc))

# Tabla de confución entre Cl_km y Cl_h
table(data.esp_z$kc,data.esp_z$hkc)

# Tamaño de los clusters
table(data.esp_z$hkc)


#Graficar el cluster por variables para ver su variablilidad entre grupos
# Primero pivotear a lo largo manteniendo el factor del cluster
df_pivot_n <- data.esp_n[,c(1:7,10)] %>% pivot_longer(cols = -hkc, names_to = "Variables", 
                                          values_to = "Valores")
df_pivot_z <- data.esp_z[,c(1:7,10)] %>% pivot_longer(cols = -hkc, names_to = "Variables", 
                                          values_to = "Valores")


# Graficar con ggplot por variables por facetas kmedias
ggplot(data = df_pivot_n) +
  geom_boxplot(aes(as.factor(hkc), Valores, color= Variables),
               position = position_dodge()) + 
  #geom_jitter(aes(Cl_km, Valores, color= Variables)) +
 # xlab("Clases 1 y 2") +
  facet_wrap(~df_pivot_n$Variables, scales = "free")
# Graficar con ggplot por variables por facetas hclust
ggplot(data = df_pivot_z) +
  geom_boxplot(aes(as.factor(hkc), Valores, color= Variables),
               position = position_dodge()) + 
  #geom_jitter(aes(Cl_km, Valores, color= Variables)) +
  # xlab("Clases 1 y 2") +
  facet_wrap(~df_pivot_z$Variables, scales = "free")

# Variables sin normalizar
df_pivot_n$hkc <- as.factor(df_pivot_n$hkc)
ggplot(df_pivot_n, aes(as.factor(x = Variables), y = Valores,group=hkc, colour = hkc)) + 
  stat_summary(fun = mean, geom="pointrange", size = 1, aes(shape = hkc))+
  stat_summary(geom="line") +  
  geom_point(aes(shape=hkc)) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) # se giran las x

# Variables sin normalizar
df_pivot_z$hkc <- as.factor(df_pivot_z$hkc)
ggplot(df_pivot_z, aes(as.factor(x = Variables), y = Valores,group=hkc, colour = hkc)) + 
  stat_summary(fun = mean, geom="pointrange", size = 1, aes(shape = hkc))+
  stat_summary(geom="line") +  
  geom_point(aes(shape=hkc)) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) # se giran las x

# regression tree using rpart
library(rpart)
m.rpart <- rpart(Vitalidad ~ Volumen + Concentr + MovRap + MovLen + MovIns + Inmoviles + Morfologia, data = E3)
# get basic information about the tree
m.rpart
# get more detailed information about the tree
summary(m.rpart)
# use the rpart.plot package to create a visualization
library(rpart.plot)
# a basic decision tree diagram
rpart.plot(m.rpart, digits = 3)

###########
#### Variables m?s importantes con Random Forest 
## Se aplica randomForest al dataset con la clasificaci?n incluida

library(randomForest)
library(ggpubr)

data.esp_z$hkc<-as.factor(data.esp_z$hkc)
df <- data.esp_z[,c(1:7,10)]

modelo_randforest <- randomForest(formula = hkc ~ . , 
                                  data = df,
                                  mtry = 5, 
                                  importance = TRUE, 
                                  na.action = na.omit,
                                  ntree = 1000)

importancia <- as.data.frame(modelo_randforest$importance) 
importancia <- rownames_to_column(importancia,var = "variable")

p1 <- ggplot(data = importancia, aes(x = reorder(variable, MeanDecreaseAccuracy), 
                               y = MeanDecreaseAccuracy, 
                               fill = MeanDecreaseAccuracy)) + 
  labs(x = "variable", title = "Reducci�n de Accuracy") + 
  geom_col() + coord_flip() + theme_bw() + 
  theme(legend.position = "bottom")
p2 <- ggplot(data = importancia, aes(x = reorder(variable, MeanDecreaseGini), 
                                     y = MeanDecreaseGini, 
                                     fill = MeanDecreaseGini)) + 
  labs(x = "variable", title = "Reducci�n de Accuracy") + 
  geom_col() + coord_flip() + theme_bw() + 
  theme(legend.position = "bottom")

ggarrange(p1, p2)

######## Comparaci?n de las variables del dataset con t-test
###     agrupadas por la claseificaci?n Cl_km

data.esp_z
t.test(Concentr ~ Cl_km, data = data.esp_z)
t.test(. ~ Cl_km, data = data.esp_z)
write.csv(data.esp_z, "D:/Databases/data.esp_z.csv", row.names = TRUE)
write.csv(data.esp_n, "D:/Databases/data.esp_n.csv", row.names = TRUE) 

########################################################################################

###############3 Comparar valores con la OMS

####################3

valores <- list()
pacientes <- 128
valores_inf <-   dim(data.esp_n %>% filter(Volumen >= 1.5))[1]
valores$vol <- (valores_inf / pacientes) * 100
valores_inf <-   dim(data.esp_n %>% filter(Concentr >= 39))[1]
valores$concentr <- (valores_inf / pacientes) * 100
valores_inf <-   dim(data.esp_n %>% filter(MovRap + MovLen >= 40))[1]
valores$movT <- (valores_inf / pacientes) * 100
valores_inf <-   dim(data.esp_n %>% filter(MovRap >= 32))[1]
valores$MovRap <- (valores_inf / pacientes) * 100
valores_inf <-   dim(data.esp_n %>% filter(Vitalidad >= 58))[1]
valores$vitalidad <- (valores_inf / pacientes) * 100
valores_inf <-   dim(data.esp_n %>% filter(Morfologia >= 4))[1]
valores$morf <- (valores_inf / pacientes) * 100
valores

### Cuantos cumplen todos los valores minimos de la OMS
valores_todos <- data.esp_n %>% filter(Volumen >= 1.5, 
                                       Concentr >= 39, 
                                       MovRap + MovLen >= 40, 
                                       MovRap >= 32, 
                                       Vitalidad >= 58, 
                                       Morfologia >= 4)
(dim(valores_todos)[1] / pacientes) * 100

# adicionar variable con la clasificaci?n segun valores OMS
data.esp_n <- data.esp_n %>% mutate(OMS = if_else(Volumen >= 1.5 & 
                       Concentr >= 39 & 
                       MovRap + MovLen >= 40 &
                       MovRap >= 32 &
                       Vitalidad >= 58 & 
                       Morfologia >= 4, 2, 1))

# comparar resultados con cluster y OMS
table(data.esp_n$OMS, data.esp_n$Cl_km)
table(data.esp_n$OMS)

#######################################################3

#### REGRESION #####################33


library(psych)
pairs.panels(E2)
pairs.panels(E2[c("Volumen", "Color", "Mov.lenta", "Viabilidad")])
plot(E3$Color,E3$Vitalidad)
table(E3$Color)

# crear modelo
E3_model <- lm(Vitalidad ~ ., data = E3)
# ver estimaciones
E3_model
# ver resumen del modelo
summary(E3_model)

# hacer predicciones con el modelo y la misma base
E3$pred <- predict(E3_model, E3)

# relación entre predicción y valor real
cor(E3$pred, E3$Vitalidad)
pairs.panels(E3[c("pred", "Vitalidad")])


###########################################
# Kmeans  con H2O

#######################3333


datos_kmeansMulti <- data.frame(E[,6],E[,10:16])
write.csv(datos_kmeansMulti, "D:/Databases/data.esp_n.csv")
library(h2o)
# Creación de un cluster local con todos los cores disponibles.
h2o.init(
  ip = "localhost",
  # -1 indica que se empleen todos los cores disponibles.
  nthreads = -1,
  # Máxima memoria disponible para el cluster.
  max_mem_size = "6g"
)
e3 <- h2o.uploadFile(path = "D:/Databases/data.esp_n.csv")
e3 <- h2o.uploadFile(path = "D:/Mabe/E3.csv")
e.km.h2o <- h2o.kmeans(training_frame = e3, k = 2)
# Train a H2O K-Means model, auto-estimate best value for k  (note: k here is maximum possible k)
fit <- h2o.kmeans(training_frame = e3, k = 20, estimate_k = TRUE)
print(e.km.h2o)



