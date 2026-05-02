#Starting

##### SCRIPT_P._percellens_SI #####

###### CDA (Canonical Discriminant Analysis) ######

# Carregar pacotes
library(tidyverse)
library(rstatix)

# Importar dados
library(readr)

dados <- read_delim(pipe("pbpaste"), delim = "\t")

head(dados)
str(dados)
names(dados)

# Ajustar fatores
dados$SEXO <- factor(dados$SEXO)
dados$MATURIDADE <- factor(dados$MATURIDADE)

# Conferir estrutura
str(dados)

# Normalidade geral
shapiro.test(dados$d13C)
shapiro.test(dados$d15N)

# Sexo
wilcox.test(d15N ~ SEXO, data = dados)
wilcox.test(d13C ~ SEXO, data = dados)

# Maturidade
wilcox.test(d15N ~ MATURIDADE, data = dados)
wilcox.test(d13C ~ MATURIDADE, data = dados)

# δ15N vs tamanho
lm_15N <- lm(d15N ~ CT, data = dados)
summary(lm_15N)

# δ13C vs tamanho
lm_13C <- lm(d13C ~ CT, data = dados)
summary(lm_13C)

# Por classe de tamanho
dados <- dados %>%
  mutate(classe_tam = ifelse(CT < median(CT), "Pequeno", "Grande"))

dados$classe_tam <- factor(dados$classe_tam)

# Testes
wilcox.test(d15N ~ classe_tam, data = dados)
wilcox.test(d13C ~ classe_tam, data = dados)

# Boxplot por sexo
ggplot(dados, aes(x = SEXO, y = d15N)) +
  geom_boxplot() +
  theme_minimal()

ggplot(dados, aes(x = MATURIDADE, y = d13C)) +
  geom_boxplot() +
  theme_minimal()

# Regressão
ggplot(dados, aes(x = CT, y = d15N)) +
  geom_point() +
  geom_smooth(method = "lm") +
  theme_minimal()

ggplot(dados, aes(x = CT, y = d13C)) +
  geom_point() +
  geom_smooth(method = "lm") +
  theme_minimal()

# Correlação
cor.test(dados$CT, dados$d15N, method = "spearman")
cor.test(dados$CT, dados$d13C, method = "spearman")

# PERMANOVA
library(vegan)

iso.scaled <- scale(dados[, c("d13C", "d15N")])

adonis2(iso.scaled ~ SEXO, data = dados, method = "euclidean")
adonis2(iso.scaled ~ MATURIDADE, data = dados, method = "euclidean")
adonis2(iso.scaled ~ classe_tam, data = dados, method = "euclidean")

# SIBER
install.packages("SIBER")
library(SIBER)
library(dplyr)

# SEXO
siber.data <- dados %>%
  mutate(
    group = ifelse(SEXO == "F", 1, 2),
    community = 1,
    iso1 = d13C,
    iso2 = d15N
  ) %>%
  select(iso1, iso2, group, community) %>%
  as.data.frame()

colnames(siber.data)

ellipses <- createSiberObject(siber.data)

plotSiberObject(ellipses)

as.data.frame()

ellipses <- createSiberObject(siber.data)
plotSiberObject(ellipses)
groupMetricsML(ellipses)

# MATURIDADE
siber.mat <- dados %>%
  mutate(
    group = ifelse(MATURIDADE == "Imaturo", 1, 2),
    community = 1,
    iso1 = d13C,
    iso2 = d15N
  ) %>%
  select(iso1, iso2, group, community) %>%
  as.data.frame()

siber.obj.mat <- createSiberObject(siber.mat)

plotSiberObject(siber.obj.mat)

groupMetricsML(siber.obj.mat)

#plot
library(ggplot2)
library(dplyr)
library(SIBER)

dados$SEXO <- trimws(dados$SEXO)

SIBER_SEX <- ggplot(dados, aes(x = d13C, y = d15N, color = SEXO, fill = SEXO)) +
  
  stat_ellipse(type = "norm", level = 0.40,
               linewidth = 0.5, alpha = 0.20,
               geom = "polygon") +
  
  geom_point(size = 5) +
  
  scale_color_manual(values = c("Female" = "#D55E00",
                                "Male" = "#0072B2")) +
  
  scale_fill_manual(values = c("Female" = "#D55E00",
                               "Male" = "#0072B2")) +
  
  labs(
    x = expression(delta^13*C~"\u2030"),
    y = expression(delta^15*N~"\u2030"),
    color = "Sex",
    fill = "Sex"
  ) +
  
  theme_classic(base_size = 20)

# mostrar
SIBER_SEX

# salvar
ggsave("SIBER_SEX.tiff",
       plot = SIBER_SEX,
       width = 8,
       height = 6,
       dpi = 600,
       compression = "lzw")




SIBER_SEX <- ggplot(dados, aes(x = d13C, y = d15N, color = SEXO, fill = SEXO)) +
  
  stat_ellipse(type = "norm", level = 0.40, linewidth = 0.5, alpha = 0.20, geom = "polygon") +
  
  geom_point(size = 5) +
  
  scale_color_manual(values = c("F" = "#D55E00", "M" = "#0072B2")) +
  scale_fill_manual(values = c("F" = "#D55E00", "M" = "#0072B2")) +
  
  labs(
    x = expression(delta^13*C~"\u2030"),
    y = expression(delta^15*N~"\u2030"),
    color = "Sex",
    fill = "Sex"
  ) +
  
  theme_classic(base_size = 20)

SIBER_SEX <- ggplot(dados, aes(x = d13C, y = d15N, color = SEXO, fill = SEXO)) +
  
  stat_ellipse(type = "norm", level = 0.40, linewidth = 0.5,
               alpha = 0.20, geom = "polygon") +
  
  geom_point(size = 5) +
  
  scale_color_manual(values = c("F" = "#D55E00", "M" = "#0072B2")) +
  scale_fill_manual(values = c("F" = "#D55E00", "M" = "#0072B2")) +
  
  labs(
    x = expression(delta^13*C~"\u2030"),
    y = expression(delta^15*N~"\u2030"),
    color = "Sex",
    fill = "Sex"
  ) +
  
  theme_classic(base_size = 20)

ggsave("SIBER_SEX.tiff",
       plot = grafico_sexo,
       width = 8,
       height = 6,
       dpi = 600)

dados$MATURIDADE <- recode(dados$MATURIDADE,
                           "Imaturo" = "Immature",
                           "Maduro" = "Mature")

ggplot(dados, aes(x = d13C, y = d15N, color = MATURIDADE, fill = MATURIDADE)) +
  
  stat_ellipse(type = "norm", level = 0.40, linewidth = 0.5, alpha = 0.20, geom = "polygon") +
  
  geom_point(size = 5) +
  
  scale_color_manual(values = c("Immature" = "#009E73", "Mature" = "#CC79A7")) +
  scale_fill_manual(values = c("Immature" = "#009E73", "Mature" = "#CC79A7")) +
  
  labs(
    x = expression(delta^13*C~"\u2030"),
    y = expression(delta^15*N~"\u2030"),
    color = "Maturity",
    fill = "Maturity"
  ) +
  
  theme_classic(base_size = 20)

SIBER_MATURITY <- ggplot(dados, aes(x = d13C, y = d15N, color = MATURIDADE, fill = MATURIDADE)) +
  
  stat_ellipse(type = "norm", level = 0.40, linewidth = 0.5,
               alpha = 0.20, geom = "polygon") +
  
  geom_point(size = 5) +
  
  scale_color_manual(values = c("Immature" = "#009E73",
                                "Mature" = "#CC79A7")) +
  
  scale_fill_manual(values = c("Immature" = "#009E73",
                               "Mature" = "#CC79A7")) +
  
  labs(
    x = expression(delta^13*C~"\u2030"),
    y = expression(delta^15*N~"\u2030"),
    color = "Maturity",
    fill = "Maturity"
  ) +
  
  theme_classic(base_size = 20)

ggsave("SIBER_MATURITY.tiff",
       plot = SIBER_MATURITY,
       width = 8,
       height = 6,
       dpi = 600)

getwd()

# NicheRover
library(nicheROVER)

sex.data <- dados[, c("d13C", "d15N")]
grupo <- dados$SEXO

X <- split(sex.data, grupo)

post <- lapply(X, function(z)
  niw.post(
    X = as.matrix(z),
    nsamples = 1000
  )
)

ov <- overlap(
  post,
  nreps = 1000,
  nprob = 1000,
  alpha = 0.95
)

ov

# F dentro de M
FM <- ov["F", "M", ]

# M dentro de F
MF <- ov["M", "F", ]

mean(FM)
quantile(FM, c(0.025, 0.975))

mean(MF)
quantile(MF, c(0.025, 0.975))

X2 <- split(dados[, c("d13C","d15N")], dados$MATURIDADE)

post2 <- lapply(X2, function(z)
  niw.post(X = as.matrix(z), nsamples = 1000)
)

ov2 <- overlap(
  post2,
  nreps = 1000,
  nprob = 1000,
  alpha = 0.95
)

ov2

IM <- ov2["Imaturo", "Maduro", ]

MI <- ov2["Maduro", "Imaturo", ]

mean(IM)
quantile(IM, c(0.025, 0.975))

mean(MI)
quantile(MI, c(0.025, 0.975))

#biplot
library(dplyr)
library(ggplot2)

#-------------------------
# 1. Padronizar variáveis
#-------------------------
dados$SEXO <- trimws(dados$SEXO)
dados$MATURIDADE <- trimws(dados$MATURIDADE)

dados$SEXO <- recode(dados$SEXO,
                     "F" = "Female",
                     "M" = "Male")

dados$MATURIDADE <- recode(dados$MATURIDADE,
                           "Imaturo" = "Immature",
                           "Maduro" = "Mature")

#-------------------------
# 2. Criar os 4 grupos
#-------------------------
dados <- dados %>%
  mutate(
    Grupo = case_when(
      SEXO == "Female" & MATURIDADE == "Immature" ~ "Immature female",
      SEXO == "Female" & MATURIDADE == "Mature"   ~ "Mature female",
      SEXO == "Male"   & MATURIDADE == "Immature" ~ "Immature male",
      SEXO == "Male"   & MATURIDADE == "Mature"   ~ "Mature male"
    )
  )

dados$Grupo <- factor(
  dados$Grupo,
  levels = c("Immature female",
             "Mature female",
             "Immature male",
             "Mature male")
)

#-------------------------
# 3. Biplot sem elipses
#-------------------------
grafico_4grupos <- ggplot(
  dados,
  aes(x = d13C,
      y = d15N,
      color = Grupo,
      shape = Grupo)
) +
  
  geom_point(size = 5, alpha = 0.95) +
  
  scale_color_manual(values = c(
    "Immature female" = "#E69F00",
    "Mature female"   = "#D55E00",
    "Immature male"   = "#009E73",
    "Mature male"     = "#0072B2"
  )) +
  
  scale_shape_manual(values = c(
    "Immature female" = 16,
    "Mature female"   = 17,
    "Immature male"   = 15,
    "Mature male"     = 18
  )) +
  
  labs(
    x = expression(delta^13*C~"\u2030"),
    y = expression(delta^15*N~"\u2030"),
    color = "Group",
    shape = "Group"
  ) +
  
  theme_classic(base_size = 20) +
  theme(
    legend.position = "right",
    legend.title = element_text(face = "bold")
  )

grafico_4grupos

ggsave("Biplot_4groups_clean.tiff",
       grafico_4grupos,
       width = 8,
       height = 6,
       dpi = 600)

##### Stomach Content Analysis

#PERMANOVA
library(vegan)

diet <- read_delim(pipe("pbpaste"), delim = "\t")

head(dados)
str(dados)
names(dados)

mat_diet <- diet[, c("DECA","MOLL","FISH","AMPH","OCRUST")]

mat_sqrt <- sqrt(mat_diet)

adonis2(mat_sqrt ~ SEXO, data = diet, method = "bray")
adonis2(mat_sqrt ~ MATUR, data = diet, method = "bray")

dist_diet <- vegdist(mat_sqrt, method = "bray")
anova(betadisper(dist_diet, diet$SEXO))
anova(betadisper(dist_diet, diet$MATUR))

#Levin
# somar peso total por categoria
totais <- colSums(diet[, c("DECA","MOLL","FISH","AMPH","OCRUST")])

# proporções
p <- totais / sum(totais)

p

B <- 1 / sum(p^2)

n <- length(p)

BA <- (B - 1) / (n - 1)

B
BA

#Costello / Amundsen
library(dplyr)
library(tidyr)
library(ggplot2)

# matriz de dieta
mat <- diet[, c("DECA","MOLL","FISH","AMPH","OCRUST")]

# total por indivíduo
diet$total <- rowSums(mat)

# função Costello
costello <- data.frame(
  Prey = colnames(mat),
  FO = NA,
  Pi = NA
)

for(i in 1:ncol(mat)){
  
  prey <- mat[,i]
  
  # indivíduos que comeram a presa
  present <- prey > 0
  
  # frequência ocorrência
  costello$FO[i] <- sum(present) / nrow(mat) * 100
  
  # abundância específica
  costello$Pi[i] <- sum(prey[present]) /
    sum(diet$total[present]) * 100
}

costello

ggplot(costello,
       aes(x = FO, y = Pi, label = Prey)) +
  
  geom_point(size = 5, color = "steelblue") +
  
  geom_text(size = 4, nudge_y = 3) +
  
  geom_hline(yintercept = 50, linetype = 2, color = "gray60") +
  geom_vline(xintercept = 50, linetype = 2, color = "gray60") +
  
  labs(
    x = "Frequency of occurrence (%)",
    y = "Prey-specific abundance (%)"
  ) +
  
  theme_classic(base_size = 14)

grafico_costello <- ggplot(costello,
                           aes(x = FO, y = Pi, label = Prey)) +
  
  geom_point(size = 5, color = "steelblue") +
  
  geom_text(size = 5, nudge_y = 5) +
  
  geom_hline(yintercept = 50, linetype = 2, color = "gray60") +
  geom_vline(xintercept = 50, linetype = 2, color = "gray60") +
  
  labs(
    x = "Frequency of occurrence (%)",
    y = "Prey-specific abundance (%)"
  ) +
  
  theme_classic(base_size = 20)

ggsave("Figure_Costello.tiff",
       plot = grafico_costello,
       width = 8,
       height = 6,
       units = "in",
       dpi = 600,
       compression = "lzw")
