
# PACOTES -----------------------------------------------------------------

pacman::p_load(pacman, ggplot2, tidyverse, magick, gganimate, optimx)


# DADOS DE ENTRADA --------------------------------------------------------

# Hidrograma
df.hidrograma <-
  data.frame(t = seq(0, 70, 5)*60,                                        # tempo [s]
             q = c(c(560, 730, 1020, 1220, 930, 730, 600, 510, 455, 410), # vazão [m³/s]
                   rep(410, 5)))

# Parâmetros do canal retangular
b <- 20        # largura do canal [m]
l <- 10000     # comprimento do canal [m]
i0 <- 0.0016   # declividade de fundo do canal [m/m]
c.chezy <- 100 # coeficiente de rugosidade de Chézy
ql <- 0        # vazão contribuição lateral [m³/s.m]
t.sim <- 70*60 # tempo total de simulação [s]
g <- 9.81      # gravidade [m/s²]


# CONDIÇÕES DE FRONTEIRA --------------------------------------------------

# Montante: hidrograma de entrada 'df.hidrograma' e característica negativa
# Jusante: curva-chave do canal 'fun.curva.chave' e característica positiva


# DISCRETIZAÇÃO -----------------------------------------------------------

# Temporal
dt <- 10
n.dt <- t.sim/dt

# Espacial
dx <- 100
n.dx <- l/dx

# Discretizar hidrograma conforme dt adotado
df.hidrograma.disc <- fun.interpolar.hidrograma(df.hidrograma, dt)

# Critério de Courant
y.pico <- uniroot(f = fun.chezy,                     # Equação de Chézy
                  interval = c(4, 10),               # Limites de variação de y
                  c.chezy = c.chezy,                 # coeficiente de Chézy
                  vazao = max(df.hidrograma.disc$q), # Vazão de pico do hidrograma
                  largura = b,                       # Largura do canal
                  decl = i0,                         # Declividade de fundo do canal
                  tol = 1e-6, maxiter = 1e6          # Tolerância e número máximo de iterações
                  )$root; y.pico

v.pico <- max(df.hidrograma.disc$q)/(b*y.pico) # velocidade no pico do hidrograma
c.pico <- sqrt(g*y.pico)                       # celeridade no pico do hidrograma
dx.courant <- dt*(abs(v.pico) + c.pico)        # passo compatível c/ critério de Courant



# Testar dx adotado
if(dx >= dx.courant){
  cat("Passo adotado é compatível com o critério de Courant\n")
}else{
  dx <- dx.courant
  cat("Δx deve ser maior que ", dx.courant)
}

dx <- 200


# INICIAR VARIÁVEIS -------------------------------------------------------

# Profundidade e velocidade
v <- matrix(0, ncol = n.dx, nrow = n.dt) %>% as.data.frame # velocidade [m/s]
y <- matrix(0, ncol = n.dx, nrow = n.dt) %>% as.data.frame # profundidade [m]


# REGIME PERMANENTE (t = 0) ------------------------------------------

# Encontrar profundidade p/ vazão no tempo = 0
y[1,] <- uniroot(f = fun.chezy,                        # Equação de Chézy
                 interval = c(4, 10),                  # Limites de variação de y
                 c.chezy = c.chezy,                    # coeficiente de Chézy
                 vazao = max(df.hidrograma.disc[1,2]), # Vazão de pico do hidrograma
                 largura = b,                          # Largura do canal
                 decl = i0,                            # Declividade de fundo do canal
                 tol = 1e-6, maxiter = 1e6             # Tolerância e número máximo de iterações
                 )$root %>% as.numeric

# Determinar velocidade a partir de y0
v[1,] <- df.hidrograma.disc[1,2]/(b*y[1,1])


# REGIME TRANSIENTE (t > 0) -----------------------------------------------

# Simulação
for(t in 2:n.dt){
  
  # Contador
  cat(paste("\nProcessamento:", round((t-1)/n.dt*100, 2)), "%\n")
  cat(paste("Contador t =", t, "\n"))
  
  # Vazão no tempo inicial
  q <- df.hidrograma.disc[t,2]
  
  # Iterações ao longo do canal
  for(i in 1:n.dx){
    
    # Velocidade e profundidade no ponto M (instante t-1, ponto i)
    vm <- v[t-1, i] # velocidade
    ym <- y[t-1, i] # profundidade
    
    # Celeridade
    cm <- sqrt(g*ym)
    
    # Declividade da linha de energia
    ifm <- ((2*ym + b)/(b*ym))*(vm/c.chezy)^2
    
    # Condições de montante (hidrograma de entrada)
    if(i == 1){
      
      # Interpolação espacial considerando regime subcrítico
      vs <- vm + (vm - v[t-1, i+1])*dt/dx*(vm - cm)
      ys <- ym + (ym - y[t-1, i+1])*dt/dx*(vm - cm)
      
      # Determinar yp (característica negativa) sabendo que vp = q(yp*b)
      y[t,i] <- uniroot(f = fun.yp.montante,
                        interval = c(0.0001, 10000),
                        vazao = q, celeridade = cm, velocidade.s = vs, prof.s = ys, prof.m = ym,
                        largura = b, velocidade.m = vm, vazao.lateral = ql, decl.energia.m = ifm,
                        decl.canal = i0, dt = dt,
                        tol = 1e-6, maxiter = 1e6)$root %>% as.numeric
      
      # Determinar vp
      v[t,i] <- q/(y[t,i]*b)
      
    }
    
    # Condições internas ao longo da adutora
    if(i > 1 & i < n.dx){
      
      va <- v[t-1, i-1]
      ya <- y[t-1, i-1]
      
      # Interpolação linear espacial em R
      vr <- vm + (v[t-1, i-1] - vm)*dt/dx*(vm + cm)
      yr <- ym + (y[t-1, i-1] - ym)*dt/dx*(vm + cm)
      
      # Interpolação linear espacial em S → analisar o regime
      if(vm < cm){ # regime subcrítico
        vs <- vm + (vm - v[t-1, i+1])*dt/dx*(vm - cm)
        ys <- ym + (ym - y[t-1, i+1])*dt/dx*(vm - cm)
      }
      if(vm > cm){ # regime supercrítico
        vs <- vm + (v[t-1, i-1] - vm)*dt/dx*(vm - cm)
        ys <- ym + (y[t-1, i-1] - ym)*dt/dx*(vm - cm)
      }
      
      # Determinar yp usando as equações de C+ e C-, igualando vp
      y[t,i] <- (yr + ys)/2 + cm/(2*g)*(vr - vs)
      
      # Determinar vp
      v[t,i] <- vr - g/cm*(y[t,i] - yr) - ql/(b*ym)*(vm - cm)*dt - g*(ifm - i0)*dt
      
    }
    
    # Condições de jusante 
    if(i == n.dx){
      
      # Interpolação linear espacial em R
      vr <- vm + (v[t-1, i-1] - vm)*dt/dx*(vm + cm)
      yr <- ym + (y[t-1, i-1] - ym)*dt/dx*(vm + cm)
      
      # Determinar yp (característica positiva) e curva-chave do canal
      y[t,i] <- uniroot(f = fun.yp.jusante,
                        interval = c(0.3848, 10000),
                        celeridade = cm, velocidade.r = vr, prof.r = yr, prof.m = ym,
                        largura = b, velocidade.m = vm, vazao.lateral = ql, decl.energia.m = ifm,
                        decl.canal = i0, dt = dt,
                        tol = 1e-6, maxiter = 1e6)$root %>% as.numeric
      
      # Determinar vp
      v[t,i] <- 111.7118*(y[t,i]-0.3848)^1.2277/(b*y[t,i])
      
    }
  }
  
}

# Adicionar coluna de tempo ao data.frame final
t <- seq(from = 0, to = t.sim, by = dt) %>% head(-1)
v <- cbind(t, v)
y <- cbind(t, y)

# VISUALIZAÇÃO ------------------------------------------------------------

# Hidrograma de entrada
plot.hidrograma.in <- ({
  
  ggplot(data = df.hidrograma.disc) +
    geom_line(aes(x = t, y = q), color = "orange", linewidth = 0.5) +
    geom_segment(x = df.hidrograma$t[which.max(df.hidrograma$q)], y = 0,
                 xend = df.hidrograma$t[which.max(df.hidrograma$q)], yend = max(df.hidrograma$q),
                 color = "red", linetype = "dashed", linewidth = 0.3) +
    labs(x = "Tempo [s]", y = "Vazão [m³/s]", title = "Hidrograma de entrada") +
    # geom_smooth(method = "gam") +
    theme_minimal() +
    theme(panel.backgroundd = NULL,
          plot.background = element_rect(color = "white", fill = "white"),
          panel.border = element_rect(fill = NA),
          text = element_text(size = 10))
  
}); plot.hidrograma.in

# Visualização animada profundidade no tempo
# Transformar em formato long
y.long <- ({
  y %>% 
    pivot_longer(cols = -t,
                 names_to = "x",
                 values_to = "y") %>% 
    mutate(x = as.numeric(gsub("V", "", x))*dx,
           t.min = t/60)
})

y.min <- min(y.long$y, na.rm = TRUE)
y.max <- max(y.long$y, na.rm = TRUE)

plot.y <- ({
  ggplot(y.long, aes(x = x, y = y, group = t.min)) +
    geom_line(color = "steelblue", linewidth = 1.2) +
    geom_hline(yintercept = y[1,2], color = "red", alpha = 0.6,
               linewidth = 0.3, linetype = "dashed") +
    scale_y_continuous(limits = c(0, y.max)) +
    labs(x = "Comprimento do canal [m]", y = "Lâmina [m]",
         title = "Tempo: {round(frame_time)} minutos") +  # Usar {frame_time} para animação
    theme_minimal() +
    theme(panel.background = element_blank(),
          plot.background = element_rect(color = "white", fill = "white"),
          panel.border = element_rect(fill = NA),
          text = element_text(size = 10)) +
    transition_time(t.min) +  # Definir a animação com tempo em minutos
    ease_aes('linear')
})

# Salvar o GIF
anim_save("profundidade.gif", plot.y,
          fps = 10, width = 16, height = 12, units = "cm", res = 100)

# Visualização animada profundidade no tempo
# Transformar em formato long
v.long <- ({
  v %>% 
    pivot_longer(cols = -t,
                 names_to = "x",
                 values_to = "v") %>% 
    mutate(x = as.numeric(gsub("V", "", x))*dx,
           t.min = t/60)
})

v.min <- min(v.long$v, na.rm = TRUE)
v.max <- max(v.long$v, na.rm = TRUE)

plot.v <- ({
  ggplot(v.long, aes(x = x, y = v, group = t.min)) +
    geom_line(color = "steelblue", linewidth = 1.2) +
    geom_hline(yintercept = v[1,2], color = "red", alpha = 0.6,
               linewidth = 0.3, linetype = "dashed") +
    scale_y_continuous(limits = c(0, v.max)) +
    labs(x = "Comprimento do canal [m]", y = "Velocidade [m/s]",
         title = "Tempo: {round(frame_time)} minutos") +  # Usar {frame_time} para animação
    theme_minimal() +
    theme(panel.background = element_blank(),
          plot.background = element_rect(color = "white", fill = "white"),
          panel.border = element_rect(fill = NA),
          text = element_text(size = 10)) +
    transition_time(t.min) +  # Definir a animação com tempo em minutos
    ease_aes('linear')
})

# Salvar o GIF
anim_save("velocidade.gif", plot.v,
          fps = 10, width = 16, height = 12, units = "cm", res = 100)
