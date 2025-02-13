
# PACOTES -----------------------------------------------------------------

pacman::p_load(pacman, ggplot2, tidyverse, magick, optimx)


# DADOS DE ENTRADA --------------------------------------------------------

# Hidrograma de cheia
df.hidrograma <- data.frame(t = seq(0, 70, 5)*60,                                                      # tempo [min]
                            q = c(c(560, 730, 1020, 1220, 930, 730, 600, 510, 455, 410), rep(410, 5))) # vazão [m³/s]

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

# Espacial
dx <- 250
n.dx <- l/dx

# Temporal
dt <- 60*1
n.dt <- t.sim/dt


# SIMULAÇÃO ---------------------------------------------------------------

# Discretizar hidrograma conforme dt adotado
df.hidrograma.disc <- fun.interpolar.hidrograma(df.hidrograma, dt)

# Assumir uma profundidade inicial y.chute
y.chute <- 4

# Critério de estabilidade de Courant
# Checar p/ pico

# Iniciar variáveis Velocidade 'v' e Profundidade 'y'
v <- matrix(0, ncol = n.dx, nrow = n.dt) %>% as.data.frame # velocidade [m/s]
y <- matrix(0, ncol = n.dx, nrow = n.dt) %>% as.data.frame # profundidade [m]

# Otimizar para encontrar y no regime permanente
suppressWarnings(
  res.optimx <-
    optimx::optimx(par = c(y.chute),                                      # parâmetro a ser otimizado - 'y'
                   fn = fun.chezy,                                        # função objetivo
                   lower = 0.0001, upper = Inf,                           # limites
                   vazao = df.hidrograma.disc[1,2], largura = b, decl = i0, c.chezy = c.chezy)) # outros parâmetros da função objetivo

# Preencher y no regime permanente
y0 <- res.optimx$p1 %>% as.numeric
y[1,] <- y0

# Calcular velocidade no regime permanente
v0 <- df.hidrograma.disc[1,2]/(y0*b)
v[1,] <- v0

# Simulação
for(t in 2:n.dt){
  
  # Contador
  cat(paste("\nProcessamento:", round((t-1)/n.dt*100, 2), "%\n"))
  cat(paste("Contador t =", t, "\n"))
  
  # Vazão no tempo inicial
  q <- df.hidrograma.disc[t,2]
  
  # Condição para regime transitório
  # Iterações ao longo da adutora
  for(i in 1:n.dx){
    
    # Velocidade e profundidade no ponto M (instante t-1, ponto i)
    vm <- v[t-1, i] # velocidade
    ym <- y[t-1, i] # profundidade
    
    # Celeridade 'cm' [m/s]
    cm <- sqrt(g*ym)
    
    # Declividade da linha de energia
    rhm <- ym*b/(2*ym + b) # raio hidráulico no ponto Ms
    ifm <- 1/rhm*(vm/c.chezy)^2
    
    # Condições de montante (hidrograma de entrada)
    if(i == 1){
      
      # Velocidade e profundidade no ponto B (instante t-1, ponto i+1)
      vb <- v[t-1, i+1]
      yb <- y[t-1, i+1]
      
      # Considerando que o regime é subcrítico
      vs <- vm + (vm - vb)*dt/dx*(vm - cm)
      ys <- ym + (ym - yb)*dt/dx*(vm - cm)
      
      # Determinar yp (característica negativa) sabendo que vp = q/(yp*b)
      vp.root.up <- uniroot(f = fun.yp.montante,         # função objetivo p/ determinar yp
                            interval = c(0.0001, 10000), # limites
                            vazao = q, celeridade = cm, velocidade.s = vs, prof.s = ys, prof.m = ym,
                            largura = b, velocidade.m = vm, vazao.lateral = ql, decl.energia.m = ifm,
                            decl.canal = i0, dt = dt,
                            tol = 1e-6, maxiter = 10e6)
      
      # Extrair yp e determinar vp
      y[t,i] <- vp.root.up$root %>% as.numeric
      v[t,i] <- q/(y[t,i]*b)
      
    }
    
    # Condições ao longo da adutora
    if(i > 1 && i < n.dx){

      # Velocidade e profundidade no ponto A (instante t-1, ponto i-1)
      va <- v[t-1, i-1]
      ya <- y[t-1, i-1]
      
      # Velocidade e profundidade no ponto B (instante t-1, ponto i+1)
      vb <- v[t-1, i+1]
      yb <- y[t-1, i+1]

      # Interpolações lineares espaciais
      vr <- vm + (va - vm)*dt/dx*(vm + cm)
      yr <- ym + (ya - ym)*dt/dx*(vm + cm)

      # Regime supercrítico → influenciado pelas condições de montante (ponto A)
      if(vm < cm){
        vs <- vm + (va - vm)*dt/dx*(vm - cm)
        ys <- ym + (ya - ym)*dt/dx*(vm - cm)
      }

      # Regime subcrítico → influenciado pelas condições de jusante (ponto B)
      if(vm > cm){
        vs <- vm + (vm - vb)*dt/dx*(vm - cm)
        ys <- ym + (ym - yb)*dt/dx*(vm - cm)
      }

      # Determinar yp subsituindo C- em C+, igualando vp
      y[t,i] <- (cm*(vr - vs)/g + yr + ys)/2

      # Determinar vp usando C+
      v[t,i] <- vs + g/cm*(y[t,i] - ys) - ql/(b*ym)*(vm - cm)*dt - g*(ifm - i0)*dt

    }

    # Condições de jusante (curva-chave)
    if(i == n.dx){
      
      # Velocidade e profundidade no ponto A (instante t-1, ponto i-1)
      va <- v[t-1, i-1]
      ya <- y[t-1, i-1]
      
      # Interpolações lineares espaciais
      vr <- vm + (va - vm)*dt/dx*(vm + cm)
      yr <- ym + (ya - ym)*dt/dx*(vm + cm)
      
      vp.root.dw <- uniroot(f = fun.yp.jusante,          # função objetivo p/ determinar yp
                            interval = c(0.3848, 10000), # limites
                            vazao = q, celeridade = cm, velocidade.r = vr, prof.r = yr, prof.m = ym,
                            largura = b, velocidade.m = vm, vazao.lateral = ql, decl.energia.m = ifm,
                            decl.canal = i0, dt = dt,
                            tol = 1e-6, maxiter = 10e6)
      
      # Extrair yp e determinar vp
      y[t,i] <- vp.root.dw$root %>% as.numeric
      v[t,i] <- 111.7118*(y[t,i] - 0.3848)^1.2277/(b*y[t,i])

    }
    
  }
  
}

# VISUALIZAÇÃO ------------------------------------------------------------

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

yp <- ym
111.7118*(yp - 0.3848)^1.2277/(b*yp)

n.dt <- 6

