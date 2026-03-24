---
title: "Informe exploratorio: GSE256394 cases"
author: "Carlos Ruiz"
date: "2026-03-23"
output:
  html_document:
    toc: true
    toc_float: true
    number_sections: true
---



# Configuracion editable


``` r
# Ruta al archivo .Rdata con el GenomicRatioSet
input_rdata <- "../results/preprocessing/GSE256394/GSE256394_cases.Rdata"

# Nombre del objeto dentro del .Rdata
object_name <- "GSE256394_cases"

# Variables principales (editar aqui en el futuro)
continuous_vars <- c()

categorical_vars <- c("Sex", "Progression", "mdsdiag", "SF3B1", "SRSF2", "U2AF1"
)

# Variables temporales candidatas (se usaran solo las que existan en colData)
temporal_vars <- c()

# PCA: numero de CpGs con mayor varianza para acelerar calculo
n_top_cpg_for_pca <- 5000

# Bins por defecto para histogramas
hist_bins <- 30
```

# Carga de datos


``` r
load(input_rdata)
```

```
## Error in `readChar()`:
## ! cannot open the connection
```

``` r
grset <- GSE105420_cases
```

```
## Error:
## ! object 'GSE105420_cases' not found
```

``` r
meta <- as.data.frame(SummarizedExperiment::colData(grset))
```

```
## Error in `h()`:
## ! error in evaluating the argument 'x' in selecting a method for function 'as.data.frame': error in evaluating the argument 'x' in selecting a method for function 'colData': object 'grset' not found
```

``` r
cat("Numero de muestras:", ncol(grset), "\n")
```

```
## Error in `h()`:
## ! error in evaluating the argument 'x' in selecting a method for function 'ncol': object 'grset' not found
```

``` r
cat("Numero de CpGs:", nrow(grset), "\n")
```

```
## Error in `h()`:
## ! error in evaluating the argument 'x' in selecting a method for function 'nrow': object 'grset' not found
```

``` r
cat("Columnas en colData:", ncol(meta), "\n")
```

```
## Error in `h()`:
## ! error in evaluating the argument 'x' in selecting a method for function 'ncol': object 'meta' not found
```


# Variables efectivas en el dataset


``` r
kable(data.frame(
  tipo = c(rep("continua", length(continuous_vars)),
           rep("categorica", length(categorical_vars)),
           rep("temporal", length(temporal_vars))),
  variable = c(continuous_vars, categorical_vars, temporal_vars)
), caption = "Variables que se analizaran")
```



Table: Variables que se analizaran

|tipo       |variable    |
|:----------|:-----------|
|categorica |Sex         |
|categorica |Progression |
|categorica |mdsdiag     |
|categorica |SF3B1       |
|categorica |SRSF2       |
|categorica |U2AF1       |

# Evaluacion univariante

## Variables continuas


``` r
for (v in continuous_vars) {
  cat("\n\n### ", v, "\n\n", sep = "")

  vec <- meta[[v]]
  df <- data.frame(value = vec)

  p <- ggplot(df, aes(x = value)) +
    geom_histogram(bins = hist_bins, fill = "#2C7FB8", color = "white", na.rm = TRUE) +
    theme_minimal(base_size = 12) +
    labs(title = paste("Histograma de", v), x = v, y = "Frecuencia")
  print(p)

  qs <- quantile(vec, probs = c(0, 0.25, 0.5, 0.75, 1), na.rm = TRUE)
  stat_tbl <- data.frame(
    estadistico = c("N total", "Min", "Q1", "Mediana", "Q3", "Max"),
    valor = c(sum(!is.na(vec)), unname(qs[1]), unname(qs[2]), unname(qs[3]), unname(qs[4]), unname(qs[5]))
  )
  print(kable(stat_tbl, caption = paste("Estadisticos de distribucion -", v)))
}
```

## Variables categoricas


``` r
for (v in categorical_vars) {
  cat("\n\n### ", v, "\n\n", sep = "")

  x <- as.factor(meta[[v]])
  tab <- as.data.frame(table(x, useNA = "ifany"), stringsAsFactors = FALSE)
  names(tab) <- c("categoria", "frecuencia")
  tab <- tab %>% mutate(proporcion = frecuencia / sum(frecuencia))

  p <- ggplot(tab, aes(x = categoria, y = frecuencia)) +
    geom_col(fill = "#E07A5F") +
    theme_minimal(base_size = 12) +
    labs(title = paste("Frecuencias de", v), x = v, y = "Frecuencia") +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
  print(p)

  print(kable(tab, caption = paste("Tabla de frecuencias -", v), digits = 3))
}
```

```
## 
## 
## ### Sex
```

```
## Error in `h()`:
## ! error in evaluating the argument 'x' in selecting a method for function 'as.factor': object 'meta' not found
```

## Variables temporales


``` r
for (v in temporal_vars) {
  cat("\n\n### ", v, "\n\n", sep = "")

  d <- meta[[v]]
  df <- data.frame(date = d)

  p <- ggplot(df, aes(x = date)) +
    geom_histogram(bins = hist_bins, fill = "#81B29A", color = "white", na.rm = TRUE) +
    theme_minimal(base_size = 12) +
    labs(title = paste("Distribucion temporal de", v), x = v, y = "Frecuencia")
  print(p)

  qs <- quantile(d, probs = c(0, 0.25, 0.5, 0.75, 1), na.rm = TRUE)
  stat_tbl <- data.frame(
    estadistico = c("N total", "Min", "Q1", "Mediana", "Q3", "Max"),
    valor = c(sum(!is.na(d)), as.character(qs[1]), as.character(qs[2]), as.character(qs[3]), as.character(qs[4]), as.character(qs[5]))
  )
  print(kable(stat_tbl, caption = paste("Estadisticos temporales -", v)))
}
```

# Evaluacion bivariante

## Continua vs continua


``` r
if (length(continuous_vars) >= 2) {
  cont_pairs <- combn(continuous_vars, 2, simplify = FALSE)

  for (pair in cont_pairs) {
    xvar <- pair[1]
    yvar <- pair[2]

    cat("\n\n### ", xvar, " vs ", yvar, "\n\n", sep = "")

    df <- data.frame(
      x = meta[[xvar]],
      y = meta[[yvar]]
    )
    df <- df %>% filter(!is.na(x), !is.na(y))

    p <- ggplot(df, aes(x = x, y = y)) +
      geom_point(alpha = 0.75, color = "#3D405B") +
      geom_smooth(method = "lm", se = TRUE, color = "#E63946") +
      theme_minimal(base_size = 12) +
      labs(title = paste("Scatterplot:", xvar, "vs", yvar), x = xvar, y = yvar)
    print(p)

    if (nrow(df) >= 3) {
      fit <- lm(y ~ x, data = df)
     fit_tab <- as.data.frame(summary(fit)$coefficients)
    colnames(fit_tab)[4] <- "p_value"
      print(kable(fit_tab, caption = "Coeficientes regresion lineal"))
    } else {
      cat("No hay suficientes datos para ajustar regresion lineal.\n")
    }
  }
}
```

## Continua vs categorica


``` r
for (xvar in continuous_vars) {
  for (gvar in categorical_vars) {
    cat("\n\n### ", xvar, " vs ", gvar, "\n\n", sep = "")

    df <- data.frame(
      x = meta[[xvar]],
      g = as.factor(meta[[gvar]])
    )
    df <- df %>% filter(!is.na(x), !is.na(g))

    p <- ggplot(df, aes(x = g, y = x, fill = g)) +
      geom_boxplot(alpha = 0.85, outlier.alpha = 0.5) +
      theme_minimal(base_size = 12) +
      labs(title = paste("Boxplot:", xvar, "por", gvar), x = gvar, y = xvar) +
      theme(axis.text.x = element_text(angle = 45, hjust = 1), legend.position = "none")
    print(p)

    df$g <- droplevels(df$g)
    if (nlevels(df$g) >= 2 && nrow(df) >= 5) {
      fit <- lm(x ~ g, data = df)
     fit_tab <- as.data.frame(summary(fit)$coefficients)
    colnames(fit_tab)[4] <- "p_value"
      print(kable(fit_tab, caption = "Coeficientes regresion lineal"))
    } else {
      cat("No hay suficientes datos/variacion para ajustar regresion lineal.\n")
    }
  }
}
```

## Categorica vs categorica


``` r
if (length(categorical_vars) >= 2) {
  cat_pairs <- combn(categorical_vars, 2, simplify = FALSE)

  for (pair in cat_pairs) {
    v1 <- pair[1]
    v2 <- pair[2]

    cat("\n\n### ", v1, " vs ", v2, "\n\n", sep = "")

    df <- data.frame(
      a = as.factor(meta[[v1]]),
      b = as.factor(meta[[v2]])
    ) %>% filter(!is.na(a), !is.na(b))

    if (nrow(df) == 0) {
      cat("No hay observaciones completas para esta comparacion.\n")
      next
    }

    tab <- table(df$a, df$b)
    tab_df <- as.data.frame(tab)
    names(tab_df) <- c(v1, v2, "frecuencia")

    p <- ggplot(tab_df, aes_string(x = v1, y = "frecuencia", fill = v2)) +
      geom_col(position = "dodge") +
      theme_minimal(base_size = 12) +
      labs(title = paste("Barras de frecuencias:", v1, "vs", v2), x = v1, y = "Frecuencia") +
      theme(axis.text.x = element_text(angle = 45, hjust = 1))
    print(p)

    print(kable(tab_df, caption = paste("Tabla de frecuencias -", v1, "vs", v2)))

    if (nrow(tab) >= 2 && ncol(tab) >= 2) {
      chisq_ok <- all(chisq.test(tab)$expected >= 5)
      chi <- suppressWarnings(chisq.test(tab))
      chi_tbl <- data.frame(
        estadistico = unname(chi$statistic),
        gl = unname(chi$parameter),
        p_value = chi$p.value,
        esperado_ge_5 = chisq_ok
      )
      print(kable(chi_tbl, caption = "Test Chi-cuadrado"))
    } else {
      cat("No hay dimensiones suficientes para test Chi-cuadrado.\n")
    }
  }
}
```

```
## 
## 
## ### Sex vs Progression
```

```
## Error in `h()`:
## ! error in evaluating the argument 'x' in selecting a method for function 'as.factor': object 'meta' not found
```

# Patrones globales de metilacion y relacion con variables principales


``` r
# Matriz beta CpGs x muestras
beta <- as.matrix(minfi::getBeta(grset))
```

```
## Error in `h()`:
## ! error in evaluating the argument 'x' in selecting a method for function 'as.matrix': error in evaluating the argument 'object' in selecting a method for function 'getBeta': object 'grset' not found
```

``` r
# Limpieza basica de filas no informativas
row_ok <- rowSums(is.na(beta)) < ncol(beta)
```

```
## Error in `base::rowSums()`:
## ! 'x' must be an array of at least two dimensions
```

``` r
beta <- beta[row_ok, , drop = FALSE]
```

```
## Error:
## ! object 'row_ok' not found
```

``` r
# Imputacion simple por media de CpG para evitar NA en PCA
if (anyNA(beta)) {
  row_means <- rowMeans(beta, na.rm = TRUE)
  idx <- which(is.na(beta), arr.ind = TRUE)
  beta[idx] <- row_means[idx[, 1]]
}
```

```
## Error:
## ! anyNA() applied to non-(list or vector) of type 'closure'
```

``` r
# Seleccion de CpGs mas variables para hacer PCA mas estable y rapido
vars <- apply(beta, 1, var)
```

```
## Error in `apply()`:
## ! dim(X) must have a positive length
```

``` r
ord <- order(vars, decreasing = TRUE)
```

```
## Error in `order()`:
## ! argument 1 is not a vector
```

``` r
ntop <- min(n_top_cpg_for_pca, length(ord))
```

```
## Error:
## ! object 'ord' not found
```

``` r
beta_top <- beta[ord[seq_len(ntop)], , drop = FALSE]
```

```
## Error:
## ! object 'ord' not found
```

``` r
# PCA sobre muestras
pca <- prcomp(t(beta_top), center = TRUE, scale. = TRUE)
```

```
## Error in `h()`:
## ! error in evaluating the argument 'x' in selecting a method for function 't': object 'beta_top' not found
```

``` r
pca_df <- as.data.frame(pca$x[, 1:10])
```

```
## Error in `h()`:
## ! error in evaluating the argument 'x' in selecting a method for function 'as.data.frame': object 'pca' not found
```

``` r
pca_df$sample <- rownames(pca_df)
```

```
## Error in `h()`:
## ! error in evaluating the argument 'x' in selecting a method for function 'rownames': object 'pca_df' not found
```

``` r
# Anadir metadata
meta$sample <- rownames(meta)
```

```
## Error in `h()`:
## ! error in evaluating the argument 'x' in selecting a method for function 'rownames': object 'meta' not found
```

``` r
pca_df <- pca_df %>% left_join(meta, by = "sample")
```

```
## Error:
## ! object 'pca_df' not found
```

``` r
var_expl <- summary(pca)$importance[2, 1:10] * 100
```

```
## Error in `h()`:
## ! error in evaluating the argument 'object' in selecting a method for function 'summary': object 'pca' not found
```

## PCA coloreado por variables principales


``` r
vars_for_pca_coloring <- unique(c(continuous_vars, categorical_vars, temporal_vars))

for (v in vars_for_pca_coloring) {
  cat("\n\n### PCA coloreado por ", v, "\n\n", sep = "")


  col_vec <- pca_df[[v]]
  p <- ggplot(pca_df, aes(x = PC1, y = PC2)) + theme_minimal(base_size = 12)

 if (v %in% categorical_vars){
    p <- p +
      geom_point(aes(color = as.factor(col_vec)), alpha = 0.85, size = 2.2) +
      labs(color = v)
  } else {
    p <- p +
      geom_point(aes(color = col_vec), alpha = 0.85, size = 2.2) +
      scale_color_viridis_c(option = "D", na.value = "grey80") +
      labs(color = v)
  }  

  p <- p +
    labs(
      title = paste("PCA de metilacion coloreado por", v),
      x = paste0("PC1 (", round(var_expl[1], 2), "%)"),
      y = paste0("PC2 (", round(var_expl[2], 2), "%)")
    )

  print(p)
}
```

```
## 
## 
## ### PCA coloreado por Sex
```

```
## Error:
## ! object 'pca_df' not found
```

## Asociacion PCA con variables principales (p-values)


``` r
vars_for_pca_assoc <- unique(c(continuous_vars, categorical_vars, temporal_vars))

assoc_results <- data.frame(
  variable = character(),
  tipo_variable = character(),
  componente = character(),
  n = integer(),
  n_niveles = integer(),
  p_value = numeric(),
  stringsAsFactors = FALSE
)

for (v in vars_for_pca_assoc) {
  if (!v %in% names(pca_df)) {
    next
  }

  vv <- pca_df[[v]]
 if (v %in% categorical_vars) {
        vv_fac <- as.factor(vv)

    for (pc in paste0("PC", 1:10)) {
      df <- data.frame(pc = pca_df[[pc]], var = vv_fac) %>%
        filter(!is.na(pc), !is.na(var))

      p_val <- NA_real_
      n_levels <- dplyr::n_distinct(df$var)
      if (nrow(df) >= 3 && n_levels >= 2) {
        fit <- lm(pc ~ var, data = df)
        an <- anova(fit)
        p_val <- an$`Pr(>F)`[1]
      }

      assoc_results <- bind_rows(
        assoc_results,
        data.frame(
          variable = v,
          tipo_variable = "categorica",
          componente = pc,
          n = nrow(df),
          n_niveles = n_levels,
          p_value = p_val,
          stringsAsFactors = FALSE
        )
      )
    }
  } else {
        for (pc in paste0("PC", 1:10)) {
      df <- data.frame(pc = pca_df[[pc]], var = vv) %>%
        filter(!is.na(pc), !is.na(var))

      p_val <- NA_real_
      if (nrow(df) >= 3 && dplyr::n_distinct(df$var) >= 2) {
        fit <- lm(pc ~ var, data = df)
        p_val <- summary(fit)$coefficients["var", "Pr(>|t|)"]
      }

      assoc_results <- bind_rows(
        assoc_results,
        data.frame(
          variable = v,
          tipo_variable = "continua",
          componente = pc,
          n = nrow(df),
          n_niveles = NA_integer_,
          p_value = p_val,
          stringsAsFactors = FALSE
        )
      )
    }
    
  }
}
```

```
## Error in `h()`:
## ! error in evaluating the argument 'table' in selecting a method for function '%in%': object 'pca_df' not found
```

``` r
if (nrow(assoc_results) == 0) {
  cat("No fue posible calcular asociaciones PCA con las variables definidas.\n")
} else {
  assoc_results <- assoc_results %>%
    mutate(
      p_adj_bh = p.adjust(p_value, method = "BH"),
      p_value = signif(p_value, 4),
      p_adj_bh = signif(p_adj_bh, 4)
    ) %>%
    arrange(componente, p_value)

  print(kable(
    assoc_results,
    caption = "Asociacion estadistica de PC1/PC2 con variables principales"
  ))
}
```

```
## No fue posible calcular asociaciones PCA con las variables definidas.
```

