#' ---------------------------
#'
#' Purpose of script:
#'
#'  Download data
#' 
#' ---------------------------

GSE105420
wget "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE105420&format=file" -O data/GSE105420_RAW.tar

# Crear la carpeta si no existe
mkdir -p data/GSE105420

# Extraer el archivo en la carpeta
tar -xvf data/GSE105420_RAW.tar -C data/GSE105420

GSE105420 <- getGEO("GSE105420", GSEMatrix = TRUE)